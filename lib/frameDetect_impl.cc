/* -*- c++ -*- */
/*
 * Copyright 2024 @dossam.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <algorithm>
#include <vector>
#include <gnuradio/io_signature.h>
#include <volk/volk_alloc.hh>
#include "frameDetect_impl.h"

namespace gr {
    namespace loraJam {

        frameDetect::sptr frameDetect::make(uint32_t center_freq,
                                            uint32_t bandwidth,
                                            uint8_t sf,
                                            std::vector<uint16_t> sync_word,
                                            uint8_t os_factor,
                                            uint16_t preamb_len)
        {
            return gnuradio::make_block_sptr<frameDetect_impl>(
                center_freq, bandwidth, sf, sync_word, os_factor, preamb_len);
        }


        /*
        * The private constructor
        */
        frameDetect_impl::frameDetect_impl(uint32_t center_freq,
                                        uint32_t bandwidth,
                                        uint8_t sf,
                                        std::vector<uint16_t> sync_word,
                                        uint8_t os_factor,
                                        uint16_t preamb_len)
            : gr::block("frameDetect",
                        gr::io_signature::make(1, 1, sizeof(gr_complex)),
                        gr::io_signature::make(0, 0, 0))
        {
            m_state = DETECT;
            m_center_freq = center_freq;
            m_bw = bandwidth;
            m_sf = sf;

            m_sync_words = sync_word;
            m_os_factor = os_factor;
            if (preamb_len < 5)
            {
                std::cerr << RED << " Preamble length should be greater than 5!" << RESET << std::endl;
            }

            m_preamb_len = preamb_len;
            net_ids.resize(2, 0);

            m_n_up_req = preamb_len - 3;
            up_symb_to_use = m_n_up_req - 1;

            m_sto_frac = 0.0;

            // Convert given sync word into the two modulated values in preamble
            if (m_sync_words.size() == 1)
            {
                uint16_t tmp = m_sync_words[0];
                m_sync_words.resize(2, 0);
                m_sync_words[0] = ((tmp & 0xF0) >> 4) << 3;
                m_sync_words[1] = (tmp & 0x0F) << 3;
            }

            m_number_of_bins = (uint32_t)(1u << m_sf);
            m_samples_per_symbol = m_number_of_bins * m_os_factor;
            additional_symbol_samp.resize(2 * m_samples_per_symbol);
            m_upchirp.resize(m_number_of_bins);
            m_downchirp.resize(m_number_of_bins);
            preamble_upchirps.resize(m_preamb_len * m_number_of_bins);
            preamble_raw_up.resize((m_preamb_len + 3) * m_samples_per_symbol);
            CFO_frac_correc.resize(m_number_of_bins);
            CFO_SFO_frac_correc.resize(m_number_of_bins);
            symb_corr.resize(m_number_of_bins);
            in_down.resize(m_number_of_bins);
            preamble_raw.resize(m_preamb_len * m_number_of_bins);
            net_id_samp.resize(m_samples_per_symbol * 2.5); // we should be able to move up to one quarter of symbol in each direction

            cor_downchirp.resize(m_number_of_bins);
            cor_upchirp.resize(m_number_of_bins);

            build_ref_chirps(&m_upchirp[0], &m_downchirp[0], m_sf);

            bin_idx = 0;
            symbol_cnt = 1;
            k_hat = 0;
            preamb_up_vals.resize(m_n_up_req, 0);
            frame_cnt = 0;

            cx_in = new kiss_fft_cpx[m_number_of_bins];
            cx_out = new kiss_fft_cpx[m_number_of_bins];
            // register message ports
            message_port_register_out(pmt::intern("detect_out"));   // detector output

        out_dict = pmt::make_dict();
        out_dict = pmt::dict_add(out_dict, pmt::intern("SF"), pmt::mp(m_sf));      // this should be done only once, normally
        out_dict = pmt::dict_add(out_dict, pmt::intern("CR"), pmt::mp(0));         // 0 is the value for 'unknown' coding rate; probably implicit mode or corrupted header
        
        samp_buf.resize(8*m_number_of_bins);
        samp_buf_offset = 0;

        }

        /*
        * Our virtual destructor.
        */
        frameDetect_impl::~frameDetect_impl() {}

        int frameDetect_impl::my_roundf(float number)
        {
            int ret_val = (int)(number > 0 ? int(number + 0.5) : std::ceil(number - 0.5));
            return ret_val;
        }
        void frameDetect_impl::forecast(int noutput_items, gr_vector_int &ninput_items_required)
        {
            ninput_items_required[0] = (m_os_factor * (m_number_of_bins + 2));
        }

        float frameDetect_impl::estimate_CFO_frac(gr_complex *samples)
        {
            int k0;
            float cfo_frac;
            double Y_1, Y0, Y1, u, v, ka, wa, k_residual;
            std::vector<gr_complex> CFO_frac_correc_aug(up_symb_to_use * m_number_of_bins); ///< CFO frac correction vector
            std::vector<gr_complex> dechirped(up_symb_to_use * m_number_of_bins);
            kiss_fft_cpx *cx_in_cfo = new kiss_fft_cpx[2 * up_symb_to_use * m_number_of_bins];
            kiss_fft_cpx *cx_out_cfo = new kiss_fft_cpx[2 * up_symb_to_use * m_number_of_bins];

            std::vector<float> fft_mag_sq(2 * up_symb_to_use * m_number_of_bins);
            kiss_fft_cfg cfg_cfo = kiss_fft_alloc(2 * up_symb_to_use * m_number_of_bins, 0, 0, 0);
            // create longer downchirp
            std::vector<gr_complex> downchirp_aug(up_symb_to_use * m_number_of_bins);
            for (int i = 0; i < up_symb_to_use; i++)
            {
                memcpy(&downchirp_aug[i * m_number_of_bins], &m_downchirp[0], m_number_of_bins * sizeof(gr_complex));
            }

            // Dechirping
            volk_32fc_x2_multiply_32fc(&dechirped[0], samples, &downchirp_aug[0], up_symb_to_use * m_number_of_bins);
            // prepare FFT
            for (int i = 0; i < 2 * up_symb_to_use * m_number_of_bins; i++)
            {
                if (i < up_symb_to_use * m_number_of_bins)
                {
                    cx_in_cfo[i].r = dechirped[i].real();
                    cx_in_cfo[i].i = dechirped[i].imag();
                }
                else
                { // add padding
                    cx_in_cfo[i].r = 0;
                    cx_in_cfo[i].i = 0;
                }
            }
            // do the FFT
            kiss_fft(cfg_cfo, cx_in_cfo, cx_out_cfo);
            // Get magnitude
            for (uint32_t i = 0u; i < 2 * up_symb_to_use * m_number_of_bins; i++)
            {
                fft_mag_sq[i] = cx_out_cfo[i].r * cx_out_cfo[i].r + cx_out_cfo[i].i * cx_out_cfo[i].i;
            }
            free(cfg_cfo);
            delete[] cx_in_cfo;
            delete[] cx_out_cfo;
            // get argmax here
            k0 = std::distance(std::begin(fft_mag_sq), std::max_element(std::begin(fft_mag_sq), std::end(fft_mag_sq)));

            // get three spectral lines
            Y_1 = fft_mag_sq[mod(k0 - 1, 2 * up_symb_to_use * m_number_of_bins)];
            Y0 = fft_mag_sq[k0];
            Y1 = fft_mag_sq[mod(k0 + 1, 2 * up_symb_to_use * m_number_of_bins)];
            // set constant coeff
            u = 64 * m_number_of_bins / 406.5506497; // from Cui yang (15)
            v = u * 2.4674;
            // RCTSL
            wa = (Y1 - Y_1) / (u * (Y1 + Y_1) + v * Y0);
            ka = wa * m_number_of_bins / M_PI;
            k_residual = fmod((k0 + ka) / 2 / up_symb_to_use, 1);
            cfo_frac = k_residual - (k_residual > 0.5 ? 1 : 0);
            // Correct CFO frac in preamble
            for (int n = 0; n < up_symb_to_use * m_number_of_bins; n++)
            {
                CFO_frac_correc_aug[n] = gr_expj(-2 * M_PI * (cfo_frac) / m_number_of_bins * n);
            }

            volk_32fc_x2_multiply_32fc(&preamble_upchirps[0], samples, &CFO_frac_correc_aug[0], up_symb_to_use * m_number_of_bins);

            return cfo_frac;
        }
        float frameDetect_impl::estimate_CFO_frac_Bernier(gr_complex *samples)
        {
            std::vector<int> k0(up_symb_to_use);
            float cfo_frac;
            std::vector<gr_complex> CFO_frac_correc_aug(up_symb_to_use * m_number_of_bins); ///< CFO frac correction vector
            std::vector<double> k0_mag(up_symb_to_use);
            std::vector<gr_complex> fft_val(up_symb_to_use * m_number_of_bins);

            std::vector<gr_complex> dechirped(m_number_of_bins);
            kiss_fft_cpx *cx_in_cfo = new kiss_fft_cpx[m_number_of_bins];
            kiss_fft_cpx *cx_out_cfo = new kiss_fft_cpx[m_number_of_bins];
            std::vector<float> fft_mag_sq(m_number_of_bins);
            for (size_t i = 0; i < m_number_of_bins; i++)
            {
                fft_mag_sq[i] = 0;
            }
            kiss_fft_cfg cfg_cfo = kiss_fft_alloc(m_number_of_bins, 0, 0, 0);
            for (int i = 0; i < up_symb_to_use; i++)
            {
                // Dechirping
                volk_32fc_x2_multiply_32fc(&dechirped[0], &samples[m_number_of_bins * i], &m_downchirp[0], m_number_of_bins);
                // prepare FFT
                for (int j = 0; j < m_number_of_bins; j++)
                {
                    cx_in_cfo[j].r = dechirped[j].real();
                    cx_in_cfo[j].i = dechirped[j].imag();
                }
                // do the FFT
                kiss_fft(cfg_cfo, cx_in_cfo, cx_out_cfo);
                // Get magnitude

                for (uint32_t j = 0u; j < m_number_of_bins; j++)
                {
                    fft_mag_sq[j] = cx_out_cfo[j].r * cx_out_cfo[j].r + cx_out_cfo[j].i * cx_out_cfo[j].i;
                    fft_val[j + i * m_number_of_bins] = gr_complex(cx_out_cfo[j].r, cx_out_cfo[j].i);
                }
                k0[i] = std::distance(std::begin(fft_mag_sq), std::max_element(std::begin(fft_mag_sq), std::end(fft_mag_sq)));

                k0_mag[i] = fft_mag_sq[k0[i]];
            }
            free(cfg_cfo);
            delete[] cx_in_cfo;
            delete[] cx_out_cfo;
            // get argmax
            int idx_max = k0[std::distance(std::begin(k0_mag), std::max_element(std::begin(k0_mag), std::end(k0_mag)))];
            gr_complex four_cum(0.0f, 0.0f);
            for (int i = 0; i < up_symb_to_use - 1; i++)
            {
                four_cum += fft_val[idx_max + m_number_of_bins * i] * std::conj(fft_val[idx_max + m_number_of_bins * (i + 1)]);
            }
            cfo_frac = -std::arg(four_cum) / 2 / M_PI;
            // Correct CFO in preamble
            for (int n = 0; n < up_symb_to_use * m_number_of_bins; n++)
            {
                CFO_frac_correc_aug[n] = gr_expj(-2 * M_PI * cfo_frac / m_number_of_bins * n);
            }
            volk_32fc_x2_multiply_32fc(&preamble_upchirps[0], samples, &CFO_frac_correc_aug[0], up_symb_to_use * m_number_of_bins);
            return cfo_frac;
        }

        float frameDetect_impl::estimate_STO_frac()
        {
            int k0;
            double Y_1, Y0, Y1, u, v, ka, wa, k_residual;
            float sto_frac = 0;

            std::vector<gr_complex> dechirped(m_number_of_bins);
            kiss_fft_cpx *cx_in_sto = new kiss_fft_cpx[2 * m_number_of_bins];
            kiss_fft_cpx *cx_out_sto = new kiss_fft_cpx[2 * m_number_of_bins];

            std::vector<float> fft_mag_sq(2 * m_number_of_bins);
            for (size_t i = 0; i < 2 * m_number_of_bins; i++)
            {
                fft_mag_sq[i] = 0;
            }
            kiss_fft_cfg cfg_sto = kiss_fft_alloc(2 * m_number_of_bins, 0, 0, 0);

            for (int i = 0; i < up_symb_to_use; i++)
            {
                // Dechirping
                volk_32fc_x2_multiply_32fc(&dechirped[0], &preamble_upchirps[m_number_of_bins * i], &m_downchirp[0], m_number_of_bins);

                // prepare FFT
                for (int j = 0; j < 2 * m_number_of_bins; j++)
                {
                    if (j < m_number_of_bins)
                    {
                        cx_in_sto[j].r = dechirped[j].real();
                        cx_in_sto[j].i = dechirped[j].imag();
                    }
                    else
                    { // add padding
                        cx_in_sto[j].r = 0;
                        cx_in_sto[j].i = 0;
                    }
                }
                // do the FFT
                kiss_fft(cfg_sto, cx_in_sto, cx_out_sto);
                // Get magnitude
                for (uint32_t j = 0u; j < 2 * m_number_of_bins; j++)
                {
                    fft_mag_sq[j] += cx_out_sto[j].r * cx_out_sto[j].r + cx_out_sto[j].i * cx_out_sto[j].i;
                }
            }
            free(cfg_sto);
            delete[] cx_in_sto;
            delete[] cx_out_sto;

            // get argmax here
            k0 = std::distance(std::begin(fft_mag_sq), std::max_element(std::begin(fft_mag_sq), std::end(fft_mag_sq)));

            // get three spectral lines
            Y_1 = fft_mag_sq[mod(k0 - 1, 2 * m_number_of_bins)];
            Y0 = fft_mag_sq[k0];
            Y1 = fft_mag_sq[mod(k0 + 1, 2 * m_number_of_bins)];

            // set constant coeff
            u = 64 * m_number_of_bins / 406.5506497; // from Cui yang (eq.15)
            v = u * 2.4674;
            // RCTSL
            wa = (Y1 - Y_1) / (u * (Y1 + Y_1) + v * Y0);
            ka = wa * m_number_of_bins / M_PI;
            k_residual = fmod((k0 + ka) / 2, 1);
            sto_frac = k_residual - (k_residual > 0.5 ? 1 : 0);

            return sto_frac;
        }

        uint32_t frameDetect_impl::get_symbol_val(const gr_complex *samples, gr_complex *ref_chirp)
        {
            double sig_en = 0;
            std::vector<float> fft_mag(m_number_of_bins);
            volk::vector<gr_complex> dechirped(m_number_of_bins);

            kiss_fft_cfg cfg = kiss_fft_alloc(m_number_of_bins, 0, 0, 0);

            // Multiply with ideal downchirp
            volk_32fc_x2_multiply_32fc(&dechirped[0], samples, ref_chirp, m_number_of_bins);

            for (int i = 0; i < m_number_of_bins; i++)
            {
                cx_in[i].r = dechirped[i].real();
                cx_in[i].i = dechirped[i].imag();
            }
            // do the FFT
            kiss_fft(cfg, cx_in, cx_out);

            // Get magnitude
            for (uint32_t i = 0u; i < m_number_of_bins; i++)
            {
                fft_mag[i] = cx_out[i].r * cx_out[i].r + cx_out[i].i * cx_out[i].i;
                sig_en += fft_mag[i];
            }
            free(cfg);
            // Return argmax here

            return sig_en ? (std::distance(std::begin(fft_mag), std::max_element(std::begin(fft_mag), std::end(fft_mag)))) : -1;
        }

        float frameDetect_impl::determine_energy(const gr_complex *samples, int length = 1)
        {
            volk::vector<float> magsq_chirp(m_number_of_bins * length);
            float energy_chirp = 0;
            volk_32fc_magnitude_squared_32f(&magsq_chirp[0], samples, m_number_of_bins * length);
            volk_32f_accumulator_s32f(&energy_chirp, &magsq_chirp[0], m_number_of_bins * length);
            return energy_chirp / m_number_of_bins / length;
        }
        float frameDetect_impl::determine_snr(const gr_complex *samples)
        {
            double tot_en = 0;
            std::vector<float> fft_mag(m_number_of_bins);
            std::vector<gr_complex> dechirped(m_number_of_bins);

            kiss_fft_cfg cfg = kiss_fft_alloc(m_number_of_bins, 0, 0, 0);

            // Multiply with ideal downchirp
            volk_32fc_x2_multiply_32fc(&dechirped[0], samples, &m_downchirp[0], m_number_of_bins);

            for (int i = 0; i < m_number_of_bins; i++)
            {
                cx_in[i].r = dechirped[i].real();
                cx_in[i].i = dechirped[i].imag();
            }
            // do the FFT
            kiss_fft(cfg, cx_in, cx_out);

            // Get magnitude
            for (uint32_t i = 0u; i < m_number_of_bins; i++)
            {
                fft_mag[i] = cx_out[i].r * cx_out[i].r + cx_out[i].i * cx_out[i].i;
                tot_en += fft_mag[i];
            }
            free(cfg);

            int max_idx = std::distance(std::begin(fft_mag), std::max_element(std::begin(fft_mag), std::end(fft_mag)));
            float sig_en = fft_mag[max_idx];
            return 10 * log10(sig_en / (tot_en - sig_en));
        }

        int frameDetect_impl::general_work(int noutput_items,
                                        gr_vector_int& ninput_items,
                                        gr_vector_const_void_star& input_items,
                                        gr_vector_void_star& output_items)
        {
            const gr_complex *in = (const gr_complex *)input_items[0];
            int items_to_output = 0;

            // check if there is enough space in the output buffer
            if (ninput_items[0] < m_os_factor * (m_number_of_bins+2))
            {
                std::cout << "frameDetect->general_work: not enough samples" << std::endl;
                return 0;
            }

            int nitems_to_process = ninput_items[0];

            // downsampling
            for (int ii = 0; ii < m_number_of_bins; ii++)
                in_down[ii] = in[(int)(m_os_factor / 2 + m_os_factor * ii - my_roundf(m_sto_frac * m_os_factor))];

            switch (m_state)
            {
            case DETECT:
            {
                bin_idx_new = get_symbol_val(&in_down[0], &m_downchirp[0]);

                if (abs(mod(abs(bin_idx_new - bin_idx) + 1, m_number_of_bins) - 1) <= 1 && bin_idx_new != -1) // look for consecutive reference upchirps(with a margin of ±1)
                {   //// 2nd reference upchirp, or 3rd or 4th, ...
                    if (symbol_cnt == 1 && bin_idx != -1)   //// 2nd upchirp
                        preamb_up_vals[0] = bin_idx;

                    preamb_up_vals[symbol_cnt] = bin_idx_new;
                    memcpy(&preamble_raw[m_number_of_bins * symbol_cnt], &in_down[0], m_number_of_bins * sizeof(gr_complex));
                    memcpy(&preamble_raw_up[m_samples_per_symbol * symbol_cnt], &in[(int)(m_os_factor / 2)], m_samples_per_symbol * sizeof(gr_complex));

                    symbol_cnt++;
                }
                else //// First reference upchirp
                {
                    memcpy(&preamble_raw[0], &in_down[0], m_number_of_bins * sizeof(gr_complex));
                    memcpy(&preamble_raw_up[0], &in[(int)(m_os_factor / 2)], m_samples_per_symbol * sizeof(gr_complex));

                    symbol_cnt = 1;
                }
                bin_idx = bin_idx_new;
                if (symbol_cnt == (int)(m_n_up_req))    //// REached the required number of reference upchirps
                {
                    additional_upchirps = 0;
                    m_state = SYNC;
                    symbol_cnt = 0;
                    cfo_frac_sto_frac_est = false;
                    k_hat = most_frequent(&preamb_up_vals[0], preamb_up_vals.size());
                    memcpy(&net_id_samp[0], &in[int(0.75 * m_samples_per_symbol - k_hat * m_os_factor)], sizeof(gr_complex) * 0.25 * m_samples_per_symbol);

                    // perform the coarse synchronization
                    items_to_consume = m_os_factor * ((int)(m_number_of_bins - k_hat));

                }
                else
                    items_to_consume = m_samples_per_symbol;
                items_to_output = 0;
                break;
            }
            case SYNC:
            {
                items_to_output = 0;
                if (!cfo_frac_sto_frac_est)
                {
                    m_cfo_frac = estimate_CFO_frac_Bernier(&preamble_raw[m_number_of_bins - k_hat]);
                    m_sto_frac = estimate_STO_frac();
                    // create correction vector
                    for (int n = 0; n < m_number_of_bins; n++)
                    {
                        CFO_frac_correc[n] = gr_expj(-2 * M_PI * m_cfo_frac / m_number_of_bins * n);
                    }
                    cfo_frac_sto_frac_est = true;
                }
                items_to_consume = m_samples_per_symbol;
                // apply cfo correction
                volk_32fc_x2_multiply_32fc(&symb_corr[0], &in_down[0], &CFO_frac_correc[0], m_number_of_bins);

                bin_idx = get_symbol_val(&symb_corr[0], &m_downchirp[0]);
                switch (symbol_cnt)
                {
                case NET_ID1:
                {
                    if (bin_idx == 0 || bin_idx == 1 || bin_idx == m_number_of_bins - 1)
                    { // look for additional upchirps. Won't work if network identifier 1 equals 2^sf-1, 0 or 1!
                        memcpy(&net_id_samp[0], &in[(int)0.75 * m_samples_per_symbol], sizeof(gr_complex) * 0.25 * m_samples_per_symbol);
                        if (additional_upchirps >= 3)
                        {
                            std::rotate(preamble_raw_up.begin(), preamble_raw_up.begin() + m_samples_per_symbol, preamble_raw_up.end());
                            memcpy(&preamble_raw_up[m_samples_per_symbol * (m_n_up_req + 3)], &in[(int)(m_os_factor / 2) + k_hat * m_os_factor], m_samples_per_symbol * sizeof(gr_complex));
                        }
                        else
                        {
                            memcpy(&preamble_raw_up[m_samples_per_symbol * (m_n_up_req + additional_upchirps)], &in[(int)(m_os_factor / 2) + k_hat * m_os_factor], m_samples_per_symbol * sizeof(gr_complex));
                            additional_upchirps++;
                        }
                    }
                    else
                    { // network identifier 1 correct or off by one
                        symbol_cnt = NET_ID2;
                        memcpy(&net_id_samp[0.25 * m_samples_per_symbol], &in[0], sizeof(gr_complex) * m_samples_per_symbol);
                        net_ids[0] = bin_idx;
                    }
                    break;
                }
                case NET_ID2:
                {
                    symbol_cnt = DOWNCHIRP1;
                    memcpy(&net_id_samp[1.25 * m_samples_per_symbol], &in[0], sizeof(gr_complex) * (m_number_of_bins + 1) * m_os_factor);
                    net_ids[1] = bin_idx;

                    break;
                }
                case DOWNCHIRP1:
                {
                    memcpy(&net_id_samp[2.25 * m_samples_per_symbol], &in[0], sizeof(gr_complex) * 0.25 * m_samples_per_symbol);
                    symbol_cnt = DOWNCHIRP2;
                    break;
                }
                case DOWNCHIRP2:
                {
                    down_val = get_symbol_val(&symb_corr[0], &m_upchirp[0]);
                    memcpy(&additional_symbol_samp[0], &in[0], sizeof(gr_complex) * m_samples_per_symbol);
                    symbol_cnt = QUARTER_DOWN;
                    break;
                }
                case QUARTER_DOWN:
                {
                    memcpy(&additional_symbol_samp[m_samples_per_symbol], &in[0], sizeof(gr_complex) * m_samples_per_symbol);
                    if (down_val < m_number_of_bins / 2)
                    {
                        m_cfo_int = floor(down_val / 2);
                    }
                    else
                    {
                        m_cfo_int = floor(double(down_val - (int)m_number_of_bins) / 2);
                    }

                    // correct STOint and CFOint in the preamble upchirps
                    std::rotate(preamble_upchirps.begin(), preamble_upchirps.begin() + mod(m_cfo_int, m_number_of_bins), preamble_upchirps.end());

                    std::vector<gr_complex> CFO_int_correc;
                    CFO_int_correc.resize((m_n_up_req + additional_upchirps) * m_number_of_bins);
                    for (int n = 0; n < (m_n_up_req + additional_upchirps) * m_number_of_bins; n++)
                    {
                        CFO_int_correc[n] = gr_expj(-2 * M_PI * (m_cfo_int) / m_number_of_bins * n);
                    }

                    volk_32fc_x2_multiply_32fc(&preamble_upchirps[0], &preamble_upchirps[0], &CFO_int_correc[0], up_symb_to_use * m_number_of_bins);

                    // correct SFO in the preamble upchirps

                    sfo_hat = float((m_cfo_int + m_cfo_frac) * m_bw) / m_center_freq;
                    double clk_off = sfo_hat / m_number_of_bins;
                    double fs = m_bw;
                    double fs_p = m_bw * (1 - clk_off);
                    int N = m_number_of_bins;
                    std::vector<gr_complex> sfo_corr_vect;
                    sfo_corr_vect.resize((m_n_up_req + additional_upchirps) * m_number_of_bins, 0);
                    for (int n = 0; n < (m_n_up_req + additional_upchirps) * m_number_of_bins; n++)
                    {
                        sfo_corr_vect[n] = gr_expj(-2 * M_PI * (pow(mod(n, N), 2) / 2 / N * (m_bw / fs_p * m_bw / fs_p - m_bw / fs * m_bw / fs) + (std::floor((float)n / N) * (m_bw / fs_p * m_bw / fs_p - m_bw / fs_p) + m_bw / 2 * (1 / fs - 1 / fs_p)) * mod(n, N)));
                    }

                    volk_32fc_x2_multiply_32fc(&preamble_upchirps[0], &preamble_upchirps[0], &sfo_corr_vect[0], up_symb_to_use * m_number_of_bins);

                    float tmp_sto_frac = estimate_STO_frac(); // better estimation of sto_frac in the beginning of the upchirps


                    float diff_sto_frac = m_sto_frac - tmp_sto_frac;

                    if (abs(diff_sto_frac) <= float(m_os_factor - 1) / m_os_factor) // avoid introducing off-by-one errors by estimating fine_sto=-0.499 , rough_sto=0.499
                        m_sto_frac = tmp_sto_frac;

                    // get SNR estimate from preamble
                    // downsample preab_raw
                    std::vector<gr_complex> corr_preamb;
                    corr_preamb.resize((m_n_up_req + additional_upchirps) * m_number_of_bins, 0);
                    
                    // apply sto correction
                    for (int i = 0; i < (m_n_up_req + additional_upchirps) * m_number_of_bins; i++)
                    {
                        corr_preamb[i] = preamble_raw_up[m_os_factor * (m_number_of_bins - k_hat + i) - int(my_roundf(m_os_factor * m_sto_frac))];
                    }
                    
                    std::rotate(corr_preamb.begin(), corr_preamb.begin() + mod(m_cfo_int, m_number_of_bins), corr_preamb.end());
                    // apply cfo correction
                    volk_32fc_x2_multiply_32fc(&corr_preamb[0], &corr_preamb[0], &CFO_int_correc[0], (m_n_up_req + additional_upchirps) * m_number_of_bins);
                    for (int i = 0; i < (m_n_up_req + additional_upchirps); i++)
                    {
                        volk_32fc_x2_multiply_32fc(&corr_preamb[m_number_of_bins * i], &corr_preamb[m_number_of_bins * i], &CFO_frac_correc[0], m_number_of_bins);
                    }
                    
                    // //apply sfo correction
                    volk_32fc_x2_multiply_32fc(&corr_preamb[0], &corr_preamb[0], &sfo_corr_vect[0], (m_n_up_req + additional_upchirps) * m_number_of_bins);

                    float snr_est = 0;
                    for (int i = 0; i < up_symb_to_use; i++)
                    {
                        snr_est += determine_snr(&corr_preamb[i * m_number_of_bins]);
                    }
                    snr_est /= up_symb_to_use;

                    // update sto_frac to its value at the beginning of the net id
                    m_sto_frac += sfo_hat * m_preamb_len;
                    // ensure that m_sto_frac is in [-0.5,0.5]
                    if (abs(m_sto_frac) > 0.5)
                    {
                        m_sto_frac = m_sto_frac + (m_sto_frac > 0 ? -1 : 1);
                    }
                    // decim net id according to new sto_frac and sto int
                    std::vector<gr_complex> net_ids_samp_dec;
                    net_ids_samp_dec.resize(2 * m_number_of_bins, 0);
                    // start_off gives the offset in the net_id_samp vector required to be aligned in time (CFOint is equivalent to STOint since upchirp_val was forced to 0)
                    int start_off = (int)m_os_factor / 2 - (my_roundf(m_sto_frac * m_os_factor)) + m_os_factor * (.25 * m_number_of_bins + m_cfo_int);
                    for (int i = 0; i < m_number_of_bins * 2; i++)
                    {
                        net_ids_samp_dec[i] = net_id_samp[start_off + i * m_os_factor];
                    }
                    volk_32fc_x2_multiply_32fc(&net_ids_samp_dec[0], &net_ids_samp_dec[0], &CFO_int_correc[0], 2 * m_number_of_bins);

                    
                    // correct CFO_frac in the network ids
                    volk_32fc_x2_multiply_32fc(&net_ids_samp_dec[0], &net_ids_samp_dec[0], &CFO_frac_correc[0], m_number_of_bins);
                    volk_32fc_x2_multiply_32fc(&net_ids_samp_dec[m_number_of_bins], &net_ids_samp_dec[m_number_of_bins], &CFO_frac_correc[0], m_number_of_bins);

                    int netid1 = get_symbol_val(&net_ids_samp_dec[0], &m_downchirp[0]);
                    int netid2 = get_symbol_val(&net_ids_samp_dec[m_number_of_bins], &m_downchirp[0]);
                    one_symbol_off = 0;

                    if (abs(netid1 - (int32_t)m_sync_words[0]) > 2) // wrong id 1, (we allow an offset of 2)
                    {
                        // check if we are in fact checking the second net ID and that the first one was considered as a preamble upchirp
                        if (abs(netid1 - (int32_t)m_sync_words[1]) <= 2)
                        {
                            net_id_off = netid1 - (int32_t)m_sync_words[1];
                            for (int i = m_preamb_len - 2; i < (m_n_up_req + additional_upchirps); i++)
                            {
                                if (get_symbol_val(&corr_preamb[i * m_number_of_bins], &m_downchirp[0]) + net_id_off == m_sync_words[0]) // found the first netID
                                {
                                    one_symbol_off = 1;
                                    if (net_id_off != 0 && abs(net_id_off) > 1)
                                        std::cout << RED << "[frame_sync_impl.cc] net id offset >1: " << net_id_off << RESET << std::endl;
             
                                    items_to_consume = -m_os_factor * net_id_off;
                                    // the first symbol was mistaken for the end of the downchirp. we should correct and output it.

                                    int start_off = (int)m_os_factor / 2 - my_roundf(m_sto_frac * m_os_factor) + m_os_factor * (0.25 * m_number_of_bins + m_cfo_int);
                                    
                                    for (int i = start_off; i < 1.25 * m_samples_per_symbol; i += m_os_factor)
                                    {
                                        samp_buf[int((i - start_off) / m_os_factor) + samp_buf_offset] = additional_symbol_samp[i];
                                    }
                                    samp_buf_offset += m_number_of_bins;
                                    m_state = SFO_COMPENSATION;
                                    symbol_cnt = 1;
                                    frame_cnt++;
                                }
                            }
                            if (!one_symbol_off)
                            {
                                m_state = DETECT;
                                symbol_cnt = 1;
                                items_to_output = 0;
                                k_hat = 0;
                                m_sto_frac = 0;
                                items_to_consume = 0;
                                samp_buf.clear();
                                samp_buf_offset = 0;
                            }
                        }
                        else
                        {
                            m_state = DETECT;
                            symbol_cnt = 1;
                            items_to_output = 0;
                            k_hat = 0;
                            m_sto_frac = 0;
                            items_to_consume = 0;
                        }
                    }
                    else // net ID 1 valid
                    {
                        net_id_off = netid1 - (int32_t)m_sync_words[0];
                        if (mod(netid2 - net_id_off, m_number_of_bins) != (int32_t)m_sync_words[1]) // wrong id 2
                        {
                            m_state = DETECT;
                            symbol_cnt = 1;
                            items_to_output = 0;
                            k_hat = 0;
                            m_sto_frac = 0;
                            items_to_consume = 0;
                        }
                        else
                        {
                            if (net_id_off != 0 && abs(net_id_off) > 1)
                                std::cout << RED << "[frame_sync_impl.cc] net id offset >1: " << net_id_off << RESET << std::endl;
                            items_to_consume = -m_os_factor * net_id_off;
                            m_state = SFO_COMPENSATION;
                            frame_cnt++;

                        }
                    }
                    if (m_state != DETECT)
                    {
                        // update sto_frac to its value at the payload beginning
                        m_sto_frac += sfo_hat * 4.25;
                        sfo_cum = ((m_sto_frac * m_os_factor) - my_roundf(m_sto_frac * m_os_factor)) / m_os_factor;

                        items_to_consume += m_samples_per_symbol / 4 + m_os_factor * m_cfo_int;
                        symbol_cnt = one_symbol_off;
                        float snr_est2 = 0;
                    }
                }
                }

                break;
            }
            case SFO_COMPENSATION:
            {
                if (symbol_cnt < 8)        /// get the header symbols (the first 8 symbols)
                {
                    memcpy(&samp_buf[0]+samp_buf_offset, &in_down[0], m_number_of_bins * sizeof(gr_complex));
                    samp_buf_offset += m_number_of_bins;
                    items_to_consume = m_samples_per_symbol;

                    //   update sfo evolution
                    if (abs(sfo_cum) > 1.0 / 2 / m_os_factor)
                    {
                        items_to_consume -= (-2 * signbit(sfo_cum) + 1);
                        sfo_cum -= (-2 * signbit(sfo_cum) + 1) * 1.0 / m_os_factor;
                    }

                    sfo_cum += sfo_hat;

                    symbol_cnt++;
                }
                
                else    /// The first 8 symbols have been received
                {   
                    //// Build uppchirp accounting for cfo_int
                    build_upchirp(&cor_upchirp[0], mod(m_cfo_int, m_samples_per_symbol), m_sf);
                    volk_32fc_conjugate_32fc(&cor_downchirp[0], &cor_upchirp[0], m_samples_per_symbol);
                    // adapt the downchirp to the cfo_frac of the frame
                    for (uint32_t n = 0; n < m_samples_per_symbol; n++)
                    {
                        cor_downchirp[n] = cor_downchirp[n] * gr_expj(-2 * M_PI * m_cfo_frac / m_samples_per_symbol * n);
                    }

                    uint16_t samp_symb[8];              //// Value of the demodulated 8 first samples
                    //// perform the FFT demod
                    for (uint8_t k(0); k<8; k++){
                        //// We assume an explicit header (case of interest), and so we use of reduced rate
                        samp_symb[k] = mod(get_symbol_val(&samp_buf[k*m_number_of_bins], &cor_downchirp[0]) - 1, (1 << m_sf))/4;
                    }

                    //// Gray demapping
                    for (uint8_t k(0); k<8; k++){
                        samp_symb[k] = (samp_symb[k] ^ (samp_symb[k] >> 1u));
                    }

                    //// Deinterleaving
                    uint8_t cw_len = 8;
                    uint8_t sf_app = m_sf - 2;
                    
                    //// Empty matrices
                    std::vector<std::vector<bool>> inter_bin(cw_len);
                    std::vector<bool> init_bit(cw_len, 0);
                    std::vector<std::vector<bool>> deinter_bin(sf_app, init_bit);
                    
                    // convert decimal vector to binary vector of vector
                    for (int i = 0; i < cw_len; i++) {
                        inter_bin[i] = int2bool(samp_symb[i], sf_app);
                    }

                    // Perform the deinterleaving
                    for (int32_t i = 0; i < cw_len; i++) {
                        for (int32_t j = 0; j < int(sf_app); j++) {
                            deinter_bin[mod((i - j - 1), sf_app)][i] = inter_bin[i][j];
                        }
                    }

                    //// Hamming decoding
                    //// each of the deinter_bin's sf_app rows is a 8 bit codeword in binary form
                    uint8_t header[sf_app];               //// contains the decoded nibbles (potential header)
                    std::vector<bool> data_nibble(4, 0);
                    bool s0, s1, s2 = 0;
                    int syndrom = 0;
                    std::vector<bool> codeword;

                    for (int i(0); i<sf_app; i++){
                        codeword = deinter_bin[i];
                        data_nibble = {codeword[3], codeword[2], codeword[1], codeword[0]};  // reorganized msb-first

                        if ((count(codeword.begin(), codeword.end(), true) % 2)){  // Don't correct if even number of errors
                            // get syndrom
                            s0 = codeword[0] ^ codeword[1] ^ codeword[2] ^ codeword[4];
                            s1 = codeword[1] ^ codeword[2] ^ codeword[3] ^ codeword[5];
                            s2 = codeword[0] ^ codeword[1] ^ codeword[3] ^ codeword[6];

                            syndrom = s0 + (s1 << 1) + (s2 << 2);

                            switch (syndrom) {
                                case 5:
                                    data_nibble[3].flip();
                                    break;
                                case 7:
                                    data_nibble[2].flip();
                                    break;
                                case 3:
                                    data_nibble[1].flip();
                                    break;
                                case 6:
                                    data_nibble[0].flip();
                                    break;
                                default:  // either parity bit wrong or no error
                                    break;
                            }
                        }

                        header[i] = bool2int(data_nibble);
                    }

                    //// Time to decode the header
                    /** [0] is payload len MSB 
                     *  [1] is payload len LSB
                     *  [2] is CR ||| + CRC_flag |
                     *  [3] and [4] are CRC
                    */

                    bool c4 = (header[0] & 0b1000) >> 3 ^ (header[0] & 0b0100) >> 2 ^ (header[0] & 0b0010) >> 1 ^ (header[0] & 0b0001);
                    bool c3 = (header[0] & 0b1000) >> 3 ^ (header[1] & 0b1000) >> 3 ^ (header[1] & 0b0100) >> 2 ^ (header[1] & 0b0010) >> 1 ^ (header[2] & 0b0001);
                    bool c2 = (header[0] & 0b0100) >> 2 ^ (header[1] & 0b1000) >> 3 ^ (header[1] & 0b0001) ^ (header[2] & 0b1000) >> 3 ^ (header[2] & 0b0010) >> 1;
                    bool c1 = (header[0] & 0b0010) >> 1 ^ (header[1] & 0b0100) >> 2 ^ (header[1] & 0b0001) ^ (header[2] & 0b0100) >> 2 ^ (header[2] & 0b0010) >> 1 ^ (header[2] & 0b0001);
                    bool c0 = (header[0] & 0b0001) ^ (header[1] & 0b0010) >> 1 ^ (header[2] & 0b1000) >> 3 ^ (header[2] & 0b0100) >> 2 ^ (header[2] & 0b0010) >> 1 ^ (header[2] & 0b0001);

                    if (header[3] == c4 && header[4] == (c3 << 3 | c2 << 2 | c1 << 1 | c0)){  //// check header checksum
                        m_cr = header[2] >> 1u;
                        out_dict = pmt::dict_add(out_dict, pmt::intern("CR"), pmt::mp(m_cr));    //// update the dictionary
                    }
                    else{
                        std::cout << "frameDetect >>> Incorrect header! " << std::endl;
                    }

                    message_port_pub(pmt::intern("detect_out"), out_dict);  //// send a notification with the SF and the detected CR


                    m_state = DETECT;
                    symbol_cnt = 1;                                 
                    items_to_consume = m_samples_per_symbol;        
                    items_to_output = 0;
                    k_hat = 0;
                    m_sto_frac = 0;
                    samp_buf.clear();
                    samp_buf_offset = 0;
                }
                break;
            }
            default:
            {
                std::cerr << "[LoRa sync] WARNING : No state! Shouldn't happen\n";
                break;
            }
            }
            consume_each(items_to_consume);
            return 0;
        }

    } /* namespace loraJam */
} /* namespace gr */
