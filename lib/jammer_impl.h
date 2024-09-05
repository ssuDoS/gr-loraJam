/* -*- c++ -*- */
/*
 * Copyright 2024 @dossam.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_LORAJAM_JAMMER_IMPL_H
#define INCLUDED_LORAJAM_JAMMER_IMPL_H

#include <gnuradio/loraJam/jammer.h>
#include <iostream>
#include <fstream>
#include <gnuradio/expj.h>
#include <volk/volk.h>



namespace gr {
namespace loraJam {

/* ********* Helper 'defines' and functions, copied from utilities (removed) ************ */

/**
 *  \brief  Return a modulated upchirp using s_f=bw
 *
 *  \param  chirp
 *          The pointer to the modulated upchirp
 *  \param  id
 *          The number used to modulate the chirp
 * \param   sf
 *          The spreading factor to use
 * \param os_factor
 *          The oversampling factor used to generate the upchirp
 */
inline void build_upchirp(gr_complex* chirp, uint32_t id, uint8_t sf, uint8_t os_factor = 1){
    double N = (1 << sf)  ;
    int n_fold = N* os_factor - id*os_factor;
    for(int n = 0; n < N* os_factor; n++){
        if(n<n_fold)
            chirp[n] = gr_complex(1.0,0.0)*gr_expj(2.0*M_PI *(n*n/(2*N)/pow(os_factor,2)+(id/N-0.5)*n/os_factor));
        else
            chirp[n] = gr_complex(1.0,0.0)*gr_expj(2.0*M_PI *(n*n/(2*N)/pow(os_factor,2)+(id/N-1.5)*n/os_factor));

    }
}

/**
 *  \brief  Return the reference chirps using s_f=bw
 *
 *  \param  upchirp
 *          The pointer to the reference upchirp
 *  \param  downchirp
 *          The pointer to the reference downchirp
 * \param   sf
 *          The spreading factor to use
 */
inline void build_ref_chirps(gr_complex* upchirp, gr_complex* downchirp, uint8_t sf, uint8_t os_factor = 1){
    double N = (1 << sf);
    build_upchirp(upchirp,0,sf,os_factor);
    volk_32fc_conjugate_32fc(&downchirp[0], &upchirp[0], N*os_factor);
}

/* *****************************************    ************************************** */


class jammer_impl : public jammer
{
private:
    uint8_t m_sf; ///< Transmission spreading factor
    uint32_t m_samp_rate; ///< Transmission sampling rate
    uint32_t m_bw; ///< Transmission bandwidth (Works only for samp_rate=bw)
    uint32_t m_number_of_bins; ///< number of bin per loar symbol
    uint32_t m_samples_per_symbol; ///< samples per symbols(Works only for 2^sf)
    
    // int m_ninput_items_required; ///< number of samples required to call this block (forecast)

    int m_os_factor; ///< ovesampling factor based on sampling rate and bandwidth

    std::vector<gr_complex> m_upchirp; ///< reference upchirp
    std::vector<gr_complex> m_downchirp; ///< reference downchirp


public:
    jammer_impl(uint8_t sf, uint32_t samp_rate, uint32_t bw);
    ~jammer_impl();

    // Where all the action really happens
    void forecast(int noutput_items, gr_vector_int& ninput_items_required);

    int general_work(int noutput_items,
                     gr_vector_int& ninput_items,
                     gr_vector_const_void_star& input_items,
                     gr_vector_void_star& output_items);
};

} // namespace loraJam
} // namespace gr

#endif /* INCLUDED_LORAJAM_JAMMER_IMPL_H */
