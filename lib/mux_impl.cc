/* -*- c++ -*- */
/*
 * Copyright 2024 @dossam.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "mux_impl.h"
#include <gnuradio/io_signature.h>

namespace gr {
namespace loraJam {

mux::sptr mux::make(uint32_t samp_rate, uint32_t bw)
{
    return gnuradio::make_block_sptr<mux_impl>(samp_rate, bw);
}


/*
 * The private constructor
 */
mux_impl::mux_impl(uint32_t samp_rate, uint32_t bw)
    : gr::block("mux",
                gr::io_signature::make(1, NB_SF, sizeof(gr_complex)),
                gr::io_signature::make(1, 1, sizeof(gr_complex))),
                m_os_factor(samp_rate/bw)
{
    for (uint8_t i(0); i<NB_SF; i++){
        m_samples_per_symbol[i] = (uint32_t) m_os_factor * (1u << (i+7));
    }

    current_in = 0;     // unsued for the moment
}

/*
 * Our virtual destructor.
 */
mux_impl::~mux_impl() {}

void mux_impl::forecast(int noutput_items, gr_vector_int& ninput_items_required)
{
    for (uint8_t i(0); i<NB_SF; i++){
        ninput_items_required[i] = (noutput_items/m_samples_per_symbol[i]) * m_samples_per_symbol[i];
    }
}

int mux_impl::general_work(int noutput_items,
                           gr_vector_int& ninput_items,
                           gr_vector_const_void_star& input_items,
                           gr_vector_void_star& output_items)
{
    const gr_complex *in[NB_SF];
    for (uint8_t i(0); i<NB_SF; i++){
        in[i] = static_cast<const gr_complex*>(input_items[i]);
    }
    
    auto out = static_cast<gr_complex*>(output_items[0]);
    m_noutput_items = noutput_items;
    uint32_t nitems_to_copy = 0;

    for (uint8_t i(0); i<NB_SF; i++){
        if (m_noutput_items < m_samples_per_symbol[0])  // if the remaining space is not enough for samples of SF=7, then it is not for higher SFs either.
            break;

        nitems_to_copy = (ninput_items[i]/m_samples_per_symbol[i])*m_samples_per_symbol[i];         // process by batch of m_samples_per_symbol[i]
        if(m_noutput_items < nitems_to_copy || nitems_to_copy < m_samples_per_symbol[i])   // ensure there the remaining output buffer space is enough for this input's samples, and there is at least a whole bacth
            continue;

        memcpy(&out[noutput_items-m_noutput_items], in[i], nitems_to_copy * sizeof(gr_complex));
        m_noutput_items -= nitems_to_copy;
        consume(i, nitems_to_copy);
    }

    return noutput_items - m_noutput_items;
}

} /* namespace loraJam */
} /* namespace gr */
