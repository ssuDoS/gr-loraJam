/* -*- c++ -*- */
/*
 * Copyright 2024 @dossam.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "jammer_impl.h"
#include <gnuradio/io_signature.h>

namespace gr {
namespace loraJam {

jammer::sptr jammer::make(uint8_t sf, uint32_t samp_rate, uint32_t bw)
{
    return gnuradio::make_block_sptr<jammer_impl>(sf, samp_rate, bw);
}


/*
 * The private constructor
 */
jammer_impl::jammer_impl(uint8_t sf, uint32_t samp_rate, uint32_t bw)
    : gr::block("jammer",
                     gr::io_signature::make(1, 1, sizeof(uint32_t)),
                     gr::io_signature::make(1, 1, sizeof(gr_complex))
                     ),
                     m_sf(sf), m_samp_rate(samp_rate), m_bw(bw)
{
    m_number_of_bins = (uint32_t)(1u << m_sf);
    m_os_factor = m_samp_rate / m_bw;
    m_samples_per_symbol = (uint32_t)(m_number_of_bins*m_os_factor);

    m_downchirp.resize(m_samples_per_symbol);
    m_upchirp.resize(m_samples_per_symbol);
    build_ref_chirps(&m_upchirp[0], &m_downchirp[0], m_sf,m_os_factor);
    
    set_output_multiple(m_samples_per_symbol);
}

/*
 * Our virtual destructor.
 */
jammer_impl::~jammer_impl() {}

void jammer_impl::forecast (int noutput_items, gr_vector_int &ninput_items_required){
    ninput_items_required[0] = noutput_items/m_samples_per_symbol;
}

int jammer_impl::general_work(int noutput_items,
                            gr_vector_int &ninput_items,
                            gr_vector_const_void_star &input_items,
                            gr_vector_void_star &output_items)
{
    const uint32_t *in = (const uint32_t *)input_items[0];
    gr_complex *out = (gr_complex *)output_items[0];

    int to_produce = 0;
    
    for (int i(0); i<ninput_items[0]; i++){     // TODO: add a condition on the output buffer, not to overflow it
        build_upchirp(&out[i*m_samples_per_symbol], in[i], m_sf, m_os_factor);
        to_produce ++;
    }

    to_produce *= m_samples_per_symbol;
    tag_t tag;
    tag.offset = nitems_written(0);
    tag.key = pmt::intern("frame_len");
    tag.value = pmt::mp(to_produce);
    add_item_tag(0, tag);
    
    consume_each(ninput_items[0]);
    return to_produce;
}

} /* namespace loraJam */
} /* namespace gr */
