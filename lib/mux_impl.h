/* -*- c++ -*- */
/*
 * Copyright 2024 @dossam.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_LORAJAM_MUX_IMPL_H
#define INCLUDED_LORAJAM_MUX_IMPL_H

#define NB_SF   6

#include <gnuradio/loraJam/mux.h>

namespace gr {
namespace loraJam {

class mux_impl : public mux
{
private:
    uint32_t m_os_factor;
    uint8_t current_in;        // index of the currently processed input (used for balancing)
    uint32_t m_samples_per_symbol[NB_SF];   // number of samples per symbols (index i corresponds to SF i+7)
    uint32_t m_noutput_items;

public:
    mux_impl(uint32_t samp_rate, uint32_t bw);
    ~mux_impl();

    // Where all the action really happens
    void forecast(int noutput_items, gr_vector_int& ninput_items_required);

    int general_work(int noutput_items,
                     gr_vector_int& ninput_items,
                     gr_vector_const_void_star& input_items,
                     gr_vector_void_star& output_items);
};

} // namespace loraJam
} // namespace gr

#endif /* INCLUDED_LORAJAM_MUX_IMPL_H */
