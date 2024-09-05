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
#include <gnuradio/loraJam/utilities.h>

namespace gr {
namespace loraJam {

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
