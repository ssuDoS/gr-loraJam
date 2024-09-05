/* -*- c++ -*- */
/*
 * Copyright 2024 @dossam.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_LORAJAM_TRIGGER_IMPL_H
#define INCLUDED_LORAJAM_TRIGGER_IMPL_H
#define NB_SF   6   // Total number of SF considered: 7->12 = 6

#include <gnuradio/loraJam/trigger.h>
#include <random>

namespace gr {
namespace loraJam {

class trigger_impl : public trigger
{
private:
    bool m_incl_cr0;        // indicates whether frame with CR=0 (undetected) should be considered valid 
    bool m_verbose;         // indicates whether capture details should be printed
    uint8_t m_cr_symb_cnt[4];   // holds the number of symbol to be generated for each CR. idx = CR - 1
    std::vector<uint16_t> pending[NB_SF];   // index 'i' holds pending data for sf 'i+7'
    uint16_t nb_gen;                    // temporarily store the generated number
    std::mt19937 rgen;                  // Mersenne Twister random number generator 
    std::uniform_int_distribution<> sf_distrib[NB_SF];     // each distrib is associated with one SF
    pmt::pmt_t m_err;       // holds the expected value to be returned in case of error (for SF and CR)
    void message_handler(pmt::pmt_t msg);

public:
    trigger_impl(bool incl_cr0,
                bool verbose,
                uint8_t cr1_symb_cnt,
                uint8_t cr2_symb_cnt,
                uint8_t cr3_symb_cnt,
                uint8_t cr4_symb_cnt);
    ~trigger_impl();

    // Where all the action really happens
    void forecast(int noutput_items, gr_vector_int& ninput_items_required);

    int general_work(int noutput_items,
                     gr_vector_int& ninput_items,
                     gr_vector_const_void_star& input_items,
                     gr_vector_void_star& output_items);

};

} // namespace loraJam
} // namespace gr

#endif /* INCLUDED_LORAJAM_TRIGGER_IMPL_H */
