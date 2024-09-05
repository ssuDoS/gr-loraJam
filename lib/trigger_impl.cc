/* -*- c++ -*- */
/*
 * Copyright 2024 @dossam.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "trigger_impl.h"
#include <gnuradio/io_signature.h>

namespace gr {
namespace loraJam {

trigger::sptr trigger::make(bool incl_cr0, bool verbose, uint8_t cr1_symb_cnt, uint8_t cr2_symb_cnt, uint8_t cr3_symb_cnt, uint8_t cr4_symb_cnt)
{
    return gnuradio::make_block_sptr<trigger_impl>(
        incl_cr0, verbose, cr1_symb_cnt, cr2_symb_cnt, cr3_symb_cnt, cr4_symb_cnt);
}

/*
 * The private constructor
 */
trigger_impl::trigger_impl(bool incl_cr0,
                        bool verbose,
                        uint8_t cr1_symb_cnt,
                        uint8_t cr2_symb_cnt,
                        uint8_t cr3_symb_cnt,
                        uint8_t cr4_symb_cnt)
    : gr::block("trigger",
                gr::io_signature::make(0, 0, 0),
                gr::io_signature::make(1, 6, sizeof(uint32_t)))  // for SF from 7 to 12; 12 is the greatest symbol size
{

    m_incl_cr0 = incl_cr0;
    m_verbose = verbose;
    m_cr_symb_cnt[0] = cr1_symb_cnt;
    m_cr_symb_cnt[1] = cr2_symb_cnt;
    m_cr_symb_cnt[2] = cr3_symb_cnt;
    m_cr_symb_cnt[3] = cr4_symb_cnt;

    rgen = std::mt19937(50885);         // seed the generator with a constant value

    // Set the distribution range according to each sf
    for (int i=0; i<NB_SF; i++){
        sf_distrib[i] = std::uniform_int_distribution<>(0, i+6);    // range is [0, SF-1]
    }


    m_err = pmt::mp(-1);

    message_port_register_in(pmt::intern("msg_in"));
    set_msg_handler(pmt::intern("msg_in"), 
        [this](pmt::pmt_t msg){ this->message_handler(msg);});
}

/*
 * Our virtual destructor.
 */
trigger_impl::~trigger_impl() {}

void trigger_impl::forecast(int noutput_items, gr_vector_int& ninput_items_required)
{
    //
}

int trigger_impl::general_work(int noutput_items,
                               gr_vector_int& ninput_items,
                               gr_vector_const_void_star& input_items,
                               gr_vector_void_star& output_items)
{
    uint32_t *out[NB_SF];
    for (int i(0); i<NB_SF; i++)
        out[i] = static_cast<uint32_t*>(output_items[i]);
    
    for (uint8_t i(0); i<NB_SF; i++){
        if (pending[i].size() && noutput_items >= pending[i].size()){
            memcpy(out[0], &pending[i][0], pending[i].size()*sizeof(uint16_t));
            produce(i, pending[i].size());
            pending[i].clear();
        }
    }

    consume_each(0);
    return WORK_CALLED_PRODUCE;
}

void trigger_impl::message_handler(pmt::pmt_t msg){
    // msg SHOULD to be a dictionary
    int sf = pmt::to_long(pmt::dict_ref(msg, pmt::intern("SF"), m_err));
    int cr = pmt::to_long(pmt::dict_ref(msg, pmt::intern("CR"), m_err));

    if (sf == -1 || cr == -1){
        std::cout << "Error on SF or CR value";
        return;
    }

    if (cr == 0 && m_incl_cr0 == false){
        return;
    }

    if (m_verbose){
        std::cout << "trigger<> SF: " << sf << ", CR: " << cr << std::endl; 
    }

    // Generate a random number and add it to the matching pending data (to be read in general_work)
    uint8_t symb_cnt = 0;   // number of symbol to generate, according to CR
    if(cr == 0){
        symb_cnt = m_cr_symb_cnt[3];  // same symb_cnt for CR 0 and 3
    }
    else{
        symb_cnt = m_cr_symb_cnt[cr - 1];
    }

    nb_gen = sf_distrib[sf-7](rgen);    // generate one random number

    for (uint8_t i=0; i<symb_cnt; i++){ // replicate the number of times configured.
        pending[sf-7].push_back(nb_gen);
    }

    
}

} /* namespace loraJam */
} /* namespace gr */
