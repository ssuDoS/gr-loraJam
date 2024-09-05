/* -*- c++ -*- */
/*
 * Copyright 2024 @dossam.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_LORAJAM_PAYLOADPROCESSING_IMPL_H
#define INCLUDED_LORAJAM_PAYLOADPROCESSING_IMPL_H

#define MAX_FRM_ID  1000            // max value expected for the frame id
#define RH_HDR_LEN  4               // length of the radiohead library header

#include <gnuradio/loraJam/payloadProcessing.h>

namespace gr {
namespace loraJam {

class payloadProcessing_impl : public payloadProcessing
{
private:
    uint8_t m_sf;
    uint8_t m_cr;
    uint8_t m_jam_sb_cnt;           // number of jamming symbols
    uint8_t m_pay_len;              // configured/expected payload length 
    uint8_t payload_len;            // payload length received in tags
    std::string m_output_path;
    std::vector<uint8_t> payload;
    bool crc_valid;
    char log_filename[255];     // let's hope it does not overflow XD
    FILE *jam_log;          // used to log frm_id in received frames
    uint16_t frm_id;

public:
    payloadProcessing_impl(uint8_t sf, uint8_t cr, uint8_t jam_sb_cnt, uint8_t pay_len, std::string output_path);
    ~payloadProcessing_impl();

    // Where all the action really happens
    void forecast(int noutput_items, gr_vector_int& ninput_items_required);

    int general_work(int noutput_items,
                     gr_vector_int& ninput_items,
                     gr_vector_const_void_star& input_items,
                     gr_vector_void_star& output_items);
};

} // namespace loraJam
} // namespace gr

#endif /* INCLUDED_LORAJAM_PAYLOADPROCESSING_IMPL_H */
