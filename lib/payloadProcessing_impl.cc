/* -*- c++ -*- */
/*
 * Copyright 2024 @dossam.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "payloadProcessing_impl.h"
#include <gnuradio/io_signature.h>

namespace gr {
namespace loraJam {

using input_type = uint8_t;
payloadProcessing::sptr
payloadProcessing::make(uint8_t sf, uint8_t cr, uint8_t jam_sb_cnt, uint8_t pay_len, std::string out_path)
{
    return gnuradio::make_block_sptr<payloadProcessing_impl>(sf, cr, jam_sb_cnt, pay_len, out_path);
}


/*
 * The private constructor
 */
payloadProcessing_impl::payloadProcessing_impl(uint8_t sf,
                                         uint8_t cr,
                                         uint8_t jam_sb_cnt,
                                         uint8_t pay_len,
                                         std::string output_path)
    : gr::block("payloadProcessing",
                gr::io_signature::make(
                    1 /* min inputs */, 1 /* max inputs */, sizeof(input_type)),
                gr::io_signature::make(
                    0 /* min outputs */, 0 /*max outputs */, 0)),
                    m_sf(sf),
                    m_cr(cr),
                    m_jam_sb_cnt(jam_sb_cnt),
                    m_pay_len(pay_len),
                    m_output_path(output_path)

{
    set_tag_propagation_policy(TPP_DONT);

    sprintf(log_filename, "%s/jam_log_sf-%u_cr-%u_jam-%u_plen-%u.txt", m_output_path.c_str(), m_sf, m_cr, m_jam_sb_cnt, m_pay_len);

    std::cout << "payloadProcessing << log_filename" << log_filename << std::endl;
    if((jam_log = fopen(log_filename, "w")) == NULL){
        perror("Error while opening log_file file");
        exit(EXIT_FAILURE);
    }

    frm_id = 0;
}

/*
 * Our virtual destructor.
 */
payloadProcessing_impl::~payloadProcessing_impl() {
    fclose(jam_log);
}

void payloadProcessing_impl::forecast(int noutput_items,
                                   gr_vector_int& ninput_items_required)
{
    ninput_items_required[0] = 1;
}

int payloadProcessing_impl::general_work(int noutput_items,
                                      gr_vector_int& ninput_items,
                                      gr_vector_const_void_star& input_items,
                                      gr_vector_void_star& output_items)
{
    auto in = static_cast<const input_type*>(input_items[0]);
    int items_to_use = ninput_items[0];
    
    std::vector<tag_t> tags;
    get_tags_in_window(tags, 0, 0, items_to_use, pmt::string_to_symbol("frame_info"));
    if (tags.size()){
        pmt::pmt_t err = pmt::string_to_symbol("error");
        payload_len = pmt::to_long(pmt::dict_ref(tags[0].value, pmt::intern("pay_len"), err));
        crc_valid = pmt::to_bool(pmt::dict_ref(tags[0].value, pmt::intern("crc_valid"), err));

        if(payload_len != m_pay_len){
            std::cout << "payloadProcessing << unexpected payload length! Expected: " << (int) m_pay_len << ", Received: " << (int) payload_len << std::endl; 
        }
    }
    
    for (int i(0); i<items_to_use; i++){
        payload.push_back(in[i]);
    }

    if(payload.size() >= payload_len){  // got a full frame
        if(crc_valid){                 // crc valid
            char tmp_msg[4];
            for (int i(0); i<4; i++)
                tmp_msg[i] = payload[RH_HDR_LEN + i];
            
            sscanf(tmp_msg, "%04hu", &frm_id);

            if(frm_id >= MAX_FRM_ID){
                fprintf(stdout, "ALERT\n");
            }

            fprintf(stdout, "frm_id: %hu\n", frm_id);
            fprintf(jam_log, "%hu\n", frm_id);
        }

        payload.erase(payload.begin(), payload.begin()+payload_len);
        
    }

    consume_each(items_to_use);

    return 0;
}

} /* namespace loraJam */
} /* namespace gr */
