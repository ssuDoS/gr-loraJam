/* -*- c++ -*- */
/*
 * Copyright 2024 @dossam.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_LORAJAM_PAYLOADPROCESSING_H
#define INCLUDED_LORAJAM_PAYLOADPROCESSING_H

#include <gnuradio/block.h>
#include <gnuradio/loraJam/api.h>

namespace gr {
namespace loraJam {

/*!
 * \brief <+description of block+>
 * \ingroup loraJam
 *
 */
class LORAJAM_API payloadProcessing : virtual public gr::block
{
public:
    typedef std::shared_ptr<payloadProcessing> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of loraJam::payloadProcessing.
     *
     * To avoid accidental use of raw pointers, loraJam::payloadProcessing's
     * constructor is in a private implementation
     * class. loraJam::payloadProcessing::make is the public interface for
     * creating new instances.
     */
    static sptr make(uint8_t sf,
                     uint8_t cr,
                     uint8_t jam_sb_cnt,
                     uint8_t pay_len,
                     std::string output_path);
};

} // namespace loraJam
} // namespace gr

#endif /* INCLUDED_LORAJAM_PAYLOADPROCESSING_H */
