/* -*- c++ -*- */
/*
 * Copyright 2024 @dossam.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_LORAJAM_FRAMEDETECT_H
#define INCLUDED_LORAJAM_FRAMEDETECT_H

#include <gnuradio/block.h>
#include <gnuradio/loraJam/api.h>

namespace gr {
namespace loraJam {

/*!
 * \brief <+description of block+>
 * \ingroup loraJam
 *
 */
class LORAJAM_API frameDetect : virtual public gr::block
{
public:
    typedef std::shared_ptr<frameDetect> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of loraJam::frameDetect.
     *
     * To avoid accidental use of raw pointers, loraJam::frameDetect's
     * constructor is in a private implementation
     * class. loraJam::frameDetect::make is the public interface for
     * creating new instances.
     */
    static sptr make(uint32_t center_freq,
                     uint32_t bandwidth,
                     uint8_t sf,
                     std::vector<uint16_t> sync_word,
                     uint8_t os_factor,
                     uint16_t preamb_len);
};

} // namespace loraJam
} // namespace gr

#endif /* INCLUDED_LORAJAM_FRAMEDETECT_H */
