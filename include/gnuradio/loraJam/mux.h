/* -*- c++ -*- */
/*
 * Copyright 2024 @dossam.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_LORAJAM_MUX_H
#define INCLUDED_LORAJAM_MUX_H

#include <gnuradio/block.h>
#include <gnuradio/loraJam/api.h>

namespace gr {
namespace loraJam {

/*!
 * \brief <+description of block+>
 * \ingroup loraJam
 *
 */
class LORAJAM_API mux : virtual public gr::block
{
public:
    typedef std::shared_ptr<mux> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of loraJam::mux.
     *
     * To avoid accidental use of raw pointers, loraJam::mux's
     * constructor is in a private implementation
     * class. loraJam::mux::make is the public interface for
     * creating new instances.
     */
    static sptr make(uint32_t samp_rate, uint32_t bw);
};

} // namespace loraJam
} // namespace gr

#endif /* INCLUDED_LORAJAM_MUX_H */
