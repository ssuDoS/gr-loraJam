/* -*- c++ -*- */
/*
 * Copyright 2024 @dossam.
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef INCLUDED_LORAJAM_TRIGGER_H
#define INCLUDED_LORAJAM_TRIGGER_H

#include <gnuradio/block.h>
#include <gnuradio/loraJam/api.h>

namespace gr {
namespace loraJam {

/*!
 * \brief <+description of block+>
 * \ingroup loraJam
 *
 */
class LORAJAM_API trigger : virtual public gr::block
{
public:
    typedef std::shared_ptr<trigger> sptr;

    /*!
     * \brief Return a shared_ptr to a new instance of loraJam::trigger.
     *
     * To avoid accidental use of raw pointers, loraJam::trigger's
     * constructor is in a private implementation
     * class. loraJam::trigger::make is the public interface for
     * creating new instances.
     */
    static sptr make(bool incl_cr0,
                     bool verbose,
                     uint8_t cr1_symb_cnt,
                     uint8_t cr2_symb_cnt,
                     uint8_t cr3_symb_cnt,
                     uint8_t cr4_symb_cnt);
};

} // namespace loraJam
} // namespace gr

#endif /* INCLUDED_LORAJAM_TRIGGER_H */
