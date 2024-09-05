#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: Not titled yet
# Author: ziree
# GNU Radio version: 3.10.3.0

from gnuradio import gr
from gnuradio.filter import firdes
from gnuradio.fft import window
import sys
import signal
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation
from gnuradio import loraJam
from gnuradio import uhd
import time




class lora_jam(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Not titled yet", catch_exceptions=True)

        ##################################################
        # Variables
        ##################################################
        self.sf = sf = 10
        self.samp_rate = samp_rate = 250000
        self.center_freq_2 = center_freq_2 = 433.1e6
        self.center_freq = center_freq = 868.1e6
        self.bw = bw = 125000
        self.Att_dB_0 = Att_dB_0 = 0
        self.Att_dB = Att_dB = 0

        ##################################################
        # Blocks
        ##################################################
        self.uhd_usrp_source_0_0 = uhd.usrp_source(
            ",".join(('', "addr=192.168.10.2")),
            uhd.stream_args(
                cpu_format="fc32",
                args='',
                channels=list(range(0,1)),
            ),
        )
        self.uhd_usrp_source_0_0.set_samp_rate(samp_rate)
        # No synchronization enforced.

        self.uhd_usrp_source_0_0.set_center_freq(center_freq, 0)
        self.uhd_usrp_source_0_0.set_antenna("RX2", 0)
        self.uhd_usrp_source_0_0.set_bandwidth(bw, 0)
        self.uhd_usrp_source_0_0.set_gain(0, 0)
        self.uhd_usrp_source_0_0.set_min_output_buffer((2**sf<<2))
        self.uhd_usrp_sink_0_0 = uhd.usrp_sink(
            ",".join(('', "addr=192.168.10.2")),
            uhd.stream_args(
                cpu_format="fc32",
                args='',
                channels=list(range(0,1)),
            ),
            "frame_len",
        )
        self.uhd_usrp_sink_0_0.set_samp_rate(samp_rate)
        # No synchronization enforced.

        self.uhd_usrp_sink_0_0.set_center_freq(center_freq, 0)
        self.uhd_usrp_sink_0_0.set_antenna('TX/RX', 0)
        self.uhd_usrp_sink_0_0.set_bandwidth(bw, 0)
        self.uhd_usrp_sink_0_0.set_gain(100, 0)
        self.loraJam_triggerPmt_0 = loraJam.triggerPmt(False, True, 2, 0, 0, 2)
        self.loraJam_jammerPmt_0 = loraJam.jammerPmt(sf, 250000, 125000, False)
        self.loraJam_frameDetect_0 = loraJam.frameDetect(int(center_freq), bw, sf, [18], (int(samp_rate/bw)), 8)


        ##################################################
        # Connections
        ##################################################
        self.msg_connect((self.loraJam_frameDetect_0, 'detect_out'), (self.loraJam_triggerPmt_0, 'msg_in'))
        self.msg_connect((self.loraJam_triggerPmt_0, 'SF9'), (self.loraJam_jammerPmt_0, 'symb_in'))
        self.connect((self.loraJam_jammerPmt_0, 0), (self.uhd_usrp_sink_0_0, 0))
        self.connect((self.uhd_usrp_source_0_0, 0), (self.loraJam_frameDetect_0, 0))


    def get_sf(self):
        return self.sf

    def set_sf(self, sf):
        self.sf = sf

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.uhd_usrp_sink_0_0.set_samp_rate(self.samp_rate)
        self.uhd_usrp_source_0_0.set_samp_rate(self.samp_rate)

    def get_center_freq_2(self):
        return self.center_freq_2

    def set_center_freq_2(self, center_freq_2):
        self.center_freq_2 = center_freq_2

    def get_center_freq(self):
        return self.center_freq

    def set_center_freq(self, center_freq):
        self.center_freq = center_freq
        self.uhd_usrp_sink_0_0.set_center_freq(self.center_freq, 0)
        self.uhd_usrp_source_0_0.set_center_freq(self.center_freq, 0)
        self.uhd_usrp_source_0_0.set_center_freq(self.center_freq, 1)

    def get_bw(self):
        return self.bw

    def set_bw(self, bw):
        self.bw = bw
        self.uhd_usrp_sink_0_0.set_bandwidth(self.bw, 0)
        self.uhd_usrp_source_0_0.set_bandwidth(self.bw, 0)
        self.uhd_usrp_source_0_0.set_bandwidth(self.bw, 1)

    def get_Att_dB_0(self):
        return self.Att_dB_0

    def set_Att_dB_0(self, Att_dB_0):
        self.Att_dB_0 = Att_dB_0

    def get_Att_dB(self):
        return self.Att_dB

    def set_Att_dB(self, Att_dB):
        self.Att_dB = Att_dB




def main(top_block_cls=lora_jam, options=None):
    if gr.enable_realtime_scheduling() != gr.RT_OK:
        print("Error: failed to enable real-time scheduling.")
    tb = top_block_cls()

    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()

        sys.exit(0)

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    tb.start()

    try:
        input('Press Enter to quit: ')
    except EOFError:
        pass
    tb.stop()
    tb.wait()


if __name__ == '__main__':
    main()
