#!/usr/bin/env python
import zmq
import random
import sys
import time
import argparse
import numpy as np

from data import ImageData
from psana import *

def parse_cmdline():
    parser = argparse.ArgumentParser(
        description='Psana plot server application'
    )

    parser.add_argument(
        '-p',
        '--port',
        metavar='PORT',
        type=int,
        default=5556,
        help='the tcp port the server listens on'
    )

    parser.add_argument(
        '-b',
        '--buffer',
        metavar='BUFFER',
        type=int,
        default=10,
        help='the size in messages of send buffer'
    )

    return parser.parse_args()


def main():
    try:
        args = parse_cmdline()

        context = zmq.Context()
        socket = context.socket(zmq.PUB)
        socket.setsockopt(zmq.SNDHWM, args.buffer)
        socket.bind("tcp://*:%d" % args.port)

        counter = 0
        events = DataSource('exp=XCS/xcs84213:run=4').events()
        input_srcs = [
            (Source('DetInfo(XcsBeamline.1:Tm6740.5)'), Camera.FrameV1, Camera.FrameV1.data16, 'yag5'),
            (Source('DetInfo(XcsEndstation.1:Opal1000.1)'), Camera.FrameV1, Camera.FrameV1.data16, 'xcs-spectrometer'),
        ]

        for evt in events:
            evt_data = evt.get(EventId)
            evt_ts, _ = evt_data.time()
            for src, data_type, data_func, topic in input_srcs:
                frame = evt.get(data_type, src)
                mydata = ImageData(evt_ts, topic, data_func(frame))
                socket.send(topic, zmq.SNDMORE)
                socket.send_pyobj(mydata)
            counter += 1
            time.sleep(.1)
    except KeyboardInterrupt:
        print '\nExitting server!'
        pass


if __name__ == '__main__':
    main()
