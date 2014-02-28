#!/usr/bin/env python
import zmq
import sys
import time
import argparse

from psdata import HistData, ImageData, XYPlotData
from psplot import Hist, Image, XYPlot, PlotInfo
import matplotlib.pyplot as plt


TypeMap = {
    HistData: Hist,
    ImageData: Image,
    XYPlotData: XYPlot,
}


def parse_cmdline():
    parser = argparse.ArgumentParser(
        description='Psana plot client application'
    )

    parser.add_argument(
        'topic',
        help='The topic from the server to suscribe too'
    )

    parser.add_argument(
        '-s',
        '--server',
        metavar='SERVER',
        default='localhost',
        help='the host name of the server'
    )

    parser.add_argument(
        '-p',
        '--port',
        metavar='PORT',
        type=int,
        default=5556,
        help='the tcp port to use on the server'
    )

    parser.add_argument(
        '-r',
        '--rate',
        metavar='RATE',
        type=float,
        default=5.0,
        help='update rate of the histogram in Hz'
    )

    parser.add_argument(
        '-b',
        '--buffer',
        metavar='BUFFER',
        type=int,
        default=10,
        help='the size in messages of recieve buffer'
    )

    parser.add_argument(
        '-x',
        '--x-range',
        metavar='X_RANGE',
        type=int,
        nargs=2,
        default=None,
        help='the fixed x range for any plots'
    )

    parser.add_argument(
        '-y',
        '--y-range',
        metavar='Y_RANGE',
        type=int,
        nargs=2,
        default=None,
        help='the fixed y range for any plots'
    )

    parser.add_argument(
        '-z',
        '--z-range',
        metavar='Z_RANGE',
        type=int,
        nargs=2,
        default=None,
        help='the fixed z range for any plots'
    )

    parser.add_argument(
        '-a',
        '--aspect',
        metavar='aspect',
        type=float,
        default=None,
        help='the aspect ratio for the plot'
    )

    return parser.parse_args()


def socket_recv(socket):
    topic = socket.recv()
    return socket.recv_pyobj()


def get_socket_gen(socket):
    poller = zmq.Poller()
    poller.register(socket, zmq.POLLIN)
    while True:
        socks = dict(poller.poll(100))
        if socks.get(socket) == zmq.POLLIN:
            yield socket_recv(socket)
        else:
            yield


def main():
    try:
        args = parse_cmdline()

        context = zmq.Context()
        socket = context.socket(zmq.SUB)
        socket.setsockopt(zmq.SUBSCRIBE, args.topic)
        socket.setsockopt(zmq.RCVHWM, args.buffer)
        socket.connect("tcp://%s:%d" % (args.server, args.port))

        # create the plot info object from cmd args
        info = PlotInfo(args.x_range, args.y_range, args.z_range, args.aspect)

        # grab an initial datagram from the server
        init_data = socket_recv(socket)
        data_type = TypeMap.get(type(init_data))

        if data_type is not None:
            plot = data_type(init_data, get_socket_gen(socket), info, rate=1.0/args.rate)
            plot_ani = plot.animate()
            try:
                plt.show()
            except:
                # sort of ugly but this can throw all kinds of errors when closing a window
                pass
        else:
            print 'Server did not return a valid datatype!'
        
    except KeyboardInterrupt:
        print '\nExitting client!'


if __name__ == '__main__':
    main()
