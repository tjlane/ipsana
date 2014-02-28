
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation

class PlotInfo(object):
    def __init__(self, xrange=None, yrange=None, zrange=None, aspect=None):
        self.xmin = None
        self.xmax = None
        self.ymin = None
        self.ymax = None
        self.zmin = None
        self.zmax = None
        self.aspect = aspect
        if xrange is not None:
            self.xmin = xrange[0]
            self.xmax = xrange[1]
        if yrange is not None:
            self.ymin = yrange[0]
            self.ymax = yrange[1]
        if zrange is not None:
            self.zmin = zrange[0]
            self.zmax = zrange[1]


class Plot(object):
    def __init__(self, init, framegen, info, rate):
        self.figure, self.ax = plt.subplots()
        self.figure.canvas.set_window_title(init.title)
        self.ax.set_title(init.ts, loc='right')
        self.framegen = framegen
        self.rate_ms = rate * 1000
        self.info = info

    def update(self, data):
        pass

    def animate(self):
        return animation.FuncAnimation(self.figure, self.update, self.ani_func, interval=self.rate_ms)

    def ani_func(self):
        yield self.framegen.next()


class Image(Plot):
    def __init__(self, init_im, framegen, info, rate=1):
        super(Image, self).__init__(init_im, framegen, info, rate)
        self.im = self.ax.imshow(init_im.image)
        self.im.set_clim(self.info.zmin, self.info.zmax)
        self.figure.colorbar(self.im)
        if self.info.aspect is not None:
            self.ax.set_aspect(self.info.aspect)

    def update(self, data):
        """
        Updates the data in the image - none means their was no update for this interval
        """
        if data is not None:
            self.ax.set_title(data.ts, loc='right')
            self.im.set_data(data.image)
        return self.im


class Hist(Plot):
    def __init__(self, init_hist, datagen, info, rate=1):
        super(Hist, self).__init__(init_hist, datagen, info, rate)
        self.hist = self.ax.hist(init_hist.data, init_hist.bins)

    def update(self, data):
        if data is not None:
            self.ax.set_title(data.ts, loc='right')
            self.hist.set_data(data.data)
        return self.hist


class XYPlot(Plot):
    def __init__(self, init_plot, datagen, info, rate=1):
        super(XYPlot, self).__init__(init_plot, datagen, info, rate)
        self.plot = self.ax.plot(init_plot.xdata, init_plot.ydata, lw=2)

    def update(self, data):
        if data is not None:
            self.ax.set_title(data.ts, loc='right')
            for i in range(data.ydata.shape[1]):
                self.plot[i].set_data(data.xdata, data.ydata[:,i])
        return self.plot


class IvsQPlot(XYPlot):
    def __init__(self, init_plot, datagen, info, rate=1):
        super(IvsQPlot, self).__init__(init_plot, datagen, info, rate)
        plt.xlabel(r'Momentum Transfer $(\AA^{-1})$')
        plt.ylabel('Intensity')

