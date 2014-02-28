class Data(object):
    def __init__(self, ts, title):
        self.ts = ts
        self.title = title


class ImageData(Data):
    """
    A data container for image data
    """

    def __init__(self, ts, title, image):
        super(ImageData, self).__init__(ts, title)
        self.image = image


class HistData(Data):
    """
    A data container for 1-d histogram data
    """

    def __init__(self, ts, title, data, bins):
        super(HistData, self).__init__(ts, title)
        self.data = data
        self.bins = bins


class XYPlotData(Data):

    def __init__(self, ts, title, xdata, ydata):
        super(XYPlotData, self).__init__(ts, title)
        self.xdata = xdata
        self.ydata = ydata
