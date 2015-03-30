from gadgetron import Gadget
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure   

#from matplotlib.axes import Subplot   
# uncomment to select /GTK/GTKAgg/GTKCairo
from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
#from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
#from matplotlib.backends.backend_gtkcairo import FigureCanvasGTKCairo as FigureCanvas

# or NavigationToolbar for classic
#from matplotlib.backends.backend_gtk import NavigationToolbar2GTK as NavigationToolbar
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar

import pygtk
pygtk.require('2.0')
import gtk

class ImageViewWindow:
    def __init__(self, img_data):
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.window.connect("delete_event", self.delete_event)
        self.window.connect('key_press_event', self.on_key_press_event)
        self.window.set_default_size(400,300)
        self.window.set_title("Gadgetron Image Viewer")

        self.vbox = gtk.VBox()
        self.window.add(self.vbox)

        self.fig = Figure(figsize=(5,4), dpi=100)
        
        plt.gray()

        self.ax = self.fig.add_subplot(111)
        self.img_ax = self.ax.imshow(np.squeeze(np.abs(img_data)))

        self.canvas = FigureCanvas(self.fig)  # a gtk.DrawingArea
        self.vbox.pack_start(self.canvas)
        self.toolbar = NavigationToolbar(self.canvas, self.window)
        self.vbox.pack_start(self.toolbar, False, False)
        self.window.show_all()

    def delete_event(self, widget, event, data=None):
        gtk.main_quit()
        return False
   
    def on_key_press_event(self, widget, event, data=None):
        keyname = gtk.gdk.keyval_name(event.keyval)
        if (keyname == "Escape"):
            self.window.destroy()
            gtk.main_quit()
            return False

    def main(self):
        gtk.main()


class ImageViewer(Gadget):
    def process_config(self, cfg):
        print "Attempting to open window"
        print "Window running"
        #Configuration Ignored

    def process(self, h,im):
        myWindow = ImageViewWindow(im)
        myWindow.main()

        self.put_next(h,im)
        return 0
