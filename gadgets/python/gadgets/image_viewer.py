from gadgetron import Gadget
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure   

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk

from matplotlib.backends.backend_gtk3 import (
    NavigationToolbar2GTK3 as NavigationToolbar)
from matplotlib.backends.backend_gtk3agg import (
    FigureCanvasGTK3Agg as FigureCanvas)

class ImageViewWindow:
    def __init__(self, img_data):
        self.window = Gtk.Window(Gtk.WindowType.TOPLEVEL)
        self.window.connect("delete_event", self.delete_event)
        self.window.connect('key_press_event', self.on_key_press_event)
        self.window.set_default_size(400,300)
        self.window.set_title("Gadgetron Image Viewer")

        self.vbox = Gtk.VBox()
        self.window.add(self.vbox)

        self.fig = Figure(figsize=(5,4), dpi=100)
        
        plt.gray()

        self.ax = self.fig.add_subplot(111)
        self.img_ax = self.ax.imshow(np.squeeze(np.abs(img_data)))

        self.canvas = FigureCanvas(self.fig)  # a gtk.DrawingArea
        self.vbox.pack_start(self.canvas, True, True, 0)
        self.toolbar = NavigationToolbar(self.canvas, self.window)
        self.vbox.pack_start(self.toolbar, False, False, 0)
        self.window.show_all()

    def delete_event(self, widget, event, data=None):
        Gtk.main_quit()
        return False
   
    def on_key_press_event(self, widget, event, data=None):
        keyname = Gdk.keyval_name(event.keyval)
        if (keyname == "Escape"):
            self.window.destroy()
            Gtk.main_quit()
            return False

    def main(self):
        Gtk.main()


class ImageViewer(Gadget):
    def process_config(self, cfg):
        print("Attempting to open window")
        print("Window running")
        #Configuration Ignored

    def process(self, h,im):
        myWindow = ImageViewWindow(im)
        myWindow.main()

        self.put_next(h,im)
        return 0
