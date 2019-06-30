#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Usage:

    Run as standalone python app

    Run as a gadget:
        1. ismrmrd_generate_cartesian_shepp_logan -r 24
        2. gadgetron_ismrmrd_client -f testdata.h5 -c python_realtime.xml -o out.h5

"""

import sys

# When load in gadget, there is no argv attr which will cause "from gi.repository import something" failed
if not hasattr(sys, 'argv'):
    sys.argv  = ['']

import time
from threading import Thread, Event

import gi
gi.require_version('Gtk', '3.0')
gi.require_version('Gdk', '3.0')
from gi.repository import GLib, Gdk, Gtk

import matplotlib
matplotlib.use('GTK3Agg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3 import NavigationToolbar2GTK3 as NavigationToolbar
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas

import numpy as np
from gadgetron import Gadget

class ImageViewWindow:

    @staticmethod
    def get_instance():
        """

        :return:
        :rtype ImageViewWindow
        """
        return ImageViewWindow._instance # type: ImageViewWindow

    def __init__(self, wait_time_before_exit=0):# event: Event):

        print('Image View Init, wait before exit is %s'%wait_time_before_exit)

        self.wait_time_before_exit=wait_time_before_exit

        self.img_index=0

        self.window = Gtk.Window(Gtk.WindowType.TOPLEVEL)
        self.window.connect("delete_event", self.delete_event)
        self.window.connect('key_press_event', self.on_key_press_event)
        self.window.set_default_size(800,600)
        self.window.set_title("Gadgetron Real Time Image Viewer")

        self.vbox = Gtk.VBox()
        self.window.add(self.vbox)

        self.fig = Figure(figsize=(5,4), dpi=100)

        #plt.gray()

        self.ax = self.fig.add_subplot(111)

        #self.img_ax = self.ax.imshow(np.squeeze(np.abs(img_data)))

        self.progress=Gtk.ProgressBar(show_text=True)
        self.vbox.pack_start(self.progress, False, False, 0)

        self.canvas = FigureCanvas(self.fig)  # a gtk.DrawingArea
        self.vbox.pack_start(self.canvas, True, True, 0)

        self.toolbar = NavigationToolbar(self.canvas, self.window)
        self.vbox.pack_start(self.toolbar, False, False, 0)
        self.window.show_all()

        ImageViewWindow._instance=self

        self.has_closed=False

        self.ui_quit_event=None # type: Event

    def wait_to_exit(self, ui_quit_event: Event):
        if not self.has_closed:
            GLib.idle_add(self.wait_to_exit_impl, ui_quit_event)
        else:
            ui_quit_event.set()
        pass

    def wait_to_exit_impl(self, ui_quit_event):
        self.ui_quit_event=ui_quit_event # mark and wait delete event callback to quit and set!
        self.window.close()
        pass

    def update_data_async(self, img_data):
        if not self.has_closed:
            GLib.idle_add(self.update_data_impl, img_data)
        else:
            print('ui has closed!')

    def update_data_impl(self, img_data):
        #print('image_data_impl %s\n'%self.img_index)
        self.img_index+=1

        self.progress.set_text(str(self.img_index))
        progress_new_value=self.progress.get_fraction()+1/61 # TODO 1/61 for what?
        self.progress.set_fraction(progress_new_value)

        self.ax.clear()
        self.img_ax = self.ax.imshow(np.squeeze(np.abs(img_data)),cmap='gray')

        # if use draw_idle the system may try to draw after window closed or when closing which will throw exception
        self.canvas.draw() # self.canvas.draw_idle()
        # prefer_width=self.window.get_preferred_width()
        pass

    def real_exit_impl(self):
        #TODO self.canvas.destroy() needed?
        Gtk.main_quit()
        pass

    def delete_event(self, widget, event, data=None):
        """
            tips: even close window by api will trigger delete_event

        :param widget:
        :param event:
        :param data:
        :return:
        """
        self.has_closed=True
        if self.ui_quit_event is not None:
            self.ui_quit_event.set()

        print('wait for %s second before exit'%self.wait_time_before_exit)
        GLib.timeout_add(int(self.wait_time_before_exit), self.real_exit_impl)

        return False # zero means succ?

    def on_key_press_event(self, widget, event, data=None):
        keyname = Gdk.keyval_name(event.keyval)
        if (keyname == "Escape"): # TODO how to deal with data arrive lately one we close main window
            self.window.close()
            return False

    @staticmethod
    def StartImpl(wait_time_before_exit, event):
        iv=ImageViewWindow(wait_time_before_exit)
        event.set()
        Gtk.main()
        pass


    @staticmethod
    def wait_to_start(wait_time_before_exit):
        ui_start_finish_event=Event()

        t=Thread(
            target=lambda: ImageViewWindow.StartImpl(wait_time_before_exit, ui_start_finish_event),
            daemon=True
        )
        t.start()

        ui_start_finish_event.wait()


    @staticmethod
    def wait_to_finish():
        print('--begin-- Wait ImageView to close and exit ui thread!')
        ui_quit_event=Event()
        ImageViewWindow.get_instance().wait_to_exit(ui_quit_event)
        ui_quit_event.wait()
        print('--end-- Wait ImageView to close and exit ui thread!')
        pass

class RealTimeImageViewerGadget(Gadget):
    ParaNameWaitBeforeExit='wait_before_exit'
    DefaultValueforWaitBeforeExit=2 # seconds
    def process_config(self, cfg):
        print("--begin-- Attempting to open window")
        ImageViewWindow.wait_to_start(
            self.params[RealTimeImageViewerGadget.ParaNameWaitBeforeExit] if RealTimeImageViewerGadget.ParaNameWaitBeforeExit in self.params else RealTimeImageViewerGadget.DefaultValueforWaitBeforeExit)
        print("--end-- Attempting to open window")


    def process(self, h,im):
        #time.sleep(0.2) # TODO
        ImageViewWindow.get_instance().update_data_async(im)
        self.put_next(h,im)
        return 0

    def end(self):
        print('--begin-- Python Gaddget end process in svc Thread!')
        ImageViewWindow.get_instance().wait_to_finish()
        print('--end-- Python Gaddget end process in svc Thread!')
        return 0
        pass


if __name__ == '__main__':
    # test by dummy data
    print('enter main thread!')

    ImageViewWindow.wait_to_start(RealTimeImageViewerGadget.DefaultValueforWaitBeforeExit+5)

    N=10
    current=0
    while current<N:
        im = np.random.random((256, 256))
        ImageViewWindow.get_instance().update_data_async(im)
        time.sleep(0.2)
        current+=1
        pass

    ImageViewWindow.wait_to_finish()

    print('exit from main thread!')

    pass
