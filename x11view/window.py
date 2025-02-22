# window.py

import sys
from Xlib import X, XK, display

# Import your two canvas implementations
from .x11_basic import X11CanvasBasic
from .x11_ss import X11CanvasSS


class X11Window:
    """
    A generic X11 window class that hosts an X11CanvasBase subclass.
    It handles the main event loop, resizing, and keyboard input,
    and calls the canvas's drawing methods as needed.
    """

    def __init__(self, width=800, height=600, title="X11 Demo", canvas_class=X11CanvasBasic):
        """
        :param width:        Initial window width (pixels)
        :param height:       Initial window height (pixels)
        :param title:        Window title
        :param canvas_class: A subclass of X11CanvasBase to use for drawing
        """
        self.display = display.Display()
        self.screen = self.display.screen()
        self.running = True

        event_mask = (X.ExposureMask |
              X.KeyPressMask |
              X.StructureNotifyMask |
              X.ButtonPressMask |
              X.ButtonReleaseMask |
              X.PointerMotionMask)

        # Create the actual X11 Window
        self.window = self.screen.root.create_window(
            x=0,
            y=0,
            width=width,
            height=height,
            border_width=0,
            depth=self.screen.root_depth,
            window_class=X.InputOutput,
            visual=self.screen.root_visual,
            colormap=self.screen.default_colormap,
            event_mask=event_mask
        )
        self.window.set_wm_name(title)
        self.window.map()

        # Instantiate the chosen canvas
        self.canvas = canvas_class(self)

    def run(self):
        """
        Main event loop.
        """
        while self.running:
            try:
                evt = self.display.next_event()
            except:
                # If the connection is closed or there's an error, exit loop.
                break
            self.handle_event(evt)

    def handle_event(self, evt):
        """
        Dispatch X11 events to appropriate methods.
        """
        if evt.type == X.Expose:
            # 'Expose' can happen multiple times; we typically redraw on the last one (count=0).
            if evt.count == 0:
                self.redraw()
        elif evt.type == X.KeyPress:
            self.handle_key(evt)
        elif evt.type == X.ButtonPress:
            self.handle_button_press(evt)
        elif evt.type == X.ButtonRelease:
            self.handle_button_release(evt)
        elif evt.type == X.MotionNotify:
            self.handle_motion(evt)
        elif evt.type == X.ConfigureNotify:
            # Window was resized
            new_width = evt.width
            new_height = evt.height
            self.canvas.create_or_resize(new_width, new_height)
            self.redraw()
        elif evt.type == X.DestroyNotify:
            # Typically triggered if user closes the window from window manager
            print("Window closed by user.")
            self.running = False

    def handle_motion(self, evt):
        pass

    def handle_button_press(self, evt):
        pass

    def handle_button_release(self, evt):
        pass

    def handle_key(self, evt):
        """
        Handle KeyPress events. Here, 'q' or 'Escape' will quit the application.
        """
        keysym = self.display.keycode_to_keysym(evt.detail, 0)
        keychar = XK.keysym_to_string(keysym)
        if keychar in ("q", "Escape"):
            print("Closing current window.")
            self.running = False
            # sys.exit(0)
        else:
            print(f"Key pressed: {keychar}")


    def redraw(self):
        self.canvas.clear()
        w = self.canvas.width
        h = self.canvas.height
    
        # Existing shapes...
        cx, cy = w // 2, h // 2
        r = min(w, h) // 6
        self.canvas.draw_filled_circle(cx, cy, r, color=(255, 0, 0))
        self.canvas.draw_circle_border(cx, cy, r + 10, color=(0, 0, 0), thickness=3)
        self.canvas.draw_dashed_line(0, cy, w, cy, thickness=2, color=(128, 128, 128))
        self.canvas.draw_arc(w // 4, h // 4, 40, ry=20, angle_start_deg=0, angle_end_deg=180,
                             thickness=2, color=(0, 0, 255))
        self.canvas.draw_filled_arc(w - 80, h - 80, 40, ry=30, angle_start_deg=180, angle_end_deg=270,
                                    color=(0, 255, 0))
    
        # Rectangles test
        self.canvas.draw_rect(10, 10, 100, 50, color=(255, 0, 255), thickness=2)
        self.canvas.draw_filled_rect(120, 10, 80, 50, color=(0, 255, 0))
    
        # ------------------------------
        # Text demonstration
        # ------------------------------
    
        # We'll define a few known "xterm" style fonts or built-in fonts
        # in ascending order of size. Some may or may not exist on your system.
        # You can check with `xlsfonts` in a terminal.
        xterm_font_list = [
            "6x13",
            "7x14",
            "8x13",
            "9x15",
            "10x20",
            # Extended XLFD forms:
            "-misc-fixed-medium-r-normal--15-140-75-75-c-90-iso8859-1",
            "-misc-fixed-medium-r-normal--20-200-75-75-c-100-iso8859-1",
        ]
    
        # We'll draw each in a new line, offset in y
        y_offset = 120
        for idx, fnt in enumerate(xterm_font_list):
            # We'll also vary color a bit
            color = (50 * idx % 256, 100 * idx % 256, 200 * idx % 256)
    
            # We'll pass the font as a single-element list for "font_candidates"
            self.canvas.draw_text(
                10, y_offset,
                text=f"This line uses '{fnt}' (if available)",
                color=color,
                font_size=12,  # We'll just ignore font_size for now in X11 core
                font_candidates=[fnt]
            )
            y_offset += 20
    
        # And let's show a line with multiple fallback fonts
        multi_fonts = [
            # Some random "rare" name that might not exist
            "-schm-rarefont-medium-r-*-13-*-*-*-*-*-iso8859-1",
            # Something else that might not exist
            "-doesnotexist-*-14-*-*-*-*-*-iso8859-1",
            # Then something that probably DOES exist
            "9x15"
        ]
    
        self.canvas.draw_text(
            10, y_offset,
            text="Trying multiple candidates, eventually falling back to '9x15'",
            color=(128, 0, 0),
            font_size=12,
            font_candidates=multi_fonts
        )
    
        self.canvas.flush()

def create_x11_window(width=800,
                      height=600,
                      title="X11 Demo",
                      ss_factor=1,
                      tile_size=128):
    """
    A factory function that returns an X11Window instance.
    If ss_factor <= 1, we use the basic canvas class.
    If ss_factor > 1, we use the supersampled canvas class.
    
    :param width:     Window width in pixels
    :param height:    Window height in pixels
    :param title:     Window title
    :param ss_factor: Supersampling factor (1 -> no SS, 2 -> 2x, etc.)
    :param tile_size: Tile size in pixels for partial put_image (only used by the SS canvas)
    :return:          An X11Window with the selected canvas
    """
    if ss_factor <= 1:
        return X11Window(width, height, title, canvas_class=X11CanvasBasic)
    else:
        # We'll define a small factory (closure) to pass ss_factor & tile_size
        def ss_canvas_factory(x11_window):
            return X11CanvasSS(
                x11_window,
                ss_factor=ss_factor,
                tile_size=tile_size
            )
        return X11Window(width, height, title, canvas_class=ss_canvas_factory)


if __name__ == "__main__":
    """
    Simple test harness. Let users run 'python -m x11view.window [ss_factor]'
    to test basic or SS mode. Example:
      python -m x11view.window 1     -> Basic
      python -m x11view.window 2     -> 2x Supersampling
    """
    import sys

    if len(sys.argv) > 1:
        try:
            sf = int(sys.argv[1])
        except ValueError:
            sf = 1
    else:
        sf = 1

    print(f"Creating window with ss_factor={sf}")
    win = create_x11_window(width=800, height=600, title="Canvas Test", ss_factor=sf)
    win.run()
