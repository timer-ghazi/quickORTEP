# window.py

import sys
import os
import platform
from Xlib import X, XK, display

# Import canvas implementation
from .x11_basic import X11CanvasBasic


class X11Window:
    """
    A generic X11 window class that hosts an X11CanvasBase subclass.
    It handles the main event loop, resizing, and keyboard input,
    and calls the canvas's drawing methods as needed.
    """

    def __init__(self, width=800, height=600, title="X11 Demo", 
                 background_color=(255, 255, 255)):
        """
        :param width:        Initial window width (pixels)
        :param height:       Initial window height (pixels)
        :param title:        Window title
        :param background_color: RGB tuple for the background color (default: white)
        """

        # Fix for macOS XQuartz DISPLAY format
        if platform.system() == "Darwin":  # Check if running on macOS
            # If DISPLAY looks like a file path and ends with ":0", fix it
            display_var = os.environ.get("DISPLAY", "")
            if display_var.startswith("/") and display_var.endswith(":0"):
                os.environ["DISPLAY"] = ":0"
                print("Detected macOS XQuartz socket path. Setting DISPLAY=:0")

        self.display = display.Display()
        self.screen = self.display.screen()
        self.running = True
        self.background_color = background_color

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

        # Instantiate the canvas with background color
        self.canvas = X11CanvasBasic(self, background_color=background_color)


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
            new_width = evt.width
            new_height = evt.height
            self.canvas.create_or_resize(new_width, new_height)
            
            # Update the view offsets if they exist (i.e., for MoleculeViewer)
            if hasattr(self, 'view_params'):
                self.view_params.x_offset = new_width // 2
                self.view_params.y_offset = new_height // 2
            
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
                      background_color=(255, 255, 255)):
    """
    A factory function that returns an X11Window instance.
    
    :param width:     Window width in pixels
    :param height:    Window height in pixels
    :param title:     Window title
    :param background_color: RGB tuple for the background color (default: white)
    :return:          An X11Window with the basic canvas
    """

    return X11Window(width, height, title, background_color=background_color)


if __name__ == "__main__":
    """
    Simple test harness. Run 'python -m x11view.window' to test the basic canvas.
    """
    import sys

    print("Creating window with basic canvas")
    win = create_x11_window(width=800, height=600, title="Canvas Test")
    win.run()
