# base.py

"""
x11view.base

Defines the abstract base class (X11CanvasBase) for all X11-based
canvas implementations. The goal is to enforce a consistent public
API (method signatures) for drawing shapes, clearing, flushing, etc.
"""

from abc import ABC, abstractmethod
from Xlib.display import Display
from Xlib import X


class X11CanvasBase(ABC):
    """
    An abstract base class defining the public drawing interface
    for X11-based canvas implementations (e.g., basic vs supersampled).

    Subclasses must implement the abstract methods in a way that
    ultimately draws onto an X11 drawable (Pixmap or similar) and
    then copies that result to the real window on flush().
    """

    def __init__(self, x11_window: "X11Window"):
        """
        Initialize the base class with references to the parent X11Window.

        :param x11_window: The X11Window instance that owns this canvas.
                           This gives us access to the display, screen, etc.
        """
        self.x11_window = x11_window
        self.display: Display = x11_window.display
        self.screen = x11_window.screen

        # Query initial geometry of the window.
        geom = x11_window.window.get_geometry()
        self.width = geom.width
        self.height = geom.height

    @abstractmethod
    def create_or_resize(self, width: int, height: int) -> None:
        """
        Handle resizing logic for this canvas implementation.

        This is typically called whenever the parent window is resized.
        Implementations should allocate or reallocate any buffers (e.g.
        Pixmaps, PIL images) that depend on the window size.

        :param width:  New width of the window
        :param height: New height of the window
        """
        pass

    @abstractmethod
    def clear(self) -> None:
        """
        Clear the drawing surface (e.g. fill with a background color).
        Subsequent drawing calls will appear on this cleared surface.
        """
        pass

    @abstractmethod
    def flush(self) -> None:
        """
        Copy or finalize the drawn content onto the actual X11 window.

        For a double-buffered approach, this typically copies
        the offscreen buffer (Pixmap or image) to the real window,
        then flushes the X11 display.
        """
        pass

    # -----------------------------------------------------------------
    # Below are the abstract drawing methods that each implementation
    # should provide. Adjust as needed for your final set of shapes.
    # -----------------------------------------------------------------

    @abstractmethod
    def draw_filled_circle(self,
                           cx: int,
                           cy: int,
                           radius: int,
                           color: tuple[int, int, int] = (255, 0, 0)
                           ) -> None:
        """
        Draw a filled circle at coordinates (cx, cy) with the specified radius and color.

        :param cx:    Center x-coordinate
        :param cy:    Center y-coordinate
        :param radius: Circle radius in pixels
        :param color: (R, G, B) tuple, each component 0..255
        """
        pass

    @abstractmethod
    def draw_circle_border(self,
                           cx: int,
                           cy: int,
                           radius: int,
                           color: tuple[int, int, int] = (0, 0, 0),
                           thickness: int = 2
                           ) -> None:
        """
        Draw an unfilled circle (border) at (cx, cy).

        :param cx:       Center x-coordinate
        :param cy:       Center y-coordinate
        :param radius:   Circle radius
        :param color:    (R, G, B) tuple, each component 0..255
        :param thickness: Line thickness in pixels
        """
        pass

    @abstractmethod
    def draw_line(self,
                  x1: int,
                  y1: int,
                  x2: int,
                  y2: int,
                  thickness: int = 4,
                  color: tuple[int, int, int] = (0, 0, 0)
                  ) -> None:
        """
        Draw a solid line from (x1, y1) to (x2, y2).

        :param x1:        Starting x-coordinate
        :param y1:        Starting y-coordinate
        :param x2:        Ending x-coordinate
        :param y2:        Ending y-coordinate
        :param thickness: Line thickness in pixels
        :param color:     (R, G, B) tuple, each component 0..255
        """
        pass

    @abstractmethod
    def draw_dashed_line(self,
                         x1: int,
                         y1: int,
                         x2: int,
                         y2: int,
                         thickness: int = 2,
                         color: tuple[int, int, int] = (128, 128, 128)
                         ) -> None:
        """
        Draw a dashed line from (x1, y1) to (x2, y2).

        :param x1:        Starting x-coordinate
        :param y1:        Starting y-coordinate
        :param x2:        Ending x-coordinate
        :param y2:        Ending y-coordinate
        :param thickness: Line thickness in pixels
        :param color:     (R, G, B) tuple, each component 0..255
        """
        pass

    @abstractmethod
    def draw_arc(self,
                 cx: int,
                 cy: int,
                 rx: int,
                 ry: int = None,
                 angle_start_deg: float = 0.0,
                 angle_end_deg: float = 360.0,
                 thickness: int = 2,
                 color: tuple[int, int, int] = (0, 0, 0)
                 ) -> None:
        """
        Draw an unfilled arc or ellipse boundary.

        :param cx:             Center x-coordinate
        :param cy:             Center y-coordinate
        :param rx:             Radius in x direction
        :param ry:             Radius in y direction (if None, use rx)
        :param angle_start_deg: Start angle (degrees)
        :param angle_end_deg:   End angle (degrees)
        :param thickness:      Line thickness in pixels
        :param color:          (R, G, B) tuple, each component 0..255
        """
        pass

    @abstractmethod
    def draw_filled_arc(self,
                        cx: int,
                        cy: int,
                        rx: int,
                        ry: int = None,
                        angle_start_deg: float = 0.0,
                        angle_end_deg: float = 360.0,
                        color: tuple[int, int, int] = (0, 0, 0)
                        ) -> None:
        """
        Draw a filled wedge or elliptical sector.

        :param cx:             Center x-coordinate
        :param cy:             Center y-coordinate
        :param rx:             Radius in x direction
        :param ry:             Radius in y direction (if None, use rx)
        :param angle_start_deg: Start angle (degrees)
        :param angle_end_deg:   End angle (degrees)
        :param color:          (R, G, B) tuple, each component 0..255
        """
        pass

    @abstractmethod
    def draw_text(self,
                  x: int,
                  y: int,
                  text: str,
                  color: tuple[int, int, int] = (0, 0, 0),
                  font_size: int = 12) -> None:
        pass
    
    @abstractmethod
    def draw_rect(self,
                  x: int,
                  y: int,
                  width: int,
                  height: int,
                  color: tuple[int, int, int] = (0, 0, 0),
                  thickness: int = 2) -> None:
        pass
    
    @abstractmethod
    def draw_filled_rect(self,
                         x: int,
                         y: int,
                         width: int,
                         height: int,
                         color: tuple[int, int, int] = (0, 0, 0)
                         ) -> None:
        pass

    @abstractmethod
    def draw_triangle(self,
                      x1: int, y1: int,
                      x2: int, y2: int,
                      x3: int, y3: int,
                      color: tuple[int, int, int] = (0, 0, 0),
                      thickness: int = 2) -> None:
        """
        Draw an unfilled triangle with vertices at (x1,y1), (x2,y2), and (x3,y3).

        :param x1, y1: First vertex coordinates
        :param x2, y2: Second vertex coordinates
        :param x3, y3: Third vertex coordinates
        :param color: (R, G, B) tuple, each component 0..255
        :param thickness: Line thickness in pixels
        """
        pass

    @abstractmethod
    def draw_filled_triangle(self,
                            x1: int, y1: int,
                            x2: int, y2: int,
                            x3: int, y3: int,
                            color: tuple[int, int, int] = (0, 0, 0)) -> None:
        """
        Draw a filled triangle with vertices at (x1,y1), (x2,y2), and (x3,y3).

        :param x1, y1: First vertex coordinates
        :param x2, y2: Second vertex coordinates
        :param x3, y3: Third vertex coordinates
        :param color: (R, G, B) tuple, each component 0..255
        """
        pass
        
    @abstractmethod
    def get_text_dimensions(self,
                           text: str,
                           font_size: int = 12) -> tuple[float, float]:
        """
        Estimate the dimensions of the given text in pixels.
        
        Parameters:
            text: The text string to measure
            font_size: The font size in pixels
            
        Returns:
            (width, height) tuple in pixels
        """
        pass