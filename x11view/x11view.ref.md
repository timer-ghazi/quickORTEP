
# `x11view` Package

A lightweight Python package for creating basic or supersampled X11 windows for 2D drawing (circles, arcs, lines, etc.). It uses python-xlib for direct X11 window and Pixmap operations, with an optional Pillow-based supersampling pipeline.

**Contents**:

1. [x11view.base](#x11viewbase)  
   - [Class: X11CanvasBase](#class-x11canvasbase)

2. [x11view.common](#x11viewcommon)  
   - [Class: X11CanvasCommon](#class-x11canvascommon)

3. [x11view.x11_basic](#x11viewx11_basic)  
   - [Class: X11CanvasBasic](#class-x11canvasbasic)

4. [x11view.x11_ss](#x11viewx11_ss)  
   - [Class: X11CanvasSS](#class-x11canvasss)

5. [x11view.window](#x11viewwindow)  
   - [Class: X11Window](#class-x11window)  
   - [Function: create_x11_window](#function-create_x11_window)

---

## x11view.base

### Class: `X11CanvasBase`

**Location**: `x11view/base.py`  
**Inheritance**: `ABC` (from Python’s `abc` module)

An **abstract base class** defining the common interface for any X11-based canvas implementation. Subclasses must implement the abstract methods below, which handle:

- Window resizing  
- Clearing the drawing surface  
- Flushing the final image to the screen  
- Drawing basic shapes (circles, arcs, lines, etc.)

**Constructor**  
```python
def __init__(self, x11_window: "X11Window"):
    """
    Initialize the canvas with a reference to the parent X11Window.
    Provides .display, .screen, .width, .height for convenience.
    """
```
- **Parameters**:  
  - `x11_window`: The `X11Window` that owns this canvas.

**Properties**  
- `display`: A reference to the python-xlib Display object.  
- `screen`: A reference to the python-xlib Screen object.  
- `width`: The current window width in pixels.  
- `height`: The current window height in pixels.

**Abstract Methods**  
All the following methods must be implemented by child classes:

1. ```python
   create_or_resize(self, width: int, height: int) -> None
   ```
   Handle window resizing. Typically reallocate any internal Pixmaps or images to match new dimensions.

2. ```python
   clear(self) -> None
   ```
   Clear or fill the drawing surface with a background color.

3. ```python
   flush(self) -> None
   ```
   Copy the (possibly double-buffered) contents to the actual window and flush the X11 display.

4. ```python
   draw_filled_circle(self, cx: int, cy: int, radius: int, color=(255, 0, 0)) -> None
   ```

5. ```python
   draw_circle_border(self, cx: int, cy: int, radius: int,
                      color=(0, 0, 0), thickness=2) -> None
   ```

6. ```python
   draw_line(self, x1: int, y1: int, x2: int, y2: int,
             thickness: int = 4, color=(0, 0, 0)) -> None
   ```

7. ```python
   draw_dashed_line(self, x1: int, y1: int, x2: int, y2: int,
                    thickness: int = 2, color=(128, 128, 128)) -> None
   ```

8. ```python
   draw_arc(self, cx: int, cy: int, rx: int, ry: int = None,
            angle_start_deg: float = 0.0, angle_end_deg: float = 360.0,
            thickness: int = 2, color=(0, 0, 0)) -> None
   ```

9. ```python
   draw_filled_arc(self, cx: int, cy: int, rx: int, ry: int = None,
                   angle_start_deg: float = 0.0, angle_end_deg: float = 360.0,
                   color=(0, 0, 0)) -> None
   ```

---

## x11view.common

### Class: `X11CanvasCommon`

**Location**: `x11view/common.py`  
**Inheritance**: Simple mixin class (no abstract base).

A small **mixin**/utility class that provides:

1. **GC caching**: A dictionary that maps `(color, thickness, line_style, fill_style)` → `GC`.  
2. **`rgb_to_pixel`**: A helper to convert `(R, G, B)` in `0..255` to a 24-bit integer.

Typically, your canvas classes inherit from both `X11CanvasBase` and `X11CanvasCommon` so they gain the `get_gc` functionality without duplicating code.

**Constructor**  
```python
def __init__(self):
    """
    Initialize the GC cache dictionary, self.gc_cache.
    """
```

**Methods**  
1. ```python
   def get_gc(self, color, thickness, line_style, fill_style=False) -> X.GC:
       """
       Return a cached GC for the given parameters if available;
       otherwise create it, store it in gc_cache, and return it.
       """
   ```
   - **Parameters**:  
     - `color`: (R, G, B) tuple.  
     - `thickness`: line width in pixels (if relevant).  
     - `line_style`: e.g., `X.LineSolid`, `X.LineOnOffDash`.  
     - `fill_style`: boolean; if True, thickness/line_style may be ignored in the GC.

2. ```python
   @staticmethod
   def rgb_to_pixel(rgb: Tuple[int, int, int]) -> int:
       """
       Convert an RGB triple into a 24-bit pixel value: 0xRRGGBB.
       """
   ```

---

## x11view.x11_basic

### Class: `X11CanvasBasic`

**Location**: `x11view/x11_basic.py`  
**Inheritance**: `X11CanvasBase`, `X11CanvasCommon`

A **direct X11** (non-supersampled) canvas that uses a Pixmap for double-buffering:

1. **Constructor** creates an offscreen Pixmap at the current window size.  
2. **All drawing** calls go to that offscreen Pixmap.  
3. **`flush()`** copies the Pixmap to the real window.  

**Constructor**  
```python
def __init__(self, x11_window: X11Window):
    """
    Creates an offscreen Pixmap of the current window size.
    Also initializes GC cache.
    """
```

**Key Methods**  
(Implementations of the abstract methods from `X11CanvasBase`):

- `create_or_resize(width, height)`: Reallocate the Pixmap at new dimensions and clear the GC cache.  
- `clear()`: Fill the Pixmap with white (or a background color).  
- `flush()`: Copy the Pixmap to the real window via `copy_area`.  
- `draw_filled_circle(...)`, `draw_circle_border(...)`, etc.: All the shape methods call `get_gc(...)` for a GC, then use python-xlib drawing calls like `fill_arc`, `poly_arc`, `poly_line`.

Example usage:
```python
canvas = X11CanvasBasic(my_window)
canvas.clear()
canvas.draw_filled_circle(100, 100, 50, color=(255,0,0))
canvas.flush()
```

---

## x11view.x11_ss

### Class: `X11CanvasSS`

**Location**: `x11view/x11_ss.py`  
**Inheritance**: `X11CanvasBase`, `X11CanvasCommon`

A **Pillow-based** supersampled canvas:

1. **Constructor** creates a high-resolution PIL `Image` (size = `width * ss_factor` by `height * ss_factor`), plus a standard Pixmap for the final displayed image.  
2. **Drawing** (circles, lines, arcs) is performed on the PIL image at scaled coordinates.  
3. **`flush()`** downscales the image back to the actual window size, then tiles it into the Pixmap via `put_image`, then copies to the window.

**Constructor**  
```python
def __init__(self,
             x11_window,
             ss_factor: int = 2,
             tile_size: int = 128):
    """
    :param x11_window: The parent X11Window
    :param ss_factor:  Supersampling factor (2 => 2x in each dimension)
    :param tile_size:  Size of tiles for put_image to balance speed/memory
    """
```

**Key Methods**  
- `create_or_resize(width, height)`: Recreate the high-res image and the final Pixmap to match new dimensions.  
- `clear()`: Fill the PIL buffer and final Pixmap with white.  
- `flush()`: Downsample the PIL image, convert to BGRA, tile-copy to the Pixmap, and then `copy_area` to the real window.  

**Supersampled Drawing**  
Methods like `draw_filled_circle(cx, cy, radius, color=...)` multiply all coordinates by `ss_factor` before calling Pillow’s `ellipse`, `line`, or `pieslice`. This yields a smoother final image.

---

## x11view.window

### Class: `X11Window`

**Location**: `x11view/window.py`  
**Purpose**: The main X11 window + event loop, which holds a reference to one of the canvas classes (`X11CanvasBasic` or `X11CanvasSS`).

**Constructor**  
```python
def __init__(self,
             width=800,
             height=600,
             title="X11 Demo",
             canvas_class=X11CanvasBasic):
    """
    Creates a new X11 window using python-xlib, sets up event masks,
    and instantiates the chosen canvas class (basic or SS).
    """
```
- **Parameters**:
  - `width`/`height`: Window dimensions.
  - `title`: A string for the window manager title bar.
  - `canvas_class`: A subclass of `X11CanvasBase` (defaults to `X11CanvasBasic`).

**Key Methods**  
1. `run()`: The main loop that calls `display.next_event()` and dispatches events.  
2. `handle_event(evt)`: Dispatch function for X11 events like `Expose`, `KeyPress`, and `ConfigureNotify`.  
3. `handle_key(evt)`: Helper to handle keyboard events; `'q'` or `'Escape'` quits by default.  
4. `redraw()`: Called after `Expose` or resize events. You can override or modify this method to define what shapes to draw. By default, draws a small demonstration of circles, lines, and arcs using `self.canvas`.

---

### Function: `create_x11_window`

**Signature**  
```python
def create_x11_window(width=800,
                      height=600,
                      title="X11 Demo",
                      ss_factor=1,
                      tile_size=128) -> X11Window:
    ...
```

**Usage**  
```python
win = create_x11_window(800, 600, title="My Window", ss_factor=2)
win.run()
```

**Behavior**  
- If `ss_factor <= 1`, returns an `X11Window` that uses `X11CanvasBasic`.  
- If `ss_factor > 1`, returns an `X11Window` that uses `X11CanvasSS`, with the specified supersampling factor.  
- `tile_size` affects how the `X11CanvasSS` does tiled `put_image`. It’s ignored by the basic canvas.

**Parameters**  
- `width`, `height`: Dimensions of the created window.  
- `title`: String for the window’s title bar.  
- `ss_factor`: Supersampling factor; `1` means no supersampling.  
- `tile_size`: The chunk size in pixels when uploading the final downsampled image to X11 (only used by the SS canvas).

**Returns**  
- An instance of `X11Window` that can be `.run()` for the event loop.

---

## Example Usage

After installing or placing this package on your PYTHONPATH, you can run:

```bash
python -m x11view.window 1
```
to test the **basic** canvas, or:
```bash
python -m x11view.window 2
```
for **2× supersampling** (using Pillow). This will open a window with some test shapes drawn in the `redraw()` method. Press `q` or `Escape` to quit.

---

## Summary

- **`X11CanvasBase`** (abstract interface for any X11 canvas)  
- **`X11CanvasCommon`** (GC caching, color utilities)  
- **`X11CanvasBasic`** (direct double-buffered Pixmap, no supersampling)  
- **`X11CanvasSS`** (Pillow-based supersampling approach)  
- **`X11Window`** (the event loop & window creation)  
- **`create_x11_window`** (simple factory to pick basic vs. SS approach)

This design offers a **modular** and **extensible** way to do 2D rendering in an X11 environment, with optional antialiasing (supersampling) provided by Pillow.
