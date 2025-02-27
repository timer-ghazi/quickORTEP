#!/usr/bin/env python3
"""
Refactored graph_viewer.py
==========================

- Minimizes repeated min/max computations by storing xmin, xmax, ymin, ymax.
- Factors out duplicate logic for drawing ticks in minimal vs full mode.
- Merges offset calculations with theme logic.
- Simplifies the test harness.
- Allows a custom theme to be passed in so you don't need to override after construction.
- Adds an update_data() method for swapping xdata/ydata on the fly.
"""

import sys
from x11view.window import create_x11_window
from trajectory import Trajectory
from config import MINIMAL_THEME, FULL_THEME

class GraphViewer:
    def __init__(
        self,
        canvas,
        xdata,
        ydata,
        mode="minimal",
        current_frame=0,
        region_x=0,
        region_y=0,
        region_width=None,
        region_height=None,
        x_axis_title="Frame",
        y_axis_title="Energy",
        title=None,
        theme=None
    ):
        """
        Initialize the graph viewer with given data and configuration.

        Parameters:
            canvas: Drawing canvas
            xdata: List of x values
            ydata: List of y values
            mode: 'minimal' or 'full' display mode
            current_frame: Index of the current frame to highlight
            region_x, region_y: Top-left corner of graph region
            region_width, region_height: Size of graph region
            x_axis_title, y_axis_title: Titles for the axes
            title: Optional title to display above the graph
            theme: Optional custom theme dict. If None, the theme is chosen
                   based on 'mode' (MINIMAL_THEME or FULL_THEME).
        """
        self.canvas = canvas
        self.xdata = xdata
        self.ydata = ydata
        self.mode = mode
        self.current_frame = current_frame
        self.title = title

        # Store min/max once to avoid repeated min/max computations
        self.xmin = min(self.xdata)
        self.xmax = max(self.xdata)
        self.ymin = min(self.ydata)
        self.ymax = max(self.ydata)

        # Graph region
        self.region_x = region_x
        self.region_y = region_y
        self.region_width = region_width if region_width is not None else canvas.width
        self.region_height = region_height if region_height is not None else canvas.height

        # Decide on a theme
        if theme is not None:
            self.theme = theme
        else:
            self.theme = MINIMAL_THEME if self.mode == "minimal" else FULL_THEME

        # Axis titles
        self.x_axis_title = x_axis_title
        self.y_axis_title = y_axis_title

        # Apply offsets/margins based on theme and data range
        self._apply_theme_offsets()

    def _apply_theme_offsets(self):
        """
        Set margins and label offsets based on data range and the chosen theme.
        This replaces the old _calculate_offsets() logic and merges theme logic.
        """
        # Unpack margins from the theme
        self.left_margin = self.theme["margin"]["left"]
        self.right_margin = self.theme["margin"]["right"]
        self.top_margin = self.theme["margin"]["top"]
        self.bottom_margin = self.theme["margin"]["bottom"]

        # For the indicator radius, set a default if not specified
        if "indicator_radius" not in self.theme:
            if self.mode == "minimal":
                self.theme["indicator_radius"] = 4
            else:
                self.theme["indicator_radius"] = 8

        # Approximate label width calculation to adjust left margin if needed
        y_min_str = f"{self.ymin:.2f}"
        y_max_str = f"{self.ymax:.2f}"
        sample_labels = [y_min_str, y_max_str]
        max_label_len = max(len(label) for label in sample_labels)
        char_width = 7  # Approx px per character

        # Adjust left margin if label is large
        if max_label_len > 8:
            self.left_margin = max(self.left_margin, 60)

    def _map_x(self, x):
        """Map data x value to pixel x coordinate."""
        if self.xmin == self.xmax:
            return int(
                self.region_x
                + self.left_margin
                + (self.region_width - self.left_margin - self.right_margin) / 2
            )
        scale = (self.region_width - self.left_margin - self.right_margin) / (self.xmax - self.xmin)
        return int(self.region_x + self.left_margin + (x - self.xmin) * scale)

    def _map_y(self, y):
        """Map data y value to pixel y coordinate."""
        if self.ymin == self.ymax:
            return int(
                self.region_y
                + self.top_margin
                + (self.region_height - self.top_margin - self.bottom_margin) / 2
            )
        scale = (self.region_height - self.top_margin - self.bottom_margin) / (self.ymax - self.ymin)
        return int(self.region_y + self.top_margin + (self.ymax - y) * scale)

    def compute_ticks(self, min_val, max_val, n_major=5, n_minor=4):
        """
        Compute positions for major and minor tick marks.
        By default we have 5 major ticks and 4 minor intervals between them.
        """
        if n_major < 2:
            return [min_val], []
        major_interval = (max_val - min_val) / (n_major - 1)
        major_ticks = [min_val + i * major_interval for i in range(n_major)]
        minor_ticks = []
        for i in range(n_major - 1):
            start = major_ticks[i]
            end = major_ticks[i + 1]
            gap = (end - start) / (n_minor + 1) if n_minor >= 1 else 0
            for j in range(1, n_minor + 1):
                minor_ticks.append(start + j * gap)
        return major_ticks, minor_ticks

    def _draw_ticks_x(self, x_axis_y):
        """
        Draws x-axis major/minor ticks and optional labels.
        Differentiates between minimal vs. full mode.
        """
        # If minimal, just do fewer major ticks, no minor
        if self.mode == "minimal":
            major_ticks, minor_ticks = self.compute_ticks(self.xmin, self.xmax, n_major=4, n_minor=0)
            show_labels = False
        else:
            major_ticks, minor_ticks = self.compute_ticks(self.xmin, self.xmax, n_major=5, n_minor=4)
            show_labels = True

        tick_len = self.theme["tick_length"]
        minor_tick_len = self.theme.get("minor_tick_length", tick_len // 2)

        # Draw major ticks (and labels if full)
        for tick in major_ticks:
            x_pixel = self._map_x(tick)
            # Major tick
            self.canvas.draw_line(
                x_pixel, x_axis_y,
                x_pixel, x_axis_y + tick_len,
                thickness=self.theme["line_thickness"]["tick"],
                color=self.theme["tick_color"]
            )
            if show_labels:
                label = f"{int(round(tick))}"
                label_width = len(label) * 7  # approximate
                self.canvas.draw_text(
                    x_pixel - (label_width // 2),
                    x_axis_y + tick_len + 12,
                    text=label,
                    color=self.theme["font_color"],
                    font_size=self.theme["font_size"]
                )

        # Draw minor ticks
        for tick in minor_ticks:
            x_pixel = self._map_x(tick)
            self.canvas.draw_line(
                x_pixel, x_axis_y,
                x_pixel, x_axis_y + minor_tick_len,
                thickness=self.theme["line_thickness"]["tick"],
                color=self.theme["tick_color"]
            )

    def _draw_ticks_y(self, y_axis_x):
        """
        Draws y-axis major/minor ticks and optional labels.
        Differentiates between minimal vs. full mode.
        """
        if self.mode == "minimal":
            major_ticks, minor_ticks = self.compute_ticks(self.ymin, self.ymax, n_major=4, n_minor=0)
            show_labels = False
        else:
            major_ticks, minor_ticks = self.compute_ticks(self.ymin, self.ymax, n_major=5, n_minor=4)
            show_labels = True

        tick_len = self.theme["tick_length"]
        minor_tick_len = self.theme.get("minor_tick_length", tick_len // 2)

        for tick in major_ticks:
            y_pixel = self._map_y(tick)
            self.canvas.draw_line(
                y_axis_x - tick_len, y_pixel,
                y_axis_x, y_pixel,
                thickness=self.theme["line_thickness"]["tick"],
                color=self.theme["tick_color"]
            )
            if show_labels:
                label = f"{tick:.2f}"
                # Approx right alignment
                label_x = y_axis_x - (tick_len + 8 + (len(label)*7))
                self.canvas.draw_text(
                    label_x, y_pixel + 4,  # small offset for vertical centering
                    text=label,
                    color=self.theme["font_color"],
                    font_size=self.theme["font_size"]
                )

        for tick in minor_ticks:
            y_pixel = self._map_y(tick)
            self.canvas.draw_line(
                y_axis_x - minor_tick_len, y_pixel,
                y_axis_x, y_pixel,
                thickness=self.theme["line_thickness"]["tick"],
                color=self.theme["tick_color"]
            )

    def draw_axes(self):
        """Draw the axes (with ticks, labels, and titles if applicable)."""
        reg_x, reg_y = self.region_x, self.region_y
        reg_w, reg_h = self.region_width, self.region_height
        L, R, T, B = self.left_margin, self.right_margin, self.top_margin, self.bottom_margin

        # Draw title
        if self.title:
            # center horizontally
            title_x = reg_x + L + (reg_w - L - R) // 2
            if self.mode == "minimal":
                title_y = reg_y + 5
                font_size = max(8, self.theme.get("font_size", 10) - 2)
            else:
                title_y = reg_y - 5
                font_size = self.theme.get("font_size", 12)

            title_width = len(self.title) * (font_size // 2)
            self.canvas.draw_text(
                title_x - (title_width // 2),
                title_y,
                text=self.title,
                color=self.theme["font_color"],
                font_size=font_size
            )

        # X axis at the bottom
        x_axis_y = reg_y + reg_h - B
        self.canvas.draw_line(
            reg_x + L, x_axis_y,
            reg_x + reg_w - R, x_axis_y,
            thickness=self.theme["line_thickness"]["axis"],
            color=self.theme["axis_color"]
        )

        # Y axis at the left
        y_axis_x = reg_x + L
        self.canvas.draw_line(
            y_axis_x, reg_y + T,
            y_axis_x, reg_y + reg_h - B,
            thickness=self.theme["line_thickness"]["axis"],
            color=self.theme["axis_color"]
        )

        # Draw ticks on both axes
        self._draw_ticks_x(x_axis_y)
        self._draw_ticks_y(y_axis_x)

        # Axis titles (full mode only)
        if self.mode == "full":
            # X-axis title
            x_title_x = reg_x + L + (reg_w - L - R) // 2
            x_title_y = x_axis_y + self.theme["font_size"] + 20
            x_title_width = len(self.x_axis_title) * 7
            self.canvas.draw_text(
                x_title_x - (x_title_width // 2), x_title_y,
                text=self.x_axis_title,
                color=self.theme["font_color"],
                font_size=self.theme["font_size"]
            )

            # Y-axis title - placed at the top left for simplicity
            y_title_x = y_axis_x - 25
            y_title_y = reg_y + T - 10
            self.canvas.draw_text(
                y_title_x, y_title_y,
                text=self.y_axis_title,
                color=self.theme["font_color"],
                font_size=self.theme["font_size"]
            )

    def plot_data(self):
        """Draw the data points and connecting lines."""
        points = [(self._map_x(x), self._map_y(y)) for x, y in zip(self.xdata, self.ydata)]
        for i in range(1, len(points)):
            x1, y1 = points[i - 1]
            x2, y2 = points[i]
            self.canvas.draw_line(
                int(x1), int(y1),
                int(x2), int(y2),
                thickness=self.theme["line_thickness"]["data"],
                color=self.theme["data_line_color"]
            )

    def draw_indicator(self):
        """Draw an indicator at the current frame position."""
        try:
            x_val = self.xdata[self.current_frame]
            y_val = self.ydata[self.current_frame]
        except IndexError:
            return
        x_pixel = self._map_x(x_val)
        y_pixel = self._map_y(y_val)

        if self.mode == "minimal":
            # scale indicator radius to region size
            plot_dimension = min(self.region_width, self.region_height)
            indicator_radius = max(3, int(plot_dimension * 0.04))
        else:
            indicator_radius = self.theme["indicator_radius"]

        self.canvas.draw_filled_circle(
            x_pixel, y_pixel, indicator_radius,
            color=self.theme["indicator_color"]
        )

    def draw_graph(self):
        """Draw the complete graph."""
        self.draw_axes()
        self.plot_data()
        self.draw_indicator()

    def update_current_frame(self, new_frame):
        """Update the current frame and redraw."""
        self.current_frame = new_frame
        self.draw_graph()

    def update_data(self, xdata, ydata, x_axis_title=None, y_axis_title=None, title=None):
        """
        Replace the current x and y data with new data, optionally update axis
        labels and/or the plot title, then redraw. This is useful for switching 
        from one data set (e.g., Energy vs. Frame) to another (Distance vs. Frame)
        without creating a new GraphViewer instance.

        Parameters:
            xdata, ydata: New lists of x and y values
            x_axis_title: Optional new X-axis label
            y_axis_title: Optional new Y-axis label
            title: Optional new graph title
        """
        self.xdata = xdata
        self.ydata = ydata

        self.xmin = min(xdata)
        self.xmax = max(xdata)
        self.ymin = min(ydata)
        self.ymax = max(ydata)

        if x_axis_title is not None:
            self.x_axis_title = x_axis_title
        if y_axis_title is not None:
            self.y_axis_title = y_axis_title
        if title is not None:
            self.title = title

        # Recompute margins & offsets in case the new data range is very different
        self._apply_theme_offsets()

        # Redraw with the new data
        self.draw_graph()


# ---------------------------
# Simplified Test Harness
# ---------------------------

def build_default_graphs(window, xdata, ydata, current_frame):
    """
    Helper function to build a main (full) graph and a thumbnail (minimal) graph.
    Returns a tuple (main_graph, thumbnail_graph).
    """
    # Main Graph
    main_width = 500
    main_height = 300
    main_x = (window.canvas.width - main_width) // 2
    main_y = (window.canvas.height - main_height) // 2

    main_graph = GraphViewer(
        canvas=window.canvas,
        xdata=xdata,
        ydata=ydata,
        mode="full",
        current_frame=current_frame,
        region_x=main_x,
        region_y=main_y,
        region_width=main_width,
        region_height=main_height,
        x_axis_title="Frame",
        y_axis_title="Energy (a.u.)",
        title="Energy Profile"  # optional
    )

    # Thumbnail (Minimal) Graph
    thumb_width = 150
    thumb_height = 80
    thumb_margin = 20
    thumb_x = window.canvas.width - thumb_width - thumb_margin
    thumb_y = window.canvas.height - thumb_height - thumb_margin

    # We can pass in a slightly tweaked minimal theme
    custom_minimal_theme = MINIMAL_THEME.copy()
    custom_minimal_theme["margin"] = {"left": 10, "right": 10, "top": 10, "bottom": 10}

    thumbnail_graph = GraphViewer(
        canvas=window.canvas,
        xdata=xdata,
        ydata=ydata,
        mode="minimal",
        current_frame=current_frame,
        region_x=thumb_x,
        region_y=thumb_y,
        region_width=thumb_width,
        region_height=thumb_height,
        x_axis_title="",
        y_axis_title="",
        title="Energy (Thumbnail)",
        theme=custom_minimal_theme
    )

    return main_graph, thumbnail_graph

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python graph_viewer.py <xyz_file> [current_frame] [ss_factor]")
        sys.exit(1)

    xyz_file = sys.argv[1]
    current_frame = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    ss_factor = int(sys.argv[3]) if len(sys.argv) > 3 else 1

    traj = Trajectory.from_xyz_file(xyz_file)
    n_frames = len(traj._raw_frames)

    x_values = list(range(n_frames))
    y_values = []
    for i in range(n_frames):
        mol = traj.get_frame(i)
        y_values.append(mol.Energy if mol.Energy is not None else 0.0)

    win = create_x11_window(width=800, height=600, title="Graph Viewer", ss_factor=ss_factor)

    # Build two graphs (main + thumbnail)
    main_graph, thumbnail_graph = build_default_graphs(win, x_values, y_values, current_frame)

    def custom_redraw():
        win.canvas.clear()
        main_graph.draw_graph()
        thumbnail_graph.draw_graph()
        win.canvas.flush()

    # Example of how you'd update both if you have a frame change:
    def update_frame(new_frame):
        main_graph.update_current_frame(new_frame)
        thumbnail_graph.update_current_frame(new_frame)

    win.redraw = custom_redraw
    win.run()
