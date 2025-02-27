#!/usr/bin/env python3
"""
graph_viewer.py
===============

Revised to avoid label overlap on the y-axis by properly positioning labels and titles.
"""

import sys
from x11view.window import create_x11_window
from trajectory import Trajectory
from config import MINIMAL_THEME, FULL_THEME

class GraphViewer:
    def __init__(self, canvas, xdata, ydata, mode="minimal", current_frame=0,
                 region_x=0, region_y=0, region_width=None, region_height=None,
                 x_axis_title="Frame", y_axis_title="Energy"):
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
        """
        self.canvas = canvas
        self.xdata = xdata
        self.ydata = ydata
        self.mode = mode
        self.current_frame = current_frame

        # Region defaults to full canvas if not provided.
        self.region_x = region_x
        self.region_y = region_y
        self.region_width = region_width if region_width is not None else canvas.width
        self.region_height = region_height if region_height is not None else canvas.height

        # Choose theme
        self.theme = MINIMAL_THEME if self.mode == "minimal" else FULL_THEME
        
        # Set default indicator sizes if not in theme
        if "indicator_radius" not in self.theme:
            if self.mode == "minimal":
                self.theme["indicator_radius"] = 4  # Default small indicator for minimal mode
            else:
                self.theme["indicator_radius"] = 8  # Default larger indicator for full mode

        # Unpack margins from the theme
        self.left_margin = self.theme["margin"]["left"]
        self.right_margin = self.theme["margin"]["right"]
        self.top_margin = self.theme["margin"]["top"]
        self.bottom_margin = self.theme["margin"]["bottom"]

        # Axis titles
        self.x_axis_title = x_axis_title
        self.y_axis_title = y_axis_title

        # Calculate appropriate offsets based on the data range
        self._calculate_offsets()

    def _calculate_offsets(self):
        """
        Calculate appropriate offsets for labels and titles based on data range.
        This ensures that labels will have enough space regardless of the data values.
        """
        # Calculate the max length of y-labels in pixels
        y_min, y_max = min(self.ydata), max(self.ydata)
        # Format labels as they would appear on the axis
        sample_labels = [f"{y_min:.2f}", f"{y_max:.2f}"]
        # Estimate max chars × approx px per char for font
        max_label_len = max(len(label) for label in sample_labels)
        # Rough estimate: 7px per character in the label at default font size
        char_width = 7  # Adjust based on your font
        
        # For scientific values like -1107.05, we need extra space
        # Check if values are large with decimals (which take a lot of space)
        large_decimal_values = (abs(y_min) > 100 and '.' in f"{y_min:.2f}") or (abs(y_max) > 100 and '.' in f"{y_max:.2f}")
        
        # Add extra space if scientific notation might be used or for large decimal values
        scientific_notation = abs(y_min) < 0.0001 or abs(y_min) > 10000 or abs(y_max) < 0.0001 or abs(y_max) > 10000
        
        # Set offsets based on calculated values
        # Reduced offset: tick length + label width + smaller safety margin
        # We're reducing the base padding and safety margin to move labels closer to axis
        self.y_label_left_offset = 10 + (max_label_len * char_width) + 5
        
        # Reduced extra padding for scientific notation or large decimal values
        if scientific_notation:
            self.y_label_left_offset += 10
        if large_decimal_values:
            self.y_label_left_offset += 8
            
        # Y-axis title doesn't need much horizontal offset anymore since we're placing it at the top
        self.y_title_left_offset = 25
        
        # For very large numbers, increase left margin if needed
        if max_label_len > 8:
            self.left_margin = max(self.left_margin, 60)

    def _map_x(self, x):
        """Map data x value to pixel x coordinate."""
        x_min, x_max = min(self.xdata), max(self.xdata)
        # Handle edge case where x_min equals x_max
        if x_min == x_max:
            return int(self.region_x + self.left_margin + (self.region_width - self.left_margin - self.right_margin) / 2)
            
        scale = (self.region_width - self.left_margin - self.right_margin) / (x_max - x_min)
        return int(self.region_x + self.left_margin + (x - x_min) * scale)

    def _map_y(self, y):
        """Map data y value to pixel y coordinate."""
        y_min, y_max = min(self.ydata), max(self.ydata)
        # Handle edge case where y_min equals y_max
        if y_min == y_max:
            return int(self.region_y + self.top_margin + (self.region_height - self.top_margin - self.bottom_margin) / 2)
            
        scale = (self.region_height - self.top_margin - self.bottom_margin) / (y_max - y_min)
        return int(self.region_y + self.top_margin + (y_max - y) * scale)

    def compute_ticks(self, min_val, max_val, n_major, n_minor):
        """Compute positions for major and minor tick marks."""
        if n_major < 2:
            return [min_val], []
        major_interval = (max_val - min_val) / (n_major - 1)
        major_ticks = [min_val + i * major_interval for i in range(n_major)]
        minor_ticks = []
        for i in range(n_major - 1):
            start = major_ticks[i]
            end = major_ticks[i + 1]
            gap = (end - start) / (n_minor + 1)
            for j in range(1, n_minor + 1):
                minor_ticks.append(start + j * gap)
        return major_ticks, minor_ticks

    def draw_axes(self):
        """Draw the x and y axes with ticks and labels."""
        reg_x, reg_y = self.region_x, self.region_y
        reg_w, reg_h = self.region_width, self.region_height
        L, R, T, B = self.left_margin, self.right_margin, self.top_margin, self.bottom_margin

        # X axis line at bottom
        x_axis_y = reg_y + reg_h - B
        self.canvas.draw_line(reg_x + L, x_axis_y, reg_x + reg_w - R, x_axis_y,
                              thickness=self.theme["line_thickness"]["axis"],
                              color=self.theme["axis_color"])
        # Y axis line at left
        self.canvas.draw_line(reg_x + L, reg_y + T, reg_x + L, reg_y + reg_h - B,
                              thickness=self.theme["line_thickness"]["axis"],
                              color=self.theme["axis_color"])

        if self.mode == "minimal":
            # For minimal mode, have 4 ticks but no labels
            x_min, x_max = min(self.xdata), max(self.xdata)
            y_min, y_max = min(self.ydata), max(self.ydata)
            
            # Create 4 evenly spaced ticks
            x_range = x_max - x_min
            y_range = y_max - y_min
            
            xticks = [
                x_min,
                x_min + x_range * 0.33,
                x_min + x_range * 0.66,
                x_max
            ]
            
            yticks = [
                y_min,
                y_min + y_range * 0.33,
                y_min + y_range * 0.66,
                y_max
            ]
            
            minor_xticks = []
            minor_yticks = []
        else:
            # 5 major ticks, 4 minor for full mode
            xticks, minor_xticks = self.compute_ticks(min(self.xdata), max(self.xdata), 5, 4)
            yticks, minor_yticks = self.compute_ticks(min(self.ydata), max(self.ydata), 5, 4)

        tick_len = self.theme["tick_length"]
        minor_tick_len = self.theme.get("minor_tick_length", tick_len // 2)

        # --- X-axis ticks ---
        for tick in xticks:
            x_pixel = self._map_x(tick)
            self.canvas.draw_line(x_pixel, x_axis_y, x_pixel, x_axis_y + tick_len,
                                  thickness=self.theme["line_thickness"]["tick"],
                                  color=self.theme["tick_color"])
            
            # Only draw labels in full mode
            if self.mode == "full":
                label = f"{int(round(tick))}"
                # Center the x-labels under their ticks
                label_width = len(label) * 7  # Approximate width in pixels
                self.canvas.draw_text(x_pixel - (label_width // 2), x_axis_y + tick_len + 12,
                                  text=label,
                                  color=self.theme["font_color"],
                                  font_size=self.theme["font_size"])
                                  
        for tick in minor_xticks:
            x_pixel = self._map_x(tick)
            self.canvas.draw_line(x_pixel, x_axis_y, x_pixel, x_axis_y + minor_tick_len,
                                  thickness=self.theme["line_thickness"]["tick"],
                                  color=self.theme["tick_color"])

        # --- Y-axis ticks ---
        # Use calculated offset for label placement
        for tick in yticks:
            y_pixel = self._map_y(tick)
            self.canvas.draw_line(reg_x + L - tick_len, y_pixel, reg_x + L, y_pixel,
                                  thickness=self.theme["line_thickness"]["tick"],
                                  color=self.theme["tick_color"])
                                  
            # Only draw labels in full mode
            if self.mode == "full":
                label = f"{tick:.2f}"
                # Right-align the y-labels
                label_x = reg_x + L - self.y_label_left_offset
                self.canvas.draw_text(label_x, y_pixel + 5,  # +5 for vertical centering
                                  text=label,
                                  color=self.theme["font_color"],
                                  font_size=self.theme["font_size"])
                                  
        for tick in minor_yticks:
            y_pixel = self._map_y(tick)
            self.canvas.draw_line(reg_x + L - minor_tick_len, y_pixel, reg_x + L, y_pixel,
                                  thickness=self.theme["line_thickness"]["tick"],
                                  color=self.theme["tick_color"])

        # Axis titles
        if self.mode == "full":
            # X-axis title
            x_title_x = reg_x + L + (reg_w - L - R) // 2
            x_title_y = reg_y + reg_h - B + self.theme["font_size"] + 20
            # Center the x-axis title
            x_title_width = len(self.x_axis_title) * 7  # Approximate width in pixels
            self.canvas.draw_text(x_title_x - (x_title_width // 2), x_title_y,
                                  text=self.x_axis_title,
                                  color=self.theme["font_color"],
                                  font_size=self.theme["font_size"])

            # Y-axis title - placed vertically along the axis
            # Since we might not have rotation support in the canvas, we'll place it at the top of the y-axis
            y_title_x = reg_x + L - 25  # Small offset from the axis
            y_title_y = reg_y + T - 10  # Positioned above the top of the axis
            
            # Draw the title above the y-axis
            self.canvas.draw_text(y_title_x, y_title_y,
                                  text=self.y_axis_title,
                                  color=self.theme["font_color"],
                                  font_size=self.theme["font_size"])

    def plot_data(self):
        """Draw the data points and connecting lines."""
        points = [(self._map_x(x), self._map_y(y)) for x, y in zip(self.xdata, self.ydata)]
        for i in range(1, len(points)):
            x1, y1 = points[i - 1]
            x2, y2 = points[i]
            self.canvas.draw_line(int(x1), int(y1), int(x2), int(y2),
                                  thickness=self.theme["line_thickness"]["data"],
                                  color=self.theme["data_line_color"])

    def draw_indicator(self):
        """Draw an indicator at the current frame position."""
        try:
            x_val = self.xdata[self.current_frame]
            y_val = self.ydata[self.current_frame]
        except IndexError:
            return
        x_pixel = self._map_x(x_val)
        y_pixel = self._map_y(y_val)
        
        # Scale indicator size based on plot size and mode
        if self.mode == "minimal":
            # For minimal mode, use a smaller indicator proportional to the plot size
            plot_dimension = min(self.region_width, self.region_height)
            # Calculate size as percentage of plot dimension, with a minimum size
            indicator_radius = max(3, int(plot_dimension * 0.04))
        else:
            # For full mode, use the theme's indicator size
            indicator_radius = self.theme["indicator_radius"]
            
        self.canvas.draw_filled_circle(x_pixel, y_pixel, indicator_radius,
                                       color=self.theme["indicator_color"])

    def draw_graph(self):
        """Draw the complete graph with axes, data, and indicator."""
        # Don't clear the canvas or flush here - this is now handled by the main redraw function
        # to allow multiple graphs on the same canvas
        self.draw_axes()
        self.plot_data()
        self.draw_indicator()

    def update_current_frame(self, new_frame):
        """Update the current frame and redraw."""
        self.current_frame = new_frame
        self.draw_graph()


# ---------------------------
# Standalone Test Harness
# ---------------------------
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python graph_viewer.py <xyz_file> [current_frame] [ss_factor]")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    current_frame = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    ss_factor = int(sys.argv[3]) if len(sys.argv) > 3 else 1

    traj = Trajectory.from_xyz_file(xyz_file)
    n_frames = len(traj._raw_frames)

    xdata = list(range(n_frames))
    ydata = []
    for i in range(n_frames):
        mol = traj.get_frame(i)
        ydata.append(mol.Energy if mol.Energy is not None else 0.0)

    win = create_x11_window(width=800, height=600, title="Graph Viewer", ss_factor=ss_factor)

    # Main full-sized graph in the center
    main_width = 500
    main_height = 300
    main_x = (win.canvas.width - main_width) // 2
    main_y = (win.canvas.height - main_height) // 2

    # Create the main graph (full mode)
    main_graph = GraphViewer(
        win.canvas,
        xdata,
        ydata,
        mode="full",
        current_frame=current_frame,
        region_x=main_x,
        region_y=main_y,
        region_width=main_width,
        region_height=main_height,
        x_axis_title="Frame",
        y_axis_title="Energy (a.u.)"
    )

    # Create the thumbnail in the lower right corner (minimal mode)
    thumb_size = 100
    thumb_margin = 20
    thumb_x = win.canvas.width - thumb_size - thumb_margin
    thumb_y = win.canvas.height - thumb_size - thumb_margin
    
    # For minimal mode, use smaller margins
    minimal_theme = MINIMAL_THEME.copy()
    minimal_theme["margin"] = {"left": 10, "right": 10, "top": 10, "bottom": 10}
    
    thumbnail_graph = GraphViewer(
        win.canvas,
        xdata,
        ydata,
        mode="minimal",
        current_frame=current_frame,
        region_x=thumb_x,
        region_y=thumb_y,
        region_width=thumb_size,
        region_height=thumb_size,
        x_axis_title="",
        y_axis_title=""
    )
    # Override the theme with smaller margins
    thumbnail_graph.theme = minimal_theme
    thumbnail_graph.left_margin = minimal_theme["margin"]["left"]
    thumbnail_graph.right_margin = minimal_theme["margin"]["right"]
    thumbnail_graph.top_margin = minimal_theme["margin"]["top"]
    thumbnail_graph.bottom_margin = minimal_theme["margin"]["bottom"]

    def custom_redraw():
        # Clear the entire canvas first
        win.canvas.clear()
        
        # Draw both graphs
        main_graph.draw_graph()
        thumbnail_graph.draw_graph()
        
        # Ensure all drawing is displayed
        win.canvas.flush()

    # Update handler to update both graphs when frame changes
    def update_frame(new_frame):
        main_graph.update_current_frame(new_frame)
        thumbnail_graph.update_current_frame(new_frame)
    
    # If your window class has a way to register frame change handlers, add it here
    # win.set_frame_change_handler(update_frame)

    win.redraw = custom_redraw
    win.run()