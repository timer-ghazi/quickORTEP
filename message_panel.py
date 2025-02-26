# message_panel.py

class MessagePanel:
    """
    MessagePanel manages the display of recent messages at the bottom of the ORTEP viewer.
    
    Responsibilities:
      - Render the most recent messages from MessageService
      - Position messages at the bottom of the canvas
      - Use appropriate styling for different message types
    
    Attributes:
      - message_service (MessageService): The service providing message data
      - x (int): X-coordinate for the message text starting point
      - y_offset (int): Vertical offset from the bottom of the canvas
      - line_spacing (int): Vertical spacing between consecutive messages
      - font_size (int): Font size for the message text
      - bg_color (tuple): Background color (R, G, B, A) for the message area
      - padding (int): Padding around the message text
    """
    def __init__(self, message_service, x=10, y_offset=10, line_spacing=16, 
                 font_size=12, bg_color=(240, 240, 240), padding=5):
        self.message_service = message_service
        self.x = x
        self.y_offset = y_offset  # Increased from default 10 to 30
        self.line_spacing = line_spacing
        self.font_size = font_size
        self.bg_color = bg_color  # Note: X11 doesn't support alpha transparency
        self.padding = padding
    
    def draw(self, canvas):
        """
        Draw the message panel on the provided canvas. The panel is rendered
        as a semi-transparent box at the bottom of the screen with the most
        recent messages displayed inside.
        
        Parameters:
          - canvas: A canvas object that must support 'draw_text' and 'draw_rectangle'
                    methods with appropriate parameters.
        """
        # Get formatted messages and their colors from the message service
        messages, colors = self.message_service.get_formatted_messages()
        
        if not messages:
            return  # No messages to display
        
        # Calculate panel dimensions
        panel_height = len(messages) * self.line_spacing + 2 * self.padding
        panel_width = canvas.width - 2 * self.x
        
        # Draw panel background
        canvas.draw_filled_rect(
            self.x, 
            canvas.height - self.y_offset - panel_height,
            panel_width, 
            panel_height,
            color=self.bg_color
        )
        
        # Draw each message with its corresponding color
        # Reverse the messages to have newest at the bottom
        for i, (message, color) in enumerate(zip(reversed(messages), reversed(colors))):
            # Calculate the y-coordinate for each line, starting from the bottom
            y = canvas.height - self.y_offset - self.padding - i * self.line_spacing
            
            # Draw the message text
            canvas.draw_text(
                self.x + self.padding, 
                y, 
                message, 
                color=color, 
                font_size=self.font_size
            )
