# message_service.py

import time

class MessageService:
    """
    MessageService manages a queue of timestamped messages for the ORTEP viewer.
    
    Responsibilities:
      - Maintain a rolling queue of the most recent messages
      - Timestamp messages for context
      - Duplicate messages to stdout
      - Support different message types (info, warning, error)
    
    Attributes:
      - max_messages (int): Maximum number of messages to keep in the queue
      - messages (list): Queue of (timestamp, type, message) tuples
      - message_types (dict): Configuration for different message types
    """
    def __init__(self, max_messages=3):
        self.max_messages = max_messages
        self.messages = []  # List of (timestamp, type, message)
        
        # Define message types with their corresponding prefixes and colors
        self.message_types = {
            "info": {
                "prefix": "INFO",
                "color": (0, 0, 0),  # Black
            },
            "warning": {
                "prefix": "WARNING",
                "color": (200, 150, 0),  # Amber
            },
            "error": {
                "prefix": "ERROR",
                "color": (200, 0, 0),  # Red
            }
        }
    
    def _add_message(self, message_type, message):
        """
        Add a message to the queue with timestamp and type.
        
        Parameters:
          - message_type (str): Type of message ('info', 'warning', 'error')
          - message (str): The message text
        """
        timestamp = time.strftime("%H:%M:%S")
        
        # Add the new message to the end of the list
        self.messages.append((timestamp, message_type, message))
        
        # Trim the list to keep only the most recent messages
        if len(self.messages) > self.max_messages:
            self.messages = self.messages[-self.max_messages:]
    
    def log_info(self, message):
        """
        Log an informational message.
        
        Parameters:
          - message (str): The information message to log
        """
        self._add_message("info", message)
        print(f"INFO: {message}")
    
    def log_warning(self, message):
        """
        Log a warning message.
        
        Parameters:
          - message (str): The warning message to log
        """
        self._add_message("warning", message)
        print(f"WARNING: {message}")
    
    def log_error(self, message):
        """
        Log an error message.
        
        Parameters:
          - message (str): The error message to log
        """
        self._add_message("error", message)
        print(f"ERROR: {message}")
    
    def get_messages(self):
        """
        Get all messages in the queue.
        
        Returns:
          - List of (timestamp, type, message) tuples
        """
        return self.messages
    
    def get_formatted_messages(self):
        """
        Get formatted strings for all messages in the queue.
        
        Returns:
          - List of formatted message strings
          - List of corresponding colors
        """
        formatted = []
        colors = []
        
        for timestamp, msg_type, message in self.messages:
            formatted.append(f"[{timestamp}] {self.message_types[msg_type]['prefix']}: {message}")
            colors.append(self.message_types[msg_type]['color'])
        
        return formatted, colors
    
    def clear(self):
        """Clear all messages from the queue."""
        self.messages = []
