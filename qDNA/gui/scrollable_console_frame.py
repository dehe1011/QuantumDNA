# pylint: skip-file

import sys
import customtkinter as ctk
from tkinter.scrolledtext import ScrolledText
import re

# ------------------------------------------------------------

import sys
import customtkinter as ctk
from tkinter.scrolledtext import ScrolledText


class RedirectText:
    def __init__(self, text_widget):
        self.output = text_widget

    def write(self, string):
        """Handles both stdout and stderr redirection."""
        self.output.configure(state="normal")
        self.output.insert(ctk.END, string.strip() + "\n")  # Ensures new line output
        self.output.see(ctk.END)
        self.output.configure(state="disabled")

    def flush(self):
        pass  # Required for compatibility with sys.stdout and sys.stderr


class ScrollableConsoleFrame(ctk.CTkFrame):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)

        # Create a ScrolledText widget for the console output
        self.console_output = ScrolledText(
            self, wrap="word", state="disabled", height=10, width=35
        )
        self.console_output.grid(row=0, column=0, sticky="nsew", padx=(10, 0), pady=10)

        scrollbar = ctk.CTkScrollbar(self, command=self.console_output.yview)
        scrollbar.grid(row=0, column=1, sticky="ns", padx=(0, 10), pady=10)
        self.console_output["yscrollcommand"] = scrollbar.set

        # Redirect both stdout and stderr to the ScrolledText widget
        self.redirected_output = RedirectText(self.console_output)
        sys.stdout = self.redirected_output
        sys.stderr = self.redirected_output

        # Override sys.excepthook to filter AssertionError traceback
        sys.excepthook = self.handle_exceptions

        # Override Tkinter's error handling to capture GUI-related assertions
        self.master.report_callback_exception = self.handle_tkinter_exceptions

    def handle_exceptions(self, exc_type, exc_value, exc_traceback):
        """Custom exception handler that filters out assertion tracebacks."""
        if issubclass(exc_type, AssertionError):
            sys.stderr.write(str(exc_value) + "\n")  # Only print the assertion message
        else:
            import traceback

            sys.stderr.write(
                "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))
            )  # Keep full error

    def handle_tkinter_exceptions(self, exc_type, exc_value, exc_traceback):
        """Handles Tkinter's internal exceptions to prevent unwanted traceback."""
        if issubclass(exc_type, AssertionError):
            sys.stderr.write(
                "Warning: " + str(exc_value) + "\n-------------------------------"
            )  # Print only assertion message
        else:
            import traceback

            sys.stderr.write(
                "".join(traceback.format_exception(exc_type, exc_value, exc_traceback))
            )  # Full traceback
