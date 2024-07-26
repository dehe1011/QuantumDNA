import sys
import customtkinter as ctk
from tkinter.scrolledtext import ScrolledText


class RedirectText:
    def __init__(self, text_widget):
        self.output = text_widget

    def write(self, string):
        self.output.configure(state="normal")
        self.output.insert(ctk.END, string)
        self.output.see(ctk.END)
        self.output.configure(state="disabled")

    def flush(self):
        pass


class ScrollableConsoleFrame(ctk.CTkFrame):
    def __init__(self, master=None, **kwargs):
        super().__init__(master, **kwargs)

        # Create a ScrolledText widget for the console output
        self.console_output = ScrolledText(
            self, wrap="word", state="disabled", height=10, width=40
        )
        self.console_output.grid(row=0, column=0, sticky="nsew", padx=(10, 0), pady=10)

        scrollbar = ctk.CTkScrollbar(self, command=self.console_output.yview)
        scrollbar.grid(row=0, column=1, sticky="ns", padx=(0, 10), pady=10)
        self.console_output["yscrollcommand"] = scrollbar.set

        # Redirect stdout to the ScrolledText widget
        self.redirected_output = RedirectText(self.console_output)
        sys.stdout = self.redirected_output
