import customtkinter as ctk


def change_state_all_widgets(frame, state):
    for widget in frame.winfo_children():
        if isinstance(
            widget, ctk.CTkBaseClass
        ):  # Check if the widget is a customtkinter widget
            widget.configure(state=state)
        elif isinstance(
            widget, ctk.CTkFrame
        ):  # Recursively disable widgets in nested frames
            change_state_all_widgets(widget, state)
