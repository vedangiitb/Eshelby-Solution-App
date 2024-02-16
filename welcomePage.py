import tkinter as tk
from homIsoPage import HomoIsoPage

class WelcomePage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Start Page")
        label.pack(pady=10,padx=10)
        
        button1 = tk.Button(self, text="Go to Page Two",
                            command=lambda: controller.show_frame(HomoIsoPage))
        button1.pack()