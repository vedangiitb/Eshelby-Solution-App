import tkinter as tk
from welcomePage import WelcomePage

class HomoIsoPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Page Two")
        label.pack(pady=10,padx=10)
        
        button1 = tk.Button(self, text="Go to Start Page",
                            command=lambda: controller.show_frame(WelcomePage))
        button1.pack()