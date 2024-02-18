import tkinter as tk

class InhomoIsoPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Isotropic Inhomogenous")
        label.pack(pady=10,padx=10)

        from welcomePage import WelcomePage
        
        button1 = tk.Button(self, text="Back",
                            command=lambda: controller.show_frame(WelcomePage))
        button1.pack()