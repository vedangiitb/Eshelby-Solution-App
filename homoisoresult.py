import tkinter as tk

class HomoIsoResult(tk.Frame):
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
        label = tk.Label(self, text="Isotropic Homogenous Result")
        label.pack(pady=10,padx=10)

        from homIsoPage import HomoIsoPage
        
        button1 = tk.Button(self, text="Back",
                            command=lambda: controller.show_frame(HomoIsoPage),)
        button1.pack()

        