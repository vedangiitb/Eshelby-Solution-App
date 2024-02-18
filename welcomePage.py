import tkinter as tk
from homIsoPage import HomoIsoPage
from homAnisoPage import HomoAnisoPage
from inhomIsoPage import InhomoIsoPage
from inhomoAniso import InhomoAnisoPage

class WelcomePage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Eshelby Solution Calculator")
        label.pack(pady=10,padx=10)
        
        button1 = tk.Button(self, text="Isotropic Homogenous",
                            command=lambda: controller.show_frame(HomoIsoPage))
        
        button2 = tk.Button(self, text="Anisotropic Homogenous",
                            command=lambda: controller.show_frame(HomoAnisoPage))
        
        button3 = tk.Button(self, text="Isotropic Inhomogenous",
                            command=lambda: controller.show_frame(InhomoIsoPage))
        
        button4 = tk.Button(self, text="Anisotropic Inhomogenous",
                            command=lambda: controller.show_frame(InhomoAnisoPage))
        
        button1.pack(side=tk.LEFT, padx=5, pady=5)
        button2.pack(side=tk.LEFT, padx=5, pady=5)
        button3.pack(side=tk.LEFT, padx=5, pady=5)
        button4.pack(side=tk.LEFT, padx=5, pady=5)
