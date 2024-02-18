import tkinter as tk
from welcomePage import WelcomePage
from homIsoPage import HomoIsoPage
from homAnisoPage import HomoAnisoPage
from inhomIsoPage import InhomoIsoPage
from inhomoAniso import InhomoAnisoPage

class SampleApp(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        
        self.frames = {}
        for F in (WelcomePage, HomoIsoPage,HomoAnisoPage,InhomoIsoPage,InhomoAnisoPage):
            frame = F(container, self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")
        
        self.show_frame(WelcomePage)
    
    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()

if __name__ == "__main__":
    app = SampleApp()
    app.mainloop()