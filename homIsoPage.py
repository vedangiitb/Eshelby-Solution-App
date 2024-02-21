import tkinter as tk
from homoisoresult import HomoIsoResult

class HomoIsoPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Isotropic Homogenous")
        label.pack(pady=10,padx=10)
        
        # Importing WelcomePage here
        from welcomePage import WelcomePage
        
        button1 = tk.Button(self, text="Back",
                            command=lambda: controller.show_frame(WelcomePage),)
        button1.pack()

        self.takeInputParams()

        # Button to process user input
        # run_button = tk.Button(self, text="Process", command=self.computeAns)
        run_button = tk.Button(self, text="Process", command=lambda:controller.show_frame(HomoIsoResult) if self.process_input() else self.process_input())
        run_button.pack()

    def process_input(self):
        # Check if the input is valid or not
        val, label = self.checkInputParams()
        
        if not val:
            output_label = tk.Label(self, text="No value entered for " + label)
            output_label.pack()

        else:
            # Replace this with actual answer compute karne wala function 
            output_label = tk.Label(self, text="All values OK!")
            output_label.pack()
            self.computeAns(self.inputParams())
            return True

    def takeInputParams(self):
        # Entry widget for user input
        self.a1 = tk.Entry(self)
        self.a1.pack(pady=10)

        self.a2 = tk.Entry(self)
        self.a2.pack(pady=10)

        self.a3 = tk.Entry(self)
        self.a3.pack(pady=10)

        self.E = tk.Entry(self)
        self.E.pack(pady=10)

        self.nu = tk.Entry(self)
        self.nu.pack(pady=10)

        self.mu = tk.Entry(self)
        self.mu.pack(pady=10)

        self.ep11 = tk.Entry(self)
        self.ep11.pack(pady=10)

        self.ep22 = tk.Entry(self)
        self.ep22.pack(pady=10)

        self.ep33 = tk.Entry(self)
        self.ep33.pack(pady=10)

        self.ep12 = tk.Entry(self)
        self.ep12.pack(pady=10)

        self.ep13 = tk.Entry(self)
        self.ep13.pack(pady=10)

        self.ep23 = tk.Entry(self)
        self.ep23.pack(pady=10)

    def inputParams(self):
        entries = {
            'a1': self.a1.get(),
            'a2': self.a2.get(),
            'a3': self.a3.get(),
            'ep11': self.ep11.get(),
            'ep22': self.ep22.get(),
            'ep33': self.ep33.get(),
            'ep12': self.ep12.get(),
            'ep13': self.ep13.get(),
            'ep23': self.ep23.get()
        }

        return entries

    def checkInputParams(self):
        
        entries = self.inputParams()

        
        for field, entry in entries.items():
            if not entry:
                return (False,field)
        
        noneCount = 0

        if self.E.get()==None:
            noneCount+=1
        if self.nu.get()==None:
            noneCount+=1
        if self.mu.get()==None:
            noneCount+=1
        
        if noneCount>0:
            return (False, "atleast 2 modulus quantities")
            
        print(entries)
        return (True,None)


    def computeAns(self,input_para):
        a1 = input_para['a1']
        a2 = input_para['a2']
        a3 = input_para['a3']
        ep11 = input_para['ep11']
        ep22 = input_para['ep22']
        ep33 = input_para['ep33']
        ep12 = input_para['ep12']
        ep13 = input_para['ep13']
        ep23 = input_para['ep23']
        E = self.E.get()
        nu = self.nu.get()
        mu = self.mu.get()

        #Call the function from from solution codes here with above parameters
        # convert the output numpy arrays into vtk format and use it to plot in ovito
        pass
            
            
        


        
       