import tkinter as tk

class HomoIsoPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text="Isotropic Homogenous")
        label.pack(pady=10,padx=10)
        
        # Importing WelcomePage here
        from welcomePage import WelcomePage
        
        button1 = tk.Button(self, text="Back",
                            command=lambda: controller.show_frame(WelcomePage))
        button1.pack()

        self.takeInputParams()

        # Button to process user input
        run_button = tk.Button(self, text="Process", command=self.process_input)
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
            self.computeAns()

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

    def checkInputParams(self):
        entries = {
            'a1': self.a1.get(),
            'a2': self.a2.get(),
            'a3': self.a3.get(),
            'E': self.E.get(),
            'nu': self.nu.get(),
            'mu': self.mu.get(),
            'ep11': self.ep11.get(),
            'ep22': self.ep22.get(),
            'ep33': self.ep33.get(),
            'ep12': self.ep12.get(),
            'ep13': self.ep13.get(),
            'ep23': self.ep23.get()
        }
        
        for field, entry in entries.items():
            if not entry:
                return (False,field)
            
        return (True,None)


    def computeAns(self):
        return True