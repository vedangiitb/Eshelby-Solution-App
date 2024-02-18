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

        # Entry widget for user input
        self.input_entry = tk.Entry(self)
        self.input_entry.pack(pady=10)
        # Button to process user input
        process_button = tk.Button(self, text="Process", command=self.process_input)
        process_button.pack()

    def process_input(self):
        # Retrieve the user input from the Entry widget
        user_input = self.input_entry.get()
        # Do something with the user input, for example, display it
        output_label = tk.Label(self, text="You entered: " + user_input)
        output_label.pack()

