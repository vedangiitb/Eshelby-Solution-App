# This is a common function that we have for input params
# Other functions will inherit this and add/remove appropriate params
class InputPagesParent:
    def __init__(self):
        # Initialize the params
        self.a1, self.a2, self.a3, self.E, self.nu, self.mu, self.ep11, self.ep22, self.ep33, self.ep12, self.ep13, self.ep23 = [None] * 12
        

    # The function currently supports only for Isotropic homogenous inside.
    # We will add other params once we are done with rest of the code
    def initInputParams(self,a1,a2,a3,ep11,ep22,ep33,ep13,ep12,ep23,E=None,nu=None,mu=None):
        # Initialize a1,a2,a3
        self.a1,self.a2,self.a3 = a1,a2,a3

        #  TODO: Add a new error function that checks if the input params are correct.
        # If the function returns false, it throws an error

        if E!=None and nu!=None:
            self.E = E
            self.nu = nu
            self.mu = E/(2*(1+nu))
        elif E!=None and mu!=None:
            self.E = E
            self.mu = self.mu
            self.nu = (E-2*mu)/(2*mu)
        elif nu!=None and mu!=None:
            self.nu = nu
            self.mu = mu
            self.E = 2*mu*(1+nu)
        else:
            pass
            # TODO: Throw an error
        
        # TODO: Check if all epij are not none
        # If the function returns false, it throws an error

        self.ep11 = ep11
        self.ep22 = ep22
        self.ep33 = ep33
        self.ep12 = ep12
        self.ep13 = ep13
        self.ep23 = ep23

        # We return True once we have confirmed that all values are OK and we have initialized all the functions 
        return True
    





