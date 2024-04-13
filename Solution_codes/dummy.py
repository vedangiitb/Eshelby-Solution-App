
import sys
import json
import numpy as np
from mayavi import mlab


def dummy_func(a,b,c,eps11,eps22,eps33,eps12,eps23,eps31,E,nu):
    a += 1
    output_data = f"a: {a}, b: {b}, c: {c}, eps11: {eps11}, eps22: {eps22}, eps33: {eps33}, eps12: {eps12}, eps23: {eps23}, eps31: {eps31}, E: {E}, nu: {nu}"
    
    return output_data

def plot():

    
    X,Y,Z = np.mgrid[-1:1:50j,-1:1:50j,-1:1:50j]
    sigma_xx = np.random.randn(50,50,50)
    mlab.points3d(X,Y,Z,sigma_xx)
    mlab.axes(line_width=1.0)
    mlab.show()
    return 

# def volume_slice():
#     mlab.volume_slice(X,Y,Z,sigma_xx)
#     mlab.axes(line_width=1.0,xlabel='x',ylabel='y',zlabel='z')
if __name__ =="__main__":

    if sys.argv[1] != 'plot':
        try:
            # Load the form data passed from Express.js
            form_data = json.loads(sys.argv[1])

            # Taking Inputs
            a = float(form_data.get('a',''))
            b = float(form_data.get('b',''))
            c = float(form_data.get('c',''))
            eps11 = float(form_data.get('eps11'))
            eps22 = float(form_data.get('eps22'))
            eps33 = float(form_data.get('eps33'))
            eps12 = float(form_data.get('eps12'))
            eps23 = float(form_data.get('eps23'))
            eps31 = float(form_data.get('eps13'))
            E = float(form_data.get('ep'))
            nu = float(form_data.get('nu'))


            axis = [a, b, c]

            # print(f"input data:{a,b,c,eps11,eps12,eps31,eps23,eps22,eps33,E,nu}")

            output_data = dummy_func(a,b,c,eps11,eps22,eps33,eps12,eps23,eps31,E,nu)

            print(output_data)

        except Exception as e:
            print("error:", e)

    else:
        plot()
