try:
    import numpy as np
    import scipy.special as sp
    import sympy as sym
    import scipy.integrate as integrate
    pi = np.pi
    import sys
    import json
    import pandas as pd
    from mayavi import mlab
except Exception as e:
    print(e)


def saveData(points,data):
    x = []
    y = []
    z = []

    for p in points:
        x.append(p[0])
        y.append(p[1])
        z.append(p[2])
    
    sigma11 = []
    sigma12 = []
    sigma13 = []
    sigma22 = []
    sigma23 = []
    sigma33 = []

    for d in data:
        sigma11.append(d[0])
        sigma22.append(d[1])
        sigma33.append(d[2])
        sigma12.append(d[3])
        sigma13.append(d[4])
        sigma23.append(d[5])

    df = pd.DataFrame(list(zip(x,y,z,sigma11,sigma22,sigma33,sigma12,sigma13,sigma23)),columns=["x","y","z",'sigma11','sigma22','sigma33','sigma12','sigma13','sigma23'])

    df.to_csv('public/temp.csv',index=False)


# Function to find Lambda (the largest root)
def find_lambda(X, axis):
    ans = sym.symbols('ans')
    eqn = sym.Eq(((X[0])**2)/((axis[0])**2 + ans) + ((X[1])**2)/((axis[1])**2 + ans) + ((X[2])**2)/((axis[2])**2 + ans), 1)
    solution = sym.solve(eqn, ans)
    sol = [complex(i) for i in solution]
    #print(sol)
    LAMBDA = 0
    for i in sol:
        LAMBDA = max(LAMBDA, i.real)
        
    return LAMBDA



def IIJ(axis, L, flag):

    theta = np.arcsin(np.sqrt((axis[0] ** 2 - axis[2] ** 2) / (axis[0] ** 2 + L)))
    k = np.sqrt((axis[0] ** 2 - axis[1] ** 2) / (axis[0] ** 2 - axis[2] ** 2))
    F = sp.ellipkinc(theta, k ** 2)
    E = sp.ellipeinc(theta, k ** 2)
    #print("theta", theta, k, F, E)
    arr = np.zeros(3)
    dels = np.sqrt((axis[0]**2+L)*(axis[1]**2+L)*(axis[2]**2+L))
    c1 = (4*pi*axis[0]*axis[1]*axis[2])/((axis[0]**2 - axis[1]**2)*(np.sqrt(axis[0]**2 - axis[2]**2)))
    arr[0] = c1*(F-E)
    c2 = (4*pi*axis[0]*axis[1]*axis[2])/((axis[1]**2 - axis[2]**2)*(np.sqrt(axis[0]**2 - axis[2]**2)))
    d1 = ((axis[1]**2 + L)*(np.sqrt(axis[0]**2 - axis[2]**2)))/dels
    arr[2] = c2*(d1-E)

    arr[1] = (4*pi*axis[0]*axis[1]*axis[2])/dels - arr[0] - arr[2]

    arr1 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i!=j:
                arr1[i][j] = (arr[j]-arr[i])/(axis[i]**2 - axis[j]**2)
    for i in range(3):
        tmp = 0
        for j in range(3):
            tmp += arr1[i][j]/3
        arr1[i][i] = (4*pi*axis[0]*axis[1]*axis[2])/(3*(axis[i]**2+L)*dels) - tmp
    
    if flag == 0:
        return arr
    else:
        return arr1

def get_Sijkl(axis, Ii, Iij):
    Sijkl = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    t1 = np.identity(3)[i, j] * np.identity(3)[k, l] * (2 * nu * Ii[i] - Ii[k] + (axis[i]**2) * Iij[k][i])
                    t2 = (np.identity(3)[i, k] * np.identity(3)[j, l] + np.identity(3)[i, l] * np.identity(3)[j, k]) * ((axis[i]**2) * Iij[i][j] - Ii[j] + (1 - nu) * (Ii[k] + Ii[l]))
                    Sijkl[i, j, k, l] = (t1 + t2) / (8 * np.pi * (1 - nu))
    return Sijkl


def lamb_der(x, a, lambda_):
    arr = []
    #Check if x in inside
    if (x[0]**2/(a[0]**2)) + (x[1]**2/(a[1]**2)) + (x[2]**2/(a[2]**2)) <= 1:
        return [0,0,0]
    denom = (x[0]**2/((a[0]**2 + lambda_)**2)) + (x[1]**2/((a[1]**2 + lambda_)**2)) + (x[2]**2/((a[2]**2 + lambda_)**2))
    for i in range(3):
        num = (2*x[i])/(a[i]**2 + lambda_)
        arr.append(num/denom)
    return arr

# Compute double derivative matrix of lambda
def lamb_der2(x, a, lambda_, lambda_der):
    arr = np.zeros((3, 3))
    if (x[0]**2/(a[0]**2)) + (x[1]**2/(a[1]**2)) + (x[2]**2/(a[2]**2)) <= 1:
        return arr
    for i in range(3):
        for j in range(3):
            arr[i, j] = ((2*a[i]*lambda_der[i]*lambda_der[j]) - 2*lambda_der[j]*(a[i]**2 + lambda_))/((a[i]**2+lambda_)*a[i])
    return arr

#Compute derivative of Ii wrt to j direction
def Ii_j_(a,lambda_,lambda_der):
    arr = np.zeros((3,3))
    c = -2*np.pi*a[0]*a[1]*a[2]
    del_l = ((a[0]**2 + lambda_)*(a[1]**2 + lambda_)*(a[2]**2 + lambda_))**(1/2)
    for i in range(3):
        for j in range(3):
            arr[i,j] = c*lambda_der[j]/((a[i]**2 + lambda_)*del_l)
    return arr


# Compute derivative of Iij wrt to k direction
def Iij_k_(a, lambda_, lambda_der):
    arr = np.zeros((3,3,3))
    c = -2*np.pi*a[0]*a[1]*a[2]
    del_l = ((a[0]**2 + lambda_)*(a[1]**2 + lambda_)*(a[2]**2 + lambda_))**(1/2)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                arr[i,j,k] = c*lambda_der[k]/((a[i]**2 + lambda_)*(a[j]**2 + lambda_)*del_l)
    return arr


# Compute double partial derivative of Iij wrt to k and l direction
def Iij_kl_(a, lambda_, lambda_der, lambda_der2):
    arr = np.zeros((3,3,3,3))
    c = -2*np.pi*a[0]*a[1]*a[2]
    del_l = ((a[0]**2 + lambda_)*(a[1]**2 + lambda_)*(a[2]**2 + lambda_))**(1/2)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    arr[i,j,k,l] = (c/((a[i]**2 + lambda_)*(a[j]**2 + lambda_)*(del_l)))*(lambda_der2[k,l] - lambda_der[k]*lambda_der[l]*(1/(a[i]**2 + lambda_) + 1/(a[j]**2 + lambda_) + 0.5*(1/(a[0]**2 + lambda_) + 1/(a[1]**2 + lambda_) + 1/(a[2]**2 + lambda_))))
    return arr


# Compute derivative of Ii wrt to j and k direction
def Ii_jk_(a, lambda_, lambda_der, lambda_der2):
    arr = np.zeros((3,3,3))
    c = -2*np.pi*a[0]*a[1]*a[2]
    del_l = ((a[0]**2 + lambda_)*(a[1]**2 + lambda_)*(a[2]**2 + lambda_))**(1/2)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                arr[i,j,k] = (c/((a[i]**2 + lambda_)*(del_l)))*(lambda_der2[j,k] - lambda_der[j]*lambda_der[k]*(1/(a[i]**2 + lambda_) + 0.5*(1/(a[0]**2 + lambda_) + 1/(a[1]**2 + lambda_) + 1/(a[2]**2 + lambda_))))
    return arr


def get_Dijkl(s, delta, x, a, IIj, IIJk, IIJkl, IIjk):
    Dijkl = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    t1 = 8*pi*(1-nu)*s[i,j,k,l] + 2*nu*delta[k,l]*x[i]*IIj[i,j]
                    t2 = (1-nu)*(delta[i,l]*x[k]*IIj[k,j] + delta[j,l]*x[k]*IIj[k,i] + delta[i,k]*x[l]*IIj[l,j] + delta[j,k]*x[l]*IIj[l,i])
                    t3 = (-1)*delta[i,j]*x[k]*(IIj[k,l] - (a[i]**2)*IIJk[k,i,l]) - (delta[i,k]*x[j] + delta[j,k]*x[i])*(IIj[j,l] - (a[i]**2)*IIJk[i,j,l])
                    t4 = (-1)*(delta[i,l]*x[j] + delta[j,l]*x[i])*(IIj[j,k]-(a[i]**2)*IIJk[i,j,k]) - x[i]*x[j]*(IIjk[j,l,k]-(a[i]**2)*IIJkl[i,j,l,k])
                    Dijkl[i, j, k, l] = (t1 + t2 + t3 + t4) / (8 * pi * (1 - nu))
    return Dijkl

def calc_interior():
    #Calculating theta and K

    theta = np.arcsin(np.sqrt((a**2 - c**2)/a**2))
    k = np.sqrt((a**2 - b**2)/(a**2 - c**2))

    F = sp.ellipkinc(theta, k**2)
    E = sp.ellipeinc(theta, k**2)
    I1 = (4*pi*a*b*c)*(F - E) / ((a**2 - b**2)*np.sqrt(a**2 - c**2))
    I3 = ((4*pi*a*b*c)*((b*np.sqrt(a**2-c**2))/(a*c) - E))/((b**2-c**2)*np.sqrt(a**2-c**2))
    I2 = 4*pi - I1 - I3


    #Solvig equations for I11 and I13
    I12 = (I2-I1)/(a**2-b**2)
    x, y = sym.symbols('x y')
    eq1 = sym.Eq(3*x+I12+y, 4*pi/(a**2))
    eq2 = sym.Eq(3*(a**2)*x+(b**2)*I12+(c**2)*y, 3*I1)
    ans = sym.solve([eq1, eq2], (x, y))
    I11 = ans[x]
    I13 = ans[y]

    #Solving equations for I21 and I22
    I23 = (I3-I2)/(b**2-c**2)
    x, y = sym.symbols('x y')
    eq1 = sym.Eq(3*x+I23+y, 4*pi/(b**2))
    eq2 = sym.Eq(3*(b**2)*x+(c**2)*I23+(a**2)*y, 3*I2)
    ans = sym.solve([eq1, eq2], (x, y))
    I21 = ans[y]
    I22 = ans[x]

    #Solving equation for I33 and I32
    I31 = (I1-I3)/(c**2-a**2)
    x, y = sym.symbols('x y')
    eq1 = sym.Eq(3*x+I31+y, 4*pi/(c**2))
    eq2 = sym.Eq(3*(c**2)*x+(a**2)*I31+(b**2)*y, 3*I3)
    ans = sym.solve([eq1, eq2], (x, y))
    I32 = ans[y]
    I33 = ans[x]
    #Solving for sigma

    #t1 t2 t3 are terms of sigma11
    t1 = (((a**2)*(3*I11 - 3*nu*I11 + nu*I21 + nu*I31))/(8*pi*(1-nu)*(1-2*nu)) + ((I1 - (nu/(1-nu))*(I2+I3))/(8*pi)) - ((1-nu)/(1-2*nu)))*eps11

    t2 = (((b**2)*(I12-I12*nu + 3*nu*I22 + nu*I32))/(8*pi*(1-nu)*(1-2*nu)) - (I1 - (nu/1-nu)*(I2-I3))/(8*pi) - (nu)/(1-2*nu))*eps22

    t3 = (((c**2)*(I13 - nu*I13 + 3*nu*I33 + nu*I23))/(8*pi*(1-nu)*(1-2*nu)) - (I1 - (nu/1-nu)*(I3-I2))/(8*pi) - (nu)/(1-2*nu))*eps33

    sigma11 = 2*mu*(t1+t2+t3)

    sigma12 = ((((a**2+b**2)*I12 + (1-2*nu)*(I1+I2))/(8*pi*(1-nu)))-1)*eps12*2*mu
    sigma23 = ((((b**2+c**2)*I12 + (1-2*nu)*(I2+I3))/(8*pi*(1-nu)))-1)*eps23*2*mu
    sigma31 = ((((a**2+c**2)*I12 + (1-2*nu)*(I1+I3))/(8*pi*(1-nu)))-1)*eps31*2*mu


    #q1 q2 q3 are for sigma22


    q1 = (((b**2)*(3*I22 - 3*nu*I22 + nu*I32 + nu*I12))/(8*pi*(1-nu)*(1-2*nu)) + ((I2 - (nu/(1-nu))*(I1+I3))/(8*pi)) - ((1-nu)/(1-2*nu)))*eps22
    q2 = (((c**2)*(I23-I23*nu + 3*nu*I33 + nu*I13))/(8*pi*(1-nu)*(1-2*nu)) - (I2 - (nu/1-nu)*(I3-I1))/(8*pi) - (nu)/(1-2*nu))*eps33
    q3 = (((a**2)*(I21 - nu*I21 + 3*nu*I11 + nu*I31))/(8*pi*(1-nu)*(1-2*nu)) - (I2 - (nu/1-nu)*(I1-I3))/(8*pi) - (nu)/(1-2*nu))*eps11

    sigma22 = 2*mu*(q1+q2+q3)


    #w1 w2 w3 are for sigma33


    w1 = (((c**2)*(3*I33 - 3*nu*I33 + nu*I13 + nu*I23))/(8*pi*(1-nu)*(1-2*nu)) + ((I3 - (nu/(1-nu))*(I2+I1))/(8*pi)) - ((1-nu)/(1-2*nu)))*eps33
    w2 = (((a**2)*(I31-I31*nu + 3*nu*I11 + nu*I21))/(8*pi*(1-nu)*(1-2*nu)) - (I3 - (nu/1-nu)*(I1-I2))/(8*pi) - (nu)/(1-2*nu))*eps11
    w3 = (((b**2)*(I32 - nu*I32 + 3*nu*I22 + nu*I12))/(8*pi*(1-nu)*(1-2*nu)) - (I3 - (nu/1-nu)*(I2-I1))/(8*pi) - (nu)/(1-2*nu))*eps22

    sigma33 = 2*mu*(w1+w2+w3)

    #stress_inside = np.array([[sigma11, sigma12, sigma31], [sigma12, sigma22, sigma23], [sigma31, sigma23, sigma33]])
    #print(stress_inside)

    
    # print(f"sigma11 = {sigma11}\nsigma22 = {sigma22}\nsigma33 = {sigma33}\nsigma12 = {sigma12}\nsigma13 = {sigma31}\nsigma23 = {sigma23}")
    # print()
    return [sigma11,sigma22,sigma33,sigma12,sigma31,sigma23]
    #return stress_inside[0][1]

def calc_exterior(X):
    epsilon_star = [[eps11, eps12, 0], [0, eps22, eps23], [eps31, 0, eps33]]
    
    
    lbd = find_lambda(X, axis)
    # print("Lambda is: ", lbd)

    #Calculating I, I1, I2, I3, I11, I22, I33, etc
        
    Ii = IIJ(axis, lbd, 0)
    Iij = IIJ(axis, lbd, 1)

    Sijkl = get_Sijkl(axis, Ii, Iij)
    
    lbd_der = lamb_der(X, axis, lbd)
    lbd_d_der = lamb_der2(X, axis, lbd, lbd_der)
    
    

    IIj = Ii_j_(axis, lbd, lbd_der)
    IIJk = Iij_k_(axis, lbd, lbd_der)
    IIJkl = Iij_kl_(axis, lbd, lbd_der, lbd_d_der)
    IIjk = Ii_jk_(axis, lbd, lbd_der, lbd_d_der)
    

    deltaij = np.identity(3)
    Dijkl = get_Dijkl(Sijkl, deltaij, X, axis, IIj, IIJk, IIJkl, IIjk)
    
    epsilon = np.zeros((3, 3))

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    epsilon[i][j] += Dijkl[i,j,k,l] * epsilon_star[k][l]


    
    sig11 = (E/((1+nu)*(1-2*nu)))*((1-nu)*epsilon[0,0] + nu*epsilon[1,1] + nu*epsilon[2,2])
    sig22 = (E/((1+nu)*(1-2*nu)))*(nu*epsilon[0,0] + (1-nu)*epsilon[1,1] + nu*epsilon[2,2])
    sig33 = (E/((1+nu)*(1-2*nu)))*(nu*epsilon[0,0] + nu*epsilon[1,1] + (1-nu)*epsilon[2,2])
    sig12 = (E/(1+nu))*epsilon[0,1]
    sig13 = (E/(1+nu))*epsilon[0,2]
    sig23 = (E/(1+nu))*epsilon[1,2]
    #stress_outside = np.array([[sig11,sig12,sig13],[sig12,sig22,sig23],[sig13,sig23,sig33]])
    #print(stress_outside)
    # print(f"For the point {X} [The point is present outside], the stress values are")
    # print(f"sigma11 = {sig11}\nsigma22 = {sig22}\nsigma33 = {sig33}\nsigma12 = {sig12}\nsigma13 = {sig13}\nsigma23 = {sig23}")
    # print()
    return [sig11,sig22,sig33,sig12,sig13,sig23]
    #return stress_outside[0][1]



# Load the form data passed from Express.js
form_data = json.loads(sys.argv[1])

# Taking Inputs
a = float(form_data.get('a'))
b = float(form_data.get('b'))
c = float(form_data.get('c'))
axis = [a,b,c]
eps11 = float(form_data.get('eps11'))
eps22 = float(form_data.get('eps22'))
eps33 = float(form_data.get('eps33'))
eps12 = float(form_data.get('eps12'))
eps23 = float(form_data.get('eps23'))
eps31 = float(form_data.get('eps13'))
E = float(form_data.get('ep'))
nu = float(form_data.get('nu'))
targets = form_data.get('targets')
plottype = form_data.get('plottype')
mu = E/(2*(1+nu))
print(form_data)

chosenDir = 0

# print("**********************************")
# print(form_data)
# print("**********************************")
if plottype == "22":
    chosenDir = 1
elif plottype=="33":
    chosenDir = 2
elif plottype=="12":
    chosenDir = 3
elif plottype=="13":
    chosenDir = 4
elif plottype == '23':
    chosenDir = 5
else:
    chosenDir = 0



intStress = calc_interior()
print(f"intStress:{intStress}")

def calcStress(x,y,z,direction):
    if((x**2/a**2) + (y**2/b**2) + (z**2/c**2) - 1 <= 0):
        if intStress[direction]/mu >0:
            print(f'stress:{intStress[direction]/mu} at cord {(x,y,z)} inside')
        return intStress[direction]/mu
    else:
        Arr = calc_exterior([x,y,z])
        if Arr[direction]/mu>0:
            print(f'stress:{Arr[direction]/mu} at cord {(x,y,z)}')
        return Arr[direction]/mu



# step = 4*a/13
print(plottype)
print(chosenDir)

x = np.linspace(-2*a,2*a,10)
y = np.linspace(-2*b,2*b,10)
z = np.linspace(-2*c,2*c,10)
B,A, C = np.meshgrid(x,y,z)
# print(X.shape)
print("plotting..")
stressArr = np.zeros(B.shape)

for i in range(A.shape[0]):
    for j in range(A.shape[1]):
        for k in range(A.shape[2]):
            stressArr[i][j][k] = calcStress(A[i][j][k],B[i][j][k],C[i][j][k],chosenDir)
            print(A[i][j][k],B[i][j][k],C[i][j][k])


# Plot the scalar field
src = mlab.pipeline.scalar_field(A,B,C, stressArr)

# Add contours to the scalar field
plot = mlab.pipeline.surface(src)
mlab.pipeline.image_plane_widget(src)



# plot = mlab.points3d(A,B,C,stressArr,scale_mode='none')
# adjust scalae_factor through mayavi gui for better visibility
# mlab.colorbar(plot)
plot.module_manager.scalar_lut_manager.show_legend = True
plot.module_manager.scalar_lut_manager.scalar_bar_representation.minimum_size = [2, 2]
plot.module_manager.scalar_lut_manager.scalar_bar_representation.position = [0.884375  , 0.41530076]
plot.module_manager.scalar_lut_manager.scalar_bar_representation.position2 = [0.075     , 0.51279791]


mlab.axes(line_width=1.0)
# mlab.volume_slice(A,B,C,stressArr,extent=[-2*a,2*a,-2*b,2*b,-2*c,2*c])
# mlab.set_picker_props(figure=None, pick_type='point_picker', tolerance=0.025, text_color=None)
# mlab.set_picker_props()
engine = mlab.get_engine()
scene = engine.scenes[0]
scene.scene.background = (0.0, 0.0, 0.0)

# labels = Labels()
# engine.add_filter(labels,plot.module_manager)
# labels.actor.mapper.label_mode = 'label_field_data'


mlab.show()
print("plotting finished")

# def volume_slice():
#     mlab.volume_slice(X,Y,Z,sigma_xx)
#     mlab.axes(line_width=1.0,xlabel='x',ylabel='y',zlabel='z')


# if __name__ =="__main__":

#     if sys.argv[1] != 'plot':
#         try:
#             # Load the form data passed from Express.js
#             form_data = json.loads(sys.argv[1])

#             # Taking Inputs
#             a = float(form_data.get('a',''))
#             b = float(form_data.get('b',''))
#             c = float(form_data.get('c',''))
#             eps11 = float(form_data.get('eps11'))
#             eps22 = float(form_data.get('eps22'))
#             eps33 = float(form_data.get('eps33'))
#             eps12 = float(form_data.get('eps12'))
#             eps23 = float(form_data.get('eps23'))
#             eps31 = float(form_data.get('eps13'))
#             E = float(form_data.get('ep'))
#             nu = float(form_data.get('nu'))


#             axis = [a, b, c]

#             # print(f"input data:{a,b,c,eps11,eps12,eps31,eps23,eps22,eps33,E,nu}")

#             output_data = dummy_func(a,b,c,eps11,eps22,eps33,eps12,eps23,eps31,E,nu)

#             print(output_data)

#         except Exception as e:
#             print("error:", e)


