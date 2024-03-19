try:
    import numpy as np
    import scipy.special as sp
    import sympy as sym
    import scipy.integrate as integrate
    pi = np.pi
    import sys
    import json
except Exception as e:
    print(e)

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


def IIJ(axis, L, i, j):
    if i==0 and j==0:
        def v(s):
            deltaS = 1 / (np.sqrt((axis[0] ** 2 + s) * (axis[1] ** 2 + s) * (axis[2] ** 2 + s)))
            return deltaS
        ans = integrate.quad(v, L, np.inf)[0]
        return ans*2*pi*axis[0]*axis[1]*axis[2]
    elif j==0:
        def v(s):
            return 1 / (np.sqrt((axis[0] ** 2 + s) * (axis[1] ** 2 + s) * (axis[2] ** 2 + s))*(axis[i]**2 + s))
        ans = integrate.quad(v, L, np.inf) [0]
        return ans*2*pi*axis[0]*axis[1]*axis[2]
    else:
        def v(s):
            return 1 / (np.sqrt((axis[0] ** 2 + s) * (axis[1] ** 2 + s) * (axis[2] ** 2 + s))*(axis[i]**2 + s)*(axis[j]**2 + s))
        ans = integrate.quad(v, L, np.inf) [0]
        return ans*2*pi*axis[0]*axis[1]*axis[2]


def get_Sijkl(axis, I, Ii, Iij):
    Sijkl = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    t1 = np.identity(3)[i, j] * np.identity(3)[k, l] * (2 * 0.3 * Ii[i] - Ii[k] + (axis[i]**2) * Iij[k][i])
                    t2 = (np.identity(3)[i, k] * np.identity(3)[j, l] + np.identity(3)[i, l] * np.identity(3)[j, k]) * ((axis[i]**2) * Iij[i][j] - Ii[j] + (1 - 0.3) * (Ii[k] + Ii[l]))
                    Sijkl[i, j, k, l] = (t1 + t2) / (8 * np.pi * (1 - 0.3))
    return Sijkl


def lamb_der(x, a, lambda_):
    arr = []
    #Check if x in inside
    if (x[0]**2/(a[0]**2)) + (x[1]**2/(a[1]**2)) + (x[2]**2/(a[2]**2)) <= 1:
        return [0,0,0]
    denom = (x[0]**2/(a[0]**2 + lambda_)) + (x[1]**2/(a[1]**2 + lambda_)) + (x[2]**2/(a[2]**2 + lambda_))
    for i in range(3):
        num = (2*x[i])/(a[i]**2 + lambda_)
        arr.append(num/denom)
    return arr

# Compute double derivative matrix of lambda
def lamb_der2(x,a,lambda_,lambda_der):
    arr = np.zeros((3,3))
    #Check if x in inside
    if (x[0]**2/(a[0]**2)) + (x[1]**2/(a[1]**2)) + (x[2]**2/(a[2]**2)) <= 1:
        return arr
    denom = (x[0]**2/(a[0]**2 + lambda_)) + (x[1]**2/(a[1]**2 + lambda_)) + (x[2]**2/(a[2]**2 + lambda_))
    for i in range(3):
        for j in range(3):
            num = 2*denom*lambda_der[i]*lambda_der[j] - 2*(x[i])*lambda_der[j]/(a[i]**2 + lambda_) - 2*(x[j])*lambda_der[i]/(a[j]**2 + lambda_)
            arr[i,j] = num/denom
    return arr


#Compute derivative of Ii wrt to j direction
def Ii_j_(x, a, lambda_, lambda_der):
    arr = np.zeros((3,3))
    c = -2*np.pi*a[0]*a[1]*a[2]
    del_l = ((a[0]**2 + lambda_)*(a[1]**2 + lambda_)*(a[2]**2 + lambda_))**(1/2)
    for i in range(3):
        for j in range(3):
            arr[i,j] = c*lambda_der[j]/((a[i]**2 + lambda_)*del_l)
    return arr


# Compute derivative of Iij wrt to k direction
def Iij_k_(x, a, lambda_, lambda_der):
    arr = np.zeros((3,3,3))
    c = -2*np.pi*a[0]*a[1]*a[2]
    del_l = ((a[0]**2 + lambda_)*(a[1]**2 + lambda_)*(a[2]**2 + lambda_))**(1/2)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                arr[i,j,k] = c*lambda_der[k]/((a[i]**2 + lambda_)*(a[j]**2 + lambda_)*del_l)
    return arr


# Compute double partial derivative of Iij wrt to k and l direction
def Iij_kl_(x, a, lambda_, lambda_der, lambda_der2):
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
def Ii_jk_(x, a, lambda_, lambda_der, lambda_der2):
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

def solve_inside(a,b,c,eps11,eps22,eps33,eps12,eps23,eps31,E,nu):

    mu = E/(2*(1+nu))

    #Calculating theta and K

    theta = np.arcsin(np.sqrt((a**2 - c**2)/a**2))
    k = np.sqrt((a**2 - b**2)/(a**2 - c**2))

    F = sp.ellipkinc(theta, k**2)
    E = sp.ellipeinc(theta, k**2)
    I1 = (4*pi*a*b*c)*(F - E) / ((a**2 - b**2)*np.sqrt(a**2 - c**2))
    I3 = (4*pi*a*b*c)*((b*np.sqrt(a**2-c**2))/(a*c) - E)
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

    t3 = (((c**2)*(I3 - nu*I3 + 3*nu*I33 + nu*I23))/(8*pi*(1-nu)*(1-2*nu)) - (I1 - (nu/1-nu)*(I3-I2))/(8*pi) - (nu)/(1-2*nu))*eps33

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

    output_data = f"sigma11: {sigma11}, sigma22: {sigma22}, sigma33: {sigma33}, sigma12: {sigma12}, sigma13: {sigma31}, sigma23: {sigma23}"

    return output_data


try:
    # Load the form data passed from Express.js
    form_data = json.loads(sys.argv[1])

    # Taking Inputs
    a = float(form_data.get('a'))
    b = float(form_data.get('b'))
    c = float(form_data.get('c'))
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

    output_data = solve_inside(a,b,c,eps11,eps22,eps33,eps12,eps23,eps31,E,nu)
    print(output_data)

except Exception as e:
    print("error:", e)


try:
    # CODE FOR OUTSIDE PART (MATRIX)
    mu = E/(2*(1+nu))

    epsilon_star = [[eps11, eps12, 0], [0, eps22, eps23], [eps31, 0, eps33]]

    X = [7, 7, 7]  # Exterior Point
    lbd = find_lambda(X, axis)
        # print("Lambda is: ", lbd)

        #Calculating I, I1, I2, I3, I11, I22, I33, etc
    I = IIJ(axis, lbd, 0, 0)
    Ii = [0, 0, 0]
    for i in range(3):
        Ii[i] = IIJ(axis, lbd, i, 0)
    Iij = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for i in range(3):
        for j in range(3):
            Iij[i][j] = IIJ(axis, lbd, i, j)

    Sijkl = get_Sijkl(axis, I, Ii, Iij)

    lbd_der = lamb_der(X, axis, lbd)
    lbd_d_der = lamb_der2(X, axis, lbd, lbd_der)

    IIj = Ii_j_(X, axis, lbd, lbd_der)
    IIJk = Iij_k_(X, axis, lbd, lbd_der)
    IIJkl = Iij_kl_(X, axis, lbd, lbd_der, lbd_d_der)
    IIjk = Ii_jk_(X, axis, lbd, lbd_der, lbd_d_der)

    deltaij = np.identity(3)
    Dijkl = get_Dijkl(Sijkl, deltaij, X, axis, IIj, IIJk, IIJkl, IIjk)

    epsilon = np.zeros((3, 3))

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    epsilon[i][j] += Dijkl[i,j,k,l] * epsilon_star[k][l]


    sigma = np.zeros((3,3))
    epsilon_kk = epsilon[0][0] + epsilon[1][1] + epsilon[2][2]
    epsilon_star_kk = epsilon_star[0][0] + epsilon_star[1][1] + epsilon_star[2][2]
    lamda = 2*mu*nu/(1-2*nu)
    for i in range(3):
        for j in range(3):
            sigma[i][j] = 2*mu*(epsilon[i][j] - epsilon_star[i][j]) + lamda*(epsilon_kk - epsilon_star_kk)

    print(sigma)
except Exception as e:
    print("error:", e)