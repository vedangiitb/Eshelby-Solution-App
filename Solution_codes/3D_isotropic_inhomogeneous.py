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
                    t1 = np.identity(3)[i, j] * np.identity(3)[k, l] * (2 * nu * Ii[i] - Ii[k] + (axis[i]**2) * Iij[k][i])
                    t2 = (np.identity(3)[i, k] * np.identity(3)[j, l] + np.identity(3)[i, l] * np.identity(3)[j, k]) * ((axis[i]**2) * Iij[i][j] - Ii[j] + (1 - nu) * (Ii[k] + Ii[l]))
                    Sijkl[i, j, k, l] = (t1 + t2) / (8 * np.pi * (1 - nu))
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


# Taking Inputs

try:
    # Load the form data passed from Express.js
    form_data = json.loads(sys.argv[1])

    # Taking Inputs
    a = float(form_data.get('a'))
    b = float(form_data.get('b'))
    c = float(form_data.get('c'))
    eps_p_11 = float(form_data.get('eps11'))
    eps_p_22 = float(form_data.get('eps22'))
    eps_p_33 = float(form_data.get('eps33'))
    eps_p_12 = float(form_data.get('eps12'))
    eps_p_23 = float(form_data.get('eps23'))
    eps_p_31 = float(form_data.get('eps13'))
    E = float(form_data.get('ep'))
    E_star = float(form_data.get('eps'))
    nu = float(form_data.get('nu'))
    nu_star = float(form_data.get('nus'))
    mu = E/(2*(1+nu))
    mu_star = E/(2*(1+nu))
    
    targets = form_data.get('targets')

    axis = [a, b, c]


    eps_0_11 = 0
    eps_0_22 = 0
    eps_0_33 = 0
    eps_0_12 = 0
    eps_0_23 = 0
    eps_0_31 = 0
# '''
#     #Calculating theta and K

#     theta = np.arcsin(np.sqrt((a**2 - c**2)/a**2))
#     k = np.sqrt((a**2 - b**2)/(a**2 - c**2))

#     F = sp.ellipkinc(theta, k**2)
#     E = sp.ellipeinc(theta, k**2)
#     I1 = (4*pi*a*b*c)*(F - E) / ((a**2 - b**2)*np.sqrt(a**2 - c**2))
#     I3 = (4*pi*a*b*c)*((b*np.sqrt(a**2-c**2))/(a*c) - E)
#     I2 = 4*pi - I1 - I3


#     #Solvig equations for I11 and I13
#     I12 = (I2-I1)/(a**2-b**2)
#     x, y = sym.symbols('x y')
#     eq1 = sym.Eq(3*x+I12+y, 4*pi/(a**2))
#     eq2 = sym.Eq(3*(a**2)*x+(b**2)*I12+(c**2)*y, 3*I1)
#     ans = sym.solve([eq1, eq2], (x, y))
#     I11 = ans[x]
#     I13 = ans[y]

#     #Solving equations for I21 and I22
#     I23 = (I3-I2)/(b**2-c**2)
#     x, y = sym.symbols('x y')
#     eq1 = sym.Eq(3*x+I23+y, 4*pi/(b**2))
#     eq2 = sym.Eq(3*(b**2)*x+(c**2)*I23+(a**2)*y, 3*I2)
#     ans = sym.solve([eq1, eq2], (x, y))
#     I21 = ans[y]
#     I22 = ans[x]

#     #Solving equation for I33 and I32
#     I31 = (I1-I3)/(c**2-a**2)
#     x, y = sym.symbols('x y')
#     eq1 = sym.Eq(3*x+I31+y, 4*pi/(c**2))
#     eq2 = sym.Eq(3*(c**2)*x+(a**2)*I31+(b**2)*y, 3*I3)
#     ans = sym.solve([eq1, eq2], (x, y))
#     I32 = ans[y]
#     I33 = ans[x]
#     #Solving for sigma

#     #t1 t2 t3 are terms of sigma11
#     t1 = (((a**2)*(3*I11 - 3*nu*I11 + nu*I21 + nu*I31))/(8*pi*(1-nu)*(1-2*nu)) + ((I1 - (nu/(1-nu))*(I2+I3))/(8*pi)) - ((1-nu)/(1-2*nu)))*eps_p_11

#     t2 = (((b**2)*(I12-I12*nu + 3*nu*I22 + nu*I32))/(8*pi*(1-nu)*(1-2*nu)) - (I1 - (nu/1-nu)*(I2-I3))/(8*pi) - (nu)/(1-2*nu))*eps_p_22

#     t3 = (((c**2)*(I3 - nu*I3 + 3*nu*I33 + nu*I23))/(8*pi*(1-nu)*(1-2*nu)) - (I1 - (nu/1-nu)*(I3-I2))/(8*pi) - (nu)/(1-2*nu))*eps_p_33

#     sigma11 = 2*mu*(t1+t2+t3)

#     sigma12 = ((((a**2+b**2)*I12 + (1-2*nu)*(I1+I2))/(8*pi*(1-nu)))-1)*eps_p_12*2*mu
#     sigma23 = ((((b**2+c**2)*I12 + (1-2*nu)*(I2+I3))/(8*pi*(1-nu)))-1)*eps_p_23*2*mu
#     sigma31 = ((((a**2+c**2)*I12 + (1-2*nu)*(I1+I3))/(8*pi*(1-nu)))-1)*eps_p_31*2*mu


#     #q1 q2 q3 are for sigma22
#     q1 = (((b**2)*(3*I22 - 3*nu*I22 + nu*I32 + nu*I12))/(8*pi*(1-nu)*(1-2*nu)) + ((I2 - (nu/(1-nu))*(I1+I3))/(8*pi)) - ((1-nu)/(1-2*nu)))*eps_p_22
#     q2 = (((c**2)*(I23-I23*nu + 3*nu*I33 + nu*I13))/(8*pi*(1-nu)*(1-2*nu)) - (I2 - (nu/1-nu)*(I3-I1))/(8*pi) - (nu)/(1-2*nu))*eps_p_33
#     q3 = (((a**2)*(I21 - nu*I21 + 3*nu*I11 + nu*I31))/(8*pi*(1-nu)*(1-2*nu)) - (I2 - (nu/1-nu)*(I1-I3))/(8*pi) - (nu)/(1-2*nu))*eps_p_11

#     sigma22 = 2*mu*(q1+q2+q3)


#     #w1 w2 w3 are for sigma33
#     w1 = (((c**2)*(3*I33 - 3*nu*I33 + nu*I13 + nu*I23))/(8*pi*(1-nu)*(1-2*nu)) + ((I3 - (nu/(1-nu))*(I2+I1))/(8*pi)) - ((1-nu)/(1-2*nu)))*eps_p_33
#     w2 = (((a**2)*(I31-I31*nu + 3*nu*I11 + nu*I21))/(8*pi*(1-nu)*(1-2*nu)) - (I3 - (nu/1-nu)*(I1-I2))/(8*pi) - (nu)/(1-2*nu))*eps_p_11
#     w3 = (((b**2)*(I32 - nu*I32 + 3*nu*I22 + nu*I12))/(8*pi*(1-nu)*(1-2*nu)) - (I3 - (nu/1-nu)*(I2-I1))/(8*pi) - (nu)/(1-2*nu))*eps_p_22

#     sigma33 = 2*mu*(w1+w2+w3)

#     print("Inside inclusion:\n" + [sigma11,sigma12,sigma22,sigma23,sigma31,sigma33] + "\n")
# '''
    # CODE FOR OUTSIDE PART (MATRIX)

    epsilon_p = [[eps_p_11, eps_p_12, eps_p_31], [eps_p_12, eps_p_22, eps_p_23], [eps_p_31, eps_p_23, eps_p_33]]
    epsilon_not = [[eps_0_11, eps_0_12, eps_0_31], [eps_0_12, eps_0_22, eps_0_23], [eps_0_31, eps_0_23, eps_0_33]]

    X = [7, 7, 7]  # Any Point
    lbd = find_lambda(X, axis)
    print("Lambda is: ", lbd)

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

    eps_star_star_12 = ( 2*(mu-mu_star)*(epsilon_not[0][1]) + 2*mu_star*epsilon_p[0][1] ) / ( 4*(mu-mu_star)*Sijkl[0][1][0][1] + 2*mu )
    eps_star_star_23 = ( 2*(mu-mu_star)*(epsilon_not[1][2]) + 2*mu_star*epsilon_p[1][2] ) / ( 4*(mu-mu_star)*Sijkl[1][2][1][2] + 2*mu )
    eps_star_star_31 = ( 2*(mu-mu_star)*(epsilon_not[2][0]) + 2*mu_star*epsilon_p[2][0] ) / ( 4*(mu-mu_star)*Sijkl[2][0][2][0] + 2*mu )

    epsilon_star_star = [[0, eps_star_star_12, eps_star_star_31], [eps_star_star_12, 0, eps_star_star_23], [eps_star_star_31, eps_star_star_23, 0]]

    # SOLVING FOR epsilon star_star 11,22,33

    K = E/(3*(1-2*nu))
    K_star = E_star/(3*(1-2*nu_star))
    t = K_star/K

    A = [ [(1-t)*Sijkl[0][0][0][0]-1, (1-t)*Sijkl[0][0][1][1], (1-t)*Sijkl[0][0][2][2]], [(1-t)*Sijkl[1][1][0][0], (1-t)*Sijkl[1][1][1][1]-1, (1-t)*Sijkl[1][1][2][2]], [(1-t)*Sijkl[2][2][0][0], (1-t)*Sijkl[2][2][1][1], (1-t)*Sijkl[2][2][2][2]-1]]

    zz1, zz2, zz3 = 0, 0, 0

    for m in range(3):
        for n in range(3):
            if m != n:
                zz1 = zz1 + Sijkl[0][0][m][n] * epsilon_star_star[m][n]
                zz2 = zz2 + Sijkl[1][1][m][n] * epsilon_star_star[m][n]
                zz3 = zz3 + Sijkl[2][2][m][n] * epsilon_star_star[m][n]


    B = [[(t-1)*epsilon_not[0][0] - t*epsilon_p[0][0] - (1-t)*zz1], [(t-1)*epsilon_not[1][1] - t*epsilon_p[1][1] - (1-t)*zz2], [(t-1)*epsilon_not[2][2] - t*epsilon_p[2][2] - (1-t)*zz3]]

    A_inv = (np.linalg.inv(A))
    epsilon_star_11_22_33 = np.dot(A_inv, B)

    epsilon_star_star[0][0] = epsilon_p[0][0] + epsilon_star_11_22_33[0]
    epsilon_star_star[1][1] = epsilon_p[1][1] + epsilon_star_11_22_33[1]
    epsilon_star_star[2][2] = epsilon_p[2][2] + epsilon_star_11_22_33[2]

    epsilon = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

    for k in range(3):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    epsilon[k][l] = Sijkl[k][l][m][n] * epsilon_star_star[m][n]

    sigmaij = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    deltaij = np.identity(3)
    epsilon_kk = epsilon[0][0] + epsilon[1][1] + epsilon[2][2]

    for i in range(3):
        for j in range(3):
            sigmaij[i][j] = (E * nu * deltaij[i][j] * epsilon_kk) / ((1+nu) * (1-2*nu)) + (E * epsilon[i][j]) / (1+nu)

except Exception as e:
    print("error:", e)
