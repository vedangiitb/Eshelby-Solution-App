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

deltaij = np.identity(3)

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
    sigma21 = []
    sigma22 = []
    sigma23 = []
    sigma31 = []
    sigma32 = []
    sigma33 = []

    for d in data:
        sigma11.append(d[0])
        sigma12.append(d[1])
        sigma13.append(d[2])
        sigma21.append(d[3])
        sigma22.append(d[4])
        sigma23.append(d[5])
        sigma31.append(d[6])
        sigma32.append(d[7])
        sigma33.append(d[8])

    df = pd.DataFrame(list(zip(x,y,z,sigma11,sigma12,sigma13,sigma21,sigma22,sigma23,sigma31,sigma32,sigma33)),columns=["x","y","z",'sigma11','sigma12','sigma13','sigma21','sigma22','sigma23','sigma31','sigma32','sigma33'])

    df.to_csv('public/temp.csv',index=False)


# Function to find Lambda (the largest root)


def find_lambda(X, axis):
    ans = sym.symbols('ans')
    eqn = sym.Eq(((X[0]) ** 2) / ((axis[0]) ** 2 + ans) + ((X[1]) ** 2) / ((axis[1]) ** 2 + ans) + ((X[2]) ** 2) / (
                (axis[2]) ** 2 + ans), 1)
    solution = sym.solve(eqn, ans)
    sol = [complex(i) for i in solution]
    LAMBDA = 0
    for i in sol:
        LAMBDA = max(LAMBDA, i.real)
    return LAMBDA


def IIJ(axis, L, flag):
    theta = np.arcsin(np.sqrt((axis[0] ** 2 - axis[2] ** 2) / (axis[0] ** 2 + L)))
    k = np.sqrt((axis[0] ** 2 - axis[1] ** 2) / (axis[0] ** 2 - axis[2] ** 2))
    F = sp.ellipkinc(theta, k ** 2)
    E = sp.ellipeinc(theta, k ** 2)
    # print("theta", theta, k, F, E)
    arr = np.zeros(3)
    dels = np.sqrt((axis[0] ** 2 + L) * (axis[1] ** 2 + L) * (axis[2] ** 2 + L))
    c1 = (4 * pi * axis[0] * axis[1] * axis[2]) / (
            (axis[0] ** 2 - axis[1] ** 2) * (np.sqrt(axis[0] ** 2 - axis[2] ** 2)))
    arr[0] = c1 * (F - E)
    c2 = (4 * pi * axis[0] * axis[1] * axis[2]) / (
            (axis[1] ** 2 - axis[2] ** 2) * (np.sqrt(axis[0] ** 2 - axis[2] ** 2)))
    d1 = ((axis[1] ** 2 + L) * (np.sqrt(axis[0] ** 2 - axis[2] ** 2))) / dels
    arr[2] = c2 * (d1 - E)

    arr[1] = (4 * pi * axis[0] * axis[1] * axis[2]) / dels - arr[0] - arr[2]

    arr1 = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            if i != j:
                arr1[i][j] = (arr[j] - arr[i]) / (axis[i] ** 2 - axis[j] ** 2)
    for i in range(3):
        tmp = 0
        for j in range(3):
            tmp += arr1[i][j] / 3
        arr1[i][i] = (4 * pi * axis[0] * axis[1] * axis[2]) / (3 * (axis[i] ** 2 + L) * dels) - tmp

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
                    t1 = np.identity(3)[i][j] * np.identity(3)[k][l] * (
                                2 * nu * Ii[i] - Ii[k] + (axis[i] ** 2) * Iij[k][i])
                    t2 = (np.identity(3)[i][k] * np.identity(3)[j][l] + np.identity(3)[i][l] * np.identity(3)[j][k]) * (
                                (axis[i] ** 2) * Iij[i][j] - Ii[j] + (1 - nu) * (Ii[k] + Ii[l]))
                    Sijkl[i][j][k][l] = (t1 + t2) / (8 * np.pi * (1 - nu))
    Sijkl[i, j, k, l] = Sijkl[j, i, k, l]
    Sijkl[i, j, k, l] = Sijkl[i, j, l, k]
    return Sijkl


def lamb_der(x, a, lambda_):
    arr = []
    # To check if x in inside the inclusion
    if (x[0] ** 2 / (a[0] ** 2)) + (x[1] ** 2 / (a[1] ** 2)) + (x[2] ** 2 / (a[2] ** 2)) <= 1:
        return [0, 0, 0]
    denom = (x[0] ** 2 / (a[0] ** 2 + lambda_) ** 2) + (x[1] ** 2 / (a[1] ** 2 + lambda_) ** 2) + (
                x[2] ** 2 / (a[2] ** 2 + lambda_) ** 2)
    for i in range(3):
        num = (2 * x[i]) / (a[i] ** 2 + lambda_)
        arr.append(num / denom)
    return arr


# Compute double derivative matrix of lambda


def lamb_der2(x, a, lambda_, lambda_der):
    arr = np.zeros((3, 3))
    if (x[0] ** 2 / (a[0] ** 2)) + (x[1] ** 2 / (a[1] ** 2)) + (x[2] ** 2 / (a[2] ** 2)) <= 1:
        return arr
    for i in range(3):
        for j in range(3):
            arr[i, j] = ((2 * a[i] * lambda_der[i] * lambda_der[j]) - 2 * lambda_der[j] * (a[i] ** 2 + lambda_)) / (
                    (a[i] ** 2 + lambda_) * a[i])
    return arr


# Compute derivative of Ii wrt to j direction


def Ii_j_(x, a, lambda_, lambda_der):
    arr = np.zeros((3, 3))
    c = -2 * np.pi * a[0] * a[1] * a[2]
    del_l = ((a[0] ** 2 + lambda_) * (a[1] ** 2 + lambda_) * (a[2] ** 2 + lambda_)) ** (0.5)
    for i in range(3):
        for j in range(3):
            arr[i, j] = c * lambda_der[j] / ((a[i] ** 2 + lambda_) * del_l)
    return arr


# Compute derivative of Iij wrt to k direction
def Iij_k_(x, a, lambda_, lambda_der):
    arr = np.zeros((3, 3, 3))
    c = -2 * np.pi * a[0] * a[1] * a[2]
    del_l = ((a[0] ** 2 + lambda_) * (a[1] ** 2 + lambda_) * (a[2] ** 2 + lambda_)) ** (1 / 2)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                arr[i, j, k] = c * lambda_der[k] / ((a[i] ** 2 + lambda_) * (a[j] ** 2 + lambda_) * del_l)
    return arr


# Compute double partial derivative of Iij wrt to k and l direction
def Iij_kl_(x, a, lambda_, lambda_der, lambda_der2):
    arr = np.zeros((3, 3, 3, 3))
    c = -2 * np.pi * a[0] * a[1] * a[2]
    del_l = ((a[0] ** 2 + lambda_) * (a[1] ** 2 + lambda_) * (a[2] ** 2 + lambda_)) ** (1 / 2)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    arr[i, j, k, l] = (c / ((a[i] ** 2 + lambda_) * (a[j] ** 2 + lambda_) * (del_l))) * (
                            lambda_der2[k, l] - lambda_der[k] * lambda_der[l] * (
                            1 / (a[i] ** 2 + lambda_) + 1 / (a[j] ** 2 + lambda_) + 0.5 * (
                            1 / (a[0] ** 2 + lambda_) + 1 / (a[1] ** 2 + lambda_) + 1 / (
                            a[2] ** 2 + lambda_))))
    return arr


# Compute derivative of Ii wrt to j and k direction
def Ii_jk_(x, a, lambda_, lambda_der, lambda_der2):
    arr = np.zeros((3, 3, 3))
    c = -2 * np.pi * a[0] * a[1] * a[2]
    del_l = ((a[0] ** 2 + lambda_) * (a[1] ** 2 + lambda_) * (a[2] ** 2 + lambda_)) ** (1 / 2)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                arr[i, j, k] = (c / ((a[i] ** 2 + lambda_) * (del_l))) * (
                        lambda_der2[j, k] - lambda_der[j] * lambda_der[k] * (1 / (a[i] ** 2 + lambda_) + 0.5 * (
                        1 / (a[0] ** 2 + lambda_) + 1 / (a[1] ** 2 + lambda_) + 1 / (a[2] ** 2 + lambda_))))
    return arr


def get_Dijkl(s, delta, x, a, IIj, IIJk, IIJkl, IIjk):
    Dijkl = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    t1 = 8 * pi * (1 - nu) * s[i, j, k, l] + 2 * nu * delta[k, l] * x[i] * IIj[i, j]
                    t2 = (1 - nu) * (
                            delta[i, l] * x[k] * IIj[k, j] + delta[j, l] * x[k] * IIj[k, i] + delta[i, k] * x[l] *
                            IIj[l, j] + delta[j, k] * x[l] * IIj[l, i])
                    t3 = delta[i, j] * x[k] * (IIj[k, l] - (a[i] ** 2) * IIJk[k, i, l]) + (
                            delta[i, k] * x[j] + delta[j, k] * x[i]) * (IIj[j, l] - (a[i] ** 2) * IIJk[i, j, l])
                    t4 = (delta[i, l] * x[j] + delta[j, l] * x[i]) * (IIj[j, k] - (a[i] ** 2) * IIJk[i, j, k]) + x[i] * \
                         x[j] * (IIjk[j, l, k] - (a[i] ** 2) * IIJkl[i, j, l, k])
                    Dijkl[i, j, k, l] = (t1 + t2 - t3 - t4) / (8 * pi * (1 - nu))
    return Dijkl






def calcStress(x1,x2,x3,chosenDir):
    eps_0_11 = 0
    eps_0_22 = 0
    eps_0_33 = 0
    eps_0_12 = 0
    eps_0_23 = 0
    eps_0_31 = 0



    X = [x1, x2, x3]

    Cijkl = np.zeros((3, 3, 3, 3))

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    Cijkl[i][j][k][l] = lamda * deltaij[i][j] * deltaij[k][l] + mu * (
                                deltaij[i][k] * deltaij[j][l] + deltaij[i][l] * deltaij[j][k])
                    Cijkl[k][l][i][j] = Cijkl[i][j][k][l]

    Cijkl_star = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    Cijkl_star[i][j][k][l] = lamda_star * deltaij[i][j] * deltaij[k][l] + mu_star * (
                                deltaij[i][k] * deltaij[j][l] + deltaij[i][l] * deltaij[j][k])
                    Cijkl_star[k][l][i][j] = Cijkl_star[i][j][k][l]

    epsilon_p = [[eps_p_11, eps_p_12, eps_p_31],
                [eps_p_12, eps_p_22, eps_p_23],
                [eps_p_31, eps_p_23, eps_p_33]]

    epsilon_not = [[eps_0_11, eps_0_12, eps_0_31],
                [eps_0_12, eps_0_22, eps_0_23],
                [eps_0_31, eps_0_23, eps_0_33]]

    sigma_not_ij = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]  # stress at infinity

    ab = np.linspace(0, 200, num=100)
    s = []
    s1 = []
    s2 = []

    lbd = find_lambda(X, axis)

    # Calculating  I1, I2, I3, I11, I22, I33, etc

    Ii = IIJ(axis, lbd, 0)
    Iij = IIJ(axis, lbd, 2)

    Sijkl = get_Sijkl(axis, Ii, Iij)

    lbd_der = lamb_der(X, axis, lbd)
    lbd_d_der = lamb_der2(X, axis, lbd, lbd_der)

    IIj = Ii_j_(X, axis, lbd, lbd_der)
    IIJk = Iij_k_(X, axis, lbd, lbd_der)
    IIJkl = Iij_kl_(X, axis, lbd, lbd_der, lbd_d_der)
    IIjk = Ii_jk_(X, axis, lbd, lbd_der, lbd_d_der)

    Dijkl = get_Dijkl(Sijkl, deltaij, X, axis, IIj, IIJk, IIJkl, IIjk)

    eps_star_star_12 = (2 * (mu - mu_star) * (epsilon_not[0][1]) + 2 * mu_star * epsilon_p[0][1]) / (
            4 * (mu - mu_star) * Sijkl[0][1][0][1] + 2 * mu)
    eps_star_star_23 = (2 * (mu - mu_star) * (epsilon_not[1][2]) + 2 * mu_star * epsilon_p[1][2]) / (
            4 * (mu - mu_star) * Sijkl[1][2][1][2] + 2 * mu)
    eps_star_star_31 = (2 * (mu - mu_star) * (epsilon_not[2][0]) + 2 * mu_star * epsilon_p[2][0]) / (
            4 * (mu - mu_star) * Sijkl[2][0][2][0] + 2 * mu)

    epsilon_star_star = [[0, eps_star_star_12, eps_star_star_31],
                        [eps_star_star_12, 0, eps_star_star_23],
                        [eps_star_star_31, eps_star_star_23, 0]]

    # SOLVING FOR epsilon star_star 11,22,33

    K = E / (3 * (1 - 2 * nu))
    K_star = E_star / (3 * (1 - 2 * nu_star))
    t = K_star / K

    A = [[(1 - t) * Sijkl[0][0][0][0] - 1, (1 - t) * Sijkl[0][0][1][1], (1 - t) * Sijkl[0][0][2][2]],
        [(1 - t) * Sijkl[1][1][0][0], (1 - t) * Sijkl[1][1][1][1] - 1, (1 - t) * Sijkl[1][1][2][2]],
        [(1 - t) * Sijkl[2][2][0][0], (1 - t) * Sijkl[2][2][1][1], (1 - t) * Sijkl[2][2][2][2] - 1]]

    zz1, zz2, zz3 = 0, 0, 0

    for m in range(3):
        for n in range(3):
            if m != n:
                zz1 = zz1 + Sijkl[0][0][m][n] * epsilon_star_star[m][n]
                zz2 = zz2 + Sijkl[1][1][m][n] * epsilon_star_star[m][n]
                zz3 = zz3 + Sijkl[2][2][m][n] * epsilon_star_star[m][n]

    B = [[(t - 1) * epsilon_not[0][0] - t * epsilon_p[0][0] - (1 - t) * zz1],
        [(t - 1) * epsilon_not[1][1] - t * epsilon_p[1][1] - (1 - t) * zz2],
        [(t - 1) * epsilon_not[2][2] - t * epsilon_p[2][2] - (1 - t) * zz3]]

    A_inv = (np.linalg.inv(A))

    epsilon_star_11_22_33 = np.dot(A_inv, B)

    epsilon_star_star[0][0] = epsilon_star_11_22_33[0]
    epsilon_star_star[1][1] = epsilon_star_11_22_33[1]
    epsilon_star_star[2][2] = epsilon_star_11_22_33[2]

    print(f"For the points {x1,x2,x3}",end=" ")
    def calc_exterior():
        epsilon_out = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        epsilon_out[i][j] += Dijkl[i][j][k][l] * epsilon_star_star[k][l]

        sigmaij_out = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

        epsilon_out_kk = epsilon_out[0][0] + epsilon_out[1][1] + epsilon_out[2][2]
        '''sigmaij_out[0][0] = (E / ((1 + nu) * (1 - 2 * nu))) * (
                    (1 - nu) * epsilon_out[0][0] + nu * epsilon_out[1][1] + nu * epsilon_out[2][2])
        sigmaij_out[1][1] = (E / ((1 + nu) * (1 - 2 * nu))) * (
                    nu * epsilon_out[0][0] + (1 - nu) * epsilon_out[1][1] + nu * epsilon_out[2][2])
        sigmaij_out[2][2] = (E / ((1 + nu) * (1 - 2 * nu))) * (
                    nu * epsilon_out[0][0] + nu * epsilon_out[1][1] + (1 - nu) * epsilon_out[2][2])
        sigmaij_out[0][1] = (E / (1 + nu)) * epsilon_out[0][1]
        sigmaij_out[0][2] = (E / (1 + nu)) * epsilon_out[0][2]
        sigmaij_out[1][2] = (E / (1 + nu)) * epsilon_out[1][2]
        sigmaij_out[1][0] = sigmaij_out[0][1]
        sigmaij_out[2][0] = sigmaij_out[0][2]
        sigmaij_out[2][1] = sigmaij_out[1][2]'''

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        for m in range(3):
                            for n in range(3):
                                sigmaij_out[i][j] += Cijkl[i][j][k][l]*Dijkl[k][l][m][n]*epsilon_star_star[m][n]

        for i in range(3):
            for j in range(3):
                sigmaij_out[i][j] += sigma_not_ij[i][j]

        return sigmaij_out
    def calc_interior():
        # calculating epsilon from Sijkl and epsilon star star
        '''
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            epsilon[k][l] += Sijkl[k][l][m][n] * epsilon_star_star[m][n]'''

        # CALCULATION OF SIGMAij

        '''
            sigmaij[0][0] = (E / ((1 + nu) * (1 - 2 * nu))) * ((1 - nu) * epsilon[0][0] + nu * epsilon[1][1] + nu * epsilon[2][2])
            sigmaij[1][1] = (E / ((1 + nu) * (1 - 2 * nu))) * (nu * epsilon[0][0] + (1 - nu) * epsilon[1][1] + nu * epsilon[2][2])
            sigmaij[2][2] = (E / ((1 + nu) * (1 - 2 * nu))) * (nu * epsilon[0][0] + nu * epsilon[1][1] + (1 - nu) * epsilon[2][2])
            sigmaij[0][1] = (E / (1 + nu)) * epsilon[0][1]
            sigmaij[0][2] = (E / (1 + nu)) * epsilon[0][2]
            sigmaij[1][2] = (E / (1 + nu)) * epsilon[1][2]
            sigmaij[1][0] = sigmaij[0][1]
            sigmaij[2][0] = sigmaij[0][2]
            sigmaij[2][1] = sigmaij[1][2]'''

        sigmaij = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]  # initialized a sigmaij array
        epsilon = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]  # initialized an epsilon array
        for i in range(3):
            for j in range(3):
                for m in range(3):
                    for n in range(3):
                        epsilon[i][j] += Sijkl[i][j][m][n]*epsilon_star_star[m][n]
        #print("EPSILON")
        #print(epsilon)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        tmp = 0
                        for m in range(3):
                            for n in range(3):
                                tmp += Sijkl[k][l][m][n]*epsilon_star_star[m][n] - epsilon_star_star[m][n]
                        sigmaij[i][j] = Cijkl[i][j][k][l]*tmp


        for i in range(3):
            for j in range(3):
                sigmaij[i][j] += sigma_not_ij[i][j]

        return sigmaij

    if (X[0]**2/axis[0]**2 + X[1]**2/axis[1]**2 + X[2]**2/axis[2]**2) <= 1:
        sigmaij = calc_interior()
        print("[the point is inside]")
    else:
        sigmaij = calc_exterior()
        print("[the point is outside]")
    
    sigma11 =sigmaij[0][0][0]
    sigma12 =sigmaij[0][1][0]
    sigma13 =sigmaij[0][2][0]
    sigma21 =sigmaij[1][0][0]
    sigma22 =sigmaij[1][1][0]
    sigma23 =sigmaij[1][2][0]
    sigma31 =sigmaij[2][0][0]
    sigma32 =sigmaij[2][1][0]
    sigma33 =sigmaij[2][2][0]
    data = [sigma11,sigma12,sigma13,sigma21,sigma22,sigma23,sigma31,sigma32,sigma33]
    
    print(f"sigma11 = {sigma11}\nsigma12 = {sigma12}\nsigma13 = {sigma13}\nsigma21 = {sigma21}\nsigma22 = {sigma22}\nsigma23 = {sigma23}\nsigma31 = {sigma31}\nsigma32 = {sigma32}\nsigma33 = {sigma33}")
    print()
    return data[chosenDir]

try:
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
    plottype = float(form_data.get('plottype'))
    mu = E/(2*(1+nu))
    mu_star = E/(2*(1+nu))
    

    axis = [a, b, c]
    lamda = 2 * mu * nu / (1 - 2 * nu)

    lamda_star = 2 * mu_star * nu_star / (1 - 2 * nu_star)
    
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
    elif plottype == '21':
        chosenDir = 6
    elif plottype == '31':
        chosenDir = 7
    elif plottype == '32':
        chosenDir = 8
    else:
        chosenDir = 0

    # step = 4*a/13
    print(plottype)
    print(chosenDir)

    x = np.linspace(-2*a,2*a,8)
    y = np.linspace(-2*b,2*b,8)
    z = np.linspace(-2*c,2*c,8)
    B,A, C = np.meshgrid(x,y,z)
    # print(X.shape)
    print("plotting..")
    stressArr = np.zeros(B.shape)

    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            for k in range(A.shape[2]):
                stressArr[i][j][k] = calcStress(A[i][j][k],B[i][j][k],C[i][j][k],chosenDir)/mu
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

except Exception as e:
    print("error:", e)
