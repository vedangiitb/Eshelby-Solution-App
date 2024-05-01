try:
    import numpy as np
    import scipy.special as sp
    import sympy as sym
    import scipy.integrate as integrate
    pi = np.pi
    import sys
    import json
    import pandas as pd
except Exception as e:
    print(e)
# Function to find Lambda (the largest root)
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


def IIJ(axis, L, i, j):
    if i == 0 and j == 0:
        def v(s):
            deltaS = 1 / (np.sqrt((axis[0] ** 2 + s) * (axis[1] ** 2 + s) * (axis[2] ** 2 + s)))
            return deltaS

        ans = integrate.quad(v, L, np.inf)[0]
        return ans * 2 * pi * axis[0] * axis[1] * axis[2]
    elif j == 0:
        def v(s):
            return 1 / (np.sqrt((axis[0] ** 2 + s) * (axis[1] ** 2 + s) * (axis[2] ** 2 + s)) * (axis[i - 1] ** 2 + s))

        ans = integrate.quad(v, L, np.inf)[0]
        return ans * 2 * pi * axis[0] * axis[1] * axis[2]
    else:
        def v(s):
            return 1 / (np.sqrt((axis[0] ** 2 + s) * (axis[1] ** 2 + s) * (axis[2] ** 2 + s)) * (
                        axis[i - 1] ** 2 + s) * (axis[j - 1] ** 2 + s))

        ans = integrate.quad(v, L, np.inf)[0]
        return ans * 2 * pi * axis[0] * axis[1] * axis[2]


def get_Sijkl(axis, I, Ii, Iij):
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
    # Check if x in inside
    if (x[0] ** 2 / (a[0] ** 2)) + (x[1] ** 2 / (a[1] ** 2)) + (x[2] ** 2 / (a[2] ** 2)) <= 1:
        return arr
    denom = (x[0] ** 2 / (a[0] ** 2 + lambda_) ** 2) + (x[1] ** 2 / (a[1] ** 2 + lambda_) ** 2) + (
                x[2] ** 2 / (a[2] ** 2 + lambda_) ** 2)
    for i in range(3):
        for j in range(3):
            num = 2 * denom * lambda_der[i] * lambda_der[j] - 2 * (x[i]) * lambda_der[j] / (a[i] ** 2 + lambda_) - 2 * (
            x[j]) * lambda_der[i] / (a[j] ** 2 + lambda_)
            arr[i, j] = num / denom
    return arr


# Compute derivative of Ii wrt to j direction
def Ii_j_(x, a, lambda_, lambda_der):
    arr = np.zeros((3, 3))
    c = -2 * np.pi * a[0] * a[1] * a[2]
    del_l = ((a[0] ** 2 + lambda_) * (a[1] ** 2 + lambda_) * (a[2] ** 2 + lambda_)) ** (1 / 2)
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
                    t2 = (1 - nu) * (delta[i, l] * x[k] * IIj[k, j] + delta[j, l] * x[k] * IIj[k, i] + delta[i, k] * x[l] * IIj[l, j] + delta[j, k] * x[l] * IIj[l, i])
                    t3 = (-1) * delta[i, j] * x[k] * (IIj[k, l] - (a[i] ** 2) * IIJk[k, i, l]) - (delta[i, k] * x[j] + delta[j, k] * x[i]) * (IIj[j, l] - (a[i] ** 2) * IIJk[i, j, l])
                    t4 = (-1) * (delta[i, l] * x[j] + delta[j, l] * x[i]) * (IIj[j, k] - (a[i] ** 2) * IIJk[i, j, k]) - x[i] * x[j] * (IIjk[j, l, k] - (a[i] ** 2) * IIJkl[i, j, l, k])
                    Dijkl[i, j, k, l] = (t1 + t2 + t3 + t4) / (8 * pi * (1 - nu))
    return Dijkl


def calc_interior():
    sigmaij = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]  # initialized a sigmaij array
    epsilon = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]  # initialized an epsilon array

    # calculating epsilon from Sijkl and epsilon star star

    for k in range(3):
        for l in range(3):
            for m in range(3):
                for n in range(3):
                    epsilon[k][l] += Sijkl[k][l][m][n] * epsilon_star_star[m][n]

    epsilon_kk = epsilon[0][0] + epsilon[1][1] + epsilon[2][2]

    # print("Epsilon is:")

    # for i in range(3):
    #     for j in range(3):
    #         print(epsilon[i][j], end=" ")
    #     print()
    # print()
    for i in range(3):
        for j in range(3):
            epsilon[i][j] -= epsilon_star_star[i][j]

    # CALCULATION OF SIGMAij

    sigmaij[0][0] = (E / ((1 + nu) * (1 - 2 * nu))) * ((1 - nu) * epsilon[0][0] + nu * epsilon[1][1] + nu * epsilon[2][2])
    sigmaij[1][1] = (E / ((1 + nu) * (1 - 2 * nu))) * (nu * epsilon[0][0] + (1 - nu) * epsilon[1][1] + nu * epsilon[2][2])
    sigmaij[2][2] = (E / ((1 + nu) * (1 - 2 * nu))) * (nu * epsilon[0][0] + nu * epsilon[1][1] + (1 - nu) * epsilon[2][2])
    sigmaij[0][1] = (E / (1 + nu)) * epsilon[0][1]
    sigmaij[0][2] = (E / (1 + nu)) * epsilon[0][2]
    sigmaij[1][2] = (E / (1 + nu)) * epsilon[1][2]
    sigmaij[1][0] = sigmaij[0][1]
    sigmaij[2][0] = sigmaij[0][2]
    sigmaij[2][1] = sigmaij[1][2]

    # print("SIGMAij for inside part by formula is:")

    sigma11 =sigmaij[0][0][0] + sigma_not_ij[0][0]
    sigma12 =sigmaij[0][1][0] + sigma_not_ij[0][1]
    sigma13 =sigmaij[0][2][0] + sigma_not_ij[0][2]
    sigma21 =sigmaij[1][0][0] + sigma_not_ij[1][0]
    sigma22 =sigmaij[1][1][0] + sigma_not_ij[1][1]
    sigma23 =sigmaij[1][2][0] + sigma_not_ij[1][2]
    sigma31 =sigmaij[2][0][0] + sigma_not_ij[2][0]
    sigma32 =sigmaij[2][1][0] + sigma_not_ij[2][1]
    sigma33 =sigmaij[2][2][0] + sigma_not_ij[2][2]
    data = [sigma11,sigma12,sigma13,sigma21,sigma22,sigma23,sigma31,sigma32,sigma33]


    print(f"sigma11 = {sigma11}\nsigma12 = {sigma12}\nsigma13 = {sigma13}\nsigma21 = {sigma21}\nsigma22 = {sigma22}\nsigma23 = {sigma23}\nsigma31 = {sigma31}\nsigma32 = {sigma32}\nsigma33 = {sigma33}")

    return data
    # for i in range(3):
    #     for j in range(3):
    #         print(sigmaij[i][j][0] + sigma_not_ij[i][j], end=" ")
    #     print()

    # print()


def calc_exterior():
    epsilon_out = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    epsilon_out[i][j] += Dijkl[i][j][k][l] * epsilon_star_star[k][l]

    sigmaij_out = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    epsilon_out_kk = epsilon_out[0][0] + epsilon_out[1][1] + epsilon_out[2][2]

    sigmaij_out[0][0] = (E / ((1 + nu) * (1 - 2 * nu))) * (
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
    sigmaij_out[2][1] = sigmaij_out[1][2]

    sigma11 =sigmaij_out[0][0][0] + sigma_not_ij[0][0]
    sigma12 =sigmaij_out[0][1][0] + sigma_not_ij[0][1]
    sigma13 =sigmaij_out[0][2][0] + sigma_not_ij[0][2]
    sigma21 =sigmaij_out[1][0][0] + sigma_not_ij[1][0]
    sigma22 =sigmaij_out[1][1][0] + sigma_not_ij[1][1]
    sigma23 =sigmaij_out[1][2][0] + sigma_not_ij[1][2]
    sigma31 =sigmaij_out[2][0][0] + sigma_not_ij[2][0]
    sigma32 =sigmaij_out[2][1][0] + sigma_not_ij[2][1]
    sigma33 =sigmaij_out[2][2][0] + sigma_not_ij[2][2]
    data = [sigma11,sigma12,sigma13,sigma21,sigma22,sigma23,sigma31,sigma32,sigma33]

    print(f"sigma11 = {sigma11}\nsigma12 = {sigma12}\nsigma13 = {sigma13}\nsigma21 = {sigma21}\nsigma22 = {sigma22}\nsigma23 = {sigma23}\nsigma31 = {sigma31}\nsigma32 = {sigma32}\nsigma33 = {sigma33}")

    return data
    # print("SIGMAij for exterior point by formula is:")

    # for i in range(3):
    #     for j in range(3):
    #         print(sigmaij_out[i][j][0]+ sigma_not_ij[i][j], end=" ")
    #     print()

    # print()

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
    lamda = 2 * mu * nu / (1 - 2 * nu)
    lamda_star = 2 * mu_star * nu_star / (1 - 2 * nu_star)

    targets = form_data.get('targets')

    axis = [a, b, c]


    eps_0_11 = 0
    eps_0_22 = 0
    eps_0_33 = 0
    eps_0_12 = 0
    eps_0_23 = 0
    eps_0_31 = 0

    targets = form_data.get('targets')

    if len(targets)==0:
        print("No points entered! Please enter the points where you want to compute the results")

    opData = []

    for x1,x2,x3 in targets:

        X = [x1, x2, x3]

        Cijkl = np.zeros((3, 3, 3, 3))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        Cijkl[i][j][k][l] = lamda * deltaij[i][j] * deltaij[k][l] + mu * (
                                    deltaij[i][k] * deltaij[j][l] + deltaij[i][l] * deltaij[j][k])
                        Cijkl[k][l][i][j] = Cijkl[i][j][k][l]

        # print("Cijkl=", Cijkl)

        Cijkl_star = np.zeros((3, 3, 3, 3))
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        Cijkl_star[i][j][k][l] = lamda_star * deltaij[i][j] * deltaij[k][l] + mu_star * (
                                    deltaij[i][k] * deltaij[j][l] + deltaij[i][l] * deltaij[j][k])
                        Cijkl_star[k][l][i][j] = Cijkl_star[i][j][k][l]

        # print("Cijkl star=", Cijkl_star)


        epsilon_p = [[eps_p_11, eps_p_12, eps_p_31], [eps_p_12, eps_p_22, eps_p_23], [eps_p_31, eps_p_23, eps_p_33]]
        epsilon_not = [[eps_0_11, eps_0_12, eps_0_31], [eps_0_12, eps_0_22, eps_0_23], [eps_0_31, eps_0_23, eps_0_33]]

        lbd = find_lambda(X, axis)
        # print("Lambda is: ", lbd)
        # print()

        # Calculating I, I1, I2, I3, I11, I22, I33, etc

        I = IIJ(axis, lbd, 0, 0)
        Ii = [0, 0, 0]
        for i in range(3):
            Ii[i] = IIJ(axis, lbd, i + 1, 0)
        Iij = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        for i in range(3):
            for j in range(3):
                Iij[i][j] = IIJ(axis, lbd, i + 1, j + 1)

        lbd_der = lamb_der(X, axis, lbd)
        lbd_d_der = lamb_der2(X, axis, lbd, lbd_der)

        IIj = Ii_j_(X, axis, lbd, lbd_der)
        IIJk = Iij_k_(X, axis, lbd, lbd_der)
        IIJkl = Iij_kl_(X, axis, lbd, lbd_der, lbd_d_der)
        IIjk = Ii_jk_(X, axis, lbd, lbd_der, lbd_d_der)

        Sijkl = get_Sijkl(axis, I, Ii, Iij)

        #print("Sijkl is dsjgfegjegfg: ")
        #print(Sijkl)
        #print("--------------------------------------------------------------------------------------------------------------")
        #print()

        Dijkl = get_Dijkl(Sijkl, deltaij, X, axis, IIj, IIJk, IIJkl, IIjk)

        eps_star_star_12 = (2 * (mu - mu_star) * (epsilon_not[0][1]) + 2 * mu_star * epsilon_p[0][1]) / (
                    4 * (mu - mu_star) * Sijkl[0][1][0][1] + 2 * mu)
        eps_star_star_23 = (2 * (mu - mu_star) * (epsilon_not[1][2]) + 2 * mu_star * epsilon_p[1][2]) / (
                    4 * (mu - mu_star) * Sijkl[1][2][1][2] + 2 * mu)
        eps_star_star_31 = (2 * (mu - mu_star) * (epsilon_not[2][0]) + 2 * mu_star * epsilon_p[2][0]) / (
                    4 * (mu - mu_star) * Sijkl[2][0][2][0] + 2 * mu)

        epsilon_star_star = [[0, eps_star_star_12, eps_star_star_31], [eps_star_star_12, 0, eps_star_star_23],
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


        sigma_not_ij = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]  # stress at infinity
        if x1 ** 2 / a ** 2 + x2 ** 2 / b ** 2 + x3 ** 2 / c ** 2 <= 1:
            print(f"For the points {x1,x2,x3}",end=" ")
            print("[the point is inside]")
            data = calc_interior()
            opData.append(data)
        else:
            print(f"For the points {x1,x2,x3}",end=" ")
            print("[the point is outside]")
            data = calc_exterior()
            opData.append(data)

    saveData(targets,opData)
except Exception as e:
    print("error:", e)
