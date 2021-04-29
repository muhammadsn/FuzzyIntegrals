import FuzzyIntegrals as FI


if __name__ == '__main__':
    # Values of g(x) , f(x) for all movies and all users
    g_x = [0.3, 0.6, 0.7, 0.3]
    f_A = [1.0, 0.8, 0.1, 0.8]
    f_B = [0.5, 0.5, 0.3, 0.7]
    f_C = [0.3, 0.3, 0.2, 0.4]
    f_D = [0.3, 0.3, 0.8, 0.8]

    # Values of g(x) , f(x) for all movies and users x1 , x2
    g_12 = [0.3, 0.6]
    f_A_12 = [1.0, 0.8]
    f_B_12 = [0.5, 0.5]
    f_C_12 = [0.3, 0.3]
    f_D_12 = [0.3, 0.3]

    # Values of g(x) , f(x) for all movies and users x3 , x4
    g_34 = [0.7, 0.3]
    f_A_34 = [0.1, 0.8]
    f_B_34 = [0.3, 0.7]
    f_C_34 = [0.2, 0.4]
    f_D_34 = [0.8, 0.8]

    # Values of g(x) , f(x) for all movies and users x2 , x3
    g_23 = [0.6, 0.7]
    f_A_23 = [0.8, 0.1]
    f_B_23 = [0.5, 0.3]
    f_C_23 = [0.3, 0.2]
    f_D_23 = [0.3, 0.8]


    A = FI.FuzzyIntegrals()

    print("Sugeno: ", A.sugeno(mu=g_23, f=f_A_23))
    print("Choquet: ", A.choquet(mu=g_23, f=f_A_23))

    print(A.get_params())
