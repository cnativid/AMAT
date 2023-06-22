import numpy as np
import matplotlib.pyplot as plt
from labellines import labelLines


def Class_Shape(weights, X, N1 = 0.5, N2 = 1, te_thickness = 0):
    C = (X**N1)*((1-X)**N2)
    n = np.size(weights) - 1
    K = np.zeros(n+1)
    for i in range(n+1):
        K[i] = np.math.factorial(n)/(np.math.factorial(n-i)*np.math.factorial(i))
    S = np.zeros(np.size(X))
    for i,x in enumerate(X):
        for j in range(n+1):
            S[i] += weights[j]*K[j]*(x**j)*((1-x)**(n-j))
    return C*S + X*te_thickness

def CST_Airfoil(lower_weights, upper_weights, N = 200, te_thickness = 0):
    i = np.linspace(0,N,N+1)
    zeta = 2*np.pi/N*i
    X = 0.5*(np.cos(zeta)+1)
    N1, N2 = 0.5, 1
    zero_index = X.argmin()
    X_lower = X[0:zero_index]
    X_upper = X[zero_index:N+1]
    Y_lower = Class_Shape(lower_weights, X_lower)
    Y_upper = Class_Shape(upper_weights, X_upper)
    Y = np.concatenate((Y_lower, Y_upper), axis=None)
    # print(Y)
    # plt.plot(X,Y)
    # # labelLines(plt.gca().get_lines(), zorder=2.5,fontsize = 7)
    X = np.round(X,6)
    Y = np.round(Y,6)
    dydx_lower = Full_Deriv(X_lower,Y_lower)
    dydx_upper = Full_Deriv(X_upper,Y_upper)
    # plt.plot(X_lower,dydx_lower, label = 'lower')
    # plt.plot(X_upper,dydx_upper, label = 'upper')
    # plt.plot()
    # labelLines(plt.gca().get_lines(), zorder=2.5,fontsize = 7)
    # plt.show()
    plt.plot(X_lower,Full_Deriv(X_lower,dydx_lower), label = 'lower')
    plt.plot(X_upper,Full_Deriv(X_upper,dydx_upper), label = 'upper')
    plt.plot()
    labelLines(plt.gca().get_lines(), zorder=2.5,fontsize = 7)
    plt.show()
    
    return X,Y

def Full_Deriv(X,Y):
    n = np.size(X)
    dX = X[1] - X[0]
    dydx = np.zeros(n)
    dydx[0] = (Y[1] - Y[0])/dX
    dydx[n-1] = (Y[n-1] - Y[n-2])/dX
    dydx[1:n-1] = (Y[2:n] - Y[0:n-2])/(2*dX)
    return dydx

X, Y = CST_Airfoil([0.17254, 0.22561, 0.20726],[-0.17254, -0.05783, -0.07251], 400)
with open('CSTfoil.txt','w') as f:
    f.write('CSTfoil')
    f.write('\n')
    for i, x in enumerate(X):
        f.write(f'{x} {Y[i]}')
        f.write('\n')
