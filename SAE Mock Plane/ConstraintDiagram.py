import numpy as np
import matplotlib.pyplot as plt
import math
import time
from ambiance import Atmosphere

g = 32.17 # ft/s^2


def Constraint_TO(WS):
    eta_p = 0.45 # Gudmunsson, 
    S_G = 100 # takeoff ground run, ft
    C_LTO =  0.6
    C_DTO =  0.02 
    C_Lmax = 1.5
    mu = 0.03
    rho = 0.0023769 # sea level, slug/ft^3
    V_LOF = 1.1*np.sqrt(2*WS/(rho*C_Lmax))
    return 1/((1.21/(g*rho*C_Lmax*S_G)*WS+0.605/C_Lmax*(C_DTO-mu*C_LTO) + mu)*V_LOF/np.sqrt(2)/(eta_p*550/745.7))

def Constraint_LDG():
    S_LGR = 200.0 # ft
    mu = 0.3
    rho = 0.0023769 # sea level, slug/ft^3
    C_LLD =  0
    C_DLD =  0.1
    C_Lmax = 1.5

    res = [1.0,1.0]
    def dx(x,V_stall):
        return np.array([x[1], -g*(mu+x[1]**2/V_stall**2*((C_DLD-C_LLD*mu)/C_Lmax))])
    
    
    dt=1e-2
    V_stall = [1,1+1e-1]

    k = 1.3 # touchdown velocity multiplier
    while np.abs(res[1]) > 1e-3:
        # x0 = np.array([0,42.4],float)
        x1 = [np.array([0,k*V_stall[0]],float)]
        x2 = [np.array([0,k*V_stall[1]],float)]
        # t = [0]
        i=0
        while x1[i][1] > 0:
            # t.append(t[i]+dt)
            x1.append(x1[i] + dt*dx(x1[i],V_stall[0]))
            i += 1 

        i=0
        while x2[i][1] > 0:
            # t.append(t[i]+dt)
            x2.append(x2[i] + dt*dx(x2[i],V_stall[1]))
            i += 1

        print(V_stall,res)
        res = [x1[-1][0]-S_LGR,x2[-1][0]-S_LGR]
        V_stall = [V_stall[1] , V_stall[1] - res[1]/((res[1]-res[0])/(V_stall[1]-V_stall[0]))]
    V_stall = V_stall[0]
    return rho*V_stall**2*C_Lmax/2
    # return S_LGR/80*C_Lmax

# def Constraint_ngTurn
def Constraint_MaxWeight():
    P_max = 1000 # our powerplant max power, Watts
    W_max = 55 # MTOW, lb
    return W_max/P_max

def Constraint_Cruise():
    eta = 0.3 # prop efficiency at trim
    P_cruise = 0.5 # frac of max power for cruise
    V_trim = 70
    LD = 4 # trim L/D
    return eta*LD/P_cruise/V_trim

if __name__ == "__main__":
    plt.figure(figsize=[5,6])
    N = 101
    WS = np.linspace(1,6,N)
    WP = np.linspace(0,.03,N)
    plt.plot(WS,Constraint_TO(WS), label = 'TO')
    plt.plot(Constraint_LDG()*np.ones(N),WP, label = 'LDG')
    plt.plot(WS,Constraint_MaxWeight()*np.ones(N), label = 'Max Weight')
    plt.plot(WS,Constraint_Cruise()*np.ones(N), label = 'Cruise')
    plt.xlabel('W/S [lbf/ft^2]')
    plt.ylabel('W/P [lbf/W]')
    plt.legend()
    plt.grid()
    plt.show()

    # Constraint_LDG()
