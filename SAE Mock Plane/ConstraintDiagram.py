import numpy as np
import matplotlib.pyplot as plt
import math
import time
from labellines import labelLines # pip install matplotlib-label-lines
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

        # print(V_stall,res)
        res = [x1[-1][0]-S_LGR,x2[-1][0]-S_LGR]
        V_stall = [V_stall[1] , V_stall[1] - res[1]/((res[1]-res[0])/(V_stall[1]-V_stall[0]))]
    V_stall = V_stall[0]
    return [rho*V_stall**2*C_Lmax/2,V_stall]
    # return S_LGR/80*C_Lmax

# def Constraint_ngTurn(WS,n):
#     C_Lmax = 1
#     rho = 0.0023769 # sea level, slug/ft^3
#     V_stall = np.sqrt(2*WS/rho/C_Lmax)
#     V_LOF = 1.3*V_stall
#     q  = rho/2*V_LOF**2
#     C_Dmin = 0.1
#     ARe = 5
    
#     k = 1/(np.pi*ARe)
#     return q*(C_Dmin/WS+k*(n/q)**2*WS)*550/745.7

def Constraint_MaxWeight():
    P_max = 1000 # our powerplant max power, Watts
    W_max = 55 # MTOW, lb
    return W_max/P_max

def Constraint_Cruise(WS,V_cruise):
    eta_p = 0.8 # prop efficiency at trim
    P_cruise = 0.5 # frac of max power for cruise
    # V_cruise = 30
    # LD = 9 # trim L/D
    C_D0 = 0.02
    rho = 0.0023769 # sea level, slug/ft^3
    q  = rho/2*V_cruise**2
    AR = 10
    e = .8
    return P_cruise/( q*V_cruise*(C_D0+WS**2/(q**2*np.pi*AR*e)) / (550/745.7*eta_p*WS))

if __name__ == "__main__":
    plt.figure(figsize=[8,6])
    N = 101
    WS = np.linspace(1e-3,6,N)
    WP = np.linspace(1e-3,Constraint_TO(WS[0]),N)

    ## Evaluate Constraints
    WP_Takeoff=Constraint_TO(WS)
    LDG_values = Constraint_LDG() # [ WS_Landing , V_stall ]
    # WS_Landing=Constraint_LDG()*np.ones(N)
    WS_Landing=LDG_values[0]*np.ones(N)
    WP_MaxWeight=Constraint_MaxWeight()*np.ones(N)
    # WP_Cruise30=Constraint_Cruise(WS,30)
    # WP_Cruise50=Constraint_Cruise(WS,50)
    V_variation = [10,20,30] # Difference between cruise velocity and stall velocity
    V_Cruise = LDG_values[1]*np.ones(len(V_variation)) + V_variation
    WP_Cruise = np.zeros([len(V_Cruise),len(WS)])
    for i in range(0,len(V_Cruise)):
        WP_Cruise[i,:] = Constraint_Cruise(WS,V_Cruise[i]) # Cruise @ V_stall + 10


    ## Plot constraints
    plt.plot(WS,Constraint_TO(WS), label = '100 ft Takeoff')
    plt.plot(WS_Landing,WP, label = '400 ft Takeoff')
    plt.plot(WS,WP_MaxWeight, label = 'Max Weight (55 lb)')
    # plt.plot(WS,WP_Cruise30, label = '30 ft/s Cruise')
    # plt.plot(WS,WP_Cruise50, label = '50 ft/s Cruise')
    for i in range(0,len(V_Cruise)):
        plt.plot(WS,WP_Cruise[i,:], label = f'V_stall + {V_variation[i]} ft/s Cruise')

    ## Plotting Options
    plt.xlabel('W/S [lbf/ft^2]')
    plt.ylabel('W/P [lbf/W]')
    # plt.legend()
    plt.ylim(0,0.1)
    plt.grid()
    plt.minorticks_on()
    labelLines(plt.gca().get_lines(), zorder=2.5,fontsize = 7)


    # plt.fill_between(
    #     x= WS, 
    #     # y1= np.minimum.reduce([WP_Takeoff,WP_MaxWeight,WP_Cruise30,WP_Cruise50]),
    #     y1= np.minimum.reduce([WP_Takeoff,WP_MaxWeight,WP_Cruise10]),
    #     where= (0 <= WS)&(WS <= Constraint_LDG()),
    #     color= "b",
    #     alpha= 0.2)
    plt.show()

    # Constraint_LDG()
