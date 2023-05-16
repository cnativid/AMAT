import numpy as np
import matplotlib.pyplot as plt
from ambiance import Atmosphere

g = 32.17 # ft/s^2


def Constraint_TO(WS):
    eta_p = 0.45 # Gudmunsson, 
    S_G = 100 # takeoff ground run, ft
    C_LTO =  0.6
    C_DTO =  0.02 
    C_Lmax = 1.5
    mu = 0.3
    rho = 0.0023769 # sea level, slug/ft^3
    V_LOF = 1.1*np.sqrt(2*WS/(rho*C_Lmax))
    return 1/((1.21/(g*rho*C_Lmax*S_G)*WS+0.605/C_Lmax*(C_DTO-mu*C_LTO) + mu)*V_LOF/np.sqrt(2)/(eta_p*550/745.7))

def Constraint_LDG():
    S_LGR = 400 # ft
    mu = 0.3
    rho = 0.0023769 # sea level, slug/ft^3

    return C_Lmax/2*rho*(2*g*mu*S_LGR)/1.3**2

if __name__ == "__main__":
    plt.figure(figsize=[5,6])
    N = 101
    WS = np.linspace(5,40,N)
    WP = np.linspace(0,.03,N)
    plt.plot(WS,Constraint_TO(WS))
    plt.plot(Constraint_LDG()*np.ones(N),WP)
    plt.xlabel('W/S [lbf/ft^2]')
    plt.ylabel('W/P [lbf/W]')
    plt.show()
