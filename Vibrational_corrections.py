from math import pi, log, exp, sqrt
from scipy.constants import speed_of_light, gas_constant, Avogadro, atm, Boltzmann, Planck
from scipy.constants import physical_constants
import numpy as np
speed_of_light=speed_of_light*100
h_to_ev = physical_constants["Hartree energy in eV"][0]
kb = physical_constants["Boltzmann constant in eV/K"][0]
hbar= physical_constants["Planck constant in eV/Hz"][0]
R = gas_constant*physical_constants["joule-electron volt relationship"][0]/(Avogadro)
T = 298.15
P= 10**5
R2= gas_constant
kb2= Boltzmann
h=Planck
sim_nu=2

def ZPE(omega: list):
    sum_start = 0.0
    for o in omega:
        sum_start += hbar*o*speed_of_light
    return 0.5*sum_start
    
def Svib(omega: list):
    sum_start=0.0
    for o in omega:
        sum_start += (hbar*o*speed_of_light)/(kb*T*(exp((hbar*o*speed_of_light)/(kb*T))-1))-log(1-exp((-hbar*o*speed_of_light)/(kb*T)))
    return R*sum_start
        
def deltaU(omega: list):
    sum_start=0.0
    for o in omega:
        sum_start += (hbar*o*speed_of_light)/(kb*(exp((hbar*o*speed_of_light)/(kb*T))-1))
    return R*sum_start

#def Strans(mass: float):
    core=(R2*log((((2.0*pi*mass*kb2*T)/(h**2))**(1.5))*((kb2*T)/(P)))+2.5*R2)
    return(core/96485)
    
#def Srot(inertia: float):
    core=(R2*log((8*pi**2*kb2*T*inertia)/(h**2*sim_nu)))+R2
    return(core/96485)

#def Srot2(inertia: float):
    #core=R2*log(((8*pi**2*kb2*T)/(h**2))**1.5*(sqrt(pi*inertia))/(sim_nu))+ 1.5*R2
    #return(core/96485)

def mu(EDFT: float, omega: list):
    #mass=(2.6567*(1/(10**26)))
    #inertia=(1.94*(1/10**46))
    corrections = ZPE(omega) - T*Svib(omega) + deltaU(omega) #- T*Strans(mass) - T*Srot(inertia) + 2.5*kb*T
   
    print(ZPE(omega))
    print(-T*Svib(omega))
    print(deltaU(omega))
    print(f"Total corrections: {corrections}")
    
    return (EDFT*h_to_ev) + corrections
 

#example
    
if __name__ == "__main__":
    #print(Strans(2.99*(1/(10**26)))*T)
    #print(T*Svib([1648.981016, 3859.398246, 3962.818706]))
    #print(Srot(5.733*(1/10**141))*T)
    #print(ZPE([1648.981016, 3859.398246, 3962.818706]))
    #print(deltaU([1648.981016, 3859.398246, 3962.818706]))
    print(mu(-5596.45436245219, [266.273123,   212.508417, 250.870699, 383.305112,  431.045829, 1273.981734, 1584.100398,  1830.148714, 2367.780582]))
    #print(3*kb*T)

#mass=3.3344*(1/(10**27))
#print(mass)
#core=(R2*log((((2.0*pi*mass*kb2*T)/(h**2))**(1.5))*((kb2*T)/(P)))+2.5*R2)
#print(core)
# core2=(core/96485)*T
# print(core2)

# H2   4.64 × 10−48 kg m2
  # O2 1.94 × 10-46 kg•m2
  # H2O 1.275365 amu^3 * Angstrom^6-> 5.833*10**-141 kg m2
  # symmetry number linear 2
  # symmetry number h2o 2
  