from math import pi, log, exp
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


def mu(EDFT: float, omega: list):
    corrections = ZPE(omega) - T*Svib(omega) + deltaU(omega)
    
    print(ZPE(omega))
    print(-T*Svib(omega))
    print(deltaU(omega))
    print(f"Total corrections: {corrections}")
    
    return (EDFT*h_to_ev) + corrections
 
def Strans(mass: float):
    core=(R2*log((((2.0*pi*mass*kb2*T)/(h**2))**(1.5))*((kb2*T)/(P)))+2.5*R2)
    return(core/96485)
    
def Srot(inertia: float):
    core=(R2*log((8*pi**2*kb2*T*inertia)/(h**2*sim_nu)))+R2
    return(core/96485)
 
#example
    
if __name__ == "__main__":
    print(Strans(3.3344*(1/(10**27)))*T)
    print(T*Svib([4732.38]))
    print(Srot(4.64*(1/10**48))*T)
    print(ZPE([4732.38]))
    print(deltaU([4732.38]))

#mass=3.3344*(1/(10**27))
#print(mass)
#core=(R2*log((((2.0*pi*mass*kb2*T)/(h**2))**(1.5))*((kb2*T)/(P)))+2.5*R2)
#print(core)
# core2=(core/96485)*T
# print(core2)

# H2   4.64 × 10−48 kg m2
  # O2 1.94 × 10-46 kg•m2
  # H2O 1.275365 amu^3 * Angstrom^6-> 2.11*10**-123 kg m2
  # symmetry number linear 2
  # symmetry number h2o 2
  
  