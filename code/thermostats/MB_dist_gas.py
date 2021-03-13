from matplotlib import pyplot as plt
import numpy as np

def MB_speed(v,m,T):
    """ Maxwell-Boltzmann speed distribution for speeds """
    kB = 1.38e-23
    v=v*100 # convert m/s to A/ps
    return (m/(2*np.pi*kB*T))**1.5 * 4*np.pi * v**2 * np.exp(-m*v**2/(2*kB*T))



fig = plt.figure(figsize=(12.5,4))
fig.tight_layout()
ax = fig.add_subplot(121)

v = np.arange(0,8,0.01)
amu = 1.66e-27
# Krypton
mass = 83.798*amu 

for T in [100,200,300,400]:
    fv = MB_speed(v,mass,T)
    ax.plot(v,fv,label=str(T)+' K',lw=2)

ax.legend(loc=0)
ax.set_xlabel('v, $[ \AA / ps ]$')
ax.set_ylabel('PDF, $f_v(v)$')
ax.set_title("Krypton, 83.798 amu")  

ax2 = fig.add_subplot(122)

v = np.arange(0,25,0.01)
amu = 1.66e-27
T=298.156
gas=['Helium', 'Oxygen', 'Argon', 'Xenon']
masses=[4.0026*amu, 15.9994*amu , 39.948*amu , 131.29*amu]
for i, mass in enumerate(masses):
    fv = MB_speed(v,mass,T)
    ax2.plot(v,fv,label=gas[i]+','+str(mass/amu),lw=2)

ax2.legend(loc=0)
ax2.set_title("T = 298.15")
ax2.set_xlabel('v, $[ \AA / ps ]$')
ax2.set_ylabel('PDF, $f_v(v)$')

plt.show()
plt.savefig("MB_ideal_gas.svg")
