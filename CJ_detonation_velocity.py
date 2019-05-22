"""
Shock and Detonation Toolbox Program

Calculates the CJ speed with various initial conditions 
using the Minimum Wave Speed Method and 
then finds the equilibrium state of the gas behind a shock wave 
traveling at the CJ speed.

"""
import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sdtoolbox.postshock import CJspeed
from sdtoolbox.postshock import PostShock_eq
from sdtoolbox.thermo import soundspeed_eq, soundspeed_fr

#Setting plotting defaults
plt.style.use('ggplot')
plt.style.use('seaborn-pastel')

# Initial state specification:
# P1 = Initial Pressure  
# T1 = Initial Temperature 
# U = Shock Speed 
# q = Initial Composition 
# mech = Cantera mechanism File name

# stoichiometric C2H6-air detonation at nominal atmospheric conditions
P1 = 101325.
T1 = 295
q ='C2H6:0.286 O2:1 N2:3.76'
mech = 'gri30_highT.cti'
gas_initial = ct.Solution(mech)
gas_initial.TPX = T1, P1, q
rho_1 = gas_initial.density
#print (gas_initial())
# compute CJ speed
[cj_speed,R2,plot_data] = CJspeed(P1, T1, q, mech, fullOutput=True)  

# compute equilibrium CJ state parameters
gas = PostShock_eq(cj_speed, P1, T1, q, mech)
ae = soundspeed_eq(gas)
af = soundspeed_fr(gas)
rho_2 = gas.density
gammae = ae**2*rho_2/gas.P
gammaf = af**2*rho_2/gas.P
w2 = cj_speed*rho_1/rho_2
u2 = cj_speed-w2
print ('CJ computation for ' + mech + ' with composition ' + q )
print ('Initial conditions: P1 = %.3e Pa & T1 = %.2f K'  % (P1,T1)  )
print ('CJ Speed   %.1f m/s' % cj_speed)
print ('CJ State')
print ('   Pressure   %.3e Pa' % gas.P)
print ('   Temperature  %.1f K' % gas.T)
print ('   Density  %.3f kg/m3' % gas.density)
print ('   Entropy  %.3e J/K' % gas.entropy_mass)
print ('   w2 (wave frame) %.1f m/s' % w2)
print ('   u2 (lab frame) %.1f m/s' % u2)
print ('   c2 frozen %.1f m/s' % af)
print ('   c2 equilbrium %.1f m/s' % ae)
print ('   gamma2 frozen %.3f ' % gammaf)
print ('   gamma2 equilbrium %.3f ' % gammae)
#print (gas.X)

"""

Calculations for different initial conditions - temperature

"""

"""
#Making a list of the temperatures to run simulations at
T = []
counter = 0
while counter<20:
    T.append(300+counter*100)
    counter+=1

#Creating a dataFrame in wchich will be stored results
CJspeeds = pd.DataFrame(data={'T1': T})
CJspeeds['CJspeed'] = np.nan
CJspeeds['pressure_ratio'] = np.nan
CJspeeds['density_ratio'] = np.nan

for i, temperature in enumerate(T):
    #Setting the gas
    T1 = temperature
    P1 = 101325.0
    mech = "gri30_highT.cti"
    q ='C2H6:0.286 O2:1 N2:3.76'
    gas_initial = ct.Solution(mech)
    gas_initial.TPX = T1, P1, q
    rho_1 = gas_initial.density
    
    # compute CJ speed
    cj_speed = CJspeed(P1, T1, q, mech)  
    
    # compute equilibrium CJ state parameters
    gas = PostShock_eq(cj_speed, P1, T1, q, mech)
    rho_2 = gas.density
    P2 = gas.P       
    
    print('Computed CJ speed: {:.1f} m/s for T={}K.'.format(cj_speed, temperature))
    
    CJspeeds.at[i, 'CJspeed'] = cj_speed
    CJspeeds.at[i, 'pressure_ratio'] = P2/P1
    CJspeeds.at[i, 'density_ratio'] = rho_2/rho_1
    
#Plotting results- CJ speed
plt.plot(CJspeeds['T1'], CJspeeds['CJspeed'], 'o-')
plt.ylabel('CJ speed (m/s)')
plt.xlabel(r'Temperature: $T(K)$')
plt.savefig('T_CJspeed.png')
plt.show()

#Plotting other properties
plt.plot(CJspeeds['T1'], CJspeeds['pressure_ratio'], 'o-')
plt.ylabel('Pressure ratio: p2/p1(-)')
plt.xlabel(r'Temperature: $T(K)$')
plt.savefig('T_pressure.png')
plt.show()

plt.plot(CJspeeds['T1'], CJspeeds['density_ratio'], 'o-')
plt.ylabel('Density ratio: rho_2/rho_1(-)')
plt.xlabel(r'Temperature: $T(K)$')
plt.savefig('T_rho.png')
plt.show()
"""

"""

Calculations for different initial conditions - pressure

"""

"""
#Making a list of the pressures to run simulations at
P = []
counter = 0
while counter<30:
    P.append(10000+5000*counter**2.5)
    counter+=1


#Creating a dataFrame in wchich will be stored results
CJspeeds = pd.DataFrame(data={'P1': P})
CJspeeds['CJspeed'] = np.nan
CJspeeds['pressure_ratio'] = np.nan
CJspeeds['density_ratio'] = np.nan

for i, pressure in enumerate(P):
    #Setting the gas
    T1 = 295
    P1 = pressure
    mech = "gri30_highT.cti"
    q ='C2H6:0.286 O2:1 N2:3.76'
    gas_initial = ct.Solution(mech)
    gas_initial.TPX = T1, P1, q
    rho_1 = gas_initial.density
    
    # compute CJ speed
    cj_speed = CJspeed(P1, T1, q, mech)  
    
    # compute equilibrium CJ state parameters
    gas = PostShock_eq(cj_speed, P1, T1, q, mech)
    rho_2 = gas.density
    P2 = gas.P
    
    print('Computed CJ speed: {:.1f} m/s for P={:3.2f}bar.'.format(cj_speed, pressure/100000))
    
    CJspeeds.at[i, 'CJspeed'] = cj_speed
    CJspeeds.at[i, 'pressure_ratio'] = P2/P1
    CJspeeds.at[i, 'density_ratio'] = rho_2/rho_1
    
#Plotting CJ speed   
fig = plt.figure()
ax = fig.add_subplot(111)
ax.semilogx(CJspeeds['P1']/100000, CJspeeds['CJspeed'], 'o-')
ax.set_ylabel('CJ speed (m/s)')
ax.set_xlabel(r'Pressure: $P(bar)$');
plt.savefig('P_CJspeed.png')
plt.show()

#Plotting other properties
fig = plt.figure()
ax = fig.add_subplot(111)
ax.semilogx(CJspeeds['P1']/100000, CJspeeds['pressure_ratio'], 'o-')
plt.ylabel('Pressure ratio: p2/p1(-)')
plt.xlabel(r'Pressure: $P(bar)$')
plt.savefig('P_pressure.png')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.semilogx(CJspeeds['P1']/100000, CJspeeds['density_ratio'], 'o-')
plt.ylabel('Density ratio: rho_2/rho_1(-)')
plt.xlabel(r'Pressure: $P(bar)$')
plt.savefig('P_rho.png')
plt.show()
"""

"""

Calculations for different initial conditions - stoichiometric coefficient

"""

"""
#Making a list of the phis to run simulations at
phis = []
counter = 0
while counter<20:
    phis.append(0.515+counter*0.1)
    counter+=1

#Creating a dataFrame in wchich will be stored results
CJspeeds = pd.DataFrame(data={'phi': phis})
CJspeeds['CJspeed'] = np.nan
CJspeeds['pressure_ratio'] = np.nan
CJspeeds['density_ratio'] = np.nan

for i, stcoef in enumerate(phis):
    #Setting the gas
    T1 = 295
    P1 = 101325.0
    mech = "gri30_highT.cti"
    q ='C2H6:'+str(0.286*stcoef)+' O2:1 N2:3.76'
    gas_initial = ct.Solution(mech)
    gas_initial.TPX = T1, P1, q
    rho_1 = gas_initial.density
    
    # compute CJ speed
    cj_speed = CJspeed(P1, T1, q, mech)  
    
    # compute equilibrium CJ state parameters
    gas = PostShock_eq(cj_speed, P1, T1, q, mech)
    rho_2 = gas.density
    P2 = gas.P
       
    print('Computed CJ speed: {:.1f} m/s for phi={:.3f}.'.format(cj_speed, stcoef))
    
    CJspeeds.at[i, 'CJspeed'] = cj_speed
    CJspeeds.at[i, 'pressure_ratio'] = P2/P1
    CJspeeds.at[i, 'density_ratio'] = rho_2/rho_1
    
#Plotting results - CJ speed
plt.plot(CJspeeds['phi'], CJspeeds['CJspeed'], 'o-')
plt.ylabel('CJ speed (m/s)')
plt.xlabel(r'Stoichiometric coefficient: $\Phi(-)$')
plt.xlim([0.515, 2.358])
plt.savefig('phi_CJspeed.png')
plt.show()

#Plotting other properties
plt.plot(CJspeeds['phi'], CJspeeds['pressure_ratio'], 'o-')
plt.ylabel('Pressure ratio: p2/p1(-)')
plt.xlabel(r'Stoichiometric coefficient: $\Phi(-)$')
plt.xlim([0.515, 2.358])
plt.savefig('phi_pressure.png')
plt.show()


plt.plot(CJspeeds['phi'], CJspeeds['density_ratio'], 'o-')
plt.ylabel('Density ratio: rho_2/rho_1(-)')
plt.xlabel(r'Stoichiometric coefficient: $\Phi(-)$')
plt.xlim([0.515, 2.358])
plt.savefig('phi_rho.png')
plt.show()
"""