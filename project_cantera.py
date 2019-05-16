import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

#Setting plotting defaults
plt.style.use('ggplot')
plt.style.use('seaborn-pastel')

#Defining the gas
gas = ct.Solution('gri30.cti')

#Defining reactor temperature and pressure
initial_state = 1200, 101325
gas.TP = initial_state

#Defining the fuel, oxidizer and set the stoichiometry
gas.set_equivalence_ratio(phi=1.0, fuel='C2H6', oxidizer={'o2':1.0, 'n2':3.76})

#Creating a batch reactor object and adding it to a reactor network
r = ct.IdealGasConstPressureReactor(gas)
reactorNetwork = ct.ReactorNet([r])

#Compiling a list of all variables for which data will be stored
stateVariableNames = [r.component_name(item) for item in range(r.n_vars)]

#Creating a DataFrame using the above list
timeHistory = pd.DataFrame(columns=stateVariableNames)

#Creating function which computes the ignition delay from the occurence
#of the peak in species' concentration.
def ignitionDelay(df, species):
    return df[species].idxmax()

estimatedIgnitionDelayTime = 1
t = 0
counter = 0
while(t < estimatedIgnitionDelayTime):
    t = reactorNetwork.step()
    #Saving only every 10th value to save time
    if (counter%10 == 0):
        timeHistory.loc[t] = reactorNetwork.get_state()
    counter+=1

#Using the 'oh' species to compute the ignition delay
tau = ignitionDelay(timeHistory, 'OH')

print('Computed Ignition Delay: {:.3e} seconds'.format(tau))

#Plotting results
plt.plot(timeHistory.index*1000, timeHistory['OH'], '-o')
plt.xlabel('Time (ms)')
plt.ylabel('$Y_{OH}$')
plt.xlim([0,10])
plt.savefig('Y_OH.png')
plt.show()

plt.plot(timeHistory.index*1000, timeHistory['O2'], '-o')
plt.plot(timeHistory.index*1000, timeHistory['C2H6'], '-o')
#plt.plot(timeHistory.index*1000, timeHistory['H2O'], '-o')
plt.plot(timeHistory.index*1000, timeHistory['CO2'], '-o')
plt.xlabel('Time (ms)')
plt.ylabel('Mole Fraction')
plt.xlim([0,10])
plt.legend(('$O2$', '$C2H6$', '$CO2$'))
plt.savefig('Y_O2.png')
plt.show()

plt.plot(timeHistory.index*1000, timeHistory['temperature'], '-o')
plt.xlabel('Time (ms)')
plt.ylabel('Temperature (K)')
plt.xlim([0,10])
plt.savefig('T.png')
plt.show()


"""

Calculations for different initial conditions - temperature

"""


"""
#Making a list of the temperatures to run simulations at
T = []
counter = 0
while counter<45:
    T.append(900+counter*25)
    counter+=1

#Creating a dataFrame in wchich will be stored results
ignitionDelays = pd.DataFrame(data={'T': T})
ignitionDelays['ignDelay'] = np.nan

for i, temperature in enumerate(T):
    #Setting the gas and reactor
    reactorTemperature = temperature
    reactorPressure = 101325.0
    gas.TP = reactorTemperature, reactorPressure
    gas.set_equivalence_ratio(phi=1.0, fuel='C2H6', oxidizer={'o2':1.0, 'n2':3.76})
    r = ct.IdealGasConstPressureReactor(gas)
    reactorNetwork = ct.ReactorNet([r])
    
    #Creating an empty data frame
    timeHistory = pd.DataFrame(columns=timeHistory.columns)
    
    t = 0
    counter = 0
    estimatedIgnitionDelayTime = 25
    while t < estimatedIgnitionDelayTime:
        t = reactorNetwork.step()
        if not counter % 20:
            timeHistory.loc[t] = r.get_state()
        counter += 1
    
    tau = ignitionDelay(timeHistory, 'OH')
    
    print('Computed Ignition Delay: {:.3e} seconds for T={}K.'.format(tau, temperature))
    
    ignitionDelays.at[i, 'ignDelay'] = tau
    
fig = plt.figure()
ax = fig.add_subplot(111)
ax.semilogy(ignitionDelays['T'], ignitionDelays['ignDelay'], 'o-')
ax.set_ylabel('Ignition Delay (s)')
ax.set_xlabel(r'Temperature: $T(K)$');
plt.savefig('tau_T.png')
"""


"""

Calculations for different initial conditions - pressure

"""


"""
#Making a list of the pressures to run simulations at
P = []
counter = 0
while counter<45:
    P.append(1000+pow(counter, 2.5)*500)
    counter+=1

#Creating a dataFrame in wchich will be stored results
ignitionDelays = pd.DataFrame(data={'P': P})
ignitionDelays['ignDelay'] = np.nan

for i, pressure in enumerate(P):
    #Setting the gas and reactor
    reactorTemperature = 1200
    reactorPressure = pressure
    gas.TP = reactorTemperature, reactorPressure
    gas.set_equivalence_ratio(phi=1.0, fuel='C2H6', oxidizer={'o2':1.0, 'n2':3.76})
    r = ct.IdealGasConstPressureReactor(gas)
    reactorNetwork = ct.ReactorNet([r])
    
    #Creating an empty data frame
    timeHistory = pd.DataFrame(columns=timeHistory.columns)
    
    t = 0
    counter = 0
    estimatedIgnitionDelayTime = 1
    while t < estimatedIgnitionDelayTime:
        t = reactorNetwork.step()
        if not counter % 20:
            timeHistory.loc[t] = r.get_state()
        counter += 1
    
    tau = ignitionDelay(timeHistory, 'OH')
    
    print('Computed Ignition Delay: {:.3e} seconds for P={}bar.'.format(tau, pressure/100000))
    
    ignitionDelays.at[i, 'ignDelay'] = tau
    
fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(ignitionDelays['P']/100000, ignitionDelays['ignDelay'], 'o-')
ax.set_ylabel('Ignition Delay (s)')
ax.set_xlabel(r'Pressure: $P(bar)$');
plt.savefig('tau_P.png')

"""


"""

Calculations for different initial conditions - stoichiometric coefficient

"""


"""
#Making a list of the phis to run simulations at
phis = []
counter = 0
while counter<40:
    phis.append(0.515+counter*0.05)
    counter+=1

#Creating a dataFrame in wchich will be stored results
ignitionDelays = pd.DataFrame(data={'phi': phis})
ignitionDelays['ignDelay'] = np.nan

for i, stcoef in enumerate(phis):
    #Setting the gas and reactor
    reactorTemperature = 1200
    reactorPressure = 101325.0
    gas.TP = reactorTemperature, reactorPressure
    gas.set_equivalence_ratio(phi=stcoef, fuel='C2H6', oxidizer={'o2':1.0, 'n2':3.76})
    r = ct.IdealGasConstPressureReactor(gas)
    reactorNetwork = ct.ReactorNet([r])
    
    #Creating an empty data frame
    timeHistory = pd.DataFrame(columns=timeHistory.columns)
    
    t = 0
    counter = 0
    estimatedIgnitionDelayTime = 1000
    while t < estimatedIgnitionDelayTime:
        t = reactorNetwork.step()
        if not counter % 20:
            timeHistory.loc[t] = r.get_state()
        counter += 1
    
    tau = ignitionDelay(timeHistory, 'OH')
    
    print('Computed Ignition Delay: {:.3e} seconds for phi={}.'.format(tau, stcoef))
    
    ignitionDelays.at[i, 'ignDelay'] = tau
    

fig = plt.figure()
ax = fig.add_subplot(111)
ax.semilogx(ignitionDelays['phi'], ignitionDelays['ignDelay']*1000, 'o-')
ax.set_ylabel('Ignition Delay (ms)')
ax.set_xlabel(r'Stoichiometric coefficient: $\Phi(-)$');
plt.xlim([0.515, 2.358])
plt.savefig('tau_phi.png')
"""