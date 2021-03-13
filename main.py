import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import sqrt, pi, exp


# #########
# Ambient pressure (P0) = 1 atm
# Ambient temperature (T0) = 22 deg C
# Linear calibration equation of the pressure transducer: delP = 25.407*V - 1.245, delP: Pressure difference in terms of mmH20, V: Voltage (V)
# #########

class CalculateVelocities():
    def __init__(self,u, x):
        self.u = u #freestream velocity m/s
        self.x = x #location from the LE direction x in cm
        self.y, self.V = self.data_read()
        self.del_pa = self.calibration_eq()
        self.rho = self.calculate_rho()
        self.vel = self.calculate_vel()
        print(self.y)
        print(self.vel)

    def data_read(self):
        dataset = pd.read_table('U'+str(self.u)+'_x'+str(self.x)+'.txt', sep='\s+')
        y = np.array(dataset[dataset.columns[1]]) # mm
        V = np.array(dataset[dataset.columns[2]]) # V
        return y, V

    def calibration_eq(self):
        delp_mmh2o = 25.407 * self.V - 1.245  # Calibration eq to mmH2o
        delp_pa = delp_mmh2o * 9.807  # Pressure diff in Pa
        return np.round(delp_pa, 6)

    def calculate_rho(self):
        T = 273.15 + 22 # Temprature in Kelvin
        R = 287 # J/kg-K
        rho = (self.del_pa+101325)/(R*T)
        return rho

    def calculate_vel(self):
        vel = []
        for i in range(len(self.rho)):
            vel.append(sqrt(2*self.del_pa[i]/self.rho[i] ))
        return np.array(vel)


b = CalculateVelocities(8, 17)
a = CalculateVelocities(8, 25)
plt.plot(b.vel, b.y)
plt.plot(a.vel, a.y)
plt.show()