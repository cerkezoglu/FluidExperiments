import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import sqrt, pow


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
        self.vel = np.round(self.calculate_vel(), 6)
        self.Re, self.theoric_dirac = self.calculate_dirac()
        self.exp_dirac = self.find_dirac()
        print("Reynolds: ", self.Re,"Theoric dirac: ", self.theoric_dirac)
        print("Experimental dirac: ", self.exp_dirac)

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

    def calculate_dirac(self):
        k_vis = 0.0000153 # interpolating from Munson @22 Celcius
        Re = self.u*self.x*1e-2/ k_vis
        if Re < 1.6e5:
            dirac = self.x * 5e1 / sqrt(Re)  # dirac mm
            print("Flow Laminar")
        elif Re > 2.5e5:
            dirac = 0.37*1e1/pow(Re, 1/5)*self.x
            print("Flow Turbulent")
        else:
            dirac = 1
        return Re, np.round(dirac, 3)

    def find_dirac(self):
        u_star = self.u*.99
        for i in range(len(self.vel)):
            if self.vel[i] > u_star:
                dirac = (u_star-self.vel[i-1])* (self.y[i]-self.y[i-1])/(self.vel[i]-self.vel[i-1])+ self.y[i-1]
                return np.round(dirac, 3)


def plot_graph(a, b):
    plt.figure(figsize=(16,9))
    plt.plot(b.vel / b.u, b.y / b.exp_dirac, linestyle='-.', label="turbulent", marker='^')
    plt.scatter(b.vel / b.u, b.y / b.exp_dirac, marker='^')
    plt.plot(a.vel / a.u, a.y / a.exp_dirac,linestyle='-.', label="laminar", marker='s')
    plt.scatter(a.vel / a.u, a.y / a.exp_dirac, marker='s')
    plt.title("Velocity profile @ x = "+str(a.x), fontsize=20)
    plt.xlabel(r"$\frac{u}{U_\infty}$", fontsize=20)
    plt.ylabel(r"$\frac{y}{\delta}$", fontsize=20)
    plt.legend(loc = 'upper left')
    plt.grid()
    plt.show()
    # plt.savefig('x'+str(a.x)+'.jpg', format='jpg')
    # plt.close()


b = CalculateVelocities(25, 17)
a = CalculateVelocities(8, 17)
c = CalculateVelocities(25, 25)
d = CalculateVelocities(8, 25)
plot_graph(a, b)
plot_graph(d, c)


