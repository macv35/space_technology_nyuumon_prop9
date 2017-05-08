# vim:fileencoding=utf-8
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

conditions = {
        "delta_v" : 7000,
        "g" : 9.8,
        "ml" : 2000, #payload
        "init_a" : 2.94,
        "rho_s" : 0.3,
        "rho_fu" : 0.07,
        "rho_ox" : 1.14,
        "f_inert_s" : 0.1,
        "H_fuel" : 0,
        "H_ox" : 0,
        "H_prod" : -241900,
        "T0" : 300,
        "yiita" : 0.96,
        "P0" : 1020000,
        "Pj" : 10000,
        "gamma" : 1.4,
        "R0" : 8314.3,
        "Cp" : 29.1
        }
class Rocket(object):
    def __init__(self, rho_fu, rho_ox, H_fu, H_ox, H_prod, ratio):
        self.rho_fu = rho_fu
                self.rho_ox = rho_ox
                self.H_fu = H_fu
                self.H_ox = H_ox
                self.H_prod = H_prod
                self.ratio = ratio
                self.delta_v = conditions["delta_v"]
                self.g = conditions["g"]
                self.ml = conditions["ml"]
                self.Cp = conditions["Cp"]
                self.init_a = conditions["init_a"]
                self.rho_s = conditions["rho_s"]
                self.f_inert_s = conditions["f_inert_s"]
                self.T0 = conditions["T0"]
                self.yiita = conditions["yiita"]
                self.P0 = conditions["P0"]
                self.Pj = conditions["Pj"]
                self.gamma = conditions["gamma"]
                self.R0 = conditions["R0"]

        def tf_equation(self, Q1, FO_ratio):
            if FO_ratio==self.ratio:
                if Q1 <= 64073.835*FO_ratio:
                    return Q1/(self.Cp*FO_ratio) + self.T0
                else:
                    return (2*Q1 + FO_ratio*self.Cp*self.T0)/(FO_ratio*self.Cp + 0.0004*Q1)
                elif FO_ratio > self.ratio:
                    if Q1 <= 64073.835 * FO_ratio:
                        return Q1/(self.Cp*FO_ratio) + self.T0
                    else:
                        return (2 * Q1 + FO_ratio * self.Cp * self.T0) / (FO_ratio * self.Cp + 0.0004 * Q1)
                elif FO_ratio < self.ratio:
                    if Q1 <= 64073.835*(1 + 0.5*FO_ratio):
                        return Q1/(self.Cp*(1 + 0.5*FO_ratio)) + self.T0
                    else:
                        return (2*Q1 + (1 + 0.5*FO_ratio)*self.Cp*self.T0)/((1 + 0.5*FO_ratio)*self.Cp + 0.0004*Q1)

        def Q1_calculate(self, fuel_ratio):
            return fuel_ratio*self.H_fu + self.H_ox - fuel_ratio*self.H_prod


if __name__ == "__main__":
    rocket = Rocket(0.07, 1.14, 0, 0, -241900, 2)
        ratio_array = np.arange(0,6.0,0.01)
        Q1_array = rocket.Q1_calculate(ratio_array)
        Tf_array = np.zeros(600)
        for i in np.arange(600):
            Tf_array = rocket.tf_equation(Q1_array[i], ratio_array[i])

        M_array = (ratio_array*2 + 16)/(ratio_array+1)
        R_array = rocket.R0/M_array
        Cp_array = (1.4/0.4)*R_array
        Vj_array = (2*rocket.yiita*Cp_array*Tf_array*(1-(rocket.Pj/rocket.P0) ** 0.2857))**0.5
        Isp_array = Vj_array/rocket.g
        rho_array = (ratio_array*rocket.rho_fu + rocket.rho_ox)/(ratio_array + 1)
        f_inert_array = 1/((1/rocket.f_inert_s-1)*(rho_array/rocket.rho_s)+1)
        payload_ratio = (np.exp(-rocket.delta_v/Vj_array)-f_inert_array)/(1-f_inert_array)
        print(Vj_array)
        plt.plot(ratio_array, payload_ratio)
        plt.show()
