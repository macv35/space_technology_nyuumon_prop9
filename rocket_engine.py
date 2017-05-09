# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import math

conditions = {
        "delta_v" : 7000, #到達デルタV(m/s)
        "g" : 9.8, #重力加速度(m/s^2)
        "m_l" : 2000, #payload(kg)
        "init_a" : 2.94, #初期加速度(m/s^2)
        "rho_s" : 0.3, #標準推進剤混合密度(g/cm^3)
        "f_inert_s" : 0.1, #標準構造質量比
        "T0" : 300, #標準温度(K)
        "eta" : 0.96, #ノズル効率
        "Pj" : 10000, #大気圧(Pa)
        "gamma" : 1.4, #比熱比
        "R0" : 8.3143, #標準気体定数(J/mol*K)
        "t_s" : 0.01, #燃焼室滞在時間(s)
        "Tt_max" : 3000 #ノズルスロート部最大耐熱温度(K)
        }

# reactants = [ [ 燃料量論モル比（化学反応式の係数）, 燃料エンタルピー, 燃料密度, 燃料分子量(kg/mol)], [ 酸化剤 ] ]
reactants = [ [2, 0, 0.07, 0.002], [1, 0, 1.14, 0.032] ]

#reactants = [ [1, 0, 0.07, 0.002], [1, 0, 1.51, 0.038] ]

# products = [ [ 生成物量論モル比（化学反応式の係数）リスト], [生成物エンタルピー(J/mol)リスト ] , [生成物分子量リスト(kg/mol)]]
products = [ [2], [-241900 ], [0.018] ]

#products = [ [2], [-271200 ], [0.02] ]

class Rocket(object):
    """ ２液式ロケットエンジンパラメータ計算クラス """

    def __init__(self, reactants, products, FO_ratio, chamber_pressure, conditions):

        self.reactants = reactants
        # 反応物密度
        self.rho_s = conditions["rho_s"]
        self.rho_fu = reactants[0][2]
        self.rho_ox = reactants[1][2]

        # 反応物生成熱
        self.H_fu = reactants[0][1]
        self.H_ox = reactants[1][1]

        #反応物量論比
        self.fu_ideal_ratio = reactants[0][0]
        self.ox_ideal_ratio = reactants[1][0]
        self.FO_ideal_ratio = self.fu_ideal_ratio/self.ox_ideal_ratio

        #反応物分子量
        self.m_fu = reactants[0][3]
        self.m_ox = reactants[1][3]

        #生成物リスト
        self.prod_ideal_ratio = products[0]
        self.prod_H = products[1]
        self.m_prod = products[2]

        self.FO_ratio = FO_ratio
        self.chamber_pressure = chamber_pressure
        self.conditions = conditions
    

    def delta_Hf(self):
        """ 温度補正なしの生成熱を計算。酸化剤を1molとして計算"""
        prod_H_total = 0

        if self.FO_ratio >= self.FO_ideal_ratio : #不完全燃焼
            for i in range( len(self.prod_ideal_ratio) ):
                prod_H_total += (self.prod_ideal_ratio[i]/self.ox_ideal_ratio) * self.prod_H[i]
        else: #完全燃焼
            for i in range( len(self.prod_ideal_ratio) ):
                prod_H_total += self.FO_ratio * (self.prod_ideal_ratio[i] / self.fu_ideal_ratio) * self.prod_H[i]

        return self.FO_ratio*self.H_fu + self.H_ox - prod_H_total


    def total_mol_after_reaction(self):
        """ 反応後の総モル数を返す。反応前の酸化剤を1molとする。"""
        total_prod_mol = 0

        if self.FO_ratio >= self.FO_ideal_ratio :
            for i in range( len(self.prod_ideal_ratio) ):
                total_prod_mol += self.prod_ideal_ratio[i]/self.ox_ideal_ratio
            return (self.FO_ratio - self.FO_ideal_ratio) + total_prod_mol
        else:
            for i in range( len(self.prod_ideal_ratio) ):
                total_prod_mol += self.FO_ratio * (self.prod_ideal_ratio[i] / self.fu_ideal_ratio)
            return (1- (self.FO_ratio * self.ox_ideal_ratio) / self.fu_ideal_ratio) + total_prod_mol


    def Tf(self):
        """ Tf の値を計算して返す """

        Cp = self.conditions["gamma"] * self.conditions["R0"] / (conditions["gamma"]-1)
        total_mol_after_reaction = self.total_mol_after_reaction()
        delta_Hf = self.delta_Hf()

        #2500K以下だと仮定
        Tf = delta_Hf / (Cp * total_mol_after_reaction) + self.conditions["T0"]
        if Tf <= 2500 :
            return Tf
        else: #高温領域で生成熱が減少することの補正
            return (Cp * total_mol_after_reaction + 2*delta_Hf) / (Cp * total_mol_after_reaction + 0.0004*delta_Hf)


    def f_inert(self):
        """ 構造質量比を計算して返す """
        mixture_ratio = self.m_ox / (self.m_fu*self.FO_ratio)
        rho = (mixture_ratio + 1) / (1/self.rho_fu + mixture_ratio/self.rho_ox)
        return 1/( (1/self.conditions["f_inert_s"] -1)*rho/self.conditions["rho_s"] + 1)


    def m_average_after_reaction(self):
        """ 反応後の平均分子量(kg/mol)を返す """
        total_tmp = 0
        if self.FO_ratio >= self.FO_ideal_ratio :
            for i in range( len(self.prod_ideal_ratio) ):
                total_tmp += (self.prod_ideal_ratio[i]/self.ox_ideal_ratio) * self.m_prod[i]
            return ((self.FO_ratio-self.FO_ideal_ratio) * self.m_fu + total_tmp) / self.total_mol_after_reaction()
        else:
            for i in range( len(self.prod_ideal_ratio) ):
                total_tmp += (self.FO_ratio * (self.prod_ideal_ratio[i] / self.fu_ideal_ratio))*self.m_prod[i]
            return ((1- self.FO_ratio*(self.ox_ideal_ratio/self.fu_ideal_ratio)) * self.m_ox + total_tmp) / self.total_mol_after_reaction()
        

    def V_j_infinity(self):
        """ 燃焼室による噴出速度損失を考慮しない排気速度 """
        Cp = self.conditions["gamma"] * self.conditions["R0"] / ((self.conditions["gamma"]-1) * self.m_average_after_reaction() )
        p_ratio_tmp = 1-math.pow(self.conditions["Pj"]/self.chamber_pressure, (self.conditions["gamma"]-1)/self.conditions["gamma"])
        return math.sqrt(2*self.conditions["eta"]*Cp*self.Tf()*p_ratio_tmp)


    def V_j(self):
        """燃焼室による噴出速度損失を考慮した排気速度 """
        V_j_infinity = self.V_j_infinity()


#
## TODO
#
        return V_j_infinity


    def Isp(self):
        return self.V_j()/self.conditions["g"]


    def payload_lambda(self):
        """ ペイロード比を計算 """
        f_inert = self.f_inert()
        return (math.exp(-1*self.conditions["delta_v"]/(self.conditions["g"]*self.Isp() ))-f_inert) / (1-f_inert)



if __name__ == "__main__":
    ratio_array = np.arange(0,6.0,0.01)
    payload_lambda_array = np.zeros(600)
    for i in np.arange(600):
        rocket = Rocket(reactants, products, ratio_array[i], 20000000, conditions)
        payload_lambda_array[i] = rocket.payload_lambda()

    plt.plot(ratio_array, payload_lambda_array)
    plt.show()
