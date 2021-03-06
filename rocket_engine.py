# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import newton
import matplotlib.pyplot as plt


class Rocket(object):
    """ ２液式ロケットエンジンパラメータ計算クラス """

    def __init__(self, reactants, products, FO_ratio, chamber_pressure, conditions):
        """
        コンストラクタ。引数は，
        reactants = [ [ 燃料量論モル比（化学反応式の係数）, 燃料エンタルピー, 燃料密度(g/cm^3), 燃料分子量(kg/mol)], [ 酸化剤についても同様に ] ]
        液酸液水のときの例：reactants = [ [2, 0, 0.07, 0.002], [1, 0, 1.14, 0.032] ]
        products = [ [ 生成物量論モル比（化学反応式の係数）リスト], [生成物エンタルピー(J/mol)リスト ] , [生成物分子量リスト(kg/mol)]]
        液酸液水のときの例：products = [ [2], [-241900 ], [0.018] ]
        FO_ratio = 酸素を1molとした時の燃料のmol比率
        chamber_pressure = 燃焼室圧力(Pa)
        conditions = 初期条件のリスト
        例：conditions = {
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
        """
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

        #初期条件
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

        Cp = self.conditions["gamma"] * self.conditions["R0"] / (self.conditions["gamma"]-1)
        total_mol_after_reaction = self.total_mol_after_reaction()
        delta_Hf = self.delta_Hf()

        #2500K以下だと仮定
        Tf = delta_Hf / (Cp * total_mol_after_reaction) + self.conditions["T0"]
        if Tf <= 2500 :
            return Tf
        else: #高温領域で生成熱が減少することの補正
            return (Cp * total_mol_after_reaction * self.conditions["T0"] + 2*delta_Hf) / (Cp * total_mol_after_reaction + 0.0004*delta_Hf)


    def f_inert(self):
        """ 構造質量比を計算して返す """
        mixture_ratio = self.mixture_ratio()
        rho = (mixture_ratio + 1) / (1/self.rho_fu + mixture_ratio/self.rho_ox)
        return 1/( (1/self.conditions["f_inert_s"] -1)*rho/self.conditions["rho_s"] + 1)


    def mixture_ratio(self):
        """ 混合質量比MRを計算して返す """
        return  self.m_ox / (self.m_fu*self.FO_ratio)

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
        Cp = self.conditions["gamma"] * self.conditions["R0"] / (self.conditions["gamma"]-1) / self.m_average_after_reaction()
        p_ratio_tmp = 1-np.power(self.conditions["Pj"]/self.chamber_pressure, (self.conditions["gamma"]-1)/self.conditions["gamma"])
        return np.sqrt(2*self.conditions["eta"]*Cp*self.Tf()*p_ratio_tmp)


    def V_j(self):
        """燃焼室による噴出速度損失を考慮した排気速度 """
        V_j_infinity = self.V_j_infinity()
## TODOだったが，推力損失めっちゃ小さいみたいだから，無視することにした
        return V_j_infinity


    def Isp(self):
        """ Ispを計算して返す """
        return self.V_j()/self.conditions["g"]


    def payload_lambda(self):
        """ ペイロード比を計算 """
        f_inert = self.f_inert()
        return (np.exp(-1*self.conditions["delta_v"]/(self.conditions["g"]*self.Isp() ))-f_inert) / (1-f_inert)

    
    def throat_temp(self):
        """スロート部の温度を計算して返す"""
        return 2*self.Tf() / (self.conditions["gamma"]+1)

    def total_rocket_mass(self):
        """ ロケット全体の質量を計算する """
        return self.conditions["m_l"] / self.payload_lambda()

    def C_F(self):
        """ 推力係数C_Fを計算する。最適膨張であることに注意。"""
        y = self.conditions["gamma"]
        return np.sqrt(2*y/(y-1)*np.power(2/(y+1), (y+1)/(y-1) ) * (1- np.power(self.conditions["Pj"]/self.chamber_pressure, (y-1)/y )))

    def A_t(self):
        """ 初期加速度からスロート面積を計算。燃焼室内の流れを無視し，燃焼室入口と出口の圧力損失を無視する。 """
        F_init = self.total_rocket_mass() * (self.conditions["g"]+self.conditions["init_a"])
        return F_init / self.C_F() / self.chamber_pressure

    
    def chamber_length(self, chamber_diameter):
        """ 燃焼室の直径(m)を受け取り，燃焼室内滞在時間から燃焼室長さを求めて返す. rho_av / rho_10 = 1と近似している。 """
        A_1 = np.pi * ( (chamber_diameter/2) ** 2)
        L_star = self.conditions["t_s"] * np.sqrt(self.conditions["gamma"]*self.conditions["R0"]*self.Tf()/self.m_average_after_reaction()) * np.power(2/(1+self.conditions["gamma"]), (self.conditions["gamma"]+1)/(2*(self.conditions["gamma"]-1)))
        return self.A_t() * L_star / A_1

    def total_prop_mass(self):
        """ 推進剤の全質量を計算する """
        return self.total_rocket_mass()*(1-self.payload_lambda())*(1-self.f_inert())

    def tank_length(self, tank_diameter):
        """ タンクの直径(m)を受け取り，全体の質量からタンクの長さを計算する。全体の密度は，燃料の平均密度で近似する。"""
        mixture_ratio = self.mixture_ratio()
        rho = (mixture_ratio + 1) / (1/self.rho_fu + mixture_ratio/self.rho_ox) *1000
        A_tank = np.pi * ((tank_diameter/2)**2)
        return self.total_prop_mass() / rho / A_tank


    def isentropic_flow_from_mach_to_A_ratio(self,M):
        """ 等エントロピー流れにおいて，マッハ数からスロート面積比を求める """
        y = self.conditions["gamma"]
        return np.power(2 * (1+(y-1)*(M**2)/2) / (y+1) , (y+1)/2 /(y-1)) / M

    def chamber_exit_mach_number(self):
        """ 燃焼室出口マッハ数M_1をニュートン・ラプソン法で求める.種は，0.1とする。"""
        return newton(self.isentropic_flow_from_mach_to_A_ratio, 0.1)

    def __throat_pressure_loss(self, M_1):
        y = self.conditions["gamma"]
        return np.power(1+(y-1)*(M_1**2)/2, y/(y-1))

    def __throat_velocity_loss(self, M_1):
        y = self.conditions["gamma"]
        p_ratio_tmp = np.power(self.conditions["Pj"]/self.chamber_pressure, (y-1)/y)
        return np.sqrt( (1- np.power(1+y*(M_1**2), (y-1)/y)*p_ratio_tmp / (1+(y-1)*(M_1**2)/2)) / (1-p_ratio_tmp))


    def throat_F_loss(self):
        """ 燃焼室による推力損失を計算する。ただ，せっかく頑張って書いたのに，どうも値が小さすぎるのか，全然計算が収束しない"""
        M_1 = self.chamber_exit_mach_number()
        return self.__throat_pressure_loss(M_1) * self.__throat_velocity_loss(M_1)


