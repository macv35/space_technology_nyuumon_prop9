# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import rocket_engine

if __name__ == "__main__":
    #初期条件
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

    #燃焼室圧力
    chamber_pressure = 20000000

    #水素・酸素
    reactants = [ [2, 0, 0.07, 0.002], [1, 0, 1.14, 0.032] ]
    products = [ [2], [-241900 ], [0.018] ]

    #水素・フッ素
#    reactants = [ [1, 0, 0.07, 0.002], [1, 0, 1.51, 0.038] ]
#    products = [ [2], [-271200 ], [0.02] ]

    #表示範囲
    range_minimum = 0.1
    range_max = 6.0

    ratio_array = np.arange(range_minimum,range_max,0.01)
    sim_times = int( (range_max-range_minimum)/0.01 )
    payload_lambda_array = np.zeros(sim_times)
    for i in np.arange(sim_times):
        rocket = rocket_engine.Rocket(reactants, products, ratio_array[i], chamber_pressure, conditions)
        payload_lambda_array[i] = rocket.payload_lambda()
        print(payload_lambda_array[i])

    plt.plot(ratio_array, payload_lambda_array)
    plt.show()
