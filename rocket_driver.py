# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import sys
import rocket_engine


reactants_name_list = ["LH2LO2", "LH2LF2", "hydrazine", "kerosene"]
reactants_dict = {
        "LH2LO2" : [ [2, 0, 0.07, 0.002], [1, 0, 1.14, 0.032] ], #水素・酸素
        "LH2LF2" : [ [1, 0, 0.07, 0.002], [1, 0, 1.51, 0.038] ],  #水素・フッ素
        "hydrazine" : [ [2, 50440, 1.00, 0.032], [1, -19.42, 1.43, 0.092] ], #ヒドラジンとN2O4
        "kerosene" : [ [1, -353500, 0.7495, 0.17034], [18.5, 0, 1.14, 0.032]] #ケロシン（ドデカン）と酸素
        }

products_dict = {
        "LH2LO2" : [ [2], [-241900 ], [0.018] ],
        "LH2LF2" : [ [2], [-271200 ], [0.02] ],
        "hydrazine" : [ [3, 4], [0, -241900], [0.028, 0.018]], #N2, H2O
        "kerosene" : [ [12, 13], [-393500, -241900], [0.044, 0.018]] #CO2, H2O
        }

#表示範囲
range_dict = {
        "LH2LO2" : [1, 8, 0.1],
        "LH2LF2" : [0.5, 8, 0.1],
        "hydrazine" : [0.5, 5, 0.5],
        "kerosene" : [0.01, 0.2, 0.001]
        }

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



#実行部分
#コマンドライン引数として，燃料の名前と燃焼室圧力をとることができる（デフォルトは"LH2LO2", 20000000）

if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] in reactants_name_list:
        reactants_name = sys.argv[1]
        chamber_pressure = 20000000
    elif len(sys.argv) == 3 and sys.argv[1] in reactants_name_list:
        reactants_name = sys.argv[1]
        chamber_pressure = int(sys.argv[2])
    else:
        reactants_name = "LH2LO2"
        chamber_pressure = 20000000

    reactants = reactants_dict[reactants_name]
    products = products_dict[reactants_name]

    range_minimum = range_dict[reactants_name][0]
    range_max = range_dict[reactants_name][1]
    range_interval = range_dict[reactants_name][2]
    ratio_array = np.arange(range_minimum,range_max,range_interval)
    sim_times = int( (range_max-range_minimum)/range_interval )
    payload_lambda_array = np.zeros(sim_times)

    for i in np.arange(sim_times):
        rocket = rocket_engine.Rocket(reactants, products, ratio_array[i], chamber_pressure, conditions)
        payload_lambda_array[i] = rocket.payload_lambda()
        print(payload_lambda_array[i])

    plt.plot(ratio_array, payload_lambda_array)
    plt.show()
