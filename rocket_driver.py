# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
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

#燃焼室・タンク直径(m)
chamber_diameter = 1
tank_diameter = 4

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
        
    Isp_array = np.zeros(sim_times)
    for i in np.arange(sim_times):
        rocket = rocket_engine.Rocket(reactants, products, ratio_array[i], chamber_pressure, conditions)
        Isp_array[i] = rocket.Isp()
        
    MR_array = np.zeros(sim_times)
    for i in np.arange(sim_times):
        rocket = rocket_engine.Rocket(reactants, products, ratio_array[i], chamber_pressure, conditions)
        MR_array[i] = rocket.mixture_ratio()
        
    Tf_array = np.zeros(sim_times)
    for i in np.arange(sim_times):
        rocket = rocket_engine.Rocket(reactants, products, ratio_array[i], chamber_pressure, conditions)
        Tf_array[i] = rocket.Tf()

    #ペイロード比最大のロケットをインスタンス化
    lambda_max_index = np.argmax(payload_lambda_array)
    lambda_max_ratio = np.arange(range_minimum,range_max,range_interval)[lambda_max_index]
    lambda_max_rocket = rocket_engine.Rocket(reactants, products, lambda_max_ratio, chamber_pressure, conditions)

    #アウトプット出力
    print( "max payload ratio: ", lambda_max_rocket.payload_lambda() )
    print( "chamber pressure: ", chamber_pressure)
    print( "propellant type: ", reactants_name)
    print( "F/O ratio in mole fraction: ", lambda_max_rocket.FO_ratio )
    print( "F/O ratio in mass fraction: ", lambda_max_rocket.mixture_ratio() )
    print( "isentropic flame temperature: ", lambda_max_rocket.Tf() )
    print( "Isp: ", lambda_max_rocket.Isp() )
    print( "chamber diameter: ", chamber_diameter)
    print( "chamber length: ", lambda_max_rocket.chamber_length(chamber_diameter) )
    print( "tank diameter: ", tank_diameter)
    print( "tank length: ", lambda_max_rocket.tank_length(tank_diameter) )
    if lambda_max_rocket.throat_temp() <= lambda_max_rocket.conditions["Tt_max"]:
        print( "note: the throat temperature is ", lambda_max_rocket.throat_temp(), " K, so it is safe. ")
    else:
        print( "WARNING: the throat temperature is ", lambda_max_rocket.throat_temp(), " K, so it is dangerous. ")

    
    G = gs.GridSpec(3,3)
    axes1 = plt.subplot(G[:,0])
    plt.plot(ratio_array, payload_lambda_array)
    plt.xlabel('fuel/oxidizer mol ratio')
    plt.ylabel('payload ratio')
    plt.style.use('ggplot')
    axes2 = plt.subplot(G[0,1:3])
    plt.plot(ratio_array, Isp_array)
    plt.xlabel('fuel/oxidizer mol ratio')
    plt.ylabel('Isp[sec]')
    plt.style.use('ggplot')
    axes3 = plt.subplot(G[1,1:3])
    plt.plot(ratio_array, MR_array)
    plt.xlabel('fuel/oxidizer mol ratio')
    plt.ylabel('mixture ratio')
    plt.style.use('ggplot')
    axes4 = plt.subplot(G[2,1:3])
    plt.plot(ratio_array, Tf_array)
    plt.xlabel('fuel/oxidizer mol ratio')
    plt.ylabel('Tf[K]')
    plt.style.use('ggplot')
    plt.subplots_adjust(wspace=1, hspace=1)
    plt.show()