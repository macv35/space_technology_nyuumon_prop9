 # -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import sys
import rocket_engine


propellant_name_list = ["LH2LO2", "LH2LF2", "hydrazine_N2O4", "kerosene", "LN2H4LHNO3", "ammonia", "hydrazine", "Al_ammonium_perchlorate",]
reactants_dict = {
        "LH2LO2" : [ [2, 0, 0.07, 0.002], [1, 0, 1.14, 0.032] ], #水素・酸素
        "LH2LF2" : [ [1, 0, 0.07, 0.002], [1, 0, 1.51, 0.038] ],  #水素・フッ素
        "hydrazine_N2O4" : [ [2, 50440, 1.00, 0.032], [1, -19.42, 1.43, 0.092] ], #ヒドラジンとN2O4
        "kerosene" : [ [1, -353500, 0.7495, 0.17034], [18.5, 0, 1.14, 0.032]], #ケロシン（ドデカン）と酸素
        "LN2H4LHNO3" : [ [5, 50440, 1.00, 0.032], [4, -174280, 1.50, 0.063] ],
        "ammonia" : [ [4, -71700, 0.60, 0.017], [3, 0, 1.14, 0.032]],
        "hydrazine": [[1, 50440, 1.00, 0.032],[1, 0, 1.14, 0.032]],
        "Al_ammonium_perchlorate" :[[10, 0, 2.7 ,0.027],[6,-29040,1.95,0.11749]],
        }
reactants_name_dict = {
        "LH2LO2" : ["H2", "O2"],
        "LH2LF2" : ["H2", "F2"],
        "hydrazine_N2O4" : ["N2H4", "N2O4"],
        "kerosene" : ["C12H26", "O2" ],
        "LN2H4LHNO3" : ["N2H4", "HNO3"],
        "ammonia": ["NH3", "O2"],
        "hydrazine":["N2H4","O2"],
        "Al_ammonium_perchlorate":["Al","NH4ClO4"],
        }

products_dict = {
        "LH2LO2" : [ [2], [-241900 ], [0.018] ],
        "LH2LF2" : [ [2], [-271200 ], [0.02] ],
        "hydrazine_N2O4" : [ [3, 4], [0, -241900], [0.028, 0.018]], #N2, H2O
        "kerosene" : [ [12, 13], [-393500, -241900], [0.044, 0.018]], #CO2, H2O
        "LN2H4LHNO3" : [ [7, 12], [0,-241900],[0.028,0.018]],
        "ammonia": [[2,6],[0,-241900],[0.028,0.018]],
        "hydrazine":[[1,2],[0,-241900],[0.028,0.018]],
        "Al_ammonium_perchlorate":[[4,2,12,3],[-1675700,-704200,-241900,0],[0.10196,0.133341,0.018,0.028]],
        }
products_name_dict = {
        "LH2LO2" : ["H2O"],
        "LH2LF2" : ["HF"],
        "hydrazine_N2O4" : ["N2", "H2O"],
        "kerosene" : ["CO2", "H2O" ],
        "LN2H4LHNO3": ["N2", "H2O"],
        "ammonia":["N2","H2O"],
        "hydrazine":["N2","H2O"],
        "Al_ammonium_perchlorate":["Al2O3","AlCl3","H20","N2"],
        }

#表示範囲
range_dict = {
        "LH2LO2" : [1, 8, 0.01],
        "LH2LF2" : [0.5, 8, 0.01],
        "hydrazine_N2O4" : [0.5, 5, 0.01],
        "kerosene" : [0.01, 0.2, 0.001],
        "LN2H4LHNO3" : [0.5, 5, 0.01],
        "ammonia":[1, 2, 0.01],
        "hydrazine" : [0.5, 2, 0.01],
        "Al_ammonium_perchlorate":[2, 11, 0.01],
        }

#燃焼室・タンク直径(m)
chamber_diameter = 0.5
tank_diameter = 3

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

def conditions_output():
    #初期条件出力
    print( "*** initial condition ***" )
    print( "target delta V in (m/s): ", conditions["delta_v"])
    print( "accelaration of gravity in (m/s^2): ", conditions["g"])
    print( "payload mass in (kg): ", conditions["m_l"])
    print( "initial accelaration in (G): ", conditions["init_a"]/conditions["g"])
    print( "standard propellant mixture density in (g/cm^3): ", conditions["rho_s"])
    print( "standard structure mass ratio: ", conditions["f_inert_s"])
    print( "Standard of temperature in (K): ", conditions["T0"])
    print( "nozzle efficiency: ", conditions["eta"])
    print( "atmosphere pressure in (Pa): ", conditions["Pj"])
    print( "specific heat ratio: ", conditions["gamma"])
    print( "standard gas constant in (J/mol*K): ", conditions["R0"])
    print( "stay time in chamber in (s): ", conditions["t_s"])
    print( "maximum heatresistant temperature at throat in (K): ", conditions["Tt_max"])
    print()


def output(lambda_max_rocket, chamber_pressure, propellant_name, chamber_diameter,tank_diameter):
    """シミュレーション結果出力"""
    print( "*** propellant information ***")
    print( "reaction")

    reaction_equation =  str(reactants_dict[propellant_name][0][0])+ " " + reactants_name_dict[propellant_name][0] + " + " + str(reactants_dict[propellant_name][1][0]) +" "+ reactants_name_dict[propellant_name][1] + "  ->  "
    for i in range( len(products_name_dict[propellant_name]) ):
        if i == 0:
            reaction_equation += str(products_dict[propellant_name][0][i]) +" "+ products_name_dict[propellant_name][i]
        else:
            reaction_equation += " + " + str(products_dict[propellant_name][0][i]) + " "+ products_name_dict[propellant_name][i]
    print( reaction_equation )
    print()

    conditions_output()

    #アウトプット出力
    print( "*** simulation output ***" )
    print( "1. ¥=ペイロード比 ~~~~ ¥= $¥Lambda = ", lambda_max_rocket.payload_lambda(),"$¥¥" )
    print( "2. ¥>燃焼室圧力 ¥> $P = ", chamber_pressure," ~¥rm{Pa}$¥¥")
    #    print( "propellant type: ", propellant_name)
    #    print( "F/O ratio in mole fraction: ", lambda_max_rocket.FO_ratio )
    print( "3. ¥>質量混合比 ¥>$¥rm{MR}=", lambda_max_rocket.mixture_ratio(),"$¥¥" )
    print( "4. ¥>断熱火炎温度 ¥>$T_f = ", lambda_max_rocket.Tf() ,"~ ¥rm{K}$¥¥")
    print( "5. ¥>比推力 ¥>$I_{¥rm{sp}}= ", lambda_max_rocket.Isp(), "~¥rm{s}$¥¥")
    print( "6. ¥>燃焼室径¥>$R_c = ", chamber_diameter," ~ ¥rm{m}$¥¥")
    print( "¥>燃焼室長さ¥>$L_c = ", lambda_max_rocket.chamber_length(chamber_diameter) ," ~¥rm{m}$¥¥")
    print( "7. ¥>タンク径¥>$R_{¥rm{tank}} = ", tank_diameter," ~ ¥rm{m}$¥¥")
    print( "¥>タンク長さ¥>$L_{¥rm{tank}} = ", lambda_max_rocket.tank_length(tank_diameter) ," ~ ¥rm{m}$¥¥")
    if lambda_max_rocket.throat_temp() <= lambda_max_rocket.conditions["Tt_max"]:
        print( "note: the throat temperature is ", lambda_max_rocket.throat_temp(), " K, so it is safe. ")
    else:
        print( "WARNING: the throat temperature is ", lambda_max_rocket.throat_temp(), " K, so it is dangerous. ")


def draw_graph_one_type(propellant_name,ratio_array, payload_lambda_array, Isp_array,MR_array, Tf_array):
    """推進剤の種類が一種類の時に，グラフを描画する"""
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
    

def show_result_for_one_type(propellant_name, chamber_pressure):
    """ 推進剤が一種類の時に結果を表示する"""
    reactants = reactants_dict[propellant_name]
    products = products_dict[propellant_name]

    range_minimum = range_dict[propellant_name][0]
    range_max = range_dict[propellant_name][1]
    range_interval = range_dict[propellant_name][2]
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
    
    #結果出力
    output(lambda_max_rocket, chamber_pressure, propellant_name, chamber_diameter, tank_diameter)
    #グラフ出力
    draw_graph_one_type(propellant_name,ratio_array, payload_lambda_array, Isp_array,MR_array, Tf_array)


def show_result_all(chamber_pressure):
    """全種類の推進剤について計算して表示する"""

    return 0


def main(argv):
    """ main関数。主にコマンドライン引数を処理する"""
    if len(argv) == 2 and argv[1] in propellant_name_list:
        propellant_name = argv[1]
        chamber_pressure = 20000000
    elif len(argv) == 3 and argv[1] in propellant_name_list:
        propellant_name = argv[1]
        chamber_pressure = int(argv[2])
    elif len(argv) == 2 and argv[1] == "all":
        propellant_name = argv[1]
        chamber_pressure = 20000000
    else:
        propellant_name = "LH2LO2"
        chamber_pressure = 20000000
    
    #結果の計算及び表示
    if propellant_name in propellant_name_list:
        show_result_for_one_type(propellant_name,chamber_pressure)
    elif propellant_name == "all":
        show_result_all(chamber_pressure)

    return 0

#実行部分
#コマンドライン引数として，燃料の名前と燃焼室圧力をとることができる（デフォルトは"LH2LO2", 20000000）
if __name__ == "__main__":
    main(sys.argv)

