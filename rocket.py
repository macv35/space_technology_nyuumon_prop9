import numpy as np
import matplotlib.pyplot as plt

def Tf_equation(Q1, FO_ratio):
	if FO_ratio==2: #化学量論比
	    return Q1/74 + 298.15
	elif FO_ratio > 2: #燃料過多
		return Q1/((FO_ratio-2)*29.1+2*FO_ratio*37) + 298.15
	elif FO_ratio < 2: #酸化剤過多
		return Q1/((1-0.5*FO_ratio)*29.1+FO_ratio*37) + 298.15

conditions = {
	"delta_v" : 7000,
	"g" : 9.8,
	"ml" : 2000, #payload
	"init_a" : 2.94, #初期加速度
	"rho_s" : 0.3, #標準推進剤混合密度
	"rho_fu" : 0.07, #燃料密度
	"rho_ox" : 1.14, #酸化剤密度
	"f_inert_s" : 0.1, #標準構造質量比
	"H_fuel" : 0,
	"H_ox" : 0,
	"H_prod" : -241900,
	"Cp" : 29.1,
	"T0" : 298.15,
	"yiita" : 0.96,
	"P0" : 10000,
	"Pj" : 15000000,
	"gamma" : 1.4
}
ratio_array = np.arange(0,6.0,0.01) #酸化剤を1とする
deltaH_array = ratio_array * conditions["H_fuel"] + conditions["H_ox"] - ratio_array * conditions["H_prod"]
Tf_array = np.zeros(600)
for i in np.arange(600):
    Tf_array = Tf_equation(deltaH_array[i], ratio_array[i])

pressure_ratio = (conditions["Pj"]/conditions["P0"])**(2/7)
Vj = (2*conditions["yiita"]*conditions["yiita"]*conditions["T0"]*(1-pressure_ratio))**0.5
Isp = Vj/conditions["g"]
rho_array = (ratio_array*conditions["rho_fu"] + conditions["rho_ox"])/(ratio_array + 1)

f_inert_array = 1/((1/conditions["f_inert_s"]-1)*(rho_array/conditions["rho_s"])+1)

payload_ratio = (np.exp(-conditions["delta_v"]/(conditions["g"]*Isp))-f_inert_array)/(1-f_inert_array)
print(payload_ratio)

plt.plot(ratio_array, f_inert_array)
plt.show()