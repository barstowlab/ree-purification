from utils.ConcentrationSolverUtils8 import InputData_3Metals_3Sites, nb_optimize_3sites_3kd_func,\
LoadState_3Metals_3Sites, BindState_3Metals_3Sites

from scipy.optimize import minimize


kdBase = 1E-9


# Original, works very well with 1 binding site
# kd1_1 = kdBase
# kd1_2 = kdBase*1.52
# kd1_3 = kdBase*3.273


invKdBound = 8
invKdStart = 4

kd1_1 = kdBase*(1/invKdStart)
kd1_2 = kdBase*1
kd1_3 = kdBase*1

kd2_1 = kdBase*1
kd2_2 = kdBase*(1/invKdStart)
kd2_3 = kdBase*1

kd3_1 = kdBase*1
kd3_2 = kdBase*1
kd3_3 = kdBase*(1/invKdStart)

# xGuess are fB1, fB2, fB3, 1/kd1_1, 1/kd2_2, 1/kd3_3
xGuess = [0.3333, 0.3333, 0.3333, 1/kd1_1, 1/kd2_2, 1/kd3_3]

# These numbers are from our low ionic strength, high REE concentration dataset. 
nMT = 72e-9
nM1T = nMT/3
nM2T = nMT/3
nM3T = nMT/3
nBT = 35.033e-9
loadVol = 400e-6


# These fractions are from our low ionic strength, high REE concentration dataset. 
# In this case, M1 is actually Eu, M2 is Yb
fMxB_target = [0.3925, 0.3538, 0.2537]


init_input_data_3M_3S = InputData_3Metals_3Sites(kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, 
kd3_1, kd3_2, kd3_3, \
nM1T, nM2T, nM3T, nBT, xGuess[0], xGuess[1], xGuess[2], loadVol)


bnds = ((0.01,0.99), (0.01,0.99), (0.01,0.99), (1, invKdBound), (1, invKdBound), (1, invKdBound))


x = minimize(nb_optimize_3sites_3kd_func, xGuess, \
args=(kdBase, fMxB_target, init_input_data_3M_3S), options={'disp': False}, bounds=bnds)

print("Fitting done")

fB1_Fit = x['x'][0]
fB2_Fit = x['x'][1]
fB3_Fit = x['x'][2]

kd1_1_Fit = (1/x['x'][3])*kdBase
kd2_2_Fit = (1/x['x'][4])*kdBase
kd3_3_Fit = (1/x['x'][5])*kdBase



init_input_data_3M_3S_verify = \
InputData_3Metals_3Sites(\
kd1_1_Fit, kd1_2, kd1_3, kd2_1, kd2_2_Fit, kd2_3, kd3_1, kd3_2, kd3_3_Fit, \
nM1T, nM2T, nM3T, nBT, \
fB1_Fit, fB2_Fit, fB3_Fit, loadVol)

load_state_verify = LoadState_3Metals_3Sites()
load_state_verify.init_with_input_data(init_input_data_3M_3S_verify)
bindState_verify = BindState_3Metals_3Sites(load_state_verify)

print('fM1b: ' + str(bindState_verify.fM1b)) 
print('fM2b: ' + str(bindState_verify.fM2b))
print('fM3b: ' + str(bindState_verify.fM3b)) 