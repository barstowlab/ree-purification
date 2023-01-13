from utils.ConcentrationSolverUtils8 import BindState_3Metals_1Site, InputData_3Metals_1Site, \
kd_optimize_1site_func, LoadState_3Metals_1Site

from scipy.optimize import minimize

kdBase = 1E-9
kd1 = 1
kd2 = 1.5
kd3 = 2

x0 = [kd2, kd3]

nMT = 72e-9
nM1T = nMT/3
nM2T = nMT/3
nM3T = nMT/3
nBT = 35.033e-9
loadVol = 400e-6

# These fractions are from our low ionic strength, high REE concentration dataset. 
# In this case, M1 is actually Eu, M2 is Yb
#fMxB_target = [0.4008, 0.3531, 0.2460]

# These are binding fractions for the true wild-type S. oneidensis
fMxB_target = [0.3925, 0.3538, 0.2537]


init_input_data_3M_1S = InputData_3Metals_1Site(kd1, kd2, kd3, nM1T, nM2T, nM3T, nBT, loadVol)

x = minimize(kd_optimize_1site_func, x0, args=(kdBase, fMxB_target, init_input_data_3M_1S), \
options={'disp': False})

print("Fitting done")

kd2Fit = x['x'][0]*kdBase
kd3Fit = x['x'][1]*kdBase

init_input_data_3M_1S_verify = InputData_3Metals_1Site(kdBase, kd2Fit, kd3Fit, nM1T, nM2T, nM3T, \
nBT, loadVol)

load_state_verify = LoadState_3Metals_1Site()
load_state_verify.init_with_input_data(init_input_data_3M_1S_verify)
bindState_verify = BindState_3Metals_1Site(load_state_verify)

print('fM1b: ' + str(bindState_verify.fM1b)) 
print('fM2b: ' + str(bindState_verify.fM2b))
print('fM3b: ' + str(bindState_verify.fM3b)) 

print('kd1: ' + str(kd1*kdBase))
print('kd2: ' + str(kd2Fit))
print('kd3: ' + str(kd3Fit))