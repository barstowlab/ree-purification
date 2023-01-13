from utils.ConcentrationSolverUtils9 import InputData_3Metals_3Sites, \
Calculate_Sub_Cycles_to_Reach_Target_Purity_and_M1_Remaining_3M_3S, \
Calculate_Binding_Site_Fraction_Arrays

from utils.vectorOutput3 import writeOutputMatrix, generateOutputMatrixWithHeaders, ensure_dir

from os.path import join

import pdb


# ------------------------------------------------------------------------------------------------ #
# Prepare output
outputDir = 'output/Fig-S9/'
fileName = 'Fig-S9E.csv'
ensure_dir(outputDir)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Calculate the number of sub-cycles to reach a target purity of M1 as a function of the fraction of 
# site 1 
# I'll start with the kd set 

kdBase = 1E-9
kd1_1 = kdBase*(1/8.0)
kd1_2 = kdBase*1
kd1_3 = kdBase*1

kd2_1 = kdBase*1
kd2_2 = kdBase*(1/7.999884764962862)
kd2_3 = kdBase*1

kd3_1 = kdBase*1
kd3_2 = kdBase*1
kd3_3 = kdBase*(1/7.99786621429086)

# These numbers are from our low ionic strength, high REE concentration dataset. 
nMT = 72e-9
nM1T = nMT/3
nM2T = nMT/3
nM3T = nMT/3
nBT = 35.033e-9
loadVol = 400e-6
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Calculate the fractions of each binding site. This is a slightly tricky calculation. 
# We start out with the original fitted fractions of each site, and then figure out how to reduce
# the fractions of binding sites 2 and 3 as the fraction of binding site 1 increases. 

fB1_Array = [0.4618352672074258, 0.475, 0.5, 0.55, 0.6, 0.65, 0.70, 0.75, 0.8, 0.85, 0.9, 0.95]
fB2_0 = 0.3703005831465813
fB3_0 = 0.16786402382522572

fB2_Array, fB3_Array = Calculate_Binding_Site_Fraction_Arrays(fB1_Array, fB2_0, fB3_0)

fB_Array = []
i = 0
while i < len(fB1_Array):
	fB_Array.append([fB1_Array[i], fB2_Array[i], fB3_Array[i]]) 
	i += 1

firstIndex = -3
fB_Array = fB_Array[firstIndex:-1]
fB1_Array = fB1_Array[firstIndex:-1]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
targetPurities = [0.99]

scenarioDict = {}

for purity in targetPurities:
	
	scenarioKey = str(purity)

	returnDict = Calculate_Sub_Cycles_to_Reach_Target_Purity_and_M1_Remaining_3M_3S(\
	kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
	nBT, fB_Array, nM1T, nM2T, nM3T, loadVol, purity)
	
	scenarioDict[scenarioKey] = returnDict
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Plot results and output to CSV

headers = ['fB1']
vectorList = [fB1_Array]

figure()

for key in scenarioDict.keys():
	subCycleNumberArray = scenarioDict[key]['subCycleNumber']
	plot(fB1_Array, subCycleNumberArray, label=key)
	
	headers.append(key + '_subCycleNumber')
	vectorList.append(subCycleNumberArray)
	

grid()
legend()
xlabel('Fraction of B1')
ylabel('Sub-cycles to Target Purity')
show()

oMatrix = generateOutputMatrixWithHeaders(vectorList, headers)
writeOutputMatrix(join(outputDir, fileName), oMatrix)
# ------------------------------------------------------------------------------------------------ #





