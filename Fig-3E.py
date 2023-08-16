# ------------------------------------------------------------------------------------------------ #
# Fig-3#.py
# Last updated: February 7th, 2023
# Buz Barstow
# Code to calculate effect REE purification with 1 microbe with 3 types of binding site. 
# Figure 3E in Genomic Characterization of Rare Earth Binding by Shewanella oneidensis by 
# Medin et al. 
# ------------------------------------------------------------------------------------------------ #


from utils.ConcentrationSolverUtils9 import InputData_3Metals_3Sites, \
Calculate_Sub_Cycles_to_Reach_Target_Purity_and_M1_Remaining_3M_3S, \
Calculate_Binding_Site_Fraction_Arrays

from utils.vectorOutput3 import writeOutputMatrix, ensure_dir, \
generateOutputMatrix_WithNonEqualVectors_WithHeaders

from os.path import join

import pdb


# ------------------------------------------------------------------------------------------------ #
# Prepare output
outputDir = 'output/Fig-3/'
fileName = 'Fig-3E.csv'
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

fB1_Array = [0.4618352672074258, 0.675, 0.70, 0.725, 0.75, 0.775, 0.8, 0.85, 0.9, 0.95]
fB2_0 = 0.3703005831465813
fB3_0 = 0.16786402382522572

fB2_Array, fB3_Array = Calculate_Binding_Site_Fraction_Arrays(fB1_Array, fB2_0, fB3_0)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
targetPurities = [0.8, 0.9, 0.95, 0.99]
firstIndexArray = [1, 4, 6, 7]

scenarioDict = {}

i = 0
for purity in targetPurities:
	
	scenarioKey = str(purity)
	
	firstIndex = firstIndexArray[i]
	
	fB1_Array_local = fB1_Array[firstIndex:]
	fB2_Array_local = fB2_Array[firstIndex:]
	fB3_Array_local = fB3_Array[firstIndex:]
	
	fB_Array_local = []
	j = 0
	while j < len(fB1_Array_local):
		fB_Array_local.append([fB1_Array_local[j], fB2_Array_local[j], fB3_Array_local[j]]) 
		j += 1

	print('Target Purity: ' + str(purity))
	
	returnDict = Calculate_Sub_Cycles_to_Reach_Target_Purity_and_M1_Remaining_3M_3S(\
	kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
	nBT, fB_Array_local, nM1T, nM2T, nM3T, loadVol, purity, 'nMT/2')
	
	
	scenarioDict[scenarioKey] = returnDict
	i += 1
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
# Plot results and output to CSV

headers = []
vectorList = []

figure()

for key in scenarioDict.keys():
	subCycleNumberArray = scenarioDict[key]['subCycleNumber']
	
	fB1_Array = scenarioDict[key]['fB1_Array']
	plot(fB1_Array, subCycleNumberArray, label=key)
	
	headers.append(key + '_fB1')
	headers.append(key + '_subCycleNumber')
	vectorList.append(fB1_Array)
	vectorList.append(subCycleNumberArray)
	
	
	

grid()
legend()
xlabel('Fraction of B1')
ylabel('Sub-cycles to Target Purity')
show()

oMatrix = generateOutputMatrix_WithNonEqualVectors_WithHeaders(vectorList, headers)
writeOutputMatrix(join(outputDir, fileName), oMatrix)
# ------------------------------------------------------------------------------------------------ #





