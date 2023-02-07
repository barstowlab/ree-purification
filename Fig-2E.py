# ------------------------------------------------------------------------------------------------ #
# Fig-S8E.py
# Last updated: January 9th, 2023
# Buz Barstow
# Code to calculate effect REE binding to a microbe with a single type of binding site. 
# Figure S8E in Genomic Characterization of Rare Earth Binding by Shewanella oneidensis by 
# Medin et al. 
# ------------------------------------------------------------------------------------------------ #

from utils.ConcentrationSolverUtils9 import InputData_3Metals_1Site, \
Calculate_Sub_Cycles_to_Reach_Target_Purity_and_M1_Remaining_3M_1S

from utils.vectorOutput3 import writeOutputMatrix, generateOutputMatrixWithHeaders, ensure_dir

from os.path import join

from scipy.interpolate import interp1d

import pdb

# ------------------------------------------------------------------------------------------------ #
# Prepare output
outputDir = 'output/Fig-S8/'
fileName = 'Fig-S8E.csv'
fileName_ex = 'Fig-S8E-examples.csv'
ensure_dir(outputDir)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
example_inv_kd1_array = \
[1.00, 1.0425974631137000, 1.0806783602831500, 1.1627844314923600, \
1.2857997031071100, 1.5519282263934800, 2.8465682481811900, 6.02694707548560]
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Calculate the number of sub-cycles to reach a target purity

kdBase = 1E-9
kd1Array = array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1], float) * kdBase

kd2 = 1.2563958631280851 * kdBase
kd3 = 2.28160979282314 * kdBase

nMT = 72e-9
nM1T = nMT/3
nM2T = nMT/3
nM3T = nMT/3
nBT = 35.033e-9
loadVol = 400e-6


targetPurities = [0.9, 0.95, 0.99]

i = 0
kdArray = []
while i < len(kd1Array):
	kd1 = kd1Array[i]
	kdArray.append([kd1, kd2, kd3])
	i += 1
	
scenarioDict = {}

for purity in targetPurities:
	
	scenarioKey = str(purity)

	returnDict = Calculate_Sub_Cycles_to_Reach_Target_Purity_and_M1_Remaining_3M_1S(\
	kdArray, nBT, nM1T, nM2T, nM3T, loadVol, purity)
	
	scenarioDict[scenarioKey] = returnDict

inv_kd_array = kdBase/kd1Array

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Figure out the sub cycles to target purity for the selected inverse kds from the earlier plots.


example_line_dict = {}

for key in scenarioDict.keys():
	example_line_dict[key] = {}
	example_line_dict[key]['example_horiz_lines'] = []
	example_line_dict[key]['example_vert_lines'] = []
	
	i = 0
	while i < len(example_inv_kd1_array):
		horizLine_x_Array = []
		horizLine_y_Array = []
		vertLine_x_Array = []
		vertLine_y_Array = []
		
		subCycleNumberArray = scenarioDict[key]['subCycleNumber']
		
		inv_kd1 = example_inv_kd1_array[i]
		invKd_to_subcycles_func = interp1d(inv_kd_array, subCycleNumberArray)
		subCycles_at_ex_invKd = float(invKd_to_subcycles_func(inv_kd1))

		vertLine_x_Array.append(inv_kd1)
		vertLine_x_Array.append(inv_kd1)
		vertLine_y_Array.append(0)
		vertLine_y_Array.append(subCycles_at_ex_invKd)
	
		horizLine_x_Array.append(0)
		horizLine_x_Array.append(inv_kd1)
		horizLine_y_Array.append(subCycles_at_ex_invKd)
		horizLine_y_Array.append(subCycles_at_ex_invKd)
		
		example_line_dict[key]['example_horiz_lines'].append([horizLine_x_Array, horizLine_y_Array])
		example_line_dict[key]['example_vert_lines'].append([vertLine_x_Array, vertLine_y_Array])
		
		i += 1
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Plot results and output to CSV

headers = ['inv_kd1']
vectorList = [inv_kd_array]

figure()

for key in scenarioDict.keys():
	subCycleNumberArray = scenarioDict[key]['subCycleNumber']
	plot(inv_kd_array, subCycleNumberArray, label=key)
	
	headers.append(key + '_subCycleNumber')
	vectorList.append(subCycleNumberArray)
	
	example_lines = example_line_dict[key]
	horiz_lines = example_lines['example_horiz_lines']
	vert_lines = example_lines['example_vert_lines']
	
	i = 0
	while i < len(horiz_lines):
		plot(horiz_lines[i][0], horiz_lines[i][1], color='black')
		plot(vert_lines[i][0], vert_lines[i][1], color='black')
		i += 1


grid()
legend()
xlabel('1E9/kd1')
ylabel('Sub-cycles to Target Purity')
show()
xlim(1, example_inv_kd1_array[-1]*1.1)

oMatrix = generateOutputMatrixWithHeaders(vectorList, headers)
writeOutputMatrix(join(outputDir, fileName), oMatrix)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
headers_examples = ['inv_kd1']
vectorList_examples = [example_inv_kd1_array]

for key in example_line_dict.keys():
	headers_examples.append(key + '_subCycleNumber')
	example_lines = example_line_dict[key]
	horiz_lines = example_lines['example_horiz_lines']
	
	tempVector = []
	
	i = 0
	while i < len(horiz_lines):
		tempVector.append(horiz_lines[i][1][0])
		i += 1
	
	vectorList_examples.append(tempVector)

oMatrixEx = generateOutputMatrixWithHeaders(vectorList_examples, headers_examples)
writeOutputMatrix(join(outputDir, fileName_ex), oMatrixEx)
# ------------------------------------------------------------------------------------------------ #
