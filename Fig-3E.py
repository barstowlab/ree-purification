# ------------------------------------------------------------------------------------------------ #
# Fig-3E.py
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

from scipy.interpolate import interp1d



# ------------------------------------------------------------------------------------------------ #
# Prepare output
outputDir = 'output/Fig-3/'
fileName = 'Fig-3E.csv'
fileName_ex = 'Fig-3E-examples.csv'
ensure_dir(outputDir)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Calculate the number of sub-cycles to reach a target purity of M1 as a function of the fraction of 
# site 1 
# I'll start with the kd set 

kdBase = 1E-6
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
# Example binding site 1 fractions from Fig-3A-B.py
example_fB1_array = [0.7, 0.8, 0.85, 0.9, 0.95]
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# Figure out the sub cycles to target purity for the selected fractions from the earlier plots.

example_line_dict = {}

for key in scenarioDict.keys():
	example_line_dict[key] = {}
	example_line_dict[key]['example_horiz_lines'] = []
	example_line_dict[key]['example_vert_lines'] = []
	
	subCycleNumberArray = scenarioDict[key]['subCycleNumber']
	fB1_Array = scenarioDict[key]['fB1_Array']
	fB1_to_subcycles_func = interp1d(fB1_Array, subCycleNumberArray)
	
	i = 0
	while i < len(example_fB1_array):

		example_fB1 = example_fB1_array[i]	
		
		if min(fB1_Array) <= example_fB1 <= max(fB1_Array):

			horizLine_x_Array = []
			horizLine_y_Array = []
			vertLine_x_Array = []
			vertLine_y_Array = []
				
			subCycles_at_ex_fB1 = float(fB1_to_subcycles_func(example_fB1))

			vertLine_x_Array.append(example_fB1)
			vertLine_x_Array.append(example_fB1)
			vertLine_y_Array.append(0)
			vertLine_y_Array.append(subCycles_at_ex_fB1)
	
			horizLine_x_Array.append(0)
			horizLine_x_Array.append(example_fB1)
			horizLine_y_Array.append(subCycles_at_ex_fB1)
			horizLine_y_Array.append(subCycles_at_ex_fB1)
		
			example_line_dict[key]['example_horiz_lines'].append(\
			[horizLine_x_Array, horizLine_y_Array])
			
			example_line_dict[key]['example_vert_lines'].append(\
			[vertLine_x_Array, vertLine_y_Array])
		
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
xlabel('Fraction of B1')
ylabel('Sub-cycles to Target Purity')
title('Fig. 3E. Effect of Increasing fB1 on Steps to Reach High Purity of Eu')
xlim(0.65, 1)
show()

oMatrix = generateOutputMatrix_WithNonEqualVectors_WithHeaders(vectorList, headers)
writeOutputMatrix(join(outputDir, fileName), oMatrix)
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# Write out the examples to a CSV table

headers_examples = []
vectorList_examples = []

for key in example_line_dict.keys():
	
	# Set up the headers
	headers_examples.append(key + '_fB1')
	headers_examples.append(key + '_subCycleNumber')
	
	example_lines = example_line_dict[key]
	horiz_lines = example_lines['example_horiz_lines']
	vert_lines = example_lines['example_vert_lines']
	
	fB1_tempVector = []
	subCycle_tempVector = []
	
	i = 0
	while i < len(vert_lines):
		fB1 = vert_lines[i][0][0]
		fB1_tempVector.append(fB1)
		
		subCycles = horiz_lines[i][1][1]
		subCycle_tempVector.append(subCycles)
		
		i += 1
		
	vectorList_examples.append(fB1_tempVector)
	vectorList_examples.append(subCycle_tempVector)
	

oMatrixEx = generateOutputMatrix_WithNonEqualVectors_WithHeaders(\
vectorList_examples, headers_examples)

writeOutputMatrix(join(outputDir, fileName_ex), oMatrixEx)
# ------------------------------------------------------------------------------------------------ #


