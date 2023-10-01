# ------------------------------------------------------------------------------------------------ #
# Fig-4B.py
# Last updated: February 7th, 2023
# Buz Barstow
# Code to calculate effect REE purification with 3 microbes with 3 types of binding site. 
# Figure S10B in Genomic Characterization of Rare Earth Binding by Shewanella oneidensis by 
# Medin et al. 
# ------------------------------------------------------------------------------------------------ #


from utils.ConcentrationSolverUtils9 import InputData_3Metals_3Sites, \
Calculate_Sub_Cycles_to_Reach_Target_Purity_and_M1_Remaining_3M_3S_3Microbes, \
Calculate_Binding_Site_Fraction_Arrays, Increment_Binding_Site_Fraction

from utils.vectorOutput3 import writeOutputMatrix, ensure_dir, \
generateOutputMatrix_WithNonEqualVectors_WithHeaders

from os.path import join

from scipy.interpolate import interp1d


import pdb


# ------------------------------------------------------------------------------------------------ #
# Prepare output
outputDir = 'output/Fig-4/'
fileName = 'Fig-4B.csv'
fileName_ex = 'Fig-4B-examples.csv'
ensure_dir(outputDir)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Example binding site 1 fractions from Fig-3A-B.py
example_fB1_increase_array = [0.015, 0.05, 0.1, 0.2, 0.3, 0.4, 0.4999, 0.5]

targetPurities = [0.8, 0.9, 0.95, 0.99]
firstIndexArray = [2, 2, 2, 4]

# I use these ones for debugging
# targetPurities = [0.8]
# firstIndexArray = [2]
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
kd2_2 = kdBase*(1/7.999895847446188)
kd2_3 = kdBase*1

kd3_1 = kdBase*1
kd3_2 = kdBase*1
kd3_3 = kdBase*(1/7.997818606180251)

# These numbers are from our low ionic strength, high REE concentration dataset. 
nMT = 72e-9
nM1T = nMT/3
nM2T = nMT/3
nM3T = nMT/3
nBT = 35.033e-9
loadVol = 400e-6
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #

increment_array = [0.015, 0.025, 0.06, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]


# Microbe 1
microbe1_fB1_0 = 0.4618352672074258
microbe1_fB2_0 = 0.3703005831465813
microbe1_fB3_0 = 0.16786402382522572

microbe1_fB1_Array, microbe1_fB1_increase_array = \
Increment_Binding_Site_Fraction(microbe1_fB1_0, increment_array)

microbe1_fB2_Array, microbe1_fB3_Array = \
Calculate_Binding_Site_Fraction_Arrays(microbe1_fB1_Array, microbe1_fB2_0, microbe1_fB3_0)


# Microbe 2
microbe2_fB2_0 = 0.3703005831465813
microbe2_fB1_0 = 0.4618352672074258
microbe2_fB3_0 = 0.16786402382522572

microbe2_fB2_Array, microbe2_fB2_increase_array = \
Increment_Binding_Site_Fraction(microbe2_fB2_0, increment_array)

microbe2_fB3_Array, microbe2_fB1_Array = \
Calculate_Binding_Site_Fraction_Arrays(microbe2_fB2_Array, microbe2_fB3_0, microbe2_fB1_0)


# Microbe 3
microbe3_fB3_0 = 0.16786402382522572
microbe3_fB1_0 = 0.4618352672074258
microbe3_fB2_0 = 0.3703005831465813

microbe3_fB3_Array, microbe3_fB3_increase_array = \
Increment_Binding_Site_Fraction(microbe3_fB3_0, increment_array)

microbe3_fB2_Array, microbe3_fB1_Array = \
Calculate_Binding_Site_Fraction_Arrays(microbe3_fB3_Array, microbe3_fB2_0, microbe3_fB1_0)

# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Start calculating


scenarioDict = {}

i = 0
for purity in targetPurities:
	
	scenarioKey = str(purity)
	
	firstIndex = firstIndexArray[i]
	
	microbe1_fB1_Array_local = microbe1_fB1_Array[firstIndex:]
	microbe1_fB2_Array_local = microbe1_fB2_Array[firstIndex:]
	microbe1_fB3_Array_local = microbe1_fB3_Array[firstIndex:]

	microbe2_fB1_Array_local = microbe2_fB1_Array[firstIndex:]
	microbe2_fB2_Array_local = microbe2_fB2_Array[firstIndex:]
	microbe2_fB3_Array_local = microbe2_fB3_Array[firstIndex:]

	microbe3_fB1_Array_local = microbe3_fB1_Array[firstIndex:]
	microbe3_fB2_Array_local = microbe3_fB2_Array[firstIndex:]
	microbe3_fB3_Array_local = microbe3_fB3_Array[firstIndex:]
	
	microbe1_fB1_increase_array_local = microbe1_fB1_increase_array[firstIndex:]
	microbe2_fB2_increase_array_local = microbe2_fB2_increase_array[firstIndex:]
	microbe3_fB3_increase_array_local = microbe3_fB3_increase_array[firstIndex:]
	
	returnDict = Calculate_Sub_Cycles_to_Reach_Target_Purity_and_M1_Remaining_3M_3S_3Microbes(\
	kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, nBT, \
	microbe1_fB1_Array_local, microbe1_fB2_Array_local, microbe1_fB3_Array_local, \
	microbe2_fB1_Array_local, microbe2_fB2_Array_local, microbe2_fB3_Array_local, \
	microbe3_fB1_Array_local, microbe3_fB2_Array_local, microbe3_fB3_Array_local, \
	nM1T, nM2T, nM3T, loadVol, purity)
	
	returnDict['microbe1_fB1_increase_array'] = microbe1_fB1_increase_array_local
	returnDict['microbe2_fB2_increase_array'] = microbe2_fB2_increase_array_local
	returnDict['microbe3_fB3_increase_array'] = microbe3_fB3_increase_array_local
	
	scenarioDict[scenarioKey] = returnDict
	i += 1
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Figure out the sub cycles to target purity for the selected fractions from the earlier plots.

example_line_dict = {}

for key in scenarioDict.keys():
	example_line_dict[key] = {}
	example_line_dict[key]['example_horiz_lines'] = []
	example_line_dict[key]['example_vert_lines'] = []
	
	subCycleNumberArray = scenarioDict[key]['subCycleNumber']
	fB1_increase_array = scenarioDict[key]['microbe1_fB1_increase_array']
	fB1_to_subcycles_func = interp1d(fB1_increase_array, subCycleNumberArray)
	
	i = 0
	while i < len(example_fB1_increase_array):

		example_fB1_increase = example_fB1_increase_array[i]	
		
		if min(fB1_increase_array) <= example_fB1_increase <= max(fB1_increase_array):

			horizLine_x_Array = []
			horizLine_y_Array = []
			vertLine_x_Array = []
			vertLine_y_Array = []
				
			subCycles_at_ex_fB1 = float(fB1_to_subcycles_func(example_fB1_increase))

			vertLine_x_Array.append(example_fB1_increase)
			vertLine_x_Array.append(example_fB1_increase)
			vertLine_y_Array.append(0)
			vertLine_y_Array.append(subCycles_at_ex_fB1)
	
			horizLine_x_Array.append(0)
			horizLine_x_Array.append(example_fB1_increase)
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
	microbe1_fB1_increase_array = scenarioDict[key]['microbe1_fB1_increase_array']
	microbe2_fB2_increase_array = scenarioDict[key]['microbe2_fB2_increase_array']
	microbe3_fB3_increase_array = scenarioDict[key]['microbe3_fB3_increase_array']
	
	
	plot(microbe1_fB1_increase_array, subCycleNumberArray, label=key)
	
	
	headers.append(key + '_m1_fb1_inc')
	headers.append(key + '_subCycleNumber')
	
	vectorList.append(microbe1_fB1_increase_array)
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
xlabel('Microbe 1, Increase in Fraction of B1')
ylabel('Sub-cycles to Target Purity')
title('Fig. 4B. Effect of Increasing fB1 on Steps to Reach High Purity of Eu')
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



