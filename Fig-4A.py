# ------------------------------------------------------------------------------------------------ #
# Fig-4A.py
# Last updated: February 7th, 2023
# Buz Barstow
# Code to calculate effect REE purification with 3 microbes with 3 types of binding site. 
# Figure 4A in Genomic Characterization of Rare Earth Binding by Shewanella oneidensis by 
# Medin et al. 
# ------------------------------------------------------------------------------------------------ #

from utils.ConcentrationSolverUtils9 import Run_Macro_Cycle_3Metals_3Sites_3Microbes, \
InputData_3Metals_3Sites, Calculate_Binding_Site_Fraction_Arrays, Increment_Binding_Site_Fraction

from utils.vectorOutput3 import writeOutputMatrix, ensure_dir, \
generateOutputMatrix_WithNonEqualVectors_WithHeaders

from os.path import join

import pdb

# ------------------------------------------------------------------------------------------------ #
# Prepare output
outputDir = 'output/Fig-4/'
fB1_vs_fMx_FileName = 'Fig-4A.csv'
ensure_dir(outputDir)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Calculate the purity of M1 and separation factor as eluant is repeatedly passed through column. 
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

increment_array = [0.015, 0.025, 0.06, 0.05, 0.05, 0.10, 0.10, 0.10]


# Microbe 1
microbe1_fB1_0 = 0.4618352672074258
microbe1_fB2_0 = 0.3703005831465813
microbe1_fB3_0 = 0.16786402382522572

microbe1_fB1_Array, fB_increase_array = \
Increment_Binding_Site_Fraction(microbe1_fB1_0, increment_array)

microbe1_fB2_Array, microbe1_fB3_Array = \
Calculate_Binding_Site_Fraction_Arrays(microbe1_fB1_Array, microbe1_fB2_0, microbe1_fB3_0)


# Microbe 2
microbe2_fB2_0 = 0.3703005831465813
microbe2_fB1_0 = 0.4618352672074258
microbe2_fB3_0 = 0.16786402382522572

microbe2_fB2_Array, fB_increase_array = \
Increment_Binding_Site_Fraction(microbe2_fB2_0, increment_array)

microbe2_fB3_Array, microbe2_fB1_Array = \
Calculate_Binding_Site_Fraction_Arrays(microbe2_fB2_Array, microbe2_fB3_0, microbe2_fB1_0)


# Microbe 3
microbe3_fB3_0 = 0.16786402382522572
microbe3_fB1_0 = 0.4618352672074258
microbe3_fB2_0 = 0.3703005831465813

microbe3_fB3_Array, fB_increase_array = \
Increment_Binding_Site_Fraction(microbe3_fB3_0, increment_array)

microbe3_fB2_Array, microbe3_fB1_Array = \
Calculate_Binding_Site_Fraction_Arrays(microbe3_fB3_Array, microbe3_fB2_0, microbe3_fB1_0)


maxSubCycles = 26
targetPurity = 0.9

# ------------------------------------------------------------------------------------------------ #
scenarioDict = {}

i = 0
while i < len(microbe1_fB1_Array):

	scenarioKey = str(fB_increase_array[i])
	
	resultsDict = {}
	
	# -------------------------------------------------------------------------------------------- #
	# Prepare Microbe 1
	microbe1_fB1 = microbe1_fB1_Array[i]
	microbe1_fB2 = microbe1_fB2_Array[i]
	microbe1_fB3 = microbe1_fB3_Array[i]
	
	inputData_microbe1 = \
	InputData_3Metals_3Sites(\
	kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
	nM1T, nM2T, nM3T, nBT, microbe1_fB1, microbe1_fB2, microbe1_fB3, \
	loadVol, adjustLoadVol=True, adjustBindingSites=True, maxSubCycles=maxSubCycles, \
	targetPurity=targetPurity, adjustBindingSiteTarget='nM1') 
	# -------------------------------------------------------------------------------------------- #
	
	# -------------------------------------------------------------------------------------------- #
	# Prepare Microbe 2
	microbe2_fB1 = microbe2_fB1_Array[i]
	microbe2_fB2 = microbe2_fB2_Array[i]
	microbe2_fB3 = microbe2_fB3_Array[i]
	
	inputData_microbe2 = \
	InputData_3Metals_3Sites(\
	kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
	nM1T, nM2T, nM3T, nBT, microbe2_fB1, microbe2_fB2, microbe2_fB3, \
	loadVol, adjustLoadVol=True, adjustBindingSites=True, maxSubCycles=maxSubCycles, \
	targetPurity=targetPurity, adjustBindingSiteTarget='nM2') 
	# -------------------------------------------------------------------------------------------- #
	
	# -------------------------------------------------------------------------------------------- #
	# Prepare Microbe 3
	microbe3_fB1 = microbe3_fB1_Array[i]
	microbe3_fB2 = microbe3_fB2_Array[i]
	microbe3_fB3 = microbe3_fB3_Array[i]
	
	inputData_microbe3 = \
	InputData_3Metals_3Sites(\
	kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
	nM1T, nM2T, nM3T, nBT, microbe3_fB1, microbe3_fB2, microbe3_fB3, \
	loadVol, adjustLoadVol=True, adjustBindingSites=True, maxSubCycles=maxSubCycles, \
	targetPurity=targetPurity, adjustBindingSiteTarget='nM3') 
	# -------------------------------------------------------------------------------------------- #
	
	
	macroCycleDiagnostic = Run_Macro_Cycle_3Metals_3Sites_3Microbes(inputData_microbe1, \
	inputData_microbe2, inputData_microbe3)
	
	macroCycleDiagnosticReturnDict = macroCycleDiagnostic.summarize_macro_cycle()

	# Parse out the output from the macro cycle
	resultsDict['cycleArray'] = macroCycleDiagnosticReturnDict['cycleArray']

	resultsDict['fM1_Microbe3_Wash_Array'] = \
	macroCycleDiagnosticReturnDict['fM1_Microbe3_Wash_Array']
	
	scenarioDict[scenarioKey] = resultsDict
	
	i += 1
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Plot results and output to csv

headers_fM1 = []
vectorList_fM1 = []

# Plot out the purity of metal 1 as a function of cycle number.
figure()
i = 0
scenarioKeys = list(scenarioDict.keys())
while i < len(scenarioKeys):
	key = scenarioKeys[i]
	cycleArray = scenarioDict[key]['cycleArray']
	fM1EluantArray = scenarioDict[key]['fM1_Microbe3_Wash_Array']
	plot(cycleArray, fM1EluantArray, label=key)
	
	headers_fM1.append(key + '_' + 'cycle')
	headers_fM1.append(key + '_' + 'fM1')
	vectorList_fM1.append(cycleArray)
	vectorList_fM1.append(fM1EluantArray)
	
	i += 1

title('Fig. 4A. 3 Sites, Purity of M1 in Microbe 3 Wash vs. Cycle Number')
legend()
grid()
xlabel('Cycle Number')
ylabel('Purity of M1 in Eluant')
show()

oMatrix_fM1 = generateOutputMatrix_WithNonEqualVectors_WithHeaders(vectorList_fM1, headers_fM1)
writeOutputMatrix(join(outputDir, fB1_vs_fMx_FileName), oMatrix_fM1)


# ------------------------------------------------------------------------------------------------ #



