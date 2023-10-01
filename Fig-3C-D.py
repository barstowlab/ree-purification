# ------------------------------------------------------------------------------------------------ #
# Fig-3C-D.py
# Last updated: February 7th, 2023
# Buz Barstow
# Code to calculate effect REE purification with 1 microbe with 3 types of binding site. 
# Figures 3C and D in Genomic Characterization of Rare Earth Binding by Shewanella oneidensis by 
# Medin et al. 
# ------------------------------------------------------------------------------------------------ #


from utils.ConcentrationSolverUtils9 import Run_Macro_Cycle_3Metals_3Sites, \
InputData_3Metals_3Sites, Calculate_Binding_Site_Fraction_Arrays

from utils.vectorOutput3 import writeOutputMatrix, ensure_dir, \
generateOutputMatrix_WithNonEqualVectors_WithHeaders

from os.path import join

import pdb

# ------------------------------------------------------------------------------------------------ #
# Prepare output
outputDir = 'output/Fig-3/'
fB1_vs_fMx_FileName = 'Fig-3C.csv'
fB1_vs_sep_factor_FileName = 'Fig-3D.csv'
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

fB1_Array = [0.4618352672074258, 0.473, 0.482, 0.502, 0.7, 0.8, 0.85, 0.9, 0.95]
fB2_0 = 0.3703005831465813
fB3_0 = 0.16786402382522572

fB2_Array, fB3_Array = Calculate_Binding_Site_Fraction_Arrays(fB1_Array, fB2_0, fB3_0)

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
scenarioDict = {}

i = 0
while i < len(fB1_Array):

	scenarioKey = str(fB1_Array[i])
	
	resultsDict = {}
	
	fB1 = fB1_Array[i]
	fB2 = fB2_Array[i]
	fB3 = fB3_Array[i]
	

	inputData = \
	InputData_3Metals_3Sites(\
	kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
	nM1T, nM2T, nM3T, nBT, fB1, fB2, fB3, loadVol, adjustLoadVol=True, adjustBindingSites=True, \
	adjustBindingSiteTarget='nMT/2') 

	macroCycleDiagnostic = Run_Macro_Cycle_3Metals_3Sites(inputData)
	macroCycleDiagnosticReturnDict = macroCycleDiagnostic.summarize_macro_cycle()

	# Parse out the output from the macro cycle
	resultsDict['cycleArray'] = macroCycleDiagnosticReturnDict['cycleArray']

	resultsDict['fM1EluantArray'] = macroCycleDiagnosticReturnDict['fM1EluantArray']
	resultsDict['fM2EluantArray'] = macroCycleDiagnosticReturnDict['fM2EluantArray']
	resultsDict['fM3EluantArray'] = macroCycleDiagnosticReturnDict['fM3EluantArray']

	resultsDict['alpha1_2_BindArray'] = macroCycleDiagnosticReturnDict['alpha1_2_BindArray']
	resultsDict['alpha1_3_BindArray'] = macroCycleDiagnosticReturnDict['alpha1_3_BindArray']
	resultsDict['alpha2_3_BindArray'] = macroCycleDiagnosticReturnDict['alpha2_3_BindArray']
	
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
	fM1EluantArray = scenarioDict[key]['fM1EluantArray']
	plot(cycleArray, fM1EluantArray, label=key)
	
	headers_fM1.append(key + '_' + 'cycle')
	headers_fM1.append(key + '_' + 'fM1')
	vectorList_fM1.append(cycleArray)
	vectorList_fM1.append(fM1EluantArray)
	
	i += 1

title('Fig. 3C. 3 Sites, Purity of M1 in Eluant vs. Cycle Number')
legend()
grid()
xlabel('Cycle Number')
ylabel('Purity of M1 in Eluant')
show()

oMatrix_fM1 = generateOutputMatrix_WithNonEqualVectors_WithHeaders(vectorList_fM1, headers_fM1)
writeOutputMatrix(join(outputDir, fB1_vs_fMx_FileName), oMatrix_fM1)


# Plot out the separation factors as a function of cycle number.
headers_sep = []
vectorList_sep = []

figure()
i = 0
while i < len(scenarioKeys):
	key = scenarioKeys[i]
	cycleArray = scenarioDict[key]['cycleArray']
	alpha1_2_BindArray = scenarioDict[key]['alpha1_2_BindArray']
	alpha1_3_BindArray = scenarioDict[key]['alpha1_3_BindArray']
	alpha2_3_BindArray = scenarioDict[key]['alpha2_3_BindArray']
	
	plot(cycleArray, alpha1_2_BindArray, label=key + ', alpha1_2')
	plot(cycleArray, alpha1_3_BindArray, label=key + ', alpha1_3')
	plot(cycleArray, alpha2_3_BindArray, label=key + ', alpha2_3')
	
	headers_sep.append(key + '_' + 'cycle')
	headers_sep.append(key + '_' + 'alpha1_2')
	headers_sep.append(key + '_' + 'alpha1_3')
	headers_sep.append(key + '_' + 'alpha2_3')
	
	vectorList_sep.append(cycleArray)
	vectorList_sep.append(alpha1_2_BindArray)
	vectorList_sep.append(alpha1_3_BindArray)
	vectorList_sep.append(alpha2_3_BindArray)
	
	i += 1

title('Fig. 3D. 3 Sites, Separation Factor vs. Cycle Number')
legend()
grid()
xlabel('Cycle Number')
ylabel('Separation Factor')
show()

oMatrix_sep = generateOutputMatrix_WithNonEqualVectors_WithHeaders(vectorList_sep, headers_sep)
writeOutputMatrix(join(outputDir, fB1_vs_sep_factor_FileName), oMatrix_sep)

# ------------------------------------------------------------------------------------------------ #



