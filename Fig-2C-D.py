# ------------------------------------------------------------------------------------------------ #
# Fig-2C-D.py
# Last updated: February 7th, 2023
# Buz Barstow
# Code to calculate effect REE binding to a microbe with a single type of binding site. 
# Figures 2C and D in Genomic Characterization of Rare Earth Binding by Shewanella oneidensis by 
# Medin et al. 
# ------------------------------------------------------------------------------------------------ #

from utils.ConcentrationSolverUtils9 import Run_Macro_Cycle_3Metals_1Site, \
InputData_3Metals_1Site

from utils.vectorOutput3 import writeOutputMatrix, ensure_dir, \
generateOutputMatrix_WithNonEqualVectors_WithHeaders

from os.path import join

# ------------------------------------------------------------------------------------------------ #
# Prepare output
outputDir = 'output/Fig-2/'
kd_vs_fMx_FileName = 'Fig-2C.csv'
kd_vs_sep_factor_FileName = 'Fig-2D.csv'
ensure_dir(outputDir)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Calculate the purity of M1 and separation factor as eluant is repeatedly passed through column.

kdBase = 1E-6

# I'm selecting some representative inverse Kds that we selected in panels A and B

# 1: 39.3% Eu binding (wild-type under LH)
# 1.0425: 0.397
# 1.0806: 0.401
# 1.1627: 0.409
# 1.2857: 0.42
# 1.5519: 0.44
# 2.8465: 0.5
# 6.0269: 0.56


inv_kd1_array = [1.00, 1.0425974631137000, 1.0806783602831500, 1.1627844314923600, \
1.2857997031071100, 1.5519282263934800, 2.8465682481811900, 6.026947075485600]

kd2 = 1.2563958631280851 * kdBase
kd3 = 2.28160979282314 * kdBase

nMT = 72e-9
nM1T = nMT/3
nM2T = nMT/3
nM3T = nMT/3
nBT = 35.033e-9
loadVol = 400e-6

scenarioDict = {}

keyFormatter = "{0:.2f}"

i = 0
while i < len(inv_kd1_array):

	kd1 = (1 / inv_kd1_array[i]) * kdBase
	
	scenarioKey = keyFormatter.format(inv_kd1_array[i])
	
	resultsDict = {}

	inputData = \
	InputData_3Metals_1Site(kd1, kd2, kd3, nM1T, nM2T, nM3T, nBT, loadVol, adjustLoadVol=True, \
	adjustBindingSites=True, maxSubCycles=100, adjustBindingSiteTarget='nMT/2') 

	macroCycleDiagnostic = Run_Macro_Cycle_3Metals_1Site(inputData)
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
title('Fig. 2C. 1 Binding Site, Purity of M1 in Eluant vs. Cycle Number')
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

legend()
grid()
xlabel('Cycle Number')
ylabel('Purity of M1 in Eluant')
show()

oMatrix_fM1 = generateOutputMatrix_WithNonEqualVectors_WithHeaders(vectorList_fM1, headers_fM1)
writeOutputMatrix(join(outputDir, kd_vs_fMx_FileName), oMatrix_fM1)


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

title('Fig. 2D. 1 Binding Site, Separation Factor vs. Cycle Number')
legend()
grid()
xlabel('Cycle Number')
ylabel('Separation Factor')
xlim(1,38)
show()

oMatrix_sep = generateOutputMatrix_WithNonEqualVectors_WithHeaders(vectorList_sep, headers_sep)
writeOutputMatrix(join(outputDir, kd_vs_sep_factor_FileName), oMatrix_sep)

# ------------------------------------------------------------------------------------------------ #



