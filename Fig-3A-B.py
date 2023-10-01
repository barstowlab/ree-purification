# ------------------------------------------------------------------------------------------------ #
# Fig-3A-B.py
# Last updated: February 7th, 2023
# Buz Barstow
# Code to calculate effect REE purification with 1 microbe with 3 types of binding site. 
# Figures 3A and B in Genomic Characterization of Rare Earth Binding by Shewanella oneidensis by 
# Medin et al. 
# ------------------------------------------------------------------------------------------------ #


from utils.ConcentrationSolverUtils9 import Calculate_Fraction_of_Each_Metal_Bound_to_Three_Sites,\
Calculate_Binding_Site_Fraction_Arrays

from utils.vectorOutput3 import writeOutputMatrix, generateOutputMatrixWithHeaders, ensure_dir

from os.path import join

from scipy.interpolate import interp1d


# ------------------------------------------------------------------------------------------------ #
# Prepare output
outputDir = 'output/Fig-3/'
fB1_vs_fMx_FileName = 'Fig-3A.csv'
fB1_vs_sep_factor_FileName = 'Fig-3B.csv'
example_FileName = 'Fig-3A-B-examples.csv'
ensure_dir(outputDir)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Calculate the fraction of each metal bound to the binding sites as a function of kd of the site
# for metal 1. I'll start with the kd set 

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

fB1_Array = [0.4618352672074258, 0.475, 0.5, 0.55, 0.6, 0.65, 0.70, 0.75, 0.8, 0.85, 0.9, 0.95]
fB2_0 = 0.3703005831465813
fB3_0 = 0.16786402382522572

fB2_Array, fB3_Array = Calculate_Binding_Site_Fraction_Arrays(fB1_Array, fB2_0, fB3_0)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Initialize arrays for storing calculation results
fM1bArray = []
fM2bArray = []
fM3bArray = []
alpha1_2Array = []
alpha1_3Array = []
alpha2_3Array = []

# Figure out the fraction of binding sites occupied by metals 1, 2, and 3
i = 0
while i < len(fB1_Array):
	
	fB1 = fB1_Array[i]
	fB2 = fB2_Array[i]
	fB3 = fB3_Array[i]
	
	fM1b, fM2b, fM3b, alpha1_2, alpha1_3, alpha2_3 = \
	Calculate_Fraction_of_Each_Metal_Bound_to_Three_Sites(\
	kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
	nM1T, nM2T, nM3T, nBT, fB1, fB2, fB3, loadVol)
	
	fM1bArray.append(fM1b)
	fM2bArray.append(fM2b)
	fM3bArray.append(fM3b)
	
	alpha1_2Array.append(alpha1_2)
	alpha1_3Array.append(alpha1_3)
	alpha2_3Array.append(alpha2_3)
	
	i += 1
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Figure out the binding site fractions for a few example binding percentages

# These examples correspond to the scenarios in Fig-2A-B.py
# fM1bExampleArray = [0.393, 0.397, 0.401, 0.409, 0.42, 0.44, 0.5, 0.53]

# These examples correspond to just the first 4 scenarios in Fig-2A-B.py (the single gene edit ones)
fM1bExampleArray = [0.393, 0.397, 0.401, 0.409]

fB1_ExampleArray = []

fM1b_to_fB1_func = interp1d(fM1bArray, fB1_Array)
fB1_to_fM1b_func = interp1d(fB1_Array, fM1bArray)


example_horiz_line_array_fM1b_plot = []
example_vert_line_array_fM1b_plot = []

# First, get fB1 examples corresponding to the fM1b examples
i = 0
while i < len(fM1bExampleArray):
	fM1bExample = fM1bExampleArray[i]
	fB1 = float(fM1b_to_fB1_func(fM1bExample))
	fB1_ExampleArray.append(fB1)
	print(fM1bExample)
	i += 1
	
# Second, do the opposite. Make some new fB1 examples and figure out the corresponding fM1b
# examples. 
fB1_ExampleArray_part2 = [0.7, 0.8, 0.85, 0.9, 0.95]
i = 0
while i < len(fB1_ExampleArray_part2):
	fB1Example = fB1_ExampleArray_part2[i]
	fM1b = float(fB1_to_fM1b_func(fB1Example))
	fM1bExampleArray.append(fM1b)
	print(fB1Example)
	i += 1
	
fB1_ExampleArray += fB1_ExampleArray_part2



# Now, make marking lines for examples
i = 0
while i < len(fM1bExampleArray):
	
	fB1 = fB1_ExampleArray[i]
	fM1bExample = fM1bExampleArray[i]
	
	horizLine_x_Array = []
	horizLine_y_Array = []
	vertLine_x_Array = []
	vertLine_y_Array = []
	
	vertLine_x_Array.append(fB1)
	vertLine_x_Array.append(fB1)
	vertLine_y_Array.append(0)
	vertLine_y_Array.append(fM1bExample)
	
	horizLine_x_Array.append(0)
	horizLine_x_Array.append(fB1)
	horizLine_y_Array.append(fM1bExample)
	horizLine_y_Array.append(fM1bExample)
	
	example_horiz_line_array_fM1b_plot.append([horizLine_x_Array, horizLine_y_Array])
	example_vert_line_array_fM1b_plot.append([vertLine_x_Array, vertLine_y_Array])
	
	i += 1
	

# Next, mark the corresponding locations on separation factor plot
fB1_to_alpha_func = interp1d(fB1_Array, alpha1_2Array)
example_horiz_line_array_alpha_plot = []
example_vert_line_array_alpha_plot = []
alpha1_2_example_array = []

i = 0 
while i < len(fB1_ExampleArray):
	horizLine_x_Array = []
	horizLine_y_Array = []
	vertLine_x_Array = []
	vertLine_y_Array = []
	
	fB1 = fB1_ExampleArray[i]
	alpha1_2_example = float(fB1_to_alpha_func(fB1))
	alpha1_2_example_array.append(alpha1_2_example)
	
	vertLine_x_Array.append(fB1)
	vertLine_x_Array.append(fB1)
	vertLine_y_Array.append(0)
	vertLine_y_Array.append(alpha1_2_example)
	
	horizLine_x_Array.append(0)
	horizLine_x_Array.append(fB1)
	horizLine_y_Array.append(alpha1_2_example)
	horizLine_y_Array.append(alpha1_2_example)
	
	example_horiz_line_array_alpha_plot.append([horizLine_x_Array, horizLine_y_Array])
	example_vert_line_array_alpha_plot.append([vertLine_x_Array, vertLine_y_Array])
	
	i += 1


# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# Plot out results of calculation on metal binding fractions
figure()
title('Fig. 3A. 3 Sites. Fraction of M1 Bound as M1-site Fraction is Increased')
plot(fB1_Array, fM1bArray, label='M1')
plot(fB1_Array, fM2bArray, label='M2')
plot(fB1_Array, fM3bArray, label='M3')
ylim(0,1)
legend()
grid()
xlabel('Fraction of B1')
ylabel('Fraction of Binding Sites Occupied by Metals')
show()

i = 0
while i < len(example_horiz_line_array_fM1b_plot):
	plot(example_horiz_line_array_fM1b_plot[i][0], example_horiz_line_array_fM1b_plot[i][1], \
	color='black')
	plot(example_vert_line_array_fM1b_plot[i][0], example_vert_line_array_fM1b_plot[i][1], \
	color='black')
	i += 1

xlim(0.46, 0.95)
ylim(0.15, 0.6)
show()






figure()
title('Fig. 3B. 3 Sites. Separation Factors as M1-site Fraction is Increased')
plot(fB1_Array, alpha1_2Array, label='alpha1_2')
plot(fB1_Array, alpha1_3Array, label='alpha1_3')
plot(fB1_Array, alpha2_3Array, label='alpha2_3')
ylim(0, 10)
legend()
grid()
xlabel('Fraction of B1')
ylabel('Separation Factor')
show()

i = 0
while i < len(example_horiz_line_array_alpha_plot):
	plot(example_horiz_line_array_alpha_plot[i][0], example_horiz_line_array_alpha_plot[i][1], \
	color='black')
	plot(example_vert_line_array_alpha_plot[i][0], example_vert_line_array_alpha_plot[i][1], \
	color='black')
	i += 1

xlim(0.46, 0.95)
ylim(0, 8)
show()

# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Output calculation results to a csv file

headers_fB1_vs_fMx = ['fB1', 'fM1', 'fM2', 'fM3']
headers_fB1_vs_sep_factor = ['fB1', 'alpha1_2', 'alpha1_3', 'alpha2_3']

vectorList_fB1_vs_fMx = [fB1_Array, fM1bArray, fM2bArray, fM3bArray]
vectorList_fB1_vs_sep_factor = [fB1_Array, alpha1_2Array, alpha1_3Array, alpha2_3Array]

oMatrix_fB1_vs_fMx = generateOutputMatrixWithHeaders(vectorList_fB1_vs_fMx, headers_fB1_vs_fMx)

oMatrix_fB1_vs_sep_factor = \
generateOutputMatrixWithHeaders(vectorList_fB1_vs_sep_factor, headers_fB1_vs_sep_factor)

writeOutputMatrix(join(outputDir, fB1_vs_fMx_FileName), oMatrix_fB1_vs_fMx)
writeOutputMatrix(join(outputDir, fB1_vs_sep_factor_FileName), oMatrix_fB1_vs_sep_factor)
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# Output examples 
headers_ex = ['fM1b', 'fB1', 'alpha1_2']
i = 0
vectorList_ex = [fM1bExampleArray, fB1_ExampleArray, alpha1_2_example_array]
oMatrix_ex = generateOutputMatrixWithHeaders(vectorList_ex, headers_ex)
writeOutputMatrix(join(outputDir, example_FileName), oMatrix_ex)

# ------------------------------------------------------------------------------------------------ #




