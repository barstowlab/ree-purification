# ------------------------------------------------------------------------------------------------ #
# Fig-2A-B.py
# Last updated: February 7th, 2023
# Buz Barstow
# Code to calculate effect REE binding to a microbe with a single type of binding site. 
# Figures S8A and B in Genomic Characterization of Rare Earth Binding by Shewanella oneidensis by 
# Medin et al. 
# ------------------------------------------------------------------------------------------------ #

from utils.ConcentrationSolverUtils9 import Calculate_Fraction_of_Each_Metal_Bound_to_Single_Site

from utils.vectorOutput3 import writeOutputMatrix, generateOutputMatrixWithHeaders, ensure_dir

from os.path import join

from scipy.interpolate import interp1d

# ------------------------------------------------------------------------------------------------ #
# Prepare output
outputDir = 'output/Fig-2/'
kd_vs_fMx_FileName = 'Fig-2A.csv'
kd_vs_sep_factor_FileName = 'Fig-2B.csv'
example_FileName = 'Fig-2A-B-examples.csv'
ensure_dir(outputDir)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Calculate the fraction of each metal bound to the binding sites as a function of kd of the site
# for metal 1. I'll start with the kd set 

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

fM1bArray = []
fM2bArray = []
fM3bArray = []
alpha1_2Array = []
alpha1_3Array = []
alpha2_3Array = []

# Figure out the fraction of binding sites occupied by metals 1, 2, and 3
i = 0
while i < len(kd1Array):
	
	kd1 = kd1Array[i]
	
	fM1b, fM2b, fM3b, alpha1_2, alpha1_3, alpha2_3 = \
	Calculate_Fraction_of_Each_Metal_Bound_to_Single_Site(\
	kd1, kd2, kd3, nM1T, nM2T, nM3T, nBT, loadVol)
	
	fM1bArray.append(fM1b)
	fM2bArray.append(fM2b)
	fM3bArray.append(fM3b)
	
	alpha1_2Array.append(alpha1_2)
	alpha1_3Array.append(alpha1_3)
	alpha2_3Array.append(alpha2_3)
	
	i += 1
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Figure out the 1/kd for a few example binding percentages
fM1bExampleArray = [0.397, 0.401, 0.409, 0.42, 0.44, 0.5, 0.56]
invKdExampleArray = []

fM1b_to_invKd_func = interp1d(fM1bArray, 1/(kd1Array/kdBase))

example_horiz_line_array_fM1b_plot = []
example_vert_line_array_fM1b_plot = []

i = 0
while i < len(fM1bExampleArray):

	horizLine_x_Array = []
	horizLine_y_Array = []
	vertLine_x_Array = []
	vertLine_y_Array = []

	fM1bExample = fM1bExampleArray[i]
	invKd = float(fM1b_to_invKd_func(fM1bExample))
	invKdExampleArray.append(invKd)
	
	vertLine_x_Array.append(invKd)
	vertLine_x_Array.append(invKd)
	vertLine_y_Array.append(0)
	vertLine_y_Array.append(fM1bExample)
	
	horizLine_x_Array.append(0)
	horizLine_x_Array.append(invKd)
	horizLine_y_Array.append(fM1bExample)
	horizLine_y_Array.append(fM1bExample)
	
	example_horiz_line_array_fM1b_plot.append([horizLine_x_Array, horizLine_y_Array])
	example_vert_line_array_fM1b_plot.append([vertLine_x_Array, vertLine_y_Array])
	
	i += 1


# Next, mark the corresponding locations on separation factor plot
invKd_to_alpha_func = interp1d(1/(kd1Array/kdBase), alpha1_2Array)
example_horiz_line_array_alpha_plot = []
example_vert_line_array_alpha_plot = []
alpha1_2_example_array = []

i = 0 
while i < len(invKdExampleArray):
	horizLine_x_Array = []
	horizLine_y_Array = []
	vertLine_x_Array = []
	vertLine_y_Array = []
	
	invKd = invKdExampleArray[i]
	alpha1_2_example = float(invKd_to_alpha_func(invKd))
	alpha1_2_example_array.append(alpha1_2_example)
	
	vertLine_x_Array.append(invKd)
	vertLine_x_Array.append(invKd)
	vertLine_y_Array.append(0)
	vertLine_y_Array.append(alpha1_2_example)
	
	horizLine_x_Array.append(0)
	horizLine_x_Array.append(invKd)
	horizLine_y_Array.append(alpha1_2_example)
	horizLine_y_Array.append(alpha1_2_example)
	
	example_horiz_line_array_alpha_plot.append([horizLine_x_Array, horizLine_y_Array])
	example_vert_line_array_alpha_plot.append([vertLine_x_Array, vertLine_y_Array])
	
	i += 1


# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
# Plot out results of calculation on metal binding fractions
figure()
title('A: Single Binding Site Type, Effect of Reducing Kd on Binding Fractions')
plot(1/(kd1Array/kdBase), fM1bArray, label='M1')
plot(1/(kd1Array/kdBase), fM2bArray, label='M2')
plot(1/(kd1Array/kdBase), fM3bArray, label='M3')
ylim(0,1)
legend()
grid()
xlabel('kd_0/kd_{new}')
ylabel('Fraction of Binding Sites Occupied by Metals')

i = 0
while i < len(example_horiz_line_array_fM1b_plot):
	plot(example_horiz_line_array_fM1b_plot[i][0], example_horiz_line_array_fM1b_plot[i][1], \
	color='black')
	plot(example_vert_line_array_fM1b_plot[i][0], example_vert_line_array_fM1b_plot[i][1], \
	color='black')
	i += 1

#xlim(1/(kd1Array/kdBase)[0], 1/(kd1Array/kdBase)[-1])
xlim(1, 6.1)
ylim(0.15, 0.6)
show()



figure()
title('B: Single Binding Site Type, Effect of Reducing Kd on Separation Factor')
plot(1/(kd1Array/kdBase), alpha1_2Array, label='alpha1_2')
plot(1/(kd1Array/kdBase), alpha1_3Array, label='alpha1_3')
plot(1/(kd1Array/kdBase), alpha2_3Array, label='alpha2_3')
legend()
grid()
xlabel('kd_0/kd_{new}')
ylabel('Separation Factor')

i = 0
while i < len(example_horiz_line_array_alpha_plot):
	plot(example_horiz_line_array_alpha_plot[i][0], example_horiz_line_array_alpha_plot[i][1], \
	color='black')
	plot(example_vert_line_array_alpha_plot[i][0], example_vert_line_array_alpha_plot[i][1], \
	color='black')
	i += 1

#xlim(1/(kd1Array/kdBase)[0], 1/(kd1Array/kdBase)[-1])
xlim(1, 6.1)
ylim(0, 25)
show()
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Output calculation results to a csv file

headers_kd_vs_fMx = ['1/kd', 'fM1', 'fM2', 'fM3']
headers_kd_vs_sep_factor = ['1/kd', 'alpha1_2', 'alpha1_3', 'alpha2_3']

vectorList_kd_vs_fMx = [1/(kd1Array/kdBase), fM1bArray, fM2bArray, fM3bArray]
vectorList_kd_vs_sep_factor = [1/(kd1Array/kdBase), alpha1_2Array, alpha1_3Array, alpha2_3Array]

oMatrix_kd_vs_fMx = generateOutputMatrixWithHeaders(vectorList_kd_vs_fMx, headers_kd_vs_fMx)

oMatrix_kd_vs_sep_factor = \
generateOutputMatrixWithHeaders(vectorList_kd_vs_sep_factor, headers_kd_vs_sep_factor)

writeOutputMatrix(join(outputDir, kd_vs_fMx_FileName), oMatrix_kd_vs_fMx)
writeOutputMatrix(join(outputDir, kd_vs_sep_factor_FileName), oMatrix_kd_vs_sep_factor)
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# Output examples 
headers_ex = ['fM1b', 'invKd', 'alpha1_2']
i = 0
vectorList_ex = [fM1bExampleArray, invKdExampleArray, alpha1_2_example_array]
oMatrix_ex = generateOutputMatrixWithHeaders(vectorList_ex, headers_ex)
writeOutputMatrix(join(outputDir, example_FileName), oMatrix_ex)

# ------------------------------------------------------------------------------------------------ #