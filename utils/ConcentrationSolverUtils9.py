# Unknown variables in competitive binding system
# cM1f = free concentration of metal 1 in Mol per Liter
# cM2f = free concentration of metal 2 in Mol per Liter
# cM1b = bound concentration of metal 1 in Mol per Liter
# cM2b = bound concentration of metal 2 in Mol per Liter
# cBf = free concentration of binding site in Mol per Liter
# cBb = bound concentration of binding site in Mol per Liter

# Known (or at least specified) variables in competitive binding system
# kd1 = dissociation constant for binding site and metal 1 in inverse Molar
# kd2 = dissociation constant for binding site and metal 2 in inverse Molar
# vol = system volume in Liter
# nM1T = total number of M1 atoms in system (in moles)
# nM2T = total number of M2 atoms in system (in moles)
# nBT  = total number of binding sites in system (in moles)



# ------------------------------------------------------------------------------------------------ #
# This calculates binding of 3 metals to 3 binding sites, each with a different dissociation 
# constant for the 3 metals. 
class BindState_3Metals_3Sites:
	def __init__(self, loadState):
		import pdb
		from sympy import nsolve, Symbol
				
		# Dissociation constant for 1st site, for 1st metal (1st number is site, 2nd is metal)
		# 1. kd1_1 = cM1f * cB1f / cM1b1
		# 2. kd1_2 = cM2f * cB1f / cM2b1		
		# 3. kd1_3 = cM3f * cB1f / cM3b1
		
		# Dissociation constant for 2nd site, for 1st metal
		# 4. kd2_1 = cM1f * cB2f / cM1b2
		# 5. kd2_2 = cM2f * cB2f / cM2b2
		# 6. kd2_3 = cM3f * cB2f / cM3b2
		# 7. kd3_1 = cM1f * cB3f / cM1b3
		# 8. kd3_2 = cM2f * cB3f / cM2b3
		# 9. kd3_3 = cM3f * cB3f / cM3b3
		# 10. nM1T = vol * (cM1f + cM1b1 + cM1b2 + cM1b3)
		# 11. nM2T = vol * (cM2f + cM2b1 + cM2b2 + cM2b3)
		# 12. nM3T = vol * (cM3f + cM3b1 + cM3b2 + cM3b3)
		# 13. nB1T = vol * (cB1f + cM1b1 + cM2b1 + cM3b1)
		# 14. nB2T = vol * (cB2f + cM1b2 + cM2b2 + cM3b2)
		# 15. nB3T = vol * (cB3f + cM1b3 + cM2b3 + cM3b3)
		
		# Concentrations of free metals
		cM1f = Symbol('cM1f', positive=True)
		cM2f = Symbol('cM2f', positive=True)
		cM3f = Symbol('cM3f', positive=True)
		
		# Concentrations of metals bound to each of the binding sites
		cM1b1 = Symbol('cM1b1', positive=True)
		cM2b1 = Symbol('cM2b1', positive=True)
		cM3b1 = Symbol('cM3b1', positive=True)
		
		cM1b2 = Symbol('cM1b2', positive=True)
		cM2b2 = Symbol('cM2b2', positive=True)
		cM3b2 = Symbol('cM3b2', positive=True)
		
		cM1b3 = Symbol('cM1b3', positive=True)
		cM2b3 = Symbol('cM2b3', positive=True)
		cM3b3 = Symbol('cM3b3', positive=True)
		
		cB1f = Symbol('cB1f', positive=True)
		cB2f = Symbol('cB2f', positive=True)
		cB3f = Symbol('cB3f', positive=True)
		
		cB1b = Symbol('cB1b', positive=True)
		cB2b = Symbol('cB2b', positive=True)
		cB3b = Symbol('cB3b', positive=True)
		
		
		
		# Dissociation constants for each of the binding sites for each of the metals
		# Site 1
		kd1_1 = loadState.kd1_1
		kd1_2 = loadState.kd1_2
		kd1_3 = loadState.kd1_3
		
		# Site 2
		kd2_1 = loadState.kd2_1
		kd2_2 = loadState.kd2_2
		kd2_3 = loadState.kd2_3
		
		# Site 3
		kd3_1 = loadState.kd3_1
		kd3_2 = loadState.kd3_2
		kd3_3 = loadState.kd3_3
		
		# Total number of each metal
		nM1T = loadState.nM1
		nM2T = loadState.nM2
		nM3T = loadState.nM3
		
		# System loading volume
		vol = loadState.vol

		# Number of binding sites
		nB1T = loadState.nB1T
		nB2T = loadState.nB2T
		nB3T = loadState.nB3T
		nBT = nB1T + nB2T + nB3T
		
		# Fractions of each metal in load solution
		fM1 = loadState.fM1
		fM2 = loadState.fM2
		fM3 = loadState.fM3
		
		# Initial guesses for values
		cM1fGuess = (nM1T - (0.99 * fM1 * nBT))/vol 
		cM2fGuess = (nM2T - (0.99 * fM2 * nBT))/vol 
		cM3fGuess = (nM3T - (0.99 * fM3 * nBT))/vol 
		
		cM1b1Guess = (0.99 * fM1 * nB1T)/vol 
		cM2b1Guess = (0.99 * fM2 * nB1T)/vol 
		cM3b1Guess = (0.99 * fM3 * nB1T)/vol 
		
		cM1b2Guess = (0.99 * fM1 * nB2T)/vol 
		cM2b2Guess = (0.99 * fM2 * nB2T)/vol 
		cM3b2Guess = (0.99 * fM3 * nB2T)/vol 
		
		cM1b3Guess = (0.99 * fM1 * nB3T)/vol 
		cM2b3Guess = (0.99 * fM2 * nB3T)/vol 
		cM3b3Guess = (0.99 * fM3 * nB3T)/vol 
		
		cB1fGuess = 0.01*nB1T/vol
		cB2fGuess = 0.01*nB2T/vol
		cB3fGuess = 0.01*nB3T/vol
		
		cB1bGuess = (0.99 * nB1T)/vol
		cB2bGuess = (0.99 * nB2T)/vol
		cB3bGuess = (0.99 * nB3T)/vol
					
			
		try:
			solutions = nsolve(\
			[\
			kd1_1 - (cM1f * cB1f / cM1b1) , \
			kd1_2 - (cM2f * cB1f / cM2b1) , \
			kd1_3 - (cM3f * cB1f / cM3b1) , \
			\
			kd2_1 - (cM1f * cB2f / cM1b2) , \
			kd2_2 - (cM2f * cB2f / cM2b2) , \
			kd2_3 - (cM3f * cB2f / cM3b2) , \
			\
			kd3_1 - (cM1f * cB3f / cM1b3) , \
			kd3_2 - (cM2f * cB3f / cM2b3) , \
			kd3_3 - (cM3f * cB3f / cM3b3) , \
			\
			nM1T - vol * (cM1f + cM1b1 + cM1b2 + cM1b3) , \
			nM2T - vol * (cM2f + cM2b1 + cM2b2 + cM2b3) , \
			nM3T - vol * (cM3f + cM3b1 + cM3b2 + cM3b3) , \
			\
			nB1T - vol * (cB1f + cM1b1 + cM2b1 + cM3b1) , \
			nB2T - vol * (cB2f + cM1b2 + cM2b2 + cM3b2) , \
			nB3T - vol * (cB3f + cM1b3 + cM2b3 + cM3b3) , \
			\
			cB1b - cM1b1 - cM2b1 - cM3b1, \
			cB2b - cM1b2 - cM2b2 - cM3b2, \
			cB3b - cM1b3 - cM2b3 - cM3b3 \
			],\
			\
			[cM1f, cM2f, cM3f, cM1b1, cM2b1, cM3b1, cM1b2, cM2b2, cM3b2, cM1b3, cM2b3, cM3b3, \
			cB1f, cB2f, cB3f, cB1b, cB2b, cB3b], \
			\
			(cM1fGuess, cM2fGuess, cM3fGuess, cM1b1Guess, cM2b1Guess, cM3b1Guess, cM1b2Guess, \
			cM2b2Guess, cM3b2Guess, cM1b3Guess, cM2b3Guess, cM3b3Guess, cB1fGuess, cB2fGuess, \
			cB3fGuess, cB1bGuess, cB2bGuess, cB3bGuess), \
			\
			dict=True, force=True, minimal=True)
		except:
			pdb.set_trace()
		
		
		self.cM1f = solutions[0][cM1f]
		self.cM2f = solutions[0][cM2f]
		self.cM3f = solutions[0][cM3f] 
		self.cM1b1 = solutions[0][cM1b1]
		self.cM2b1 = solutions[0][cM2b1]
		self.cM3b1 = solutions[0][cM3b1]
		self.cM1b2 = solutions[0][cM1b2]
		self.cM2b2 = solutions[0][cM2b2]
		self.cM3b2 = solutions[0][cM3b2]
		self.cM1b3 = solutions[0][cM1b3]
		self.cM2b3 = solutions[0][cM2b3]
		self.cM3b3 = solutions[0][cM3b3]
		self.cB1f = solutions[0][cB1f]
		self.cB2f = solutions[0][cB2f]
		self.cB3f = solutions[0][cB3f]
		self.cB1b = solutions[0][cB1b]
		self.cB2b = solutions[0][cB2b]
		self.cB3b = solutions[0][cB3b]
			
		self.vol = vol
		
		self.nM1f = self.cM1f * self.vol
		self.nM2f = self.cM2f * self.vol
		self.nM3f = self.cM3f * self.vol
		
		self.cM1b = (self.cM1b1 + self.cM1b2 + self.cM1b3)
		self.cM2b = (self.cM2b1 + self.cM2b2 + self.cM2b3)
		self.cM3b = (self.cM3b1 + self.cM3b2 + self.cM3b3)
		
		self.nM1b = self.cM1b * vol
		self.nM2b = self.cM2b * vol
		self.nM3b = self.cM3b * vol
		
		totalMetalsF = (self.nM1f + self.nM2f + self.nM3f)
		totalMetalsB = (self.nM1b + self.nM2b + self.nM3b)

		self.fM1f = self.nM1f / totalMetalsF
		self.fM2f = self.nM2f / totalMetalsF
		self.fM3f = self.nM3f / totalMetalsF
		
		self.fM1b = self.nM1b / totalMetalsB
		self.fM2b = self.nM2b / totalMetalsB
		self.fM3b = self.nM3b / totalMetalsB
		
		# Calculate the distribution coefficients
		# Distribution coefficient for M1: the free concentration of M1 / bound concentration
		self.DM1 = self.cM1f / self.cM1b
		self.DM2 = self.cM2f / self.cM2b
		self.DM3 = self.cM3f / self.cM3b
		
		# Calculate separation factors
		# If they are less than 1, invert them
		self.alpha1_2 = self.DM1 / self.DM2
		self.alpha1_3 = self.DM1 / self.DM3
		self.alpha2_3 = self.DM2 / self.DM3
		
		if self.alpha1_2 < 1:
			self.alpha1_2 = 1/self.alpha1_2
		
		if self.alpha1_3 < 1:
			self.alpha1_3 = 1/self.alpha1_3
			
		if self.alpha2_3 < 1:
			self.alpha2_3 = 1/self.alpha2_3
		
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# This calculates binding of 3 metals to 3 binding sites, each with a different dissociation 
# constant for the 3 metals. 
class BindState_3Metals_2Sites:
	def __init__(self, loadState):
		import pdb
		from sympy import nsolve, Symbol
				
		# Concentrations of free metals
		cM1f = Symbol('cM1f', positive=True)
		cM2f = Symbol('cM2f', positive=True)
		cM3f = Symbol('cM3f', positive=True)
		
		# Concentrations of metals bound to each of the binding sites
		cM1b1 = Symbol('cM1b1', positive=True)
		cM2b1 = Symbol('cM2b1', positive=True)
		cM3b1 = Symbol('cM3b1', positive=True)
		cM1b2 = Symbol('cM1b2', positive=True)
		cM2b2 = Symbol('cM2b2', positive=True)
		cM3b2 = Symbol('cM3b2', positive=True)
		
		cB1f = Symbol('cB1f', positive=True)
		cB2f = Symbol('cB2f', positive=True)
		cB1b = Symbol('cB1b', positive=True)
		cB2b = Symbol('cB2b', positive=True)
		
		# Dissociation constants for each of the binding sites for each of the metals
		# Site 1
		kd1_1 = loadState.kd1_1
		kd1_2 = loadState.kd1_2
		kd1_3 = loadState.kd1_3
		
		# Site 2
		kd2_1 = loadState.kd2_1
		kd2_2 = loadState.kd2_2
		kd2_3 = loadState.kd2_3
		
		# Total number of each metal
		nM1T = loadState.nM1
		nM2T = loadState.nM2
		nM3T = loadState.nM3
		
		# System loading volume
		vol = loadState.vol

		# Number of binding sites
		nB1T = loadState.nB1T
		nB2T = loadState.nB2T
		nBT = nB1T + nB2T
		
		# Fractions of each metal in load solution
		fM1 = loadState.fM1
		fM2 = loadState.fM2
		fM3 = loadState.fM3
		
		# Initial guesses for values
		cM1fGuess = (nM1T - (0.99 * fM1 * nBT))/vol 
		cM2fGuess = (nM2T - (0.99 * fM2 * nBT))/vol 
		cM3fGuess = (nM3T - (0.99 * fM3 * nBT))/vol 
		
		cM1b1Guess = (0.99 * fM1 * nB1T)/vol 
		cM2b1Guess = (0.99 * fM2 * nB1T)/vol 
		cM3b1Guess = (0.99 * fM3 * nB1T)/vol 
		
		cM1b2Guess = (0.99 * fM1 * nB2T)/vol 
		cM2b2Guess = (0.99 * fM2 * nB2T)/vol 
		cM3b2Guess = (0.99 * fM3 * nB2T)/vol 
		
		cB1fGuess = 0.01*nB1T/vol
		cB2fGuess = 0.01*nB2T/vol
		
		cB1bGuess = 0.99*nB1T/vol
		cB2bGuess = 0.99*nB2T/vol
		
					
			
		try:
			solutions = nsolve(\
			[\
			kd1_1 - (cM1f * cB1f / cM1b1), \
			kd1_2 - (cM2f * cB1f / cM2b1), \
			kd1_3 - (cM3f * cB1f / cM3b1), \
			\
			kd2_1 - (cM1f * cB2f / cM1b2), \
			kd2_2 - (cM2f * cB2f / cM2b2), \
			kd2_3 - (cM3f * cB2f / cM3b2), \
			\
			nM1T - vol * (cM1f + cM1b1 + cM1b2), \
			nM2T - vol * (cM2f + cM2b1 + cM2b2), \
			nM3T - vol * (cM3f + cM3b1 + cM3b2), \
			\
			nB1T - vol * (cB1f + cM1b1 + cM2b1 + cM3b1), \
			nB2T - vol * (cB2f + cM1b2 + cM2b2 + cM3b2), \
			\
			cB1b - cM1b1 - cM2b1 - cM3b1, \
			cB2b - cM1b2 - cM2b2 - cM3b2 \
			],\
			\
			[cM1f, cM2f, cM3f, \
			cM1b1, cM2b1, cM3b1, cM1b2, cM2b2, cM3b2, \
			cB1f, cB2f, cB1b, cB2b], \
			\
			(cM1fGuess, cM2fGuess, cM3fGuess, \
			cM1b1Guess, cM2b1Guess, cM3b1Guess, cM1b2Guess, cM2b2Guess, cM3b2Guess, \
			cB1fGuess, cB2fGuess, cB1bGuess, cB2bGuess), \
			\
			dict=True)
		except:
			pdb.set_trace()
		
		
		self.cM1f = solutions[0][cM1f]
		self.cM2f = solutions[0][cM2f]
		self.cM3f = solutions[0][cM3f]
		
		self.cM1b1 = solutions[0][cM1b1]
		self.cM2b1 = solutions[0][cM2b1]
		self.cM3b1 = solutions[0][cM3b1]
		self.cM1b2 = solutions[0][cM1b2]
		self.cM2b2 = solutions[0][cM2b2]
		self.cM3b2 = solutions[0][cM3b2]
		
		self.cB1f = solutions[0][cB1f]
		self.cB2f = solutions[0][cB2f]
		self.cB1b = solutions[0][cB1b]
		self.cB2b = solutions[0][cB2b]
			
		self.vol = vol
		
		self.nM1f = self.cM1f * self.vol
		self.nM2f = self.cM2f * self.vol
		self.nM3f = self.cM3f * self.vol
		
		self.cM1b = (self.cM1b1 + self.cM1b2)
		self.cM2b = (self.cM2b1 + self.cM2b2)
		self.cM3b = (self.cM3b1 + self.cM3b2)
		
		self.nM1b = self.cM1b * vol
		self.nM2b = self.cM2b * vol
		self.nM3b = self.cM3b * vol
		
		totalMetalsF = (self.nM1f + self.nM2f + self.nM3f)
		totalMetalsB = (self.nM1b + self.nM2b + self.nM3b)

		self.fM1f = self.nM1f / totalMetalsF
		self.fM2f = self.nM2f / totalMetalsF
		self.fM3f = self.nM3f / totalMetalsF
		
		self.fM1b = self.nM1b / totalMetalsB
		self.fM2b = self.nM2b / totalMetalsB
		self.fM3b = self.nM3b / totalMetalsB
		
		# Calculate the distribution coefficients
		# Distribution coefficient for M1: the free concentration of M1 / bound concentration
		self.DM1 = self.cM1f / self.cM1b
		self.DM2 = self.cM2f / self.cM2b
		self.DM3 = self.cM3f / self.cM3b
		
		# Calculate separation factors
		# If they are less than 1, invert them
		self.alpha1_2 = self.DM1 / self.DM2
		self.alpha1_3 = self.DM1 / self.DM3
		self.alpha2_3 = self.DM2 / self.DM3
		
		if self.alpha1_2 < 1:
			self.alpha1_2 = 1/self.alpha1_2
		
		if self.alpha1_3 < 1:
			self.alpha1_3 = 1/self.alpha1_3
			
		if self.alpha2_3 < 1:
			self.alpha2_3 = 1/self.alpha2_3
		
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
class BindState_3Metals_1Site:
	def __init__(self, loadState):
		import pdb
		from sympy import nsolve, Symbol

		cM1f = Symbol('cM1f', positive=True)
		cM2f = Symbol('cM2f', positive=True)
		cM3f = Symbol('cM3f', positive=True)	
		cM1b = Symbol('cM1b', positive=True)
		cM2b = Symbol('cM2b', positive=True)
		cM3b = Symbol('cM3b', positive=True)
		cBf = Symbol('cBf', positive=True)
		cBb = Symbol('cBb', positive=True)
	
		kd1 = loadState.kd1
		kd2 = loadState.kd2
		kd3 = loadState.kd3
		nM1T = loadState.nM1
		nM2T = loadState.nM2
		nM3T = loadState.nM3
		nBT = loadState.nBT
		vol = loadState.vol
		
		cM1fGuess = 0.5*nM1T/vol
		cBfGuess = 0.01*nBT/vol
		cM1bGuess = 0.5*nM1T/vol
		cM2fGuess = 0.5*nM2T/vol
		cM2bGuess = 0.5*nM2T/vol
		cBbGuess = 0.99*nBT/vol
		cM3fGuess = 0.5*nM3T/vol
		cM3bGuess = 0.5*nM3T/vol

		try:
			solutions = nsolve(\
			[\
			kd1 - ((cM1f * cBf) / cM1b), \
			kd2 - ((cM2f * cBf) / cM2b), \
			kd3 - ((cM3f * cBf) / cM3b), \
			\
			nM1T - (vol * (cM1f + cM1b)), \
			nM2T - (vol * (cM2f + cM2b)), \
			nM3T - (vol * (cM3f + cM3b)), \
			\
			nBT - (vol * (cBf + cBb)), \
			cBb - cM1b - cM2b - cM3b \
			],\
			[cM1f, cBf, cM1b, cM2f, cM2b, cBb, cM3f, cM3b], \
			(cM1fGuess, cBfGuess, cM1bGuess, cM2fGuess, cM2bGuess, cBbGuess, cM3fGuess, \
			cM3bGuess),\
			dict=True)
		except:
			pdb.set_trace()
		
		self.cM1f = solutions[0][cM1f]
		self.cBf = solutions[0][cBf]
		self.cM1b = solutions[0][cM1b]
		self.cM2f = solutions[0][cM2f]
		self.cM2b = solutions[0][cM2b]
		self.cBb = solutions[0][cBb]
		self.cM3f = solutions[0][cM3f]
		self.cM3b = solutions[0][cM3b]
		
		self.vol = vol
		
		self.nM1f = self.cM1f * self.vol
		self.nM2f = self.cM2f * self.vol
		self.nM3f = self.cM3f * self.vol
		
		self.nM1b = self.cM1b * self.vol
		self.nM2b = self.cM2b * self.vol
		self.nM3b = self.cM3b * self.vol
		
		totalMetalsF = (self.nM1f + self.nM2f + self.nM3f)
		totalMetalsB = (self.nM1b + self.nM2b + self.nM3b)

		self.fM1f = self.nM1f / totalMetalsF
		self.fM2f = self.nM2f / totalMetalsF
		self.fM3f = self.nM3f / totalMetalsF
		
		self.fM1b = self.nM1b / totalMetalsB
		self.fM2b = self.nM2b / totalMetalsB
		self.fM3b = self.nM3b / totalMetalsB
		
		# Calculate the distribution coefficients
		# Distribution coefficient for M1: the free concentration of M1 / bound concentration
		self.DM1 = self.cM1f / self.cM1b
		self.DM2 = self.cM2f / self.cM2b
		self.DM3 = self.cM3f / self.cM3b
		
		# Calculate separation factors
		# If they are less than 1, invert them
		self.alpha1_2 = self.DM1 / self.DM2
		self.alpha1_3 = self.DM1 / self.DM3
		self.alpha2_3 = self.DM2 / self.DM3
		
		if self.alpha1_2 < 1:
			self.alpha1_2 = 1/self.alpha1_2
		
		if self.alpha1_3 < 1:
			self.alpha1_3 = 1/self.alpha1_3
			
		if self.alpha2_3 < 1:
			self.alpha2_3 = 1/self.alpha2_3
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
class EluteSolution_3Metals:
	def __init__(self):
		self.cM1 = 0
		self.cM2 = 0
		self.cM3 = 0
		
		self.fM1 = 0
		self.fM2 = 0
		self.fM3 = 0
		
		self.vol = 0
		
		self.nM1 = 0
		self.nM2 = 0
		self.nM3 = 0
		
		self.fullyInitialized = False
		
	def init_with_bind_state(self, bindState, eluteVol):
		self.cM1 = bindState.cM1b * (bindState.vol / eluteVol)
		self.cM2 = bindState.cM2b * (bindState.vol / eluteVol)
		self.cM3 = bindState.cM3b * (bindState.vol / eluteVol)
		
		self.fM1 = bindState.fM1b
		self.fM2 = bindState.fM2b
		self.fM3 = bindState.fM3b
		
		self.vol = eluteVol
		
		self.nM1 = bindState.nM1b
		self.nM2 = bindState.nM2b
		self.nM3 = bindState.nM3b
		
		self.fullyInitialized = True
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
class WashSolution_3Metals:
	def __init__(self):
		self.cM1 = 0
		self.cM2 = 0
		self.cM3 = 0
		
		self.fM1 = 0
		self.fM2 = 0
		self.fM3 = 0
		
		self.vol = 0
		
		self.nM1 = 0
		self.nM2 = 0
		self.nM3 = 0
		
		self.fullyInitialized = False
	
	def init_with_bind_state(self, bindState):
		self.cM1 = bindState.cM1f
		self.cM2 = bindState.cM2f
		self.cM3 = bindState.cM3f
		
		self.fM1 = bindState.fM1f
		self.fM2 = bindState.fM2f
		self.fM3 = bindState.fM3f		
		
		self.vol = bindState.vol
		
		self.nM1 = bindState.nM1f
		self.nM2 = bindState.nM2f
		self.nM3 = bindState.nM3f
	
		self.fullyInitialized = True
	
	def update_concentrations_fractions_with_metal_numbers_and_vol(self):
		self.cM1 = self.nM1 / self.vol
		self.cM2 = self.nM2 / self.vol
		self.cM3 = self.nM3 / self.vol
		
		totalMetals = self.nM1 + self.nM2 + self.nM3
		
		self.fM1 = self.nM1 / totalMetals
		self.fM2 = self.nM2 / totalMetals
		self.fM3 = self.nM3 / totalMetals
		
		self.fullyInitialized = True
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
class LoadState_3Metals_3Sites:
	def __init__(self):
		self.vol = 0
		self.nM1 = 0
		self.nM2 = 0
		self.nM3 = 0
		
		self.kd1_1 = 0
		self.kd1_2 = 0
		self.kd1_3 = 0
		
		self.kd2_1 = 0
		self.kd2_2 = 0
		self.kd2_3 = 0
		
		self.kd3_1 = 0
		self.kd3_2 = 0
		self.kd3_3 = 0
		
		self.nBT = 0
		
		self.fM1 = 0
		self.fM2 = 0
		self.fM3 = 0
		self.fullyInitialized = False

	def init_with_elute_or_wash_solution(self, eluteOrWashSolution, \
	kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
	nB1T, nB2T, nB3T, loadVol):
		
		self.vol = loadVol
		self.nM1 = eluteOrWashSolution.nM1
		self.nM2 = eluteOrWashSolution.nM2
		self.nM3 = eluteOrWashSolution.nM3
		
		# Dissociation constants for each of the binding sites for each of the metals
		# Site 1
		self.kd1_1 = kd1_1
		self.kd1_2 = kd1_2
		self.kd1_3 = kd1_3
		
		# Site 2
		self.kd2_1 = kd2_1
		self.kd2_2 = kd2_2
		self.kd2_3 = kd2_3
		
		# Site 3
		self.kd3_1 = kd3_1
		self.kd3_2 = kd3_2
		self.kd3_3 = kd3_3

		self.nB1T = nB1T
		self.nB2T = nB2T
		self.nB3T = nB3T
		
		self.nBT = nB1T + nB2T + nB3T
		
		totalMetals = self.nM1 + self.nM2 + self.nM3
		
		self.fM1 = self.nM1 / totalMetals
		self.fM2 = self.nM2 / totalMetals
		self.fM3 = self.nM3 / totalMetals
		
		self.fullyInitialized = True
		
	def init_with_elute_or_wash_solution_v2(self, eluteOrWashSolution, \
	inputData, nB1T, nB2T, nB3T, loadVol):
		
		self.vol = loadVol
		self.nM1 = eluteOrWashSolution.nM1
		self.nM2 = eluteOrWashSolution.nM2
		self.nM3 = eluteOrWashSolution.nM3
		
		# Dissociation constants for each of the binding sites for each of the metals
		# Site 1
		self.kd1_1 = inputData.kd1_1
		self.kd1_2 = inputData.kd1_2
		self.kd1_3 = inputData.kd1_3
		
		# Site 2
		self.kd2_1 = inputData.kd2_1
		self.kd2_2 = inputData.kd2_2
		self.kd2_3 = inputData.kd2_3
		
		# Site 3
		self.kd3_1 = inputData.kd3_1
		self.kd3_2 = inputData.kd3_2
		self.kd3_3 = inputData.kd3_3

		self.nB1T = nB1T
		self.nB2T = nB2T
		self.nB3T = nB3T
		
		self.nBT = nB1T + nB2T + nB3T
		
		totalMetals = self.nM1 + self.nM2 + self.nM3
		
		self.fM1 = self.nM1 / totalMetals
		self.fM2 = self.nM2 / totalMetals
		self.fM3 = self.nM3 / totalMetals
		
		self.fullyInitialized = True
		
		
	def init_with_input_data(self, inputData):
		
		self.vol = inputData.initial_loadVol
		
		self.nM1 = inputData.initial_nM1T
		self.nM2 = inputData.initial_nM2T
		self.nM3 = inputData.initial_nM3T
		
		self.kd1_1 = inputData.kd1_1
		self.kd1_2 = inputData.kd1_2
		self.kd1_3 = inputData.kd1_3
		
		self.kd2_1 = inputData.kd2_1
		self.kd2_2 = inputData.kd2_2
		self.kd2_3 = inputData.kd2_3
		
		self.kd3_1 = inputData.kd3_1
		self.kd3_2 = inputData.kd3_2
		self.kd3_3 = inputData.kd3_3
		
		self.nB1T = inputData.initial_nB1T
		self.nB2T = inputData.initial_nB2T
		self.nB3T = inputData.initial_nB3T
		
		self.nBT = inputData.initial_nBT
		
		totalMetals = self.nM1 + self.nM2 + self.nM3
		
		self.fM1 = self.nM1 / totalMetals
		self.fM2 = self.nM2 / totalMetals
		self.fM3 = self.nM3 / totalMetals
		self.fullyInitialized = True
	
		
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
class LoadState_3Metals_2Sites:
	def __init__(self):
		self.vol = 0
		self.nM1 = 0
		self.nM2 = 0
		self.nM3 = 0
		
		self.kd1_1 = 0
		self.kd1_2 = 0
		self.kd1_3 = 0
		
		self.kd2_1 = 0
		self.kd2_2 = 0
		self.kd2_3 = 0
		
		self.nBT = 0
		
		self.fM1 = 0
		self.fM2 = 0
		self.fM3 = 0
		self.fullyInitialized = False

	def init_with_elute_or_wash_solution(self, eluteOrWashSolution, \
	kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, nB1T, nB2T, loadVol):
		
		self.vol = loadVol
		self.nM1 = eluteOrWashSolution.nM1
		self.nM2 = eluteOrWashSolution.nM2
		self.nM3 = eluteOrWashSolution.nM3
		
		# Dissociation constants for each of the binding sites for each of the metals
		# Site 1
		self.kd1_1 = kd1_1
		self.kd1_2 = kd1_2
		self.kd1_3 = kd1_3
		
		# Site 2
		self.kd2_1 = kd2_1
		self.kd2_2 = kd2_2
		self.kd2_3 = kd2_3

		self.nB1T = nB1T
		self.nB2T = nB2T
		
		self.nBT = nB1T + nB2T
		
		totalMetals = self.nM1 + self.nM2 + self.nM3
		
		self.fM1 = self.nM1 / totalMetals
		self.fM2 = self.nM2 / totalMetals
		self.fM3 = self.nM3 / totalMetals
		
		self.fullyInitialized = True
		
	def init_with_input_data(self, inputData):
		
		self.vol = inputData.initial_loadVol
		
		self.nM1 = inputData.initial_nM1T
		self.nM2 = inputData.initial_nM2T
		self.nM3 = inputData.initial_nM3T
		
		self.kd1_1 = inputData.kd1_1
		self.kd1_2 = inputData.kd1_2
		self.kd1_3 = inputData.kd1_3
		
		self.kd2_1 = inputData.kd2_1
		self.kd2_2 = inputData.kd2_2
		self.kd2_3 = inputData.kd2_3
		
		self.nB1T = inputData.initial_nB1T
		self.nB2T = inputData.initial_nB2T
		
		self.nBT = inputData.initial_nBT
		
		totalMetals = self.nM1 + self.nM2 + self.nM3
		
		self.fM1 = self.nM1 / totalMetals
		self.fM2 = self.nM2 / totalMetals
		self.fM3 = self.nM3 / totalMetals
		self.fullyInitialized = True
	
		
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
class LoadState_3Metals_1Site:
	def __init__(self):
		self.vol = 0
		self.nM1 = 0
		self.nM2 = 0
		self.nM3 = 0
		self.kd1 = 0
		self.kd2 = 0
		self.kd3 = 0
		self.fM1 = 0
		self.fM2 = 0
		self.fM3 = 0
		self.fullyInitialized = False

	def init_with_elute_or_wash_solution(self, eluteOrWashSolution, kd1, kd2, kd3, nBT, loadVol):
		self.vol = loadVol
		self.nM1 = eluteOrWashSolution.nM1
		self.nM2 = eluteOrWashSolution.nM2
		self.nM3 = eluteOrWashSolution.nM3
		self.kd1 = kd1
		self.kd2 = kd2
		self.kd3 = kd3
		self.nBT = nBT
		
		totalMetals = self.nM1 + self.nM2 + self.nM3
		
		self.fM1 = self.nM1 / totalMetals
		self.fM2 = self.nM2 / totalMetals
		self.fM3 = self.nM3 / totalMetals
		
		self.fullyInitialized = True
		
	def init_with_input_data(self, inputData):
		self.vol = inputData.initial_loadVol
		self.nM1 = inputData.initial_nM1T
		self.nM2 = inputData.initial_nM2T
		self.nM3 = inputData.initial_nM3T
		self.kd1 = inputData.kd1
		self.kd2 = inputData.kd2
		self.kd3 = inputData.kd3
		self.nBT = inputData.initial_nBT
		
		totalMetals = self.nM1 + self.nM2 + self.nM3
		
		self.fM1 = self.nM1 / totalMetals
		self.fM2 = self.nM2 / totalMetals
		self.fM3 = self.nM3 / totalMetals
		self.fullyInitialized = True
	
		
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
class InputData_3Metals_3Sites:
	def __init__(self, \
	kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
	initial_nM1T, initial_nM2T, initial_nM3T, \
	initial_nBT, fB1, fB2, fB3, \
	initial_loadVol, \
	initial_eluteVol = 1.0, \
	targetPurity=0.9, preferred_operating_concentration=0.1, \
	lowest_allowed_vol=0.001, \
	printFigures=False, printDiagnostics=True, maxSubCycles=20, maxMacroCycles=100, \
	printPurity=False, adjustLoadVol=True, adjustBindingSites=True):
		
		import pdb
		
		self.printFigures = printFigures
		self.printDiagnostics = printDiagnostics
		self.maxSubCycles = maxSubCycles
		self.maxMacroCycles = maxMacroCycles
		self.printPurity = printPurity
		self.adjustLoadVol = adjustLoadVol
		self.adjustBindingSites = adjustBindingSites
		
		# Dissociation constants for each of the binding sites for each of the metals
		# Site 1
		self.kd1_1 = kd1_1
		self.kd1_2 = kd1_2
		self.kd1_3 = kd1_3
		
		# Site 2
		self.kd2_1 = kd2_1
		self.kd2_2 = kd2_2
		self.kd2_3 = kd2_3
		
		# Site 3
		self.kd3_1 = kd3_1
		self.kd3_2 = kd3_2
		self.kd3_3 = kd3_3

		
		self.initial_nM1T = initial_nM1T
		self.initial_nM2T = initial_nM2T
		self.initial_nM3T = initial_nM3T
		
		self.initial_nBT = initial_nBT
		self.fB1 = fB1
		self.fB2 = fB2
		self.fB3 = fB3
		
		self.initial_nB1T = self.initial_nBT * self.fB1
		self.initial_nB2T = self.initial_nBT * self.fB2
		self.initial_nB3T = self.initial_nBT * self.fB3
		
		self.initial_loadVol = initial_loadVol
		self.initial_eluteVol = initial_eluteVol
		
		self.lowest_allowed_vol = lowest_allowed_vol
		
		self.targetPurity = targetPurity
		
		self.preferred_operating_concentration = preferred_operating_concentration
		
		
	
	
	def update_with_elute_or_wash_solution(self, eluteOrWashSolution):
		
		import pdb
		
		self.initial_nM1T = eluteOrWashSolution.nM1
		self.initial_nM2T = eluteOrWashSolution.nM2
		self.initial_nM3T = eluteOrWashSolution.nM3
		
		self.initial_nBT, self.initial_nB1T, self.initial_nB2T, self.initial_nB3T, \
		self.initial_loadVol = \
		Calculate_New_Binding_Site_Number_and_Load_Volume(eluteOrWashSolution, self)
		
		
		
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
class InputData_3Metals_2Sites:
	def __init__(self, \
	kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, \
	initial_nM1T, initial_nM2T, initial_nM3T, \
	initial_nBT, fB1, fB2, \
	initial_loadVol, 
	initial_eluteVol = 1.0, \
	targetPurity=0.9, preferred_operating_concentration=0.1, \
	lowest_allowed_vol=0.001, \
	printFigures=False, printDiagnostics=True, maxSubCycles=20, maxMacroCycles=100, \
	printPurity=False, adjustLoadVol=True, adjustBindingSites=True):
		
		import pdb
		
		self.printFigures = printFigures
		self.printDiagnostics = printDiagnostics
		self.maxSubCycles = maxSubCycles
		self.maxMacroCycles = maxMacroCycles
		self.printPurity = printPurity
		self.adjustLoadVol = adjustLoadVol
		self.adjustBindingSites = adjustBindingSites
		
		# Dissociation constants for each of the binding sites for each of the metals
		# Site 1
		self.kd1_1 = kd1_1
		self.kd1_2 = kd1_2
		self.kd1_3 = kd1_3
		
		# Site 2
		self.kd2_1 = kd2_1
		self.kd2_2 = kd2_2
		self.kd2_3 = kd2_3
		
		self.initial_nM1T = initial_nM1T
		self.initial_nM2T = initial_nM2T
		self.initial_nM3T = initial_nM3T
		
		self.initial_nBT = initial_nBT
		self.fB1 = fB1
		self.fB2 = fB2
		
		self.initial_nB1T = self.initial_nBT * self.fB1
		self.initial_nB2T = self.initial_nBT * self.fB2
		
		self.initial_loadVol = initial_loadVol
		self.initial_eluteVol = initial_eluteVol
		
		self.lowest_allowed_vol = lowest_allowed_vol
		self.targetPurity = targetPurity
		self.preferred_operating_concentration = preferred_operating_concentration
		
		
	
	
	def update_with_elute_or_wash_solution(self, eluteOrWashSolution):
		
		import pdb
		
		self.initial_nM1T = eluteOrWashSolution.nM1
		self.initial_nM2T = eluteOrWashSolution.nM2
		self.initial_nM3T = eluteOrWashSolution.nM3
		
		self.initial_nBT, self.initial_nB1T, self.initial_nB2T, self.initial_loadVol = \
		Calculate_New_2_Binding_Site_Numbers_and_Load_Volume(eluteOrWashSolution, self)
		
		
		
# ------------------------------------------------------------------------------------------------ #







# ------------------------------------------------------------------------------------------------ #
class InputData_3Metals_1Site:
	def __init__(self, kd1, kd2, kd3, initial_nM1T, initial_nM2T, initial_nM3T, initial_nBT, \
	initial_loadVol, 
	initial_eluteVol = 1.0, \
	targetPurity=0.9, preferred_operating_concentration=0.1, \
	lowest_allowed_vol=0.001, \
	printFigures=False, printDiagnostics=True, maxSubCycles=20, maxMacroCycles=100, \
	printPurity=False, adjustLoadVol=True, adjustBindingSites=True):
		
		import pdb
		
		self.printFigures = printFigures
		self.printDiagnostics = printDiagnostics
		self.maxSubCycles = maxSubCycles
		self.maxMacroCycles = maxMacroCycles
		self.printPurity = printPurity
		self.adjustLoadVol = adjustLoadVol
		self.adjustBindingSites = adjustBindingSites
		
		self.kd1 = kd1
		self.kd2 = kd2
		self.kd3 = kd3
		
		self.initial_nM1T = initial_nM1T
		self.initial_nM2T = initial_nM2T
		self.initial_nM3T = initial_nM3T
		
		self.initial_nBT = initial_nBT
		self.initial_loadVol = initial_loadVol
		self.initial_eluteVol = initial_eluteVol
		
		self.lowest_allowed_vol = lowest_allowed_vol
		
		self.targetPurity = targetPurity
		
		self.preferred_operating_concentration = preferred_operating_concentration
		
# 		pdb.set_trace()
		
		
	
	
	def update_with_elute_or_wash_solution(self, eluteOrWashSolution):
		
		import pdb
		
		self.initial_nM1T = eluteOrWashSolution.nM1
		self.initial_nM2T = eluteOrWashSolution.nM2
		self.initial_nM3T = eluteOrWashSolution.nM3
		
		self.initial_nBT, self.initial_loadVol = \
		Calculate_New_Binding_Site_Number_and_Load_Volume(eluteOrWashSolution, self)
		
		#pdb.set_trace()





# ------------------------------------------------------------------------------------------------ #
class SubCycleDiagnostic:
	def __init__(self, loadState, bindState, washSolution, eluteSolution):
		
		self.loadState = loadState
		self.bindState = bindState
		self.washSolution = washSolution
		self.eluteSolution = eluteSolution
		
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
class SubCycleDiagnostic_3M_3S_3Microbes:
	def __init__(self, loadState_microbe1, bindState_microbe1, washSolution_microbe1, \
	eluteSolution_microbe1, \
	loadState_microbe2, bindState_microbe2, washSolution_microbe2, eluteSolution_microbe2, \
	loadState_microbe3, bindState_microbe3, washSolution_microbe3, eluteSolution_microbe3 ):
		
		self.loadState_microbe1 = loadState_microbe1
		self.bindState_microbe1 = bindState_microbe1
		self.washSolution_microbe1 = washSolution_microbe1
		self.eluteSolution_microbe1 = eluteSolution_microbe1
		
		self.loadState_microbe2 = loadState_microbe2
		self.bindState_microbe2 = bindState_microbe2
		self.washSolution_microbe2 = washSolution_microbe2
		self.eluteSolution_microbe2 = eluteSolution_microbe2
		
		self.loadState_microbe3 = loadState_microbe3
		self.bindState_microbe3 = bindState_microbe3
		self.washSolution_microbe3 = washSolution_microbe3
		self.eluteSolution_microbe3 = eluteSolution_microbe3
		
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
class MacroCycleDiagnostic_3Metals_1Site:
	def __init__(self, inputData):
		from copy import deepcopy
		self.copiedInputData = deepcopy(inputData)
		self.inputData = inputData
		self.subCycleDiagnosticArray = []
		
	def summarize_macro_cycle(self):
		fM1EluantArray = []
		fM2EluantArray = []
		fM3EluantArray = []
		
		nM1EluantArray = []
		nM2EluantArray = []
		nM3EluantArray = []
		
		DM1BindArray = []
		DM2BindArray = []
		DM3BindArray = []
		
		alpha1_2_BindArray = []
		alpha1_3_BindArray = []
		alpha2_3_BindArray = []
		
		nBTArray = []

		volArray = []
		
		cycleArray = []
		
		fM1EluantArray.append(self.subCycleDiagnosticArray[0].loadState.fM1)
		fM2EluantArray.append(self.subCycleDiagnosticArray[0].loadState.fM2)
		fM3EluantArray.append(self.subCycleDiagnosticArray[0].loadState.fM3)
		
		nM1EluantArray.append(self.subCycleDiagnosticArray[0].loadState.nM1)
		nM2EluantArray.append(self.subCycleDiagnosticArray[0].loadState.nM2)
		nM3EluantArray.append(self.subCycleDiagnosticArray[0].loadState.nM3)
		
		nBTArray.append(self.subCycleDiagnosticArray[0].loadState.nBT)
		
		volArray.append(self.subCycleDiagnosticArray[0].loadState.vol)
		
		# How do we define the DM and alpha arrays for the zeroth step? These are a bit of a kludge
		
		DM1BindArray.append(1)
		DM2BindArray.append(1)
		DM3BindArray.append(1)
		
		nM1_load = self.subCycleDiagnosticArray[0].loadState.nM1
		nM2_load = self.subCycleDiagnosticArray[0].loadState.nM2
		nM3_load = self.subCycleDiagnosticArray[0].loadState.nM3
		
		alpha1_2_load = nM1_load / nM2_load
		alpha1_3_load = nM1_load / nM3_load
		alpha2_3_load = nM2_load / nM3_load

		alpha1_2_BindArray.append(nM1_load / nM2_load)
		alpha1_3_BindArray.append(nM1_load / nM3_load)
		alpha2_3_BindArray.append(nM2_load / nM3_load)
		


		cycleArray.append(0)

		i = 0
		while i < len(self.subCycleDiagnosticArray):
			cycleArray.append(i + 1)
			
			fM1EluantArray.append(self.subCycleDiagnosticArray[i].eluteSolution.fM1)
			fM2EluantArray.append(self.subCycleDiagnosticArray[i].eluteSolution.fM2)
			fM3EluantArray.append(self.subCycleDiagnosticArray[i].eluteSolution.fM3)
			
			nM1EluantArray.append(self.subCycleDiagnosticArray[i].eluteSolution.nM1)
			nM2EluantArray.append(self.subCycleDiagnosticArray[i].eluteSolution.nM2)
			nM3EluantArray.append(self.subCycleDiagnosticArray[i].eluteSolution.nM3)
			
			DM1BindArray.append(self.subCycleDiagnosticArray[i].bindState.DM1)
			DM2BindArray.append(self.subCycleDiagnosticArray[i].bindState.DM2)
			DM3BindArray.append(self.subCycleDiagnosticArray[i].bindState.DM3)

			alpha1_2_BindArray.append(self.subCycleDiagnosticArray[i].bindState.alpha1_2)
			alpha1_3_BindArray.append(self.subCycleDiagnosticArray[i].bindState.alpha1_3)
			alpha2_3_BindArray.append(self.subCycleDiagnosticArray[i].bindState.alpha2_3)
			
			nBTArray.append(self.subCycleDiagnosticArray[i].loadState.nBT)
			
			volArray.append(self.subCycleDiagnosticArray[i].loadState.vol)

			i += 1
		
		returnDict = {}
		returnDict['cycleArray'] = cycleArray
		returnDict['fM1EluantArray'] = fM1EluantArray
		returnDict['fM2EluantArray'] = fM2EluantArray
		returnDict['fM3EluantArray'] = fM3EluantArray
		returnDict['nM1EluantArray'] = nM1EluantArray
		returnDict['nM2EluantArray'] = nM2EluantArray
		returnDict['nM3EluantArray'] = nM3EluantArray
		returnDict['DM1BindArray'] = DM1BindArray
		returnDict['DM2BindArray'] = DM2BindArray
		returnDict['DM3BindArray'] = DM3BindArray
		returnDict['alpha1_2_BindArray'] = alpha1_2_BindArray
		returnDict['alpha1_3_BindArray'] = alpha1_3_BindArray
		returnDict['alpha2_3_BindArray'] = alpha2_3_BindArray
		returnDict['volArray'] = volArray
		returnDict['nBTArray'] = nBTArray
		
		
		return returnDict
		
# ------------------------------------------------------------------------------------------------ #






# ------------------------------------------------------------------------------------------------ #
class MacroCycleDiagnostic_3Metals_2Sites(MacroCycleDiagnostic_3Metals_1Site):
	def __init__(self, inputData):
		super().__init__(inputData)

		
	def summarize_macro_cycle(self):
		returnDict = super().summarize_macro_cycle()
		
		nB1Array = []
		nB2Array = []

		
		nB1Array.append(self.subCycleDiagnosticArray[0].loadState.nB1T)
		nB2Array.append(self.subCycleDiagnosticArray[0].loadState.nB2T)
		

		i = 0
		while i < len(self.subCycleDiagnosticArray):
			
			nB1Array.append(self.subCycleDiagnosticArray[i].loadState.nB1T)
			nB2Array.append(self.subCycleDiagnosticArray[i].loadState.nB2T)
			
			i += 1
		
		returnDict['nB1Array'] = nB1Array
		returnDict['nB2Array'] = nB2Array
		
		
		return returnDict
		
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
class MacroCycleDiagnostic_3Metals_3Sites(MacroCycleDiagnostic_3Metals_1Site):
	def __init__(self, inputData):
		super().__init__(inputData)

		
	def summarize_macro_cycle(self):
		returnDict = super().summarize_macro_cycle()
		
		nB1Array = []
		nB2Array = []
		nB3Array = []
		
		nB1Array.append(self.subCycleDiagnosticArray[0].loadState.nB1T)
		nB2Array.append(self.subCycleDiagnosticArray[0].loadState.nB2T)
		nB3Array.append(self.subCycleDiagnosticArray[0].loadState.nB3T)
		

		i = 0
		while i < len(self.subCycleDiagnosticArray):
			
			nB1Array.append(self.subCycleDiagnosticArray[i].loadState.nB1T)
			nB2Array.append(self.subCycleDiagnosticArray[i].loadState.nB2T)
			nB3Array.append(self.subCycleDiagnosticArray[i].loadState.nB3T)
			
			i += 1
		
		returnDict['nB1Array'] = nB1Array
		returnDict['nB2Array'] = nB2Array
		returnDict['nB3Array'] = nB3Array
		
		
		return returnDict
		
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
class MacroCycleDiagnostic_3Metals_3Sites_3Microbes(MacroCycleDiagnostic_3Metals_1Site):
	def __init__(self, inputData):
		super().__init__(inputData)

		
	def summarize_macro_cycle(self):
		
		fM1_Microbe3_Wash_Array = []
		fM2_Microbe3_Wash_Array = []
		fM3_Microbe3_Wash_Array = []
		
		cycleArray = []
		
		fM1_Microbe3_Wash_Array.append(self.subCycleDiagnosticArray[0].loadState_microbe1.fM1)
		fM2_Microbe3_Wash_Array.append(self.subCycleDiagnosticArray[0].loadState_microbe1.fM2)
		fM3_Microbe3_Wash_Array.append(self.subCycleDiagnosticArray[0].loadState_microbe1.fM3)
		
		
		cycleArray.append(0)

		i = 0
		while i < len(self.subCycleDiagnosticArray):
			cycleArray.append(i + 1)
			
			fM1 = self.subCycleDiagnosticArray[i].washSolution_microbe3.fM1
			fM2 = self.subCycleDiagnosticArray[i].washSolution_microbe3.fM2
			fM3 = self.subCycleDiagnosticArray[i].washSolution_microbe3.fM3
			
			fM1_Microbe3_Wash_Array.append(fM1)
			fM2_Microbe3_Wash_Array.append(fM2)
			fM3_Microbe3_Wash_Array.append(fM3)
			
			i += 1
		
		returnDict = {}
		returnDict['cycleArray'] = cycleArray
		returnDict['fM1_Microbe3_Wash_Array'] = fM1_Microbe3_Wash_Array
		returnDict['fM2_Microbe3_Wash_Array'] = fM2_Microbe3_Wash_Array
		returnDict['fM3_Microbe3_Wash_Array'] = fM3_Microbe3_Wash_Array
		
		
		return returnDict	
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Calculate_New_Binding_Site_Numbers_and_Load_Volume(eluteOrWashSolution, inputData):
	
	import pdb
	
	lowest_allowed_vol = inputData.lowest_allowed_vol
	preferred_operating_concentration = inputData.preferred_operating_concentration
	

	if inputData.adjustBindingSites == True:
		# Reduce the number of binding sites to the number of M1 atoms
		nM1 = eluteOrWashSolution.nM1
		newNBT = nM1 * 1
		newNB1 = newNBT * inputData.fB1
		newNB2 = newNBT * inputData.fB2
		newNB3 = newNBT * inputData.fB3
		
	else:
		newNBT = inputData.initial_nBT 
		newNB1 = newNBT * inputData.fB1
		newNB2 = newNBT * inputData.fB2
		newNB3 = newNBT * inputData.fB3
		

	
	if inputData.adjustLoadVol == True:
		# Adjust the loading volume to keep the concentration high enough
		newLoadVol = eluteOrWashSolution.nM1 / preferred_operating_concentration
		
	else:
		newLoadVol = inputData.initial_loadVol

	#pdb.set_trace()
	
	return newNBT, newNB1, newNB2, newNB3, newLoadVol
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Calculate_New_Binding_Site_Numbers_and_Load_Volume_2(eluteOrWashSolution, \
inputData, target='nM1'):
	
	import pdb
	
	lowest_allowed_vol = inputData.lowest_allowed_vol
	preferred_operating_concentration = inputData.preferred_operating_concentration
	

	if inputData.adjustBindingSites == True:
		# Reduce the number of binding sites to the number of M1 atoms
		
		if target == 'nM1':
			nM1 = eluteOrWashSolution.nM1
			newNBT = nM1 * 1
		elif target == 'nM2':
			nM2 = eluteOrWashSolution.nM2
			newNBT = nM2 * 1
		elif target == 'nM3':
			nM3 = eluteOrWashSolution.nM3
			newNBT = nM3 * 1
		
		newNB1 = newNBT * inputData.fB1
		newNB2 = newNBT * inputData.fB2
		newNB3 = newNBT * inputData.fB3
		
	else:
		newNBT = inputData.initial_nBT 
		newNB1 = newNBT * inputData.fB1
		newNB2 = newNBT * inputData.fB2
		newNB3 = newNBT * inputData.fB3
		

	
	if inputData.adjustLoadVol == True:
		# Adjust the loading volume to keep the concentration high enough
		newLoadVol = eluteOrWashSolution.nM1 / preferred_operating_concentration
		
	else:
		newLoadVol = inputData.initial_loadVol

	#pdb.set_trace()
	
	return newNBT, newNB1, newNB2, newNB3, newLoadVol
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def Calculate_New_2_Binding_Site_Numbers_and_Load_Volume(eluteOrWashSolution, inputData):
	
	import pdb
	
	lowest_allowed_vol = inputData.lowest_allowed_vol
	preferred_operating_concentration = inputData.preferred_operating_concentration
	

	if inputData.adjustBindingSites == True:
		# Reduce the number of binding sites to the number of M1 atoms
		nM1 = eluteOrWashSolution.nM1
		newNBT = nM1 * 1
		newNB1 = newNBT * inputData.fB1
		newNB2 = newNBT * inputData.fB2
		
	else:
		newNBT = inputData.initial_nBT 
		newNB1 = newNBT * inputData.fB1
		newNB2 = newNBT * inputData.fB2
		

	
	if inputData.adjustLoadVol == True:
		# Adjust the loading volume to keep the concentration high enough
		newLoadVol = eluteOrWashSolution.nM1 / preferred_operating_concentration
		
	else:
		newLoadVol = inputData.initial_loadVol


	return newNBT, newNB1, newNB2, newLoadVol
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Calculate_New_Binding_Site_Number_and_Load_Volume(eluteOrWashSolution, inputData):
	
	import pdb
	
	lowest_allowed_vol = inputData.lowest_allowed_vol
	preferred_operating_concentration = inputData.preferred_operating_concentration
	

	if inputData.adjustBindingSites == True:
		# Reduce the number of binding sites to the number of M1 atoms
		nM1 = eluteOrWashSolution.nM1
		newNBT = nM1 * 1
	else:
		newNBT = inputData.initial_nBT 

	
	if inputData.adjustLoadVol == True:
		# Adjust the loading volume to keep the concentration high enough
		newLoadVol = eluteOrWashSolution.nM1 / preferred_operating_concentration
		
	else:
		newLoadVol = inputData.initial_loadVol


	return newNBT, newLoadVol
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
# This function performs sub-cycles until a target purity of metal 1 is reached

def Run_Macro_Cycle_3Metals_3Sites(inputData):
	
	import pdb

	initialLoadState = LoadState_3Metals_3Sites()
	initialLoadState.init_with_input_data(inputData)
		
	i = 1
	loadState = initialLoadState
	currentPurity = loadState.fM1
	
	macroCycleDiagnostic = MacroCycleDiagnostic_3Metals_3Sites(inputData)

	# Iterative sub-cycles for purification starts here. 
	# The idea is to repeatedly pass the eluant from the previous cycle back through the column
	# until the target purity is achieved
	while i <= inputData.maxSubCycles and currentPurity < inputData.targetPurity:
		
		print(i)
		
		# This is the most important calculation. Calculate the concentrations of metals that 
		# are bound and free in the binding column
		bindState = BindState_3Metals_3Sites(loadState)
		
		
		# Calculate the contents of the wash and eluant
		washSolution = WashSolution_3Metals()
		washSolution.init_with_bind_state(bindState)
		eluteSolution = EluteSolution_3Metals()
		eluteSolution.init_with_bind_state(bindState, inputData.initial_eluteVol)

		# Prepare report out on cycle		
		subCycleDiagnostic = SubCycleDiagnostic(loadState, bindState, washSolution, eluteSolution)
		macroCycleDiagnostic.subCycleDiagnosticArray.append(subCycleDiagnostic)
		
		currentPurity = eluteSolution.fM1

		# Get ready for next cycle
		# Reduce the loading volume to keep the concentrations of metals in the column reasonable
		# In this case where we have three binding sites, we are going to maintain equal numbers
		# of each type. 
		newNBT, newNB1, newNB2, newNB3, newLoadVol = \
		Calculate_New_Binding_Site_Numbers_and_Load_Volume(eluteSolution, inputData)
		
		#pdb.set_trace()


		# Create a description of the solution to be loaded into the column on the next cycle
		loadState = LoadState_3Metals_3Sites()
		loadState.init_with_elute_or_wash_solution(eluteSolution, \
		inputData.kd1_1, inputData.kd1_2, inputData.kd1_3, \
		inputData.kd2_1, inputData.kd2_2, inputData.kd2_3, \
		inputData.kd3_1, inputData.kd3_2, inputData.kd3_3, \
		newNB1, newNB2, newNB3, newLoadVol)
		
		
		
		i += 1

			
	return macroCycleDiagnostic
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
# This function performs sub-cycles until a target purity of metal 1 is reached

def Run_Macro_Cycle_3Metals_2Sites(inputData):
	
	import pdb

	initialLoadState = LoadState_3Metals_2Sites()
	initialLoadState.init_with_input_data(inputData)
		
	i = 1
	loadState = initialLoadState
	currentPurity = loadState.fM1
	
	macroCycleDiagnostic = MacroCycleDiagnostic_3Metals_2Sites(inputData)

	# Iterative sub-cycles for purification starts here. 
	# The idea is to repeatedly pass the eluant from the previous cycle back through the column
	# until the target purity is achieved
	while i <= inputData.maxSubCycles and currentPurity < inputData.targetPurity:
		
		print(i)
		
		# This is the most important calculation. Calculate the concentrations of metals that 
		# are bound and free in the binding column
		bindState = BindState_3Metals_2Sites(loadState)
		
		# Calculate the contents of the wash and eluant
		washSolution = WashSolution_3Metals()
		washSolution.init_with_bind_state(bindState)
		eluteSolution = EluteSolution_3Metals()
		eluteSolution.init_with_bind_state(bindState, inputData.initial_eluteVol)

		# Prepare report out on cycle		
		subCycleDiagnostic = SubCycleDiagnostic(loadState, bindState, washSolution, eluteSolution)
		macroCycleDiagnostic.subCycleDiagnosticArray.append(subCycleDiagnostic)
		
		currentPurity = eluteSolution.fM1

		# Get ready for next cycle
		# Reduce the loading volume to keep the concentrations of metals in the column reasonable
		# In this case where we have three binding sites, we are going to maintain equal numbers
		# of each type. 
		newNBT, newNB1, newNB2, newLoadVol = \
		Calculate_New_2_Binding_Site_Numbers_and_Load_Volume(eluteSolution, inputData)


		# Create a description of the solution to be loaded into the column on the next cycle
		loadState = LoadState_3Metals_2Sites()
		loadState.init_with_elute_or_wash_solution(eluteSolution, \
		inputData.kd1_1, inputData.kd1_2, inputData.kd1_3, \
		inputData.kd2_1, inputData.kd2_2, inputData.kd2_3, \
		newNB1, newNB2, newLoadVol)
		
		i += 1

			
	return macroCycleDiagnostic
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
# This function performs sub-cycles until a target purity of metal 1 is reached

def Run_Macro_Cycle_3Metals_1Site(inputData):
	
	import pdb

	initialLoadState = LoadState_3Metals_1Site()
	initialLoadState.init_with_input_data(inputData)
		
	i = 1
	loadState = initialLoadState
	currentPurity = loadState.fM1
	
	macroCycleDiagnostic = MacroCycleDiagnostic_3Metals_1Site(inputData)

	# Iterative sub-cycles for purification starts here. 
	# The idea is to repeatedly pass the eluant from the previous cycle back through the column
	# until the target purity is achieved
	while i <= inputData.maxSubCycles and currentPurity < inputData.targetPurity:
		
		print(i)
		
		# This is the most important calculation. Calculate the concentrations of metals that 
		# are bound and free in the binding column
		bindState = BindState_3Metals_1Site(loadState)
		
		
		# Calculate the contents of the wash and eluant
		washSolution = WashSolution_3Metals()
		washSolution.init_with_bind_state(bindState)
		eluteSolution = EluteSolution_3Metals()
		eluteSolution.init_with_bind_state(bindState, inputData.initial_eluteVol)

		# Prepare report out on cycle		
		subCycleDiagnostic = SubCycleDiagnostic(loadState, bindState, washSolution, eluteSolution)
		macroCycleDiagnostic.subCycleDiagnosticArray.append(subCycleDiagnostic)
		
		currentPurity = eluteSolution.fM1

		# Get ready for next cycle
		# Reduce the loading volume to keep the concentrations of metals in the column reasonable
		# In this case where we have three binding sites, we are going to maintain equal numbers
		# of each type. 
		newNBT, newLoadVol = \
		Calculate_New_Binding_Site_Number_and_Load_Volume(eluteSolution, inputData)
		
		#pdb.set_trace()


		# Create a description of the solution to be loaded into the column on the next cycle
		loadState = LoadState_3Metals_1Site()
		loadState.init_with_elute_or_wash_solution(eluteSolution, \
		inputData.kd1, inputData.kd2, inputData.kd3, \
		newNBT, newLoadVol)
		
		#pdb.set_trace()
		
		
		
		i += 1

			
	return macroCycleDiagnostic
# ------------------------------------------------------------------------------------------------ #





# ------------------------------------------------------------------------------------------------ #
def nb_optimize_3sites_3kd_func(x, kdBase, fMxB_target, init_input_data_3M_3S):

# x will contain the fraction of metal binding sites

	from numpy import zeros, dot
	
	tempLoadState = LoadState_3Metals_3Sites()
	tempLoadState.init_with_input_data(init_input_data_3M_3S)
	
	nBT = tempLoadState.nBT
	fB1 = x[0]
	fB2 = x[1]
	fB3 = x[2]
	
	kd1_1 = (1/x[3]) * kdBase
	kd2_2 = (1/x[4]) * kdBase
	kd3_3 = (1/x[5]) * kdBase
	
	tempLoadState.nB1T = nBT*fB1
	tempLoadState.nB2T = nBT*fB2
	tempLoadState.nB3T = nBT*fB3
	
	tempLoadState.kd1_1 = kd1_1
	tempLoadState.kd2_2 = kd2_2
	tempLoadState.kd3_3 = kd3_3
	
	print(str(fB1) + ', ' + str(fB2) + ', ' + str(fB3) + ', ' \
	+ str(x[3]) + ', ' + str(x[4]) + ', ' + str(x[5]))

	bindState_3M_3S = BindState_3Metals_3Sites(tempLoadState)
	
	fM1b = bindState_3M_3S.fM1b
	fM2b = bindState_3M_3S.fM2b
	fM3b = bindState_3M_3S.fM3b

	y = zeros(4)
	
	y[0] = fM1b - fMxB_target[0]
	y[1] = fM2b - fMxB_target[1]
	y[2] = fM3b - fMxB_target[2]
	y[3] = 1 - (fB1 + fB2 + fB3)
	
	return dot(y,y)
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def kd_optimize_1site_func(x, kdBase, fMxB_target, init_input_data_3M_1S):

# x will contain the kd

	from numpy import zeros, dot
	
	kd1 = kdBase 
	kd2 = x[0]*kdBase
	kd3 = x[1]*kdBase
	
	print(str(kd1) + ', ' + str(kd2) + ', ' + str(kd3))

	tempLoadState = LoadState_3Metals_1Site()
	tempLoadState.init_with_input_data(init_input_data_3M_1S)
	
	tempLoadState.kd1 = kd1
	tempLoadState.kd2 = kd2
	tempLoadState.kd3 = kd3
	
	bindState_3M_1S = BindState_3Metals_1Site(tempLoadState)
	
	fM1b = bindState_3M_1S.fM1b
	fM2b = bindState_3M_1S.fM2b
	fM3b = bindState_3M_1S.fM3b

	y = zeros(3)
	
	y[0] = fM1b - fMxB_target[0]
	y[1] = fM2b - fMxB_target[1]
	y[2] = fM3b - fMxB_target[2]
	
	return dot(y,y)
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Calculate_Fraction_of_Each_Metal_Bound_to_Single_Site(kd1, kd2, kd3, \
initial_nM1T, initial_nM2T, initial_nM3T, initial_nBT, initial_loadVol):
	
	import pdb
	
	inputData = InputData_3Metals_1Site(kd1, kd2, kd3, initial_nM1T, initial_nM2T, initial_nM3T, \
	initial_nBT, initial_loadVol)

	initialLoadState = LoadState_3Metals_1Site()
	
	initialLoadState.init_with_input_data(inputData)

	bindState = BindState_3Metals_1Site(initialLoadState)
	
	fM1b = bindState.fM1b
	fM2b = bindState.fM2b
	fM3b = bindState.fM3b
	
	alpha1_2 = bindState.alpha1_2
	alpha1_3 = bindState.alpha1_3
	alpha2_3 = bindState.alpha2_3

	return fM1b, fM2b, fM3b, alpha1_2, alpha1_3, alpha2_3
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Calculate_Fraction_of_Each_Metal_Bound_to_Three_Sites(\
kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
initial_nM1T, initial_nM2T, initial_nM3T, initial_nBT, fB1, fB2, fB3, initial_loadVol):
	
	import pdb
	
	inputData = InputData_3Metals_3Sites(\
	kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
	initial_nM1T, initial_nM2T, initial_nM3T, initial_nBT, fB1, fB2, fB3, initial_loadVol)

	initialLoadState = LoadState_3Metals_3Sites()	
	initialLoadState.init_with_input_data(inputData)

	bindState = BindState_3Metals_3Sites(initialLoadState)
	
	fM1b = bindState.fM1b
	fM2b = bindState.fM2b
	fM3b = bindState.fM3b
	
	alpha1_2 = bindState.alpha1_2
	alpha1_3 = bindState.alpha1_3
	alpha2_3 = bindState.alpha2_3

	return fM1b, fM2b, fM3b, alpha1_2, alpha1_3, alpha2_3
# ------------------------------------------------------------------------------------------------ #




# ------------------------------------------------------------------------------------------------ #
def Calculate_Sub_Cycles_to_Reach_Target_Purity_and_M1_Remaining_3M_1S(\
kdArray, initial_nBT, initial_nM1T,initial_nM2T, initial_nM3T, initial_loadVol, targetPurity):

	import pdb

	i = 0
	subCycleNumberArray = []
	M1RemainingArray = []

	while i < len(kdArray):

		kdSet = kdArray[i]
		
		kd1 = kdArray[i][0]
		kd2 = kdArray[i][1]
		kd3 = kdArray[i][2]
			
		inputData = \
		InputData_3Metals_1Site(\
		kd1, kd2, kd3, initial_nM1T, initial_nM2T, initial_nM3T, initial_nBT, \
		initial_loadVol, 
		initial_eluteVol=1.0, \
		targetPurity=targetPurity, preferred_operating_concentration=0.1, \
		lowest_allowed_vol=0.001, \
		printFigures=False, printDiagnostics=False, maxSubCycles=1000, maxMacroCycles=100, \
		printPurity=False, adjustLoadVol=True, adjustBindingSites=True)
	
		macroCycleDiagnostic = Run_Macro_Cycle_3Metals_1Site(inputData)
	
		macroCycleDiagnosticReturnDict = macroCycleDiagnostic.summarize_macro_cycle()
		
# 		pdb.set_trace()
		
		nM1EluantArray = macroCycleDiagnosticReturnDict['nM1EluantArray']
		cycleArray = macroCycleDiagnosticReturnDict['cycleArray']
	
		m1Remaining = nM1EluantArray[-1]
		M1RemainingArray.append(m1Remaining)
	
		numberSubCycles = len(cycleArray) - 1
		subCycleNumberArray.append(numberSubCycles)
	
		i += 1
	
	returnDict = {}	
	returnDict['kdArray'] = kdArray
	returnDict['subCycleNumber'] = subCycleNumberArray
	returnDict['M1RemainingArray'] = M1RemainingArray
		
	return returnDict
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Calculate_Binding_Site_Fraction_Arrays(fB1_Array, fB2_0, fB3_0):
# Calculate the fractions of each binding site. This is a slightly tricky calculation. 
# We start out with the original fitted fractions of each site, and then figure out how to reduce
# the fractions of binding sites 2 and 3 as the fraction of binding site 1 increases. 

	# Calculate the sum fraction of sites 2 and 3
	fB_2_and_3_total = 1 - fB1_Array[0]

	# Calculate the fraction of site 2 out of sites 2 and 3
	fB2_other = fB2_0 / fB_2_and_3_total
	fB3_other = fB3_0 / fB_2_and_3_total

	# Now, calculate the fractions of sites 2 and 3 as the fraction of site 1 increases
	fB2_Array = [fB2_0]
	fB3_Array = [fB3_0]
	i = 1
	while i < len(fB1_Array):
		fB_2_and_3_total_new = 1 - fB1_Array[i]
		fB2_new = fB2_other * fB_2_and_3_total_new
		fB3_new = fB3_other * fB_2_and_3_total_new
	
		fB2_Array.append(fB2_new)
		fB3_Array.append(fB3_new)
	
		fB_total = fB1_Array[i] + fB2_new + fB3_new
		
		print(str(fB1_Array[i]) + ', ' + str(fB2_new) + ', ' + str(fB3_new) + ', ' + str(fB_total))
	
		i += 1
	
	return fB2_Array, fB3_Array
# ------------------------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------------------------ #
def Increment_Binding_Site_Fraction(microbe_fB_0, increment_array):
# Produce an array of incremented binding site fractions

	microbe_fB_array = [microbe_fB_0]
	microbe_fB = microbe_fB_0

	fB_increase_array = [0]
	fB_increase = 0

	i = 0
	while i < len(increment_array):
	
		increment = increment_array[i]
	
		microbe_fB += increment 
		fB_increase += increment
	
		microbe_fB_array.append(microbe_fB)
		fB_increase_array.append(fB_increase)
	
		i += 1
	
	return microbe_fB_array, fB_increase_array
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Calculate_Sub_Cycles_to_Reach_Target_Purity_and_M1_Remaining_3M_3S(\
kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
initial_nBT, fB_Array, initial_nM1T,initial_nM2T, initial_nM3T, initial_loadVol, targetPurity):

	import pdb

	i = 0
	subCycleNumberArray = []
	M1RemainingArray = []
	fB1_Array = []

	while i < len(fB_Array):
	
		fB1 = fB_Array[i][0]
		fB2 = fB_Array[i][1]
		fB3 = fB_Array[i][2]
		
		fB1_Array.append(fB1)
		
		inputData = \
		InputData_3Metals_3Sites(\
		kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
		initial_nM1T, initial_nM2T, initial_nM3T, \
		initial_nBT, fB1, fB2, fB3, \
		initial_loadVol, \
		initial_eluteVol=1.0, \
		targetPurity=targetPurity, preferred_operating_concentration=0.1, \
		lowest_allowed_vol=0.001, \
		printFigures=False, printDiagnostics=False, maxSubCycles=1000, maxMacroCycles=100, \
		printPurity=False, adjustLoadVol=True, adjustBindingSites=True)
	
		macroCycleDiagnostic = Run_Macro_Cycle_3Metals_3Sites(inputData)
	
		macroCycleDiagnosticReturnDict = macroCycleDiagnostic.summarize_macro_cycle()
		
# 		pdb.set_trace()
		
		nM1EluantArray = macroCycleDiagnosticReturnDict['nM1EluantArray']
		cycleArray = macroCycleDiagnosticReturnDict['cycleArray']
	
		m1Remaining = nM1EluantArray[-1]
		M1RemainingArray.append(m1Remaining)
	
		numberSubCycles = len(cycleArray) - 1
		subCycleNumberArray.append(numberSubCycles)
	
		i += 1
	
	returnDict = {}	
	returnDict['fB1_Array'] = fB1_Array
	returnDict['subCycleNumber'] = subCycleNumberArray
	returnDict['M1RemainingArray'] = M1RemainingArray
		
	return returnDict
# ------------------------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------------------------ #
def Calculate_Sub_Cycles_to_Reach_Target_Purity_and_M1_Remaining_3M_3S_3Microbes(\
kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, nBT, \
microbe1_fB1_Array, microbe1_fB2_Array, microbe1_fB3_Array, \
microbe2_fB1_Array, microbe2_fB2_Array, microbe2_fB3_Array, \
microbe3_fB1_Array, microbe3_fB2_Array, microbe3_fB3_Array, \
nM1T, nM2T, nM3T, loadVol, targetPurity):

	import pdb

	i = 0
	subCycleNumberArray = []
		
	maxSubCycles = 1000
	maxMacroCycles = 100

	while i < len(microbe1_fB1_Array):
			
		microbe1_fB1 = microbe1_fB1_Array[i]
		microbe1_fB2 = microbe1_fB2_Array[i]
		microbe1_fB3 = microbe1_fB3_Array[i]
	
		inputData_microbe1 = \
		InputData_3Metals_3Sites(\
		kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
		nM1T, nM2T, nM3T, nBT, microbe1_fB1, microbe1_fB2, microbe1_fB3, \
		loadVol, adjustLoadVol=True, adjustBindingSites=True, maxSubCycles=maxSubCycles, \
		targetPurity=targetPurity) 
		# ---------------------------------------------------------------------------------------- #
	
		# ---------------------------------------------------------------------------------------- #
		# Prepare Microbe 2
		microbe2_fB1 = microbe2_fB1_Array[i]
		microbe2_fB2 = microbe2_fB2_Array[i]
		microbe2_fB3 = microbe2_fB3_Array[i]
	
		inputData_microbe2 = \
		InputData_3Metals_3Sites(\
		kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
		nM1T, nM2T, nM3T, nBT, microbe2_fB1, microbe2_fB2, microbe2_fB3, \
		loadVol, adjustLoadVol=True, adjustBindingSites=True, maxSubCycles=maxSubCycles, \
		targetPurity=targetPurity) 
		# ---------------------------------------------------------------------------------------- #
	
		# ---------------------------------------------------------------------------------------- #
		# Prepare Microbe 3
		microbe3_fB1 = microbe3_fB1_Array[i]
		microbe3_fB2 = microbe3_fB2_Array[i]
		microbe3_fB3 = microbe3_fB3_Array[i]
	
		inputData_microbe3 = \
		InputData_3Metals_3Sites(\
		kd1_1, kd1_2, kd1_3, kd2_1, kd2_2, kd2_3, kd3_1, kd3_2, kd3_3, \
		nM1T, nM2T, nM3T, nBT, microbe3_fB1, microbe3_fB2, microbe3_fB3, \
		loadVol, adjustLoadVol=True, adjustBindingSites=True, maxSubCycles=maxSubCycles, \
		targetPurity=targetPurity) 
		# ---------------------------------------------------------------------------------------- #
	
		# ---------------------------------------------------------------------------------------- #
		# Run the macro cycle	
		macroCycleDiagnostic = Run_Macro_Cycle_3Metals_3Sites_3Microbes(\
		inputData_microbe1, inputData_microbe2, inputData_microbe3)
	
		macroCycleDiagnosticReturnDict = macroCycleDiagnostic.summarize_macro_cycle()
		
		cycleArray = macroCycleDiagnosticReturnDict['cycleArray']
	
		numberSubCycles = len(cycleArray) - 1
		subCycleNumberArray.append(numberSubCycles)
	
		i += 1
	
	returnDict = {}	
	returnDict['microbe1_fB1_Array'] = microbe1_fB1_Array
	returnDict['subCycleNumber'] = subCycleNumberArray
		
	return returnDict
# ------------------------------------------------------------------------------------------------ #







# ------------------------------------------------------------------------------------------------ #
# This function performs sub-cycles until a target purity of metal 1 is reached

def Run_Macro_Cycle_3Metals_3Sites_3Microbes(inputData_microbe1, inputData_microbe2, \
inputData_microbe3):
	
	import pdb

	initialLoadState = LoadState_3Metals_3Sites()
	initialLoadState.init_with_input_data(inputData_microbe1)
		
	i = 1
	loadState_microbe1 = initialLoadState
	currentPurity = loadState_microbe1.fM1
	
	macroCycleDiagnostic = MacroCycleDiagnostic_3Metals_3Sites_3Microbes(inputData_microbe1)

	# Iterative sub-cycles for purification starts here. 
		
	# I'm going to cycle the solution through three microbes instead of just one. 
	# The idea is to use the first microbe to enrich for the metal I want, use the second to deplete
	# the second metal, and the third microbe to deplete the third metal. 
	
	while i <= inputData_microbe1.maxSubCycles and currentPurity < inputData_microbe1.targetPurity:
		
		print(i)
		
		# ---------------------------------------------------------------------------------------- #
		# Microbe 1	to Enrich for Metal 1
		bindState_microbe1 = BindState_3Metals_3Sites(loadState_microbe1)		
		
		# Calculate the contents of the eluant of the microbe 1 column
		eluteSolution_microbe1 = EluteSolution_3Metals()
		eluteSolution_microbe1.init_with_bind_state(bindState_microbe1, \
		inputData_microbe1.initial_eluteVol)
		
		washSolution_microbe1 = WashSolution_3Metals()
		washSolution_microbe1.init_with_bind_state(bindState_microbe1)
		print('Microbe 1 Eluant: ' + str(eluteSolution_microbe1.fM1))
		# ---------------------------------------------------------------------------------------- #

		# ---------------------------------------------------------------------------------------- #
		# Microbe 2 to Deplete Metal 2

		# Adjust the total number of binding sites in the microbe 2 column to match the number 
		# of M2 ions in the eluant from the microbe 1 column
		newNBT, newNB1, newNB2, newNB3, newLoadVol = \
		Calculate_New_Binding_Site_Numbers_and_Load_Volume_2(eluteSolution_microbe1,\
		inputData_microbe2, target='nM2')
		
		# Prepare load solution for microbe 2 and calculate wash solution
		loadState_microbe2 = LoadState_3Metals_3Sites()
		loadState_microbe2.init_with_elute_or_wash_solution_v2(eluteSolution_microbe1, \
		inputData_microbe2, newNB1, newNB2, newNB3, newLoadVol)

		bindState_microbe2 = BindState_3Metals_3Sites(loadState_microbe2)
		washSolution_microbe2 = WashSolution_3Metals()
		washSolution_microbe2.init_with_bind_state(bindState_microbe2)
		
		eluteSolution_microbe2 = EluteSolution_3Metals()
		eluteSolution_microbe2.init_with_bind_state(bindState_microbe2, \
		inputData_microbe2.initial_eluteVol)
		print('Microbe 2 Wash: ' + str(washSolution_microbe2.fM1))
		# ---------------------------------------------------------------------------------------- #
		
		# ---------------------------------------------------------------------------------------- #
		# Microbe 3 to Deplete Metal 3
		
		# Adjust the total number of binding sites in the microbe 3 column to match the number 
		# of M3 ions in the wash from the microbe 2 column
		newNBT, newNB1, newNB2, newNB3, newLoadVol = \
		Calculate_New_Binding_Site_Numbers_and_Load_Volume_2(washSolution_microbe2,\
		inputData_microbe3, target='nM3')
		
		# Prepare load solution for microbe 3 and calculate wash solution
		loadState_microbe3 = LoadState_3Metals_3Sites()
		loadState_microbe3.init_with_elute_or_wash_solution_v2(washSolution_microbe2, \
		inputData_microbe3, newNB1, newNB2, newNB3, newLoadVol)

		bindState_microbe3 = BindState_3Metals_3Sites(loadState_microbe3)
		washSolution_microbe3 = WashSolution_3Metals()
		washSolution_microbe3.init_with_bind_state(bindState_microbe3)
		
		eluteSolution_microbe3 = EluteSolution_3Metals()
		eluteSolution_microbe3.init_with_bind_state(bindState_microbe3, \
		inputData_microbe3.initial_eluteVol)
		print('Microbe 3 Wash: ' + str(washSolution_microbe3.fM1))
		# ---------------------------------------------------------------------------------------- #
		
		
		# ---------------------------------------------------------------------------------------- #
		# Prepare report out on cycle		
		subCycleDiagnostic = SubCycleDiagnostic_3M_3S_3Microbes(\
		loadState_microbe1, bindState_microbe1, washSolution_microbe1, eluteSolution_microbe1, \
		loadState_microbe2, bindState_microbe2, washSolution_microbe2, eluteSolution_microbe2, \
		loadState_microbe3, bindState_microbe3, washSolution_microbe3, eluteSolution_microbe3)

		macroCycleDiagnostic.subCycleDiagnosticArray.append(subCycleDiagnostic)
		# ---------------------------------------------------------------------------------------- #
				
		# ---------------------------------------------------------------------------------------- #
		# Prepare Microbe 1 for next cycle
		loadState_microbe1 = LoadState_3Metals_3Sites()
		
		newNBT, newNB1, newNB2, newNB3, newLoadVol = \
		Calculate_New_Binding_Site_Numbers_and_Load_Volume_2(washSolution_microbe3,\
		inputData_microbe1, target='nM1')
		
		loadState_microbe1.init_with_elute_or_wash_solution_v2(washSolution_microbe3, \
		inputData_microbe1, newNB1, newNB2, newNB3, newLoadVol)
		# ---------------------------------------------------------------------------------------- #
		
		currentPurity = washSolution_microbe3.fM1
		i += 1

			
	return macroCycleDiagnostic
# ------------------------------------------------------------------------------------------------ #




