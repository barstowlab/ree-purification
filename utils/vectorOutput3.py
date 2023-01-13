# ---------------------------------------------------------------------------- #
def generateOutputMatrix(vectorList, delimeter='\t'):
	
	outputMatrix = []
	outputString = ''
	i = 0
	
	
	
	while i < len(vectorList[0]):
		
		outputString = ''
		j = 0
		
		while j < len(vectorList):
			outputString += str(vectorList[j][i])
			if j < len(vectorList) - 1:
				outputString += delimeter
			j += 1
		
		outputString += '\n'
		outputMatrix.append(outputString)
		i += 1
	
	return outputMatrix
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
def generateOutputMatrixWithHeaders(vectorList, headers, delimeter=','):
	
	if len(headers) != len(vectorList):
		print("Header array length != vector list length")
		return
	
	outputMatrix = []
	
	k = 0
	headerString = ''
	while k < len(headers):
		headerString += str(headers[k]) + delimeter
		k += 1
	
	headerString += '\n'
	outputMatrix.append(headerString)
	
	outputString = ''
	i = 0
	
	while i < len(vectorList[0]):
		
		outputString = ''
		j = 0
		
		while j < len(vectorList):
			outputString += str(vectorList[j][i])
			if j < len(vectorList) - 1:
				outputString += delimeter
			
			j += 1
		
		outputString += '\n'
		outputMatrix.append(outputString)
		i += 1
	
	return outputMatrix
# ---------------------------------------------------------------------------- #


# ----------------------------------------------------------------------------------------------- #
def generateOutputMatrix_WithNonEqualVectors_WithHeaders(vectorList, headers, delimeter=','):
	
	import pdb
	
	if len(headers) != len(vectorList):
		print("Header array length != vector list length")
		return
	
	outputMatrix = []
	
	k = 0
	headerString = ''
	while k < len(headers):
		headerString += str(headers[k])
		
		if k < len(vectorList) - 1:
			headerString += delimeter
		
		k += 1
	
	headerString += '\n'
	outputMatrix.append(headerString)
	
	outputString = ''
	i = 0
	
	# Figure out the longest vector in the vector list
	longestVectorLength = 0
	i = 0
	while i < len(vectorList):
		lengthVector = len(vectorList[i])
		if lengthVector > longestVectorLength:
			longestVectorLength = lengthVector
		i += 1
	
	
	i = 0
	while i < longestVectorLength:
		
		outputString = ''
		j = 0
		
		
		while j < len(vectorList):
			lineLength = len(vectorList[j])
			
			try:
				outputString += str(vectorList[j][i])
			except:
				outputString += ''
			
			if j < len(vectorList) - 1:
				outputString += delimeter
			
			j += 1
		
		outputString += '\n'
		outputMatrix.append(outputString)
		i += 1
	
	return outputMatrix
# ----------------------------------------------------------------------------------------------- #





# ---------------------------------------------------------------------------- #
def writeOutputMatrix(fileName, outputMatrix):
	outputFile = open(fileName, 'w')
	i = 0
	while i < len(outputMatrix):
		outputFile.write(outputMatrix[i])
		i+=1
	return	
# ---------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------------------ #
def ensure_dir(f):
    import os
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
# ------------------------------------------------------------------------------------------------ #
