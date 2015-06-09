'''
Isaac Gibbs
5/9/2014
FASTA State_Path_Probability Generator
'''

#dnaAnalyzer
import sys
import csv
import math

emissions = [] #emissions data
del emissions [:]
gcrN = 0.87308085977482 #probability of GC rich non-centromeric (calculated with Emissions_Calculations.py)
ngcrN = 0.126919140225179 #probability of not GC rich non-centromeric (calculated with Emissions_Calculations.py)
gcrC = 0.00291545189504373 #probability of GC rich centromeric (calculated with Emissions_Calculations.py)
ngcrC = 0.997084548104956 #probability of not GC rich centromeric (calculated with Emissions_Calculations.py)
ln0 = (2 * math.log(1/9)) + (7 * math.log(8/9)) 
progress = 0 #used for printing the starting location of the current set of state paths being calculated
header = [['C1', 'C2', 'Result', "GCR N: " + repr(gcrN), "!GCR N: " + repr(ngcrN), "GCR C: " + repr(gcrC), "!GCR C: " + repr(ngcrC)]] #header for state_path_probabilities.csv
writer = csv.writer(open('c:\\State_Path_Probabilities.csv', 'w', newline=''), delimiter = ',') #creates state_path_probabilities.csv
writer.writerows(header) #appends header to state_path_probabilities.csv
writer = csv.writer(open('c:\\State_Path_Probabilities.csv', 'a', newline=''), delimiter = ',') #opens a writer for state_path_probabilities.csv (must be the same path as the file creation above)
	
with open('c:\\Emissions.csv') as csvfile: #reads the emissions file from the path
    a = [ row.strip().split(',') for row in csvfile]
    die = dict((data[0],data[1:]) for data in a)

for i in a:
	emissions.append(int(i[6])) #adds emissions from file to the emissions variable
print(len(a)) #prints the length of the emissions found

def calculateState(data, gcrN, ngcrN, gcrC, ngcrC): #statepath is two dimensional array of states (non-centromeric, centromeric) and emissions(gcr, !gcr) 

	ntype0 = 0
	ntype1 = 0
	ctype0 = 0
	ctype1 = 0
	global ln0
	for i in data:	#counts the number of different state path types
		if i[0] is 'n':
			if i[1] is 0:
				ntype0 = ntype0 + 1
			if i[1] is 1:
				ntype1 = ntype1 + 1
		if i[0] is 'c':
			if i[1] is 0:
				ctype0 = ctype0 + 1
			if i[1] is 1:
				ctype1 = ctype1 + 1
	ln1 = math.log(ngcrN) * ntype0
	ln2 = math.log(gcrN) * ntype1
	ln3 = math.log(ngcrC) * ctype0
	ln4 = math.log(gcrC) * ctype1
	result = ln0 + ln1 + ln2 + ln3 + ln4

	return result

def stateProbabilities(emissions):
	pathlength = len(emissions)
	states = [] #statepath tuple, (statepath chars, emissions)
	pathtuples = []
	results = []
	statepath = []
	startindex = 1
	global writer
	global progress
	i = startindex
	while i < pathlength: #creates tuples representing all possible paths given a pathLength
		j = i + 1
		while j < pathlength:
			pathtuples.append((i, j))
			j = j + 1
		i = i + 1

	for i in pathtuples: #creates all possible state paths given a pathLength
		if i[0] != progress:
			print(i[0])
			progress = i[0]
		statepath = []
		for x in range(pathlength):
			if x >= i[0]:
				if x < i[1]:
					statepath.append('c')
				else:
					statepath.append('n')
			else:
				statepath.append('n')
		
		result = calculateState(zip(statepath,emissions), gcrN, ngcrN, gcrC, ngcrC) #calculates the probability of a given state path
		states.append((i, result))
		if len(states) > 500: #writes to state_path_probabilities.csv after 500 calculations (prevents running out of ram)
			lines = []
			del lines[:]
			for i in states:
				lines.append([i[0][0], i[0][1], i[1]])
			writer.writerows(lines)
			states = []
	lines = []
	del lines[:]
	for i in states: #writes any remaining data to state_path_probabilities.csv after all states are calculated
		lines.append([i[0][0], i[0][1], i[1]])
	writer.writerows(lines)
	states = []


stateProbabilities(emissions) #Calculates the state path probabilities given emission data and GCR data provided at the top
