'''
Isaac Gibbs
5/9/2014
FASTA Emissions Generator (Requires BioPython to be installed)
'''

from Bio import SeqIO
import sys
import csv
import math

input_file = open('lgvi.fasta', 'r') #FASTA file to be analyzed
gcr = [['Frame', 'GCR Content']]
NonCentromeric = [['Global Window', 'GC Content', 'G Count', 'C Count', 'A Count', 'T Count', 'GCR', 'gcrN', 'ngcrN']] #header for the non-centromeric emissions calculations
NonCentromeric1 = [[]]
del NonCentromeric1[:]
NonCentromeric2 = [[]]
del NonCentromeric2[:]
Centromeric = [['Global Window', 'GC Content', 'G Count', 'C Count', 'A Count', 'T Count', 'GCR', 'gcrC', 'ngcrC']] #header for the centromeric emissions calculations
Centromeric1 = [[]]
del Centromeric1[:]
Emissions = [[]]
del Emissions[:]
rounding = 5
line = 0
gcrN = 0.0
ngcrN = 0.0
gcrC = 0.0
ngcrC = 0.0
GCRcount = 0
windowCount = 0
progress = 0

def loadSequence(ifile): #reads data from the fasta file
	for seq_record in SeqIO.parse(ifile,"fasta"):
		record = seq_record
	print("___FASTA LOADED___")
	print("{}: {}".format("Sequence ID", record.id))
	print("{}: {}".format("Sequence Length", len(record)))
	print("\n")
	return record

def gcContentWindowed(sequence, data, bound1, bound2, window, tolerance, gcrthreshold): #Calculates the GCR emissions from bound1 to bound2 with a set window size a window overlap tolerance between (0, 1) and a set GCR threshold
	place = bound1 #the lower bound of the current window
	gcr = 0
	gcrcount = 0
	framesize = 0
	windowcount = 0
	global GCRcount
	global windowCount
	while (place + (window * (1 - tolerance))) < bound2:	#tolerence is how much of the window can push beyond bound2
		windowcount = windowcount + 1
		gcontent = sequence.seq[place: place + window].count('G')
		ccontent = sequence.seq[place: place + window].count('C')
		acontent = sequence.seq[place: place + window].count('A')
		tcontent = sequence.seq[place: place + window].count('T')
		contentpercent = ((gcontent + ccontent) / (window))
		if contentpercent > gcrthreshold:
			gcr = 1
			gcrcount = gcrcount + 1
		else:
			gcr = 0
		data.append(["{} {}".format(place, place + window - 1),"{}".format(round(contentpercent, rounding)), gcontent, ccontent, acontent, tcontent, gcr])
		place = place + window
	GCRcount = GCRcount + gcrcount
	windowCount = windowCount + windowcount

def centromereProb(sequence, bound1, bound2):

	for seq_record in SeqIO.parse(sequence,"fasta"):
		record = seq_record
	data.append([(bound2 - bound1) / len(record)])

#CALCULATES THE gcrN and ngcrN probabilities and writes out the emissions for the non-centromeric regions specified
'''gcContentWindowed(sequence, NonCentromeric1, 1, 2135113, 1000 , 0.75, 0.42)
gcContentWindowed(sequence, NonCentromeric2, 2482904, 4255303, 1000, 0.75, 0.42)
NonCentromeric = NonCentromeric + NonCentromeric1 + NonCentromeric2
gcrN = GCRcount/windowCount
ngcrN = 1 - (GCRcount/windowCount)
NonCentromeric.append(['','','','','','','',gcrN, ngcrN])
with open('c:\\NonCentromeric.csv', 'w', newline='') as fp:
    a = csv.writer(fp, delimiter=',')
    data = NonCentromeric
    a.writerows(data)'''

#CALCULATES THE gcrC and ngcrC probabilitiesand writes out the emissions for the centromeric region specified
'''GCRcount = 0
windowCount = 0
gcContentWindowed(sequence, Centromeric1, 2135983, 2478592, 1000 , 0.75, 0.42)
Centromeric = Centromeric + Centromeric1
gcrC = GCRcount/windowCount
ngcrC = 1 - (GCRcount/windowCount)
Centromeric.append(['','','','','','','',gcrC, ngcrC])
with open('c:\\Centromeric.csv', 'w', newline='') as fp:
    a = csv.writer(fp, delimiter=',')
    data = Centromeric
    a.writerows(data)'''

#Reads in the FASTA sequence
sequence = loadSequence(input_file) 

#Calculates emissions for the given sequence with the given arguments (fasta_sequence, variable_to_write_to, start_boundary, end_boundary, window_size, window_tolerance, gcr_threshold)
#Currently the Emissions.csv columns are [Current_Window, GC_Content, G_Count, C_Count, A_Count, T_Count, GCR] 
gcContentWindowed(sequence, Emissions, 1, len(sequence), 1000, 0.75, 0.42)
data = Emissions
with open('c:\\Emissions.csv', 'w', newline='') as fp: #writes emissions to Emissions.csv
    a = csv.writer(fp, delimiter=',')
    data = data
    a.writerows(data)



