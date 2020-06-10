from __future__ import division
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import IUPACData

import itertools
dipeptides = []
for i, j in itertools.product(list(IUPACData.protein_letters), list(IUPACData.protein_letters)) :
    dipeptides.append(i+j)

matrix = [list(IUPACData.protein_letters) + dipeptides + ['Label']]

for record in SeqIO.parse('na_thermo_test.fasta', 'fasta'):
    X = ProteinAnalysis(str(record.seq))
    line = [] 
    
    for letter in IUPACData.protein_letters :
    	line.append(X.get_amino_acids_percent()[letter])

    seq = str(record.seq)
    dps = []
    for i in range(0, len(seq)-1) :
    	dps.append(seq[i:i+2])
    for dp in dipeptides:
        line.append(dps.count(dp)/len(dps))

    matrix.append(line + ['thermo'])

for record in SeqIO.parse('na_meso_test.fasta', 'fasta'):
    X = ProteinAnalysis(str(record.seq))
    line = [] 
    
    for letter in IUPACData.protein_letters :
    	line.append(X.get_amino_acids_percent()[letter])

    seq = str(record.seq)
    dps = []
    for i in range(0, len(seq)-1) :
    	dps.append(seq[i:i+2])
    for dp in dipeptides:
        line.append(dps.count(dp)/len(dps))
    
    matrix.append(line + ['meso'])

import csv

with open("AAc_DPc_Test.csv","w+") as my_csv:
    csvWriter = csv.writer(my_csv, delimiter=',')
    csvWriter.writerows(matrix)
