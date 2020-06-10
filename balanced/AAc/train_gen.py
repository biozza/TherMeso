from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import IUPACData

matrix = [list(IUPACData.protein_letters) + ['Label']]

for record in SeqIO.parse('train_thermo.fasta', 'fasta'):
    X = ProteinAnalysis(str(record.seq))
    line = [] 
    for letter in IUPACData.protein_letters :
    	line.append(X.get_amino_acids_percent()[letter])
    matrix.append(line + ['thermo'])

for record in SeqIO.parse('train_meso.fasta', 'fasta'):
    X = ProteinAnalysis(str(record.seq))
    line = [] 
    for letter in IUPACData.protein_letters :
    	line.append(X.get_amino_acids_percent()[letter])
    matrix.append(line + ['meso'])

import csv

with open("AAc_Train.csv","w+") as my_csv:
    csvWriter = csv.writer(my_csv, delimiter=',')
    csvWriter.writerows(matrix)
