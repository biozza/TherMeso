from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import random

frag_train_thermo = {'Hi'}
frag_train_meso = {'Hi'}

for record in SeqIO.parse('train_thermo_single.fasta', 'fasta'):
    #TP_seq.append(str(record.seq))
    seq = str(record.seq)
    for i in range(0, len(seq)-4) :
    	frag_train_thermo.add(seq[i:i+5])

for record in SeqIO.parse('train_meso_single.fasta', 'fasta'):
    #MP_seq.append(str(record.seq))
    seq = str(record.seq)
    for i in range(0, len(seq)-4) :
    	frag_train_meso.add(seq[i:i+5])

train_thermo_prune = frag_train_thermo.difference(frag_train_meso)
train_meso_prune = frag_train_meso.difference(frag_train_thermo)

print(len(frag_train_thermo), len(frag_train_meso))
print(len(train_thermo_prune), len(train_meso_prune))

thermo_thermo = 0
thermo_meso = 0
thermo_na = 0
meso_meso = 0
meso_thermo = 0
meso_na = 0

na_ter = []

for record in SeqIO.parse('test_thermo_single.fasta', 'fasta'):
    #TP_seq.append(str(record.seq))
    frag_test_thermo = {'Hi'}
    seq = str(record.seq)
    for i in range(0, len(seq)-4) :
    	frag_test_thermo.add(seq[i:i+5])
    
    tp_score = len(frag_test_thermo.intersection(train_thermo_prune))
    mp_score = len(frag_test_thermo.intersection(train_meso_prune))

    if tp_score >= 29 :
        thermo_thermo += 1
    elif mp_score >= 28 :
        thermo_meso += 1
    else :
        thermo_na += 1
        na_ter.append(record)

with open("na_thermo_test.fasta", "w") as output_handle:
    SeqIO.write(na_ter, output_handle, "fasta")
    
na_mes = []

for record in SeqIO.parse('test_meso_single.fasta', 'fasta'):
    #TP_seq.append(str(record.seq))
    frag_test_meso = {'Hi'}
    seq = str(record.seq)
    for i in range(0, len(seq)-4) :
    	frag_test_meso.add(seq[i:i+5])
    
    tp_score = len(frag_test_meso.intersection(train_thermo_prune))
    mp_score = len(frag_test_meso.intersection(train_meso_prune))

    if tp_score >= 29 :
        meso_thermo += 1
    elif mp_score >= 28 :
        meso_meso += 1
    else :
        meso_na += 1
        na_mes.append(record)

with open("na_meso_test.fasta", "w") as output_handle:
    SeqIO.write(na_mes, output_handle, "fasta")

print(thermo_thermo, thermo_meso, thermo_na, len(na_ter))
print(meso_thermo, meso_meso, meso_na, len(na_mes))

