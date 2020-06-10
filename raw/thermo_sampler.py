from Bio import SeqIO
import random

sequences = list(SeqIO.parse(open('cdhit_output_thermo_merged_40.fasta'),'fasta'))
print(len(sequences))
random.shuffle(sequences)

samples = random.sample(sequences, 1407)

random.shuffle(samples)

test = samples[:469]
train = samples[469:]

with open('test_thermo.fasta', "w") as output_handle:
    SeqIO.write(test, output_handle, "fasta")
with open('train_thermo.fasta', "w") as output_handle:
    SeqIO.write(train, output_handle, "fasta")
