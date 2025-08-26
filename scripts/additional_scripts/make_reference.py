import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

#### load data
ref = sys.argv[1]
fasta = list(SeqIO.parse(ref, "fasta"))

#### get contig lengths, sequences, ids and descriptions
lens = []; seqs = []; ids = []; desc = []
for i in fasta:
    lens.append(len(i))
    seqs.append(str(i.seq))
    ids.append(i.id)
    desc.append(i.description)

#### sort everything by length
sortedSeqs = [x for _, x in sorted(zip(lens, seqs))]
sortedIDs = [x for _, x in sorted(zip(lens, ids))]
sortedDesc = [x for _, x in sorted(zip(lens, desc))]
sortedLens = sorted(lens)

#### determine where is the first contig larger than 1000 bp
for i in range(len(sortedLens)):
    if sortedLens[i] > 1000:
        break

#### take only sequences larger than 1000 bp and reverse order such that it goes from largest to smallest contig
truncSeqs = (sortedSeqs[i:])[::-1]
truncIDs = (sortedIDs[i:])[::-1]
truncDesc = (sortedDesc[i:])[::-1]

#### write to new fasta file; write names and lengths of chromosomes to file
newFasta = []
o = open(ref.split(".f")[0] + "_LargerThan1000bp.chrs", "w")
for i in range(len(truncSeqs)):
    newFasta.append(SeqRecord(Seq(truncSeqs[i]), id=truncIDs[i], description=truncDesc[i]))
    o.write(truncIDs[i] + "\t" + str(len(truncSeqs[i])) + "\n")

SeqIO.write(newFasta, ref.split(".f")[0] + "_LargerThan1000bp.fasta", "fasta")

