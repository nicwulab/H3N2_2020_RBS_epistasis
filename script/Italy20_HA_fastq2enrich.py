#!/usr/bin/python
import os
import sys
import string
import operator
from Bio import SeqIO
from collections import Counter

def hamming(str1, str2):
    assert len(str1) == len(str2)
    return sum(map(operator.ne, str1, str2))

#def rc(seq):
  #seq = str(seq)
  #complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
  #rcseq = seq.translate(complements)[::-1]
  #return rcseq

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def ProcessMultilib(merged_file, trim_5, trim_3):
  print ("Reading %s" % merged_file)
  records = SeqIO.parse(merged_file,"fastq")
  variants = [] 
  read_count = 0
  for record in records:
    read_count += 1
    seq  = record.seq
    if len(seq) != 369: continue
    trim_seq  = seq[trim_5:-trim_3]
    if 'N' in trim_seq: continue
    pep = translation(trim_seq)
    variants.append(pep)
    #if read_count == 100: break
  return Counter(variants)

def Ita20HAmut2ID(mut,refseq):
  shift = 120
  haplo = []
  assert(len(mut)==len(refseq))
  for n in range(len(mut)):
    pos = n+shift
    if refseq[n]!=mut[n]:
       haplo.append(refseq[n]+str(pos)+mut[n])
  return '-'.join(haplo)

def Output(InputDict, Rep1Dict, Rep2Dict, outfile, refseq):
  print ("Compiling results into %s" % outfile)
  outfile = open(outfile,'w')
  muts = list(set(list(InputDict.keys())+list(Rep1Dict.keys())+list(Rep2Dict.keys())))
  Ipt_total_count  = sum(InputDict.values())
  Rep1_total_count = sum(Rep1Dict.values())
  Rep2_total_count = sum(Rep2Dict.values())
  outfile.write("\t".join(['mut','mutclass','InputCount',
                           'Rep1Count','Rep2Count',
                           'Rep1Freq','Rep2Freq',
                           'Rep1Enrich','Rep2Enrich'])+"\n")
  for mut in muts:
    IptFreq  = float(InputDict[mut]+1)/float(Ipt_total_count)
    Rep1Freq = float(Rep1Dict[mut]+1)/float(Rep1_total_count)
    Rep2Freq = float(Rep2Dict[mut]+1)/float(Rep2_total_count)
    Rep1Enrich = Rep1Freq/IptFreq
    Rep2Enrich = Rep2Freq/IptFreq
    mutclass    = 0 if refseq == mut else hamming(refseq,mut)
    ID          = 'WT' if refseq == mut else Ita20HAmut2ID(mut,refseq)
    outfile.write("\t".join(map(str,[ID, mutclass, InputDict[mut],
                                     Rep1Dict[mut], Rep2Dict[mut],
                                     Rep1Freq, Rep2Freq,
                                     Rep1Enrich, Rep2Enrich]))+"\n")
  outfile.close()

def main():
  reffile = 'Italy20HA_mutlib_ref.fasta' 
  outfile = '../results/Ita20HA_MultiMutLib.tsv'
  trim_5  = 25
  trim_3  = 23
  refseq  = next(SeqIO.parse(reffile,"fasta")).seq
  Ita20HA_InputDict   = ProcessMultilib('../fastq_merged/Italy20_DMSlib_input_AGTTCC_L001_merged_001.assembled.fastq', trim_5, trim_3)
  Ita20HA_Rep1Dict    = ProcessMultilib('../fastq_merged/Italy20_DMSlib_1_ATGTCA_L001_merged_001.assembled.fastq', trim_5, trim_3)
  Ita20HA_Rep2Dict    = ProcessMultilib('../fastq_merged/Italy20_DMSlib_2_CCGTCC_L001_merged_001.assembled.fastq', trim_5, trim_3)
  Output(Ita20HA_InputDict,Ita20HA_Rep1Dict,Ita20HA_Rep2Dict,outfile,refseq)

if __name__ == "__main__":
  main()
