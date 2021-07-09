# asb_assess
This is a script to assess assemble with reference sequence.

for Assessment.py:
usage: Assessment.py [-h] -i INPUT -r REFERENCE [-k KMER_LENGTH] -o OUT_PREFIX

  -i INPUT, --input INPUT                             assemble result with fasta format.
  -r REFERENCE, --reference REFERENCE                 reference sequence with fasta format.
  -k KMER_LENGTH, --kmer-length KMER_LENGTH           the kmer length used in assessment, default=21
   -o OUT_PREFIX, --out-prefix OUT_PREFIX             prefix of transforming assemble result and reference sequence into fasta format

The result will be writen in OUT_PREFIX.result file locally.
