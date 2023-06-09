# inSTRbility

Contributers: Harriet Dashnow, Michael Goldberg

inSTRbility: STR instability

Takes aligned cram input. Extracts the portions of reads that overlap a given locus. Writes these as a fasta file. These will be used to measure the degree of STR instability at a given locus. This is being developed with ONT data.

Usage:
`python extract-region.py --out sample.fasta sample.cram`
