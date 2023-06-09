# Parse TRF ngs dat file into tab delimited file
# Assumes TRF was run on aligned contigs.
# Each @ is a contig labeled with it's mapping position,
# each line is an STR

import pdb
import argparse
import os
import sys
from ref2read import translate
import itertools
import pysam
import statistics as stats

def parse_args():
    # top-level parser
    parser = argparse.ArgumentParser()

    # subcommands
    parser.add_argument('BAM', type=str,
                        help='bam or cram file')
    parser.add_argument('--out', type=str, default = '',
                        help='output file of reads segments in fasta format')

    return parser.parse_args()

def get_ref_pos(pos, intervals):
    results = list(intervals.find((pos, pos)))
    if len(results) == 0:
        return
    all_results = []
    for result in results:
        result_contig = (result[0], result[1])
        result_ref = result[2]
        # Result reference is 1 bp
        if result_ref[0] == result_ref[1]:
            all_results.append(result_ref[0])
        # Result is a range so need to calculate position in range
        else:
            # If both ranges have a non-zero lenth, check the are equal
            if result_contig[0] != result_contig[1]:
                try:
                    assert result_ref[1] - result_ref[0] == result_contig[1] - result_contig[0]
                except AssertionError:
                    sys.stderr.write('WARNING, inconsistent range lengths:\n')
                    sys.stderr.write(f'result_ref {result_ref} len: {result_ref[1] - result_ref[0]}\n')
                    sys.stderr.write(f'result_contig {result_contig} len: {result_contig[1] - result_contig[0]}\n')
                    sys.stderr.write(f'result would be {result_ref[0] + pos - result_contig[0]}\n')
                    assert False
            all_results.append(result_ref[0] + pos - result_contig[0])
    # Check all results are equal
    if len(set(all_results)) == 1:
        return all_results[0]
    else:
        try:
            return stats.mode(all_results)
        except stats.StatisticsError:
            return round(stats.mean(all_results))

def write_read(outfile, sequence, SeqID):
    """Add read segment to fasta file"""
    outfile.write(f'>{SeqID}\n{sequence}\n')

def extract_read_chunks(samfilename, outfilename = ''):
    """
    """
    if outfilename != '':
        outfile = open(outfilename, 'w')
    #ref2readobj = Reference2Read(samfile)
    samfile = pysam.AlignmentFile(samfilename, "rc")

    locus = ('chr4', 3074876, 3074940) #HD locus hardcoded for now
    for read in samfile.fetch(*locus):
        if read.is_supplementary:
            continue
        #pdb.set_trace()
        # start is the position in the read corresponding to the start of the locus in the reference
        start = translate(read, locus[1], True)
        end = translate(read, locus[2], False)

        # Skip reads that end within the locus
        if start == None or end == None:
            continue
        # Include the portions of reads that end within the locus
#        if start == None:
#            start = 0
#        if end == None:
#            end = read.query_length # Need to +1?

        # extract sequence of read between those coordinates
        read_chunk = read.query_sequence[start:end]

        if outfilename != '':
            write_read(outfile, read_chunk, read.query_name)

def main():
    args = parse_args()

    extract_read_chunks(args.BAM, args.out)

if __name__ == '__main__':
    main()
