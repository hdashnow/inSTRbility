import pdb
import os
import collections
import pysam
from interlap import InterLap

#               consumes_qry    consumes_ref
MATCH  = 0  # M yes             yes
INS    = 1  # I yes             no
DEL    = 2  # D no              yes
SKIP   = 3  # N no              yes
SOFT   = 4  # S yes             no
HARD   = 5  # H no              no
PAD    = 6  # P no              no
EQUAL  = 7  # = yes             yes
DIFF   = 8  # X yes             yes

def consumes_query(cigartuple):
    """If cigar operation consumes query return True, else return False"""
    return cigartuple[0] in [MATCH, INS, SOFT, EQUAL, DIFF]

def consumes_reference(cigartuple):
    """If cigar operation consumes reference return True, else return False"""
    return cigartuple[0] in [MATCH, DEL, SKIP, EQUAL, DIFF]

def translate(aln, pos, start = True):
        """translate a genomic reference coordinate to a read (aln) coordinate.
        if start is True, return leftmost, otherwise return rightmost"""
        inter = collections.defaultdict(InterLap)
        I = inter[aln.query_name]
        ref_name = aln.reference_name

        aln_off = 0
        ref_off = aln.reference_start
        for c in aln.cigartuples:
            op_len = c[1]
            cons_q = consumes_query(c)
            cons_r = consumes_reference(c)

            I.add((ref_off, ref_off + op_len * int(cons_r),
                dict(start=aln_off, stop=aln_off + op_len *
                    int(cons_q))))

            aln_off += int(cons_q) * op_len
            ref_off += int(cons_r) * op_len
        
        results = [(start, stop, read_dict) for (start, stop, read_dict) in I.find((pos, pos))]

        results.sort()

        # No match. This may indicate a non-spanning read
        if len(results) == 0:
            return None

        all_pos = []

        for start, stop, read_dict in results:
            over = pos - start # how far past start
            # print(f"start: {start}")
            # print(f"stop: {stop}")
            # print(f"over: {over}")
            # print(f"read_dict: {read_dict}")
            if read_dict["start"] + over > 0: # Not sure if this is correct
                all_pos.append(read_dict["start"] + over)


        if start:
            return all_pos[0]
        else:
            return all_pos[-1]


if __name__ == "__main__":
    import sys
    #c2r = Reference2Read(sys.argv[1])

