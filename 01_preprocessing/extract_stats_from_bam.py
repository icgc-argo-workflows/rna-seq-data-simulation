import sys
import pysam
import re
import numpy as np
import numpy.random as npr
import pickle
npr.seed(23)

transitions = {('A', 'C'):0,
               ('A', 'G'):1,
               ('A', 'T'):2,
               ('C', 'A'):3,
               ('C', 'G'):4,
               ('C', 'T'):5,
               ('G', 'A'):6,
               ('G', 'C'):7,
               ('G', 'T'):8,
               ('T', 'A'):9,
               ('T', 'C'):10,
               ('T', 'G'):11,
               ('A', 'N'):12,
               ('C', 'N'):13,
               ('G', 'N'):14,
               ('T', 'N'):15}

def main():
    
    ### parse command line options
    if len(sys.argv) < 4:
        sys.stderr.write('Usage: %s <alignment.bam> <variants.vcf> <out_fname> [--dump-blocks]\n' % sys.argv[0])
        sys.exit(1)
    fname_bam = sys.argv[1]
    fname_vcf = sys.argv[2]
    outfname = sys.argv[3]
    outfname_pickle = sys.argv[3] + '.pickle'
    dump_blocks = '--dump-blocks' in sys.argv
    if dump_blocks:
        out = open(outfname, 'w')

    qmatrix = np.zeros((40, 101), dtype='int')
    insertions_len = np.zeros((101, ), dtype='int')
    insertions_pos = np.zeros((101, ), dtype='int')
    deletions_len = np.zeros((101, ), dtype='int')
    deletions_pos = np.zeros((101, ), dtype='int')
    mismatches = np.zeros((16, 40, 101), dtype='int')
    total_bases = 0

    ### load likely variant positions (to be ignored in the mismatch stats)
    variants = dict()
    for line in open(fname_vcf, 'r'):
        if line[0] == '#':
            continue
        sl = line.strip().split('\t')
        chrm = sl[0]
        pos = int(sl[1]) - 1 ### vcf is 1-based
        if not chrm in variants:
            variants[chrm] = set()
        for i in range(pos, pos + len(sl[3])):
            variants[chrm].add(i)

    ### process alignments
    with pysam.AlignmentFile(fname_bam, 'rb') as bamfile:
        linecounter = 0
        for read in bamfile:

            linecounter += 1
            if linecounter % 10000 == 0:
                sys.stdout.write('.')
                if linecounter % 100000 == 0:
                    sys.stdout.write('%i\n' % linecounter)
                sys.stdout.flush()

            if read.is_unmapped:
                continue
            if read.is_secondary:
                continue

            total_bases += len(read.seq)

            if npr.randint(0, 100) < 2:
                ### build mask for updating qmatrix
                qs = np.diag(np.ones((40,), dtype='int'))[:, [ord(_) - 66 for _ in read.qual]]
                ### update qmatrix, adjusting for read length differences
                qmatrix += np.c_[qs, np.zeros((40, 101 - len(read.qual)), dtype='int')]

            ### get blocks
            gen_pos = read.pos ### first position in alignment consuming the reference
            read_pos = 0
            insertion_starts = [(0, 0)]
            exon_starts = [(0, 0)]
            blocks = [[gen_pos, gen_pos]]
            chrm = read.reference_name
            for op in read.cigar:
                ### N --> 3
                if op[0] == 3:
                    exon_starts.append((gen_pos, op[1]))
                    blocks[-1][1] = gen_pos
                    blocks.append([gen_pos + op[1], gen_pos + op[1]])
                ### I --> 1
                elif op[0] == 1:
                    ### INSERTION
                    if not chrm in variants or not gen_pos in variants[chrm]:
                        ### report 0-based genome position of last matching base
                        insertions_len[op[1]] += 1
                        insertions_pos[read_pos] += 1
                    insertion_starts.append((read_pos, op[1]))

                if not op[0] in [1, 4, 5, 6]:
                    gen_pos += op[1]
                if not op[0] in [2, 3, 5, 6]:
                    read_pos += op[1]
            blocks[-1][1] = gen_pos
            
            if dump_blocks:
                out.write('\t'.join([read.qname, str(int(read.is_reverse)), read.reference_name, ':'.join([','.join(_) for _ in np.array(blocks, dtype='str')])]) + '\n')

            ### collect deletion and mismatch info
            md = read.get_tag('MD')

            if md.isdigit():
                continue

            if not md[-1].isdigit():
                md += '0'
            (match, var) = ([int(i) for i in re.split('[^0-9]+', md)], re.split('[0-9]+', md))
            var = var[1:-1]
            assert(len(match) >= len(var))
     
            read_pos = 0
            if read.cigar[0][0] == 4:
                read_pos += read.cigar[0][1]
            gen_pos = read.pos
            e_idx = 1
            i_idx = 1
            for o in range(len(match)):
                gen_pos += match[o]
                read_pos += match[o]
                if o < len(var):
                    ### get intron offset
                    if len(exon_starts) > 1:
                        while e_idx < len(exon_starts) and gen_pos >= exon_starts[e_idx][0]:
                            gen_pos += exon_starts[e_idx][1]
                            e_idx += 1
                    ### get insertion offset
                    if len(insertion_starts) > 1:
                        while i_idx < len(insertion_starts) and read_pos >= insertion_starts[i_idx][0]:
                            read_pos += insertion_starts[i_idx][1]
                            i_idx += 1
                        
                    if var[o][0] == '^':
                        ### DELETION
                        ### check if the deletion spans an intron
                        ### if yes, split it up accordingly
                        del_len = (len(var[o]) - 1)
                        if not chrm in variants or not gen_pos in variants[chrm]:
                            deletions_len[del_len] += 1
                        deletions_pos[read_pos] += 1
                        if del_len > 1 and len(exon_starts) > 1 and e_idx < len(exon_starts) and gen_pos + del_len > exon_starts[e_idx][0]:
                            before_len = exon_starts[e_idx][0] - gen_pos
                            gen_pos += before_len
                            gen_pos += exon_starts[e_idx][1]
                            e_idx += 1
                            gen_pos += (del_len - before_len)
                        else:
                            gen_pos += del_len
                    else:
                        ### SNP (REF --> Read):
                        if var[o] != 'N' and (not chrm in variants or not gen_pos in variants[chrm]):
                            mismatches[transitions[(var[o], read.seq[read_pos])], ord(read.qual[read_pos]) - 66, read_pos] += 1
                        read_pos += 1
                        gen_pos += 1
            
        pickle.dump((qmatrix,  mismatches,  insertions_len,  insertions_pos,  deletions_len, deletions_pos, total_bases), open(outfname_pickle, 'wb'), -1) 

        if dump_blocks:
            out.close()

if __name__ == "__main__":
    main()
