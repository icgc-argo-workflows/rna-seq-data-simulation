import sys
import pickle
import gzip
import numpy as np
import numpy.random as npr
npr.seed(23)
import re

REV_MAP = {'A':'T',
           'T':'A',
           'G':'C',
           'C':'G',
           'N':'N',
           'a':'t',
           't':'a',
           'g':'c',
           'c':'g',
           'n':'n'}
ALPH = 'ACGT'
TRANS = [('A', 'C'),
         ('A', 'G'),
         ('A', 'T'),
         ('C', 'A'),
         ('C', 'G'),
         ('C', 'T'),
         ('G', 'A'),
         ('G', 'C'),
         ('G', 'T'),
         ('T', 'A'),
         ('T', 'C'),
         ('T', 'G'),
         ('A', 'N'),
         ('C', 'N'),
         ('G', 'N'),
         ('T', 'N')]
MUT_MAP = {'A': [0, 1, 2, 12],
           'C': [3, 4, 5, 13],
           'G': [6, 7, 8, 14],
           'T': [9, 10, 11, 15]}

def _rev_comp(s):
    """computre reverse complement of s"""

    return ''.join([REV_MAP[_] for _ in s][::-1])


def _gen_qvalues(l, n, qprobs, qrange, offset=66):
    qv = np.zeros((n, l), dtype='|S1')
    for i in range(l):
        qv[:, i] = [chr(_ + offset) for _ in npr.choice(qrange, n, p=qprobs[:, i])]
    return [''.join(_) for _ in qv.astype('str')]


def _get_random_pos(probs, length):
    try:
        return npr.choice(np.arange(length), 1, p=probs[:length])[0]
    except:
        ### choose uniformly amongst positions if probablity vector does not sum to one (should rarely happen)
        return npr.choice(np.arange(length), 1)[0]


def _get_random_len(probs):
    return npr.choice(np.arange(len(probs)) + 1, 1, p=probs)[0]


def _get_random_seq(length):
    return ''.join([ALPH[npr.randint(4)] for _ in range(length)]) 
    

def _get_alternate(q, r, pos, mismatches):
    
    probs = mismatches[:, ord(q) - 66, pos][MUT_MAP[r]]
    if np.sum(probs) == 0:
        idx = npr.choice(MUT_MAP[r], 1)[0]
    else:
        probs = probs / np.sum(probs).astype('float')
        idx = npr.choice(MUT_MAP[r], 1, p=probs)[0]
    return TRANS[idx][1]


def _get_expected_readlen(fname, thresh=10000):
    
    lens = []
    for line in gzip.open(fname, 'rt', encoding='utf-8'):
        if line.startswith('>'):
            continue
        lens.append(len(line.strip()))
        if len(lens) >= thresh:
            break
    return int(np.median(lens))


def _mutate_read(read_seq, qs, mposprobs, mismatches):

    ### get mutated position
    ### (in the offchance that the read_seq is shorter, we just mutate the last position; 
    ### this should only happen, if a prior deletion has shortened our read due to being close
    ### to the transcript's start or end - this is very rare)
    mpos = min(_get_random_pos(mposprobs, len(read_seq)), len(read_seq) - 1)

    if read_seq[mpos] == 'N':
        return read_seq

    ### sample mutation according to quality value and current base:
    alternate = _get_alternate(qs[mpos], read_seq[mpos], mpos, mismatches)

    ### return new seq
    if mpos == 0:
        return alternate +  read_seq[mpos + 1:]
    elif mpos == len(read_seq) - 1:
        return read_seq[:mpos] + alternate
    else:
        return read_seq[:mpos] + alternate + read_seq[mpos + 1:]

def _get_del_compensation(read_header, read_seq, dlen, txs):
    ### get header info
    tx_id, mate1, mate2, this = read_header.split(';')
    ### get transcript sequence
    tx_seq = txs[tx_id.split('/')[1]]
    if this.split(':')[0] == 'mate1':
        start, end = mate1.split(':')[1].split('-')
    else:
        start, end = mate2.split(':')[1].split('-')
    start = int(start) - 1
    end = int(end)

    ### check that we are were we are supposed to be
    ### and generate compensation sequences
    if this.split(':')[1] == 'fwd':
        assert tx_seq[start:end] == read_seq
        return tx_seq[end:end + dlen]
    else:
        assert tx_seq[start:end] == _rev_comp(read_seq)
        return _rev_comp(tx_seq[max(0, start - dlen):start])


def _handle_read(read_seq, read_header, txs, dprob, dposprobs, dlenprobs, iprob, iposprobs, ilenprobs, qs, mposprobs, mismatches):

    readlen = len(read_seq)

    ### introduce deletions
    if npr.random() < dprob:
        dpos = min(_get_random_pos(dposprobs, dposprobs.shape[0]), readlen - 1)
        dlen = min(_get_random_len(dlenprobs), readlen - dpos)
        comp = _get_del_compensation(read_header, read_seq, dlen, txs)
        if read_header.split(':')[-1] == 'fwd':
            read_seq = read_seq[:dpos] + read_seq[dpos+dlen:] + comp
        else:
            read_seq = comp + read_seq[:dpos] + read_seq[dpos+dlen:]

    ### introduce insertions
    if npr.random() < iprob:
        ipos = min(_get_random_pos(iposprobs, iposprobs.shape[0]), readlen - 1)
        ilen = _get_random_len(ilenprobs)
        iseq = _get_random_seq(ilen)
        read_seq = read_seq[:ipos] + iseq + read_seq[ipos:]
        if read_header.split(':')[-1] == 'fwd':
            read_seq = read_seq[:-ilen]
        else:
            read_seq = read_seq[ilen:]

    ### introduce mutations
    while npr.random() < mutrate:
        read_seq = _mutate_read(read_seq, qs, mposprobs, mismatches)

    return (read_seq, qs)


def _print_record(readid, read_seq, qs, outstream=sys.stdout):
    outstream.write('@%s\n%s\n+\n%s\n' % (readid, read_seq, qs))

def _load_transcripts(tx_fa):

    txs = dict()
    seq = []
    for line in open(tx_fa, 'r'):
        if line.startswith('>'):
            if len(seq) > 0:
                txs[header] = ''.join(seq)
            seq = []
            header = line.strip()[1:]
            continue
        seq.append(line.strip())
    if len(seq) > 0:
        txs[header] = ''.join(seq)
    return txs


if len(sys.argv) < 4:
    sys.stderr.write('Usage: %s <fasta> <quality_stats> <transcripts.fa>\n' % sys.argv[0])
    sys.exit(1)
fasta = fasta = sys.argv[1]
qs_pickle = sys.argv[2]
tx_fasta = sys.argv[3]

### get transcript sequences
txs = _load_transcripts(tx_fasta)

### load quality information
(qmatrix,  mismatches,  insertions_len,  insertions_pos,  deletions_len, deletions_pos, total_bases) = pickle.load(open(qs_pickle, 'rb'))
qsum = np.sum(qmatrix, axis=0).astype('float')
qsum[qsum == 0] = 1
qprobs = qmatrix/qsum
qrange = np.arange(qmatrix.shape[0])
iposprobs = insertions_pos / np.sum(insertions_pos).astype('float')
ilenprobs = insertions_len / np.sum(insertions_len)
iprob = np.sum(insertions_len) / np.sum(total_bases)
dposprobs = deletions_pos / np.sum(deletions_pos).astype('float')
dlenprobs = deletions_len / np.sum(deletions_len)
dprob = np.sum(deletions_len) / np.sum(total_bases)
mutrate = np.sum(mismatches) / float(total_bases)
mposprobs = np.sum(np.sum(mismatches, axis=1), axis=0) / np.sum(mismatches).astype('float')

### estimate expected read length
expected_readlen = _get_expected_readlen(fasta)

### iterate through summarized, sorted bam file
### iterate through input fasta file
reads = []
qvalues = []
out = gzip.open(re.sub(r'.fasta.gz$', '', fasta) + '.fastq.gz', 'wt', encoding='utf-8')
for l, line in enumerate(gzip.open(fasta, 'rt', encoding='utf-8')):

    if l > 0 and l % 10000 == 0:
        sys.stdout.write('.')
        if l % 100000 == 0:
            sys.stdout.write('%i\n' % (l))
        sys.stdout.flush()
    if line.startswith('>'):
        read_header = line[1:].strip()
        continue
    
    read_seq = line.strip()

    ### check whether we have enough quality values left
    if len(qvalues) < 2:
        qvalues = _gen_qvalues(expected_readlen, 10000, qprobs, qrange)

    ### introduce sequencing errors into the read
    (read_seq, qs) = _handle_read(read_seq,
                                 read_header,
                                 txs,
                                 dprob, 
                                 dposprobs,
                                 dlenprobs,
                                 iprob,
                                 iposprobs,
                                 ilenprobs,
                                 qvalues.pop(),
                                 mposprobs,
                                 mismatches)
    qs = qs[:len(read_seq)]
    assert len(read_seq) == len(qs)
    assert(len(read_seq) <= expected_readlen)
    _print_record(read_header, read_seq, qs, out)
        
