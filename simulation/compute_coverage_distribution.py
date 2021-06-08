import sys
import os
import pandas as pd
import numpy as np
import numpy.random as npr
import scipy.stats as spst
import matplotlib
import matplotlib.pyplot as plt
import pickle
import re

if len(sys.argv) < 9:
    sys.stderr.write('Usage: %s <gene_counts.tsv> <icgc_sample_id> <gene_lengths.tsv> <annotation.gtf> <transcripts.fa> <ID_tag> <outdir> <random_seed>\n' % sys.argv[0])
    sys.exit(1)

fname_exp = sys.argv[1]
icgc_sample = sys.argv[2]
fname_lens = sys.argv[3]
fname_gtf_sim = sys.argv[4]
fname_fa_sim = sys.argv[5]
idtag = sys.argv[6]
outdir = sys.argv[7]
random_seed = int(sys.argv[8])

REPLICATES = 3
NUM_SOM_ASE = 50
NUM_GERM_ASE = 500
NUM_DGE = 1500

outbase = os.path.join(outdir, re.sub(r'.fa$', '', os.path.basename(fname_fa_sim))) + '.' + idtag

### these files come from pre-processing
fname_germline_pickle = outbase + '.tx_hets_germline.pickle'
fname_somatic_pickle = outbase + '.tx_hets_somatic.pickle'
### output filenames
fname_out_readcounts = outbase + '.read_count.csv'
fname_out_factors_normal_ht1 = outbase + '.factors_normal_ht1.csv'
fname_out_factors_normal_ht2 = outbase + '.factors_normal_ht2.csv'
fname_out_factors_tumor_ht1 = outbase + '.factors_tumor_ht1.csv'
fname_out_factors_tumor_ht2 = outbase + '.factors_tumor_ht2.csv'
fname_out_metadata_all = outbase + '.metadata_complete.tsv'

### get data from expression
gene_lens = pd.read_csv(fname_lens, index_col=0, sep='\t')
gene_exp = pd.read_csv(fname_exp, index_col=0, usecols=["feature", icgc_sample], sep='\t')

### remove uncounted features
kidx = np.where([not _.startswith('__') for _ in gene_exp.index])[0]
gene_exp = gene_exp.iloc[kidx]

### remove zero counts
kidx = np.where(gene_exp.iloc[:, 0] > 0)[0]
gene_exp = gene_exp.iloc[kidx]

### merge
gene_exp = gene_exp.join(gene_lens)
del gene_lens

### parse simulated transcript file
tx_len = 0
header = None
tx_names = []
ge_names = []
tx_lens = []
for line in open(fname_fa_sim, 'r'):
    if line[0] == '>':
        if tx_len > 0:
            tx_names.append(header.split('|')[0])
            ge_names.append(header.split('|')[1])
            tx_lens.append(tx_len)
        tx_len = 0
        header = line.strip()[1:]
        continue
    tx_len += len(line.strip())
simulated = pd.DataFrame(index=tx_names, data={'gene_id':ge_names, 'tx_length':tx_lens})

### normalize to number of fragments
gene_exp['norm_exp'] = gene_exp.iloc[:, 0] / gene_exp['length'] * 100
gene_exp['log_norm_exp'] = np.log10(gene_exp.iloc[:, 0] / gene_exp['length'] * 100)

### estimate empirical distribution
empirical_pdf = spst.gaussian_kde(gene_exp['log_norm_exp'])

### generate sample in the size of number of simulated genes
npr.seed(random_seed)
num_genes = simulated['gene_id'].unique().shape[0]
emp_sample = 10**empirical_pdf.resample(num_genes)[0]

### plot sampled expression distribution
fig = plt.figure(figsize=(15, 10))
ax = fig.add_subplot(111)
_, bins, _ = ax.hist(gene_exp['log_norm_exp'], label='ICGC expression', bins=200, histtype='step', linewidth=1.5, density=True)
ax.hist(emp_sample, label='Simulation', bins=bins, histtype='step', linewidth=1.5, density=True)
ax.plot(bins, empirical_pdf(bins), label='KDE estimate')
plt.legend()
plt.savefig('expression_dist.pdf', fmt='pdf', bbox_inches='tight')

### we take the expression weights in multi-transcript genes from Figure 2 in https://www.nature.com/articles/s41598-020-73081-5
### we express at most 6 transcripts per gene here
weights = {1:[1.00],
           2:[0.83, 0.17],
           3:[0.75, 0.20, 0.05],
           4:[0.70, 0.20, 0.08, 0.02],
           5:[0.65, 0.20, 0.08, 0.05, 0.02],
           6:[0.62, 0.20, 0.07, 0.05, 0.04, 0.02],
          }

### distribute empirical expression across transcripts for each of the genes
gene_exp_sim = np.zeros((simulated.shape[0],), dtype='int')
tx_exp = np.zeros((simulated.shape[0],), dtype='int')
for i, gid in enumerate(simulated['gene_id'].unique()):
    if i > 0 and i % 100 == 0:
        sys.stderr.write('.')
        if i % 1000 == 0:
            sys.stderr.write('%i\n' % i)
        sys.stderr.flush()
    gidx = np.where(simulated['gene_id'] == gid)[0]
    gene_exp_sim[gidx] = int(emp_sample[i] / 2) + 1 # we divide by two, as we are simulating each haplotype independently
    num_tx = min(6, gidx.shape[0])
    for k, j in enumerate(npr.permutation(num_tx)):
        tx_exp[gidx[j]] = np.floor(weights[num_tx][k] * (int(emp_sample[i]) + 1)).astype('int')
simulated['gene_exp'] = gene_exp_sim
simulated['tx_exp'] = tx_exp
simulated['reads'] = np.ceil(simulated['tx_length'] * simulated['tx_exp'] / 100).astype('int')

def filter_var_pos(var_set, simulated):
    keep = {}
    for k in var_set:
        txs = []
        for tx in var_set[k]:
            txid = tx.split('|')[0]
            ### only keep variants that affect transcripts we have in our list and that are expressed
            if txid in simulated.index and simulated.loc[txid]['tx_exp'] > 0:
                txs.append(tx)
        if len(txs) > 0:
            keep[k] = txs
    return keep

### loading variants
hets_germline = pickle.load(open(fname_germline_pickle, 'rb')) # dictionary with key=(chrm, pos) and value=[tx_id1, tx_id2, tx_id3, ...]
hets_germline = filter_var_pos(hets_germline, simulated)
hets_somatic = pickle.load(open(fname_somatic_pickle, 'rb')) # dictionary with key=(chrm, pos) and value=[tx_id1, tx_id2, tx_id3, ...]
hets_somatic = filter_var_pos(hets_somatic, simulated)

### compute factors for transcript expression in each sample
### normal sample factors: germline ASE
### tumor sample factors: germline ASE, somatic ASE, DGE (driven by and independent of somatic variant), DTE
### we will only modulate the factors in *ht2* and achieve increase/decrease by factors >/< than 1.0
simulated['factor_ase_ht1_g'] = np.ones((simulated.shape[0],), dtype='int') 
simulated['factor_ase_ht2_g'] = np.ones((simulated.shape[0],), dtype='int') 
simulated['factor_ase_ht1_s'] = np.ones((simulated.shape[0],), dtype='int') 
simulated['factor_ase_ht2_s'] = np.ones((simulated.shape[0],), dtype='int') 
simulated['factor_dge_ht1'] = np.ones((simulated.shape[0],), dtype='int') 
simulated['factor_dge_ht2'] = np.ones((simulated.shape[0],), dtype='int') 

def select_and_remove(var_set, N):

    keys = [_ for _ in var_set]
    ### select
    idx = npr.choice(range(len(keys)), size=N, replace=False)
    mut_set = [var_set[keys[i]] for i in idx]
    ### remove
    for i in idx:
        del var_set[keys[i]]
    ### return
    return mut_set

### germline ASE -- affecting both tumor and normal samples
### strategy:
###  - select a subset of N germline variant positions
###  - for each position retrieve the affected transcripts
###  - for each transcript select an effect size of the variant in form of factor that is multiplied with tx_expression
for p in select_and_remove(hets_germline, NUM_GERM_ASE):
    txs = [_.split('|')[0] for _ in p]
    ### choose factor
    factor = npr.randint(2, 16)
    down = npr.random() < 0.5
    for tx in txs:
        if down:
            simulated.at[tx, 'factor_ase_ht1_g'] = factor * simulated.loc[tx, 'factor_ase_ht1_g']
        else:
            simulated.at[tx, 'factor_ase_ht2_g'] = factor * simulated.loc[tx, 'factor_ase_ht2_g']

### somatic ASE -- affecting only tumor samples
### strategy:
###  - select a subset of N somatic variant positions
###  - for each position retrieve the affected transcripts
###  - for each transcript select an effect size of the variant in form of factor that is multiplied with tx_expression
for p in select_and_remove(hets_germline, NUM_SOM_ASE):
    txs = [_.split('|')[0] for _ in p]
    ### choose factor, we use a 10% chance that there is NMD and the factor becomes 0
    if npr.random() < 0.1:
        factor = 0
    else:
        factor = npr.randint(2, 16)
    down = npr.random() < 0.5
    for tx in txs:
        if down:
            simulated.at[tx, 'factor_ase_ht1_s'] = factor * simulated.loc[tx, 'factor_ase_ht1_s']
        else:
            simulated.at[tx, 'factor_ase_ht2_s'] = factor * simulated.loc[tx, 'factor_ase_ht2_s']

### tumor DGE -- affecting only tumor samples
### strategy:
###  - select a set of genes, irrespective of somatic variants
###  - for each transcript select an effect size of the variant in form of factor that is multiplied with expression of all transcripts
gene_ids_u = np.unique(simulated['gene_id'])
idx = npr.choice(range(gene_ids_u.shape[0]), size=NUM_DGE, replace=False)
for gid in gene_ids_u[sorted(idx)]:
    gidx = simulated.index[np.where(simulated['gene_id'] == gid)[0]]
    factor = 2**((npr.random() + 1) * 2)
    down = npr.random() < 0.5
    if down:
        simulated.at[gidx, 'factor_dge_ht1'] = factor * simulated.loc[gidx, 'factor_dge_ht1'] 
    else:
        simulated.at[gidx, 'factor_dge_ht2'] = factor * simulated.loc[gidx, 'factor_dge_ht2'] 

### writing outputs
### Note: We will only write factors and read counts for transcripts that have expression > 0
###       as otherwise Polyester will error out
out_idx = np.where(simulated['reads'] > 0)[0]

### writing base read counts to output
simulated[['reads' for _ in range(REPLICATES)]].iloc[out_idx].to_csv(fname_out_readcounts, sep=',', header=False, index=False)

### writing tumor and normal factors to output (one output file per haplotype)
simulated[['factor_ase_ht1_g' for _ in range(REPLICATES)]].iloc[out_idx].to_csv(fname_out_factors_normal_ht1, sep=',', header=False, index=False)
simulated[['factor_ase_ht2_g' for _ in range(REPLICATES)]].iloc[out_idx].to_csv(fname_out_factors_normal_ht2, sep=',', header=False, index=False)
simulated['factor_tumor_ht1'] = simulated['factor_ase_ht1_g'] * simulated['factor_ase_ht1_s'] * simulated['factor_dge_ht1']
simulated[['factor_tumor_ht1' for _ in range(REPLICATES)]].iloc[out_idx].to_csv(fname_out_factors_tumor_ht1, sep=',', header=False, index=False)
simulated['factor_tumor_ht2'] = simulated['factor_ase_ht2_g'] * simulated['factor_ase_ht2_s'] * simulated['factor_dge_ht2']
simulated[['factor_tumor_ht2' for _ in range(REPLICATES)]].iloc[out_idx].to_csv(fname_out_factors_tumor_ht2, sep=',', header=False, index=False)

### output all metadata
simulated.to_csv(fname_out_metadata_all, sep='\t')
