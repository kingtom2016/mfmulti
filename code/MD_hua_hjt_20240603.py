# Modified by Zhongyi Hua, njbxhzy@hotmail.com

# -------------------FROM HERE IS THE ORIGINAL DESCRIPTION---------------------------------------------

# Metagenomic diversity (MD) calculator 
# Measures size and (dis)similarity of protein clusters from 
# MMSeqs2 cluster and alignment outputs

# Note: this code assumes you have installed MMseqs2 

# Usage without random subsampling: python MD.py -i <predicted_protein-encoding_genes_file.faa> -r F 
# With random subsampling: python MD.py -i <predicted_protein-encoding_genes_file.faa> -r T -N xxx
# where xxx is the total number of observations to randomly subsample
# Where users may wish to define their own MMSeqs2 E value and/or sequence identity cut-offs: python MD.py -i <predicted_protein-encoding_genes_file.faa> -r F -e xxx -I yyy
# where xxx and yyy are the user degined inputs for these cut-offs

# Written by Damien Finn, damien.finn@thuenen.de



import os
import argparse
from collections import defaultdict
import datetime
from pathlib import Path
import random
import shutil
import subprocess as sp

import numpy as np
from itertools import product
from pathos.pools import ProcessPool as Pool


def run_mmseqs2(args):
    # Run MMseqs2 
    sp.run(['mmseqs', 'createdb', args.input_file, args.tmp / 'DB', '-v', '1'])
    _cmd_lst = [args.bin, 'cluster', args.tmp / 'DB', args.tmp / 'Clu', args.tmp, '--cluster-mode', '1', '-v', '1']
    if args.Eval:
        _cmd_lst.append('-e', args.Eval[0])
    if args.identity:
        _cmd_lst.append('--min-seq-id', args.identity[0])
    sp.run(_cmd_lst, check=True)
    sp.run([args.bin, 'createtsv', args.tmp / 'DB', args.tmp / 'DB', args.tmp / 'Clu', args.tmp / 'Clu.tsv', '-v', '1'])
    sp.run([args.bin, 'search', args.tmp / 'DB', args.tmp / 'DB', args.tmp / 'res', args.tmp, '--search-type','3', '-v', '1'])
    sp.run([args.bin, 'convertalis', args.tmp / 'DB', args.tmp / 'DB', args.tmp / 'res', args.tmp / 'res.m8', '-v', '1'])


def gen_cluster(input_clust: Path, resample: bool, seq_obs: int):
    # Generate a dictionary of clusters 
    print(' -- Making Protein dictionary -- ')
    clustlist = [_line for _line in input_clust.read_text().splitlines()]
    if resample:
        clustlist = random.sample(clustlist, seq_obs)
    clusters = defaultdict(list)
    for _line in clustlist:
        _col1, _col2 = _line.split('\t')
        clusters[_col1].append(_col2)
    return clusters


# Generate a dictionary of pairwise distances
def gen_sim(input_aln: Path):
    print(' -- Collating pair-wise dissimilarities -- ')
    pairsim = defaultdict(float)

    for _line in input_aln.read_text().splitlines():
        _col1, _col2, _col3, *_ = _line.split('\t')
        pairsim[tuple(sorted((_col1, _col2)))] = float(_col3) # sorted to accelrate get results
    return pairsim

# Now to bring them together
def summary_res(clusters: dict, pairsim: dict):
    print(' -- Bringing things together -- ')
    clustdist = defaultdict(list)

    O = 0
    for _key, _val in clusters.items():
        query_keys = [tuple(sorted(_)) for _ in product([_key], _val)]
        O += len(query_keys)
        clustdist[_key] = [(1- pairsim.get(_, 0)) for _ in query_keys]
    
    # Derive indices from the cluster dictionary
    print(' -- Deriving diversity indices -- ')

    # Total Contigs?
    print('Total protein-encoding genes: ', O)

    # Protein Richness
    P = len(clustdist)
    print('Protein Richness: ', P)

    CDlens = [ len(k) for k in clustdist.values()] 
    tmp = np.sum(CDlens)

    # Shannon Diversity
    H = [ (len(k)/tmp)*np.log(len(k)/tmp) for k in clustdist.values() ]
    print('Shannon Diversity: ', np.sum(H)*-1) 

    # Simpson Diversity
    S = [ np.square(len(k)/tmp) for k in clustdist.values() ]
    print('Simpson Evenness: ', 1 - np.sum(S))

    # Metagenomic diversity
    MD = [ (1 + (np.sum(k))/len(k)) for k in clustdist.values()]
    print('Log10 Protein Dissimilarity:', np.log10(np.sum(MD)))
    print('Metagenomic Diversity index:', (1/O)*np.sum(MD))


def parse_args():
    # Define input data
    parser = argparse.ArgumentParser(description= 'Derive a Metagenomic Diversity index from clustered amino acid sequences')
    parser.add_argument('-i', '--input-file', required = True, 
                        help = 'Amino acid sequences as .faa')
    parser.add_argument('-r', '--random-sampling', required = False, action='store_true',
                        help = 'To randomly subsample observations or not')
    parser.add_argument('-N', '--seq-obs', required = False, type = int, 
                        help = 'If randomly subsampling, N denotes the number of observations to subsample without replacement')
    parser.add_argument('-e', '--Eval', required = False, nargs = 1, type = str, 
                        help = 'Redefine the MMSeqs2 minimum Evalue cutoff score')
    parser.add_argument('-I', '--identity', required = False, nargs = 1, type = str, 
                        help = 'Redefine the MMSeqs2 minimum sequence identity cutoff score')
    parser.add_argument('--bin', required=False, type=Path, default='mmseqs',
                        help='MMSeqs2 bin path. Default: mmseqs')
    parser.add_argument('--threads', required=False, default=os.cpu_count(), type=int,
                        help='Parallel cores. Default: all cores detected')
    parser.add_argument('--tmp', required=False, type=Path,
                        help='The temporarty directory. Default: The tmp directory under the working directory')
    args = parser.parse_args()
    return args


def main():
    """main interface
    """
    args = parse_args()
    path_cur_dir = Path(__file__).absolute().parent
    # Check and create temporary directory
    if not args.tmp:
        args.tmp = path_cur_dir / 'tmp'
    if not args.tmp.exists():
        args.tmp.mkdir()
    start_time = datetime.datetime.now()
    print(' -- Running MMseqs2 -- ')
    run_mmseqs2(args)
    #start_time = datetime.datetime.now()
    if args.random_sampling:
        sample_N = int(args.seq_obs)
    else:
        sample_N = -1 # just occupy a seat
    clusters = gen_cluster(args.tmp / 'Clu.tsv',
                           args.random_sampling,
                           sample_N)
    pairsim = gen_sim(input_aln = args.tmp / 'res.m8')
    summary_res(clusters, pairsim)
    excute_time = datetime.datetime.now() - start_time
    print(f'Time: {excute_time} Seconds')
    shutil.rmtree(args.tmp)


if __name__ == '__main__':
    main()