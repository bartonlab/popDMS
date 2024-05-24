import sys
import os
import re
import copy

import numpy as np
import scipy as sp
import scipy.stats as st

import itertools
from itertools import combinations

import pandas as pd
pd.set_option('future.no_silent_downcasting', True)

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from matplotlib.font_manager import FontProperties
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker

from colorsys import hls_to_rgb

import mplot as mp


# GLOBAL VARIABLES

## MaveDB column identifiers
MAVEDB_NT     = 'hgvs_nt'
MAVEDB_PRO    = 'hgvs_pro'
MAVEDB_SPLICE = 'hgvs_splice'
MAVEDB_WT     = '_wt'
MAVEDB_ACC    = 'accession'

## popDMS column identifiers
COL_SITE = 'site'
COL_AA   = 'amino_acid'
COL_WT   = 'WT_indicator'
COL_S    = 'joint'

## Sequence conventions
NUC    = ['A', 'C', 'G', 'T']                                      # Ordered list of nucleotides
CODONS = ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT',  # Ordered list of codons
          'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
          'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT',
          'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
          'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
          'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
          'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT',
          'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
AA     = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L',   # Ordered list of amino acids
          'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']
CODON2AA = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',            # Map from codons to amino acids
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W' }
AA2CODON = {'I': ['ATA', 'ATC', 'ATT'],                            # Map from amino acids to codons
            'M': ['ATG'], 
            'T': ['ACA', 'ACC', 'ACG', 'ACT'], 
            'N': ['AAC', 'AAT'], 
            'K': ['AAA', 'AAG'], 
            'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'], 
            'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'], 
            'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], 
            'P': ['CCA', 'CCC', 'CCG', 'CCT'], 
            'H': ['CAC', 'CAT'], 
            'Q': ['CAA', 'CAG'], 
            'V': ['GTA', 'GTC', 'GTG', 'GTT'], 
            'A': ['GCA', 'GCC', 'GCG', 'GCT'], 
            'D': ['GAC', 'GAT'], 
            'E': ['GAA', 'GAG'], 
            'G': ['GGA', 'GGC', 'GGG', 'GGT'], 
            'F': ['TTC', 'TTT'], 
            'Y': ['TAC', 'TAT'], 
            'C': ['TGC', 'TGT'], 
            'W': ['TGG'], 
            '*': ['TAA', 'TAG', 'TGA']}
CODON2AANUM = dict(zip(CODONS, [AA.index(CODON2AA[c]) for c in CODONS]))

MU  = { 'GC': 1.0e-7,
        'AT': 7.0e-7,
        'CG': 5.0e-7,
        'AC': 9.0e-7,
        'GT': 2.0e-6,
        'TA': 3.0e-6,
        'TG': 3.0e-6,
        'CA': 5.0e-6,
        'AG': 6.0e-6,
        'TC': 1.0e-5,
        'CT': 1.2e-5,
        'GA': 1.6e-5  }

# Figure variables

cm2inch = lambda x: x/2.54
SINGLE_COLUMN   = cm2inch(8.8)
ONE_FIVE_COLUMN = cm2inch(11.4)
DOUBLE_COLUMN   = cm2inch(18.0)
SLIDE_WIDTH     = 10.5
GOLDR           = (1.0 + np.sqrt(5)) / 2.0

FONTFAMILY    = 'Arial'
SIZESUBLABEL  = 8
SIZELABEL     = 6
SIZETICK      = 6
SMALLSIZEDOT  = 6.
SIZELINE      = 0.6
AXES_FONTSIZE = 6
AXWIDTH       = 0.4

BKCOLOR = '#252525'

FIGPROPS = {
    'transparent' : True,
    'bbox_inches' : 'tight'
}


# Built-in plot

def fig_dms(pop_file, fig_file):
    ''' Show performance comparison across data sets '''
    
    # Read in data
    df_pop = pd.read_csv(pop_file)
     
    # Plot variables
    w = DOUBLE_COLUMN

    box_t = 0.95
    box_b = 0.05
    box_l = 0.05
    box_r = 0.95

    sites_per_line = 100
    n_rows = 21
    h_per_line = w * (box_r - box_l) * (n_rows / sites_per_line)
    dh_per_line = 0.1 * h_per_line

    n_sites = len(np.unique(df_pop[COL_SITE]))
    n_lines = int(np.ceil(n_sites / sites_per_line))

    h = ((n_lines * h_per_line) + ((n_lines - 1) * dh_per_line)) / (box_t - box_b)

    fig = plt.figure(figsize = (w, h))
    
    box_y  = (box_t - box_b) / (n_lines + 0.1 * (n_lines - 1))
    box_dy = 0.1 * box_y

    cur_sites = n_sites
    box = []
    for i in range(n_lines):
        if cur_sites > sites_per_line:
            box.append(dict(left=box_l, right=box_r, bottom=box_t - (i + 1) * box_y - i * box_dy, top=box_t - i * box_y - i * box_dy))
        else:
            box.append(dict(left=box_l, right=box_l + (cur_sites / sites_per_line) * (box_r - box_l), bottom=box_t - (i + 1) * box_y - i * box_dy, top=box_t - i * box_y - i * box_dy))
        cur_sites -= sites_per_line
    
    gs = [gridspec.GridSpec(1, 1, **box[i]) for i in range(n_lines)]
    ax = [plt.subplot(gs[i][0, 0]) for i in range(n_lines)]
    
    # Plot selection heatmaps

    sites  = np.unique(df_pop[COL_SITE])
    df_WT  = df_pop[df_pop[COL_WT]==True]
    WT     = [df_WT[df_WT[COL_SITE]==s].iloc[0][COL_AA] for s in sites]
    s_WT   = [df_WT[df_WT[COL_SITE]==s].iloc[0][COL_S]  for s in sites]
    s_vec  = [[df_pop[(df_pop[COL_SITE]==sites[i]) & (df_pop[COL_AA]==aa)].iloc[0][COL_S] - s_WT[i] for aa in AA] for i in range(len(sites))]
    s_norm = np.max(np.fabs(s_vec))
    
    for i in range(n_lines):
        sub_sites = sites[i * sites_per_line:(i + 1) * sites_per_line]
        df_pop_sub = df_pop[df_pop[COL_SITE].isin(sub_sites)]
        plot_selection(ax[i], df_pop_sub, s_norm=s_norm, legend=(i==n_lines-1))
    
    # Save figure
    
    plt.show()
    fig.savefig(fig_file, **FIGPROPS)


def plot_selection(ax, df_pop, s_norm=1, legend=False):
    """ Plot selection heatmap. """

    # process stored data

    sites = np.unique(df_pop[COL_SITE])
    df_WT = df_pop[df_pop[COL_WT]==True]
    WT    = [df_WT[df_WT[COL_SITE]==s].iloc[0][COL_AA] for s in sites]
    s_WT  = [df_WT[df_WT[COL_SITE]==s].iloc[0][COL_S]  for s in sites]
    s_vec = [[df_pop[(df_pop[COL_SITE]==sites[i]) & (df_pop[COL_AA]==aa)].iloc[0][COL_S] - s_WT[i] for aa in AA] for i in range(len(sites))]

    # plot selection across the protein, normalizing WT residues to zero

    site_rec_props = dict(height=1, width=1, ec=None, lw=AXWIDTH/2, clip_on=False)
    prot_rec_props = dict(height=len(AA), width=len(sites), ec=BKCOLOR, fc='none', lw=AXWIDTH/2, clip_on=False)
    cBG             = '#F5F5F5'
    rec_patches     = []
    WT_dots_x       = []
    WT_dots_y       = []
    
    for i in range(len(sites)):
        WT_dots_x.append(i + 0.5)
        WT_dots_y.append(len(AA)-AA.index(WT[i])-0.5)
        for j in range(len(AA)):
            
            # skip WT
            if AA[j]==WT[i]:
                continue

            # fill BG for unobserved
            c = cBG
            if s_vec[i][j]!=-s_WT[i]:
                t = s_vec[i][j] / s_norm
                if np.fabs(t)>1:
                    t /= np.fabs(t)
                if t>0:
                    c = hls_to_rgb(0.02, 0.53 * t + 1. * (1 - t), 0.83)
                else:
                    c = hls_to_rgb(0.58, 0.53 * np.fabs(t) + 1. * (1 - np.fabs(t)), 0.60)
            
            rec = matplotlib.patches.Rectangle(xy=(i, len(AA)-1-j), fc=c, **site_rec_props)
            rec_patches.append(rec)
    
    rec = matplotlib.patches.Rectangle(xy=(0, 0), **prot_rec_props)
    rec_patches.append(rec)
    
    # add patches and plot
    
    for patch in rec_patches:
        ax.add_artist(patch)

    pprops = { 'colors':    [BKCOLOR],
               'xlim':      [0, len(sites) + 1],
               'ylim':      [0, len(AA) + 0.5],
               'xticks':    [],
               'yticks':    [],
               'plotprops': dict(lw=0, s=0.2*SMALLSIZEDOT, marker='o', clip_on=False),
               #'xlabel':    'Sites',
               #'ylabel':    'Amino acids',
               'theme':     'open',
               'axoffset':  0,
               'hide' :     ['top', 'bottom', 'left', 'right'] }

    mp.plot(type='scatter', ax=ax, x=[WT_dots_x], y=[WT_dots_y], **pprops)
    
    # legend
    
    rec_patches = []

    if legend:
    
        invt = ax.transData.inverted()
        xy1  = invt.transform((0,0))
        xy2  = invt.transform((0,9))
        legend_dy = (xy1[1]-xy2[1])/3 # multiply by 3 for slides/poster

        xloc = len(sites) + 2 + 6  # len(sites) - 6.5
        yloc = 17  # -4
        for i in range(-5, 5+1, 1):
            c = cBG
            t = i/5
            if t>0:
                c = hls_to_rgb(0.02, 0.53 * t + 1. * (1 - t), 0.83)
            else:
                c = hls_to_rgb(0.58, 0.53 * np.fabs(t) + 1. * (1 - np.fabs(t)), 0.60)
            rec = matplotlib.patches.Rectangle(xy=(xloc + i, yloc), fc=c, **site_rec_props)
            rec_patches.append(rec)
            
        txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
        ax.text(xloc+0.5, yloc-4*legend_dy, 'Mutation effects', clip_on=False, **txtprops)
        ax.text(xloc-4.5, yloc+2*legend_dy, '-', clip_on=False, **txtprops)
        ax.text(xloc+5.5, yloc+2*legend_dy, '+', clip_on=False, **txtprops)

        xloc = len(sites) + 2  # 0
        yloc = 11  # -4
        c    = cBG
        rec  = matplotlib.patches.Rectangle(xy=(xloc, yloc + 4*legend_dy), fc=c, **site_rec_props) # paper
    #    rec  = matplotlib.patches.Rectangle(xy=(xloc, yloc + 6*legend_dy), fc=c, **site_rec_props) # slides
        rec_patches.append(rec)

        for patch in rec_patches:
            ax.add_artist(patch)
            
        mp.scatter(ax=ax, x=[[xloc + 0.5]], y=[[yloc + 0.5]], **pprops)

        txtprops['ha'] = 'left'
        ax.text(xloc + 1.5, yloc + 0.5, 'WT amino acid', clip_on=False, **txtprops)
        ax.text(xloc + 1.5, yloc + 0.5 + 4*legend_dy, 'Not observed', clip_on=False, **txtprops) # paper
    #    ax.text(xloc + 1.3, yloc + 0.5 + 6*legend_dy, 'Not observed', clip_on=False, **txtprops) # slides


# FUNCTIONS

def get_variant_sites_nucs(mavedb_nt, shift_by_one=True):
    '''
    Return a list of variants and their locations from an 'hvgs_nt' column in MaveDB
    format
    '''
    ################################################################################
    
    # MaveDB variant format c.[XA>B;YC>D]
    # X, Y: nucleotide numbers
    # A, C: reference nucleotides
    # B, D: mutant nucleotides
    # ;: separator between substitutions
    
    if mavedb_nt==MAVEDB_WT:
        return [], []
    
    var_str = mavedb_nt.replace('[', '').replace(']', '').split('.')[1].split(';')
    
    sites = []
    nucs  = []
    for var in var_str:
        sites.append(int(''.join([x for x in var.split('>')[0] if not x.isalpha()])))
        nucs.append(var.split('>')[-1])
    
    ord = np.argsort(sites)
    sites = [sites[i] for i in ord]
    nucs = [nucs[i] for i in ord]
    
    if shift_by_one:
        sites = np.array(sites) - 1
    
    return sites, nucs
    

def get_selection_file(dir, name, file_ext='.csv.gz'):
    ''' Return the path to a file recording selection coefficients '''
    return os.path.join(dir, '_'.join([name, 'selection_coefficients']) + file_ext)


def get_reads_file(dir, name, replicate, file_ext='.csv'):
    ''' Return the path to a file recording the number of reads for each time point '''
    return os.path.join(dir, '_'.join([name, 'reads_rep', str(replicate)]) + file_ext)


def get_codon_count_file(dir, name, type, replicate, file_ext='.csv.gz'):
    ''' Return the path to a codon count file (temporary, to be converted to frequencies) '''
    return os.path.join(dir, '_'.join([name, type, 'codon_counts_rep', str(replicate)]) + file_ext)


def get_codon_freq_file(dir, name, type, replicate, file_ext='.csv.gz'):
    ''' Return the path to a codon frequency file '''
    return os.path.join(dir, '_'.join([name, type, 'codon_rep', str(replicate)]) + file_ext)


def get_aa_freq_file(dir, name, type, replicate, file_ext='.csv.gz'):
    ''' Return the path to an amino acid frequency file '''
    return os.path.join(dir, '_'.join([name, type, 'aa_rep', str(replicate)]) + file_ext)
    
    
def read_fancy_comments(file, comment_char=None, sep=','):
    '''
    In count files stored in MaveDB format, '#' characters are sometimes used as
    both comment lines (at the beginning of the file) and as part of accession
    numbers. read_csv() in pandas removes all instances of the comment character.
    Instead, this function will remove only lines that start with the comment
    character. h/t Thymen on stack overflow.
    '''
    
    if not isinstance(comment_char, str):
        return pd.read_csv(file)
    
    headers = None
    results = []
    with open(file, 'r') as f:
        first_line = True
        for line in f.readlines():
            if line.startswith(comment_char):
                continue

            if line.strip():
                if first_line:
                    headers = [word for word in line.strip().split(sep) if word]
                    headers = list(map(str.strip, headers))
                    first_line = False
                else:
                    results.append([word for word in line.strip().split(sep) if word])

    return pd.DataFrame(results, columns=headers)


def convert_codon_counts_to_variant_frequencies(dir, name, types, n_replicates, reference_sequence_file, comment_char=None):
    '''
    Convert codon counts (popDMS format) to variant frequencies
    '''
    
    # Read in codon counts
    replicates = [i+1 for i in range(n_replicates)]
    codon_counts = {}
    for type in types:
        codon_list = []
        for rep in replicates:
            codon_list.append(pd.read_csv(get_codon_count_file(dir, name, type, rep)))
        codon_counts[type] = codon_list
    
    # Read in reference sequence
    ref_seq = ''.join(open(reference_sequence_file).read().split())
    
    # Extract numbered sites from data and construct the reference sequence
    codon_length = 3
    site_list = sorted(codon_counts['single'][0]['site'].astype('str').str.extractall('(\d+)')[0].astype(int).unique())
    ref_codons = np.array([''.join(ref_seq[i:i+codon_length]) for i in range(0, len(ref_seq), codon_length)])
    ref_codons_by_site = dict(zip(site_list, ref_codons))
    ref_aas = np.array([CODON2AA[c] for c in ref_codons])
    ref_aas_by_site = dict(zip(site_list, ref_aas))

    # Compute total reads and save to file
    reads_cols = ['generation', 'reads']
    time_points = [sorted(list(codon_counts['single'][r_idx]['generation'].unique())) for r_idx in range(len(replicates))]
    total_count = []
    for r_idx in range(len(replicates)):
        total_count.append([])
        for t in time_points[r_idx]:
            total_count[r_idx].append(np.sum(codon_counts['single'][r_idx][(codon_counts['single'][r_idx]['generation']==t) & (codon_counts['single'][r_idx]['site']==site_list[0])]['counts']))
        df_temp = pd.DataFrame(data=[[time_points[r_idx][i], total_count[r_idx][i]] for i in range(len(time_points[r_idx]))], columns=reads_cols)
        df_temp.to_csv(get_reads_file(dir, name, replicates[r_idx], file_ext='.csv'), index=False)

    # Map dataframe row data to dictionaries
    def get_key(type, row):
        if type=='single': return (row['site'], CODON2AA[row['codon']])
        elif type=='double': return (row['site_1'], CODON2AA[row['codon_1']], row['site_2'], CODON2AA[row['codon_2']])
    
    # Initialize dictionary of codon/aa counts and fill in counts
    aa_counts = dict(zip(types, [[[{} for t in time_points[r_idx]] for r_idx in range(len(replicates))] for type in types]))
    for r_idx in range(len(replicates)):
        for i in range(len(time_points[r_idx])):
            for type in types:
                for df_iter, row in codon_counts[type][r_idx][codon_counts[type][r_idx]['generation']==time_points[r_idx][i]].iterrows():
                    key = get_key(type, row)
                    if key not in aa_counts[type][r_idx][i]:
                        aa_counts[type][r_idx][i][key] = row['counts']
                    else:
                        aa_counts[type][r_idx][i][key] += row['counts']
    
    # Save single aa frequencies to file, including WT indicator
    single_aa_columns = ['generation', 'site', 'aa', 'frequency', 'WT_indicator']
    for r_idx in range(len(replicates)):
        counts_list = []
        for i in range(len(time_points[r_idx])):
            t = time_points[r_idx][i]
            for aas, counts in aa_counts['single'][r_idx][i].items():
                if counts>0:
                    counts_list.append([t] + list(aas) + [counts/total_count[r_idx][i], aas[1]==ref_aas_by_site[aas[0]]])
        df_temp = pd.DataFrame(data=counts_list, columns=single_aa_columns)
        df_temp.to_csv(get_aa_freq_file(dir, name, 'single', replicates[r_idx]), index=False, compression='gzip')

    # Save double aa frequencies to file
    if 'double' in types:
        double_aa_columns = ['generation', 'site_1', 'aa_1', 'site_2', 'aa_2', 'frequency']
        for r_idx in range(len(replicates)):
            counts_list = []
            for i in range(len(time_points[r_idx])):
                t = time_points[r_idx][i]
                for aas, counts in aa_counts['double'][r_idx][i].items():
                    if counts>0:
                        counts_list.append([t] + list(aas) + [counts/total_count[r_idx][i]])

        df_temp = pd.DataFrame(data=counts_list, columns=double_aa_columns)
        df_temp.to_csv(get_aa_freq_file(dir, name, 'double', replicates[r_idx]), index=False, compression='gzip')

    # Convert codon counts to frequencies and save to file
    df_cols = dict(single=['generation', 'site', 'codon', 'frequency'],
                   double=['generation', 'site_1', 'codon_1', 'site_2', 'codon_2', 'frequency'])
    id_cols = dict(single=['site', 'codon'], 
                   double=['site_1', 'codon_1', 'site_2', 'codon_2'])
    for type in types:
        for r_idx in range(len(replicates)):
            counts_list = []
            for i in range(len(time_points[r_idx])):
                t = time_points[r_idx][i]
                for df_iter, row in codon_counts[type][r_idx][codon_counts[type][r_idx]['generation']==t].iterrows():
                    counts_list.append([t] + list(row[id_cols[type]]) + [row['counts']/total_count[r_idx][i]])
            df_temp = pd.DataFrame(data=counts_list, columns=df_cols[type])
            df_temp.to_csv(get_codon_freq_file(dir, name, type, replicates[r_idx]), index=False, compression='gzip')


def compute_freq_MaveDB(name, rep, types, haplotype_counts_file, reference_sequence, time_points, time_cols, freq_dir, comment_char=None):
    '''
    Compute variant frequencies (single and double codons) from an input haplotype
    counts file in MaveDB format
    '''
    ################################################################################
    
    # Read in haplotype counts, fill NA counts, and drop unneeded columns
    df_data = read_fancy_comments(haplotype_counts_file, comment_char=comment_char).replace('NA', np.nan).fillna(0).drop([MAVEDB_ACC], axis=1)
    if MAVEDB_SPLICE in df_data.columns:
        df_data = df_data.drop([MAVEDB_SPLICE], axis=1)
    df_data.reset_index(drop=True, inplace=True)
    
    # Remove rows with ambiguous nucleotides
    df_data = df_data[~df_data[MAVEDB_NT].astype('str').str.contains('X', regex=False)]
    
    # Extract numbered sites from data and construct the reference sequence
    codon_length = 3
    site_list = sorted(df_data[MAVEDB_PRO].astype('str').str.extractall('(\d+)')[0].astype(int).unique())
    ref_codons = np.array([''.join(reference_sequence[i:i+codon_length]) for i in range(0, len(reference_sequence), codon_length)])
    ref_codons_by_site = dict(zip(site_list, ref_codons))
    ref_aas = np.array([CODON2AA[c] for c in ref_codons])
    ref_aas_by_site = dict(zip(site_list, ref_aas))
    
    # Explicitly cast numerical columns
    col_types = dict(zip(time_cols, ['float64' for c in time_cols]))
    df_data = df_data.astype(col_types)
    
    # Compute total counts
    total_count = []
    time_sort = np.argsort(time_points)
    time_points = list(np.array(time_points)[time_sort])
    time_cols = np.array(time_cols)[time_sort]
    for i in range(len(time_points)):
        total_count.append(np.sum(df_data[time_cols[i]]))
    
    # Initialize dictionary of codon/aa counts and fill in WT counts
    codon_counts = dict(zip(['single', 'double'], [[{} for t in time_points] for i in range(2)]))
    aa_counts = dict(zip(types, [[{} for t in time_points] for type in types]))
    for i in range(len(time_points)):
        for seq_i in range(len(site_list)):
            idx_i = site_list[seq_i]
            
            for c in CODONS:
                codon_counts['single'][i][(idx_i, c)] = 0
            codon_counts['single'][i][(idx_i, ref_codons[seq_i])] = total_count[i]
            
            for a in AA:
                aa_counts['single'][i][(idx_i, a)] = 0
            aa_counts['single'][i][(idx_i, CODON2AA[ref_codons[seq_i]])] = total_count[i]
                
            for seq_j in range(seq_i+1, len(site_list)):
                idx_j = site_list[seq_j]
                codon_counts['double'][i][(idx_i, ref_codons[seq_i], idx_j, ref_codons[seq_j])] = total_count[i]
                aa_counts['double'][i][(idx_i, CODON2AA[ref_codons[seq_i]], idx_j, CODON2AA[ref_codons[seq_j]])] = total_count[i]
    
    # Iterate through haplotypes and compute variant counts
    reference_sequence = list(reference_sequence)
    for df_iter, row in df_data.iterrows():
    
        # Progress bar
        if int(df_iter)%1000==0:
            print('{} replicate {}, computing codon frequencies... {:2.1%}'.format(name, rep, int(df_iter) / len(df_data)), end="\r", flush=True)
        
        # Extract list of mutations (skip WT)
        variants = row[MAVEDB_NT]
        if variants==MAVEDB_WT:
            continue
        variant_sites, variant_nucs = get_variant_sites_nucs(variants, shift_by_one=True)
        variant_sequence = reference_sequence.copy()
        for v_site, v_nuc in zip(variant_sites, variant_nucs):
            variant_sequence[v_site] = v_nuc
        variant_codons = np.array([''.join(variant_sequence[i:i+codon_length]) for i in range(0, len(variant_sequence), codon_length)])
        variant_aas = np.array([CODON2AA[c] for c in variant_codons])
        
        # Get mismatched codons/aas and sparsely fill counts
        for var_state, ref_state, single_state_counts, double_state_counts, states in [
            [variant_codons, ref_codons, codon_counts['single'], codon_counts['double'], CODONS],
            [variant_aas, ref_aas, aa_counts['single'], aa_counts['double'], AA]]:
            
            mismatches = [i for i in range(len(var_state)) if var_state[i]!=ref_state[i]]
            
            for ii in range(len(mismatches)):
                seq_i = mismatches[ii]
                idx_i = site_list[seq_i]
                st_i  = var_state[seq_i]
                
                # Single codon/aa counts
                for i in range(len(time_points)):
                    single_state_counts[i][(idx_i, st_i)] += row[time_cols[i]]
                
                # Double codon/aa counts
                for jj in range(ii+1, len(mismatches)):
                    seq_j = mismatches[jj]
                    idx_j = site_list[seq_j]
                    st_j  = var_state[seq_j]
                    
                    for i in range(len(time_points)):
                        if (idx_i, st_i, idx_j, st_j) not in double_state_counts[i]:
                            double_state_counts[i][(idx_i, st_i, idx_j, st_j)] = row[time_cols[i]]
                        else:
                            double_state_counts[i][(idx_i, st_i, idx_j, st_j)] += row[time_cols[i]]
                        
    # Add counts for one mutant and one reference codon/aa
    for i in range(len(time_points)):
        for ref_state, single_state_counts, double_state_counts, states in [
            [ref_codons, codon_counts['single'], codon_counts['double'], CODONS],
            [ref_aas, aa_counts['single'], aa_counts['double'], AA]]:
            for seq_i, st_i in itertools.product(range(len(site_list)), states):
                idx_i = site_list[seq_i]
                if st_i!=ref_state[seq_i]:
                    total_single = single_state_counts[i][(idx_i, st_i)]
                    for seq_j in range(seq_i):
                        idx_j = site_list[seq_j]
                        total_double = 0
                        for st_j in states:
                            if st_j!=ref_state[seq_j] and (idx_j, st_j, idx_i, st_i) in double_state_counts[i]:
                                total_double += double_state_counts[i][(idx_j, st_j, idx_i, st_i)]
                        double_state_counts[i][(idx_j, ref_state[seq_j], idx_i, st_i)] = total_single - total_double
                    for seq_j in range(seq_i+1, len(site_list)):
                        idx_j = site_list[seq_j]
                        total_double = 0
                        for st_j in states:
                            if st_j!=ref_state[seq_j] and (idx_i, st_i, idx_j, st_j) in double_state_counts[i]:
                                total_double += double_state_counts[i][(idx_i, st_i, idx_j, st_j)]
                        double_state_counts[i][(idx_i, st_i, idx_j, ref_state[seq_j])] = total_single - total_double
                    
    # Correct reference/WT counts
    # Future: add function to get key for each data type (single/double, codon/aa) and iterate generically through all data structures
    for i in range(len(time_points)):
        for codons, counts in codon_counts['single'][i].items():
            if codons[1]!=ref_codons_by_site[codons[0]]:
                codon_counts['single'][i][(codons[0], ref_codons_by_site[codons[0]])] -= counts
                
        for codons, counts in codon_counts['double'][i].items():
            if codons[1]!=ref_codons_by_site[codons[0]] or codons[3]!=ref_codons_by_site[codons[2]]:
                codon_counts['double'][i][(codons[0], ref_codons_by_site[codons[0]], codons[2], ref_codons_by_site[codons[2]])] -= counts
                
        for aas, counts in aa_counts['single'][i].items():
            if aas[1]!=ref_aas_by_site[aas[0]]:
                aa_counts['single'][i][(aas[0], ref_aas_by_site[aas[0]])] -= counts
                
        for aas, counts in aa_counts['double'][i].items():
            if aas[1]!=ref_aas_by_site[aas[0]] or aas[3]!=ref_aas_by_site[aas[2]]:
                aa_counts['double'][i][(aas[0], ref_aas_by_site[aas[0]], aas[2], ref_aas_by_site[aas[2]])] -= counts

    # Save total reads to file
    reads_cols = ['generation', 'reads']
    df_temp = pd.DataFrame(data=[[time_points[i], total_count[i]] for i in range(len(time_points))], columns=reads_cols)
    df_temp.to_csv(get_reads_file(freq_dir, name, rep, file_ext='.csv'), index=False)
    
    # Save single aa frequencies to file, including WT indicator
    single_aa_columns = ['generation', 'site', 'aa', 'frequency', 'WT_indicator']
    counts_list = []
    for i in range(len(time_points)):
        t = time_points[i]
        for aas, counts in aa_counts['single'][i].items():
            if counts>0:
                counts_list.append([t] + list(aas) + [counts/total_count[i], aas[1]==ref_aas_by_site[aas[0]]])

    path = get_aa_freq_file(freq_dir, name, 'single', rep)
    df_temp = pd.DataFrame(data=counts_list, columns=single_aa_columns)
    df_temp.to_csv(path, index=False, compression='gzip')

    # Save other frequencies to file
    single_codon_columns = ['generation', 'site', 'codon', 'frequency']
    double_codon_columns = ['generation', 'site_1', 'codon_1', 'site_2', 'codon_2', 'frequency']
    double_aa_columns = ['generation', 'site_1', 'aa_1', 'site_2', 'aa_2', 'frequency']
    for counts_dict, columns, path in [
        [codon_counts['single'], single_codon_columns, get_codon_freq_file(freq_dir, name, 'single', rep)],
        [codon_counts['double'], double_codon_columns, get_codon_freq_file(freq_dir, name, 'double', rep)],
        [aa_counts['double'], double_aa_columns, get_aa_freq_file(freq_dir, name, 'double', rep)]]:
        
        counts_list = []
        for i in range(len(time_points)):
            t = time_points[i]
            for pairs, counts in counts_dict[i].items():
                if counts>0:
                    counts_list.append([t] + list(pairs) + [counts/total_count[i]])

        df_temp = pd.DataFrame(data=counts_list, columns=columns)
        df_temp.to_csv(path, index=False, compression='gzip')
    
    
def compute_variant_frequencies_Bloom(name, codon_counts_files, replicates, times, freq_dir='.', error_file=None, comment_char=None):
    '''
    Compute variant frequencies from codon counts in Bloom lab format. If an error
    file is provided, corrects for sequencing errors.

    Required arguments:
        - name: Name of the experiment
        - codon_counts_files: List of paths to codon counts files
        - replicates: List of replicate numbers
        - times: List of time points
        - freq_dir: Path to directory where frequencies will be saved

    Optional arguments:
        - error_file (default: None): Path to a file recording sequence counts from
            sequencing the reference sequence only, used to estimate sequencing
            error rates and correct variant frequencies
        - comment_char (default: None): Character used to indicate comments in the
            codon counts files
    '''
    ################################################################################
    
    # Collect basic information
    unique_reps = np.sort(np.unique(replicates))
    n_replicates = len(unique_reps)
    replicates = np.array(replicates)
    times = np.array(times)

    # Read in sequencing error frequencies, if provided
    error_freqs = None
    if error_file is not None:
        error_freqs = get_codon_error_freqs_Bloom(error_file, comment_char=comment_char)

    # Read in codon counts for each replicate and save to file
    codon_counts_files = np.array(codon_counts_files)
    for r_idx in range(n_replicates):
        
        # Read in codon counts and time order them
        n_files = len(codon_counts_files[replicates==unique_reps[r_idx]])
        codon_counts = [pd.read_csv(codon_counts_files[replicates==unique_reps[r_idx]][i]) for i in range(n_files)]
        codon_times = [times[replicates==unique_reps[r_idx]][i] for i in range(n_files)]
        time_sort = np.argsort(codon_times)
        codon_counts = [codon_counts[i] for i in time_sort]
        codon_times = [codon_times[i] for i in time_sort]

        # Get reference sequence
        sites = list(np.sort(np.unique(codon_counts[0]['site'])))
        ref_codons = np.array([str(codon_counts[0][codon_counts[0]['site']==s].iloc[0]['wildtype']) for s in sites])
        ref_aas = np.array([CODON2AA[c] for c in ref_codons])
        ref_aas_by_site = dict(zip(sites, ref_aas))

        # Compute amino acid frequencies and save to file
        max_reads = [0 for t in codon_times]
        aa_freqs = [{} for t in codon_times]
        for i in range(len(codon_times)):
            for df_iter, row in codon_counts[i].iterrows():
                total_counts = np.sum(row[CODONS])
                if total_counts>max_reads[i]:
                    max_reads[i] = total_counts

                if error_freqs is not None:
                    f_vec_measured = np.array(row[CODONS])/total_counts
                    err_matrix = np.eye(len(CODONS))
                    wt_idx = CODONS.index(ref_codons[sites.index(row['site'])])
                    for c_i in range(len(CODONS)):
                        if c_i!=wt_idx:
                            err_matrix[c_i, wt_idx] = error_freqs[row['site']][c_i]
                        else:
                            err_matrix[c_i, c_i] = 1 - np.sum(error_freqs[row['site']])
                    f_vec_corrected = np.matmul(np.linalg.inv(err_matrix), f_vec_measured)

                    for aa in AA:
                        aa_idxs = np.array([CODONS.index(c) for c in AA2CODON[aa]])
                        aa_freqs[i][(row['site'], aa)] = np.sum(f_vec_corrected[aa_idxs])

                else:
                    for aa in AA:
                        aa_cols = AA2CODON[aa]
                        aa_counts = np.sum(row[aa_cols])
                        aa_freqs[i][(row['site'], aa)] = aa_counts/total_counts

        # Save reads to file
        reads_cols = ['generation', 'reads']
        df_temp = pd.DataFrame(data=[[codon_times[i], max_reads[i]] for i in range(len(codon_times))], columns=reads_cols)
        df_temp.to_csv(get_reads_file(freq_dir, name, unique_reps[r_idx], file_ext='.csv'), index=False)

        # Save amino acid frequencies to file, with WT indicator
        aa_freqs_list = []
        for i in range(len(codon_times)):
            t = codon_times[i]
            for aas, freq in aa_freqs[i].items():
                aa_freqs_list.append([t] + list(aas) + [freq, aas[1]==ref_aas_by_site[aas[0]]])

        aa_freqs_cols = ['generation', 'site', 'aa', 'frequency', 'WT_indicator']
        df_temp = pd.DataFrame(data=aa_freqs_list, columns=aa_freqs_cols)
        df_temp.to_csv(get_aa_freq_file(freq_dir, name, 'single', unique_reps[r_idx]), index=False, compression='gzip')


def get_codon_error_freqs_Bloom(path, comment_char=None):
    ''' Get sequencing error frequencies from data in Bloom lab format. '''

    # Read in dataframe and get list of sites
    df = pd.read_csv(path, comment=comment_char)
    sites = np.sort(np.unique(df['site']))

    # Get wildtype codons
    wt_codons = np.array([str(df[df['site']==s].iloc[0]['wildtype']) for s in sites])

    # Get total number of reads for each site
    total_reads = np.zeros(len(sites))
    for i in range(len(sites)):
        total_reads[i] = np.sum(df[df['site']==sites[i]].iloc[0][CODONS])

    # Get error frequencies
    error_freqs = {}
    for s in sites:
        error_freqs[s] = np.zeros(len(CODONS))
    for seq_i in range(len(sites)):
        df_site = df[df['site']==sites[seq_i]].iloc[0]
        for c in CODONS:
            error_freqs[sites[seq_i]][CODONS.index(c)] = df_site[c]/total_reads[seq_i]
        
        # Set error for WT codon to zero
        error_freqs[sites[seq_i]][CODONS.index(wt_codons[seq_i])] = 0

    return error_freqs


def compute_variant_frequencies(name, reference_sequence_file, haplotype_counts_file, n_replicates, time_points, time_cols, freq_dir='.', with_epistasis=False, comment_char=None, format='MaveDB'):
    '''
    Computes variant (allele) frequencies, used by popDMS, from haplotype counts. By
    default, we assume that counts are saved in MaveDB format. For more information
    on this format, see https://www.mavedb.org/docs/mavedb/data_formats.html

    Required arguments:
        - name: String that will be associated with saved files, serving as an
            identifier for the data set
        - reference_sequence_file: File path to the reference sequence, which is
            needed to specify 'WT' variants
        - haplotype_counts_file: Path to an existing file recording haplotype
            frequencies
        - n_replicates: Number of experimental replicates
        - time_points: Time points in the data; for normalization across data sets,
            this should be roughly the number of 'generations' between sequencing
            samples in the experiment
        - time_cols: A dictionary mapping between the number for each experimental
            replicate and the column in the haplotype counts file that contains
            count information

    Optional arguments:
        - freq_dir (default: '.'): Path to directory where frequencies will be saved
        - with_epistasis (default: False): If True, will compute 3- and 4-amino acid
            frequencies for epistasis analysis
        - comment_char (default: None): Character used as a comment in the haplotype
            counts file; in some MaveDB files, this is '#'
        - format (default: 'MaveDB'): Format of the input data files (options:
            'MaveDB', 'Bloom')
    '''
    ################################################################################
    
    # Read in reference sequence
    ref_seq = ''.join(open(reference_sequence_file).read().split())
    
    # Compute and save counts for each experimental replicate
    replicates = [i+1 for i in range(n_replicates)]
    for rep in replicates:
        types = ['single', 'double']
        if with_epistasis:
            types += ['triple', 'quadruple']
        if format=='MaveDB':
            compute_freq_MaveDB(name, rep, types, haplotype_counts_file, ref_seq, time_points, time_cols[rep], freq_dir, comment_char)
        # elif format=='Bloom':
        #     compute_freq_Bloom(haplotype_counts_file, ref_seq, time_points, time_cols[rep], path_single, path_double, comment_char)
        print('\r' + ' ' * 80, end='\r')
        print(name + ' replicate '+str(rep)+' complete')
    

def compute_dx_covariance(aa_freqs):
    '''
    Compute net change in amino acid frequency and integrated covariance matrix
    from amino acid frequency data.
    '''
    
    # Gather information
    reps  = len(aa_freqs['single'])
    sites = list(np.unique(aa_freqs['single'][0]['site']))
    sites.sort()

    q = len(AA)
    L = len(sites)
    bigL = q*L
    
    p2i = {}
    for seq_i in range(len(sites)):
        for a in AA:
            p2i[(sites[seq_i], a)] = (q * seq_i) + AA.index(a)
    
    # Shape dx vector (reps x [q*L]) and covariance matrix (reps x [q*L, q*L]), compute for each replicate
    dx   = [np.zeros(bigL) for i in range(reps)]
    icov = [np.zeros((bigL, bigL)) for i in range(reps)]

    for r_idx in range(reps):
        # Get times
        times = list(np.unique(aa_freqs['single'][r_idx]['generation']))
        times.sort()

        dtsum = np.array([times[1]-times[0]] + [times[i+1]-times[i-1] for i in range(1, len(times)-1)] + [times[-1]-times[-2]])

        # Compute dense frequency vector to speed calculations
        x = np.array([np.zeros(bigL) for i in range(len(times))])
        for i in range(len(times)):
            t = times[i]
            df_t = aa_freqs['single'][r_idx][aa_freqs['single'][r_idx]['generation']==t]
            for df_iter, row in df_t.iterrows():
                x[i][p2i[(row['site'], row['aa'])]] += row['frequency']
            
        # Compute dx (final - initial frequency)
        dx[r_idx] = x[-1] - x[0]

        # Compute integrated covariance
        ## Off-diagonal terms, same time (note: diagonals will temporarily be incorrect)
        #icov[r_idx] = np.einsum('i,ijk->jk', dtsum, -np.array([np.outer(x[i], x[i])/3 for i in range(len(times))]), optimize=True)
        icov[r_idx] = np.tensordot(dtsum, -np.array([np.outer(x[i], x[i])/3 for i in range(len(times))]), axes=1)

        ## Off-diagonal terms, cross times
        for i in range(1, len(times)):
            dt = times[i] - times[i-1]
            icov[r_idx] -= dt * (np.outer(x[i-1], x[i]) + np.outer(x[i], x[i-1]))/6

        ## Off-diagonal terms, sparse pair correlations
        for i in range(len(times)):
            df_t = aa_freqs['double'][r_idx][aa_freqs['double'][r_idx]['generation']==times[i]]
            for df_iter, row in df_t.iterrows():
                icov[r_idx][p2i[(row['site_1'], row['aa_1'])], p2i[(row['site_2'], row['aa_2'])]] += dtsum[i] * row['frequency']/2
                icov[r_idx][p2i[(row['site_2'], row['aa_2'])], p2i[(row['site_1'], row['aa_1'])]] += dtsum[i] * row['frequency']/2

        ## Diagonal terms, same time (overwrite previous diagonals)
        icov[r_idx][np.diag_indices_from(icov[r_idx])] = dtsum.dot((x/2) - (x**2/3))

        ## Diagonal terms, cross times
        for i in range(1, len(times)):
            dt = times[i] - times[i-1]
            icov[r_idx][np.diag_indices_from(icov[r_idx])] -= dt * (x[i] * x[i-1])/3

    return dx, icov, p2i


def compute_dx_covariance_independent(aa_freqs):
    '''
    Compute net change in amino acid frequency and integrated covariance matrix
    from amino acid frequency data. This version assumes that each site is
    independent, so that the covariance matrix is block diagonal.
    '''
    
    # Gather information
    reps  = len(aa_freqs['single'])
    sites = list(np.unique(aa_freqs['single'][0]['site']))
    sites.sort()

    q = len(AA)
    L = len(sites)
    bigL = q*L
    
    aa2i = {}
    for a in AA:
        aa2i[a] = AA.index(a)
    
    # Shape dx vector (reps x [q*L]) and covariance matrix (reps x [q*L, q*L]), compute for each replicate
    dx   = [[np.zeros(q) for i in range(L)] for j in range(reps)]
    icov = [[np.zeros((q, q)) for i in range(L)] for j in range(reps)]

    for r_idx in range(reps):
        # Get times
        times = list(np.unique(aa_freqs['single'][r_idx]['generation']))
        times.sort()

        dtsum = np.array([times[1]-times[0]] + [times[i+1]-times[i-1] for i in range(1, len(times)-1)] + [times[-1]-times[-2]])

        # Iterate over individual sites
        for seq_i in range(L):
            df_site = aa_freqs['single'][r_idx][aa_freqs['single'][r_idx]['site']==sites[seq_i]]

            # Compute dense frequency vector to speed calculations
            x = np.array([np.zeros(q) for i in range(len(times))])
            for i in range(len(times)):
                t = times[i]
                df_t = df_site[df_site['generation']==t]
                for df_iter, row in df_t.iterrows():
                    x[i][aa2i[row['aa']]] += row['frequency']
                
            # Compute dx (final - initial frequency)
            dx[r_idx][seq_i] = x[-1] - x[0]

            # Compute integrated covariance
            ## Off-diagonal terms, same time (note: diagonals will temporarily be incorrect)
            #icov[r_idx][seq_i] = np.einsum('i,ijk->jk', dtsum, -np.array([np.outer(x[i], x[i])/3 for i in range(len(times))]), optimize=True)
            icov[r_idx][seq_i] = np.tensordot(dtsum, -np.array([np.outer(x[i], x[i])/3 for i in range(len(times))]), axes=1)

            ## Off-diagonal terms, cross times
            for i in range(1, len(times)):
                dt = times[i] - times[i-1]
                icov[r_idx][seq_i] -= dt * (np.outer(x[i-1], x[i]) + np.outer(x[i], x[i-1]))/6

            ## Diagonal terms, same time (overwrite previous diagonals)
            icov[r_idx][seq_i][np.diag_indices_from(icov[r_idx][seq_i])] = dtsum.dot((x/2) - (x**2/3))

            ## Diagonal terms, cross times
            for i in range(1, len(times)):
                dt = times[i] - times[i-1]
                icov[r_idx][seq_i][np.diag_indices_from(icov[r_idx][seq_i])] -= dt * (x[i] * x[i-1])/3

    return dx, icov, aa2i


def get_aa_freqs(freq_dir, name, freq_types, replicates):
    '''
    Load amino acid frequencies from file.
    '''
    
    aa_freqs = {}
    for type in freq_types:
        freq_list = []
        for rep in replicates:
            freq_list.append(pd.read_csv(get_aa_freq_file(freq_dir, name, type, rep)))
        aa_freqs[type] = freq_list

    return aa_freqs


def get_max_reads(freq_dir, name, replicates):
    '''
    Get maximum number of reads across all replicates.
    '''
    
    max_reads = 1e5
    for rep in replicates:
        if os.path.isfile(get_reads_file(freq_dir, name, rep)):
            df = pd.read_csv(get_reads_file(freq_dir, name, rep))
            if np.max(df['reads'])>max_reads:
                max_reads = np.max(df['reads'])
        else:
            print('No file %s, setting max_reads=%.2e' % (get_reads_file(freq_dir, name, rep), max_reads))
            
    return max_reads


def get_best_regularization(corrs, gamma_values, corr_cutoff_pct=0.05):
    '''
    Compute best regularization strength from correlation data.
    '''

    corr_thresh = (np.max(corrs)**2 - corrs[0]**2)*corr_cutoff_pct
    gamma_opt = 0.1
    if np.fabs(np.max(corrs)**2-corrs[0]**2)<0.01:
        gamma_opt = gamma_values[0]
    else:
        gamma_set = False
        for i in range(np.argmax(corrs), 0, -1):
            if np.fabs((corrs[i]**2-corrs[i-1]**2)/(np.log10(gamma_values[i])-np.log10(gamma_values[i-1]))) >= corr_thresh:
                gamma_opt = gamma_values[i]
                gamma_set = True
                break
        if not gamma_set:
            gamma_opt = gamma_values[np.argmax(corrs)]
    
    return gamma_opt


def plot_regularization(corrs, gamma_values):
    '''
    Plot correlation as a function of regularization strength.
    '''
    
    plt.plot(gamma_values, corrs)
    plt.xscale('log')
    plt.xlabel('Regularization strength (gamma)')
    plt.ylabel('Average correlation between replicates')
    plt.show()


def infer_correlated(name, n_replicates, corr_cutoff_pct, freq_dir='.', output_dir='.', with_epistasis=False, plot_gamma=True):
    '''
    Infer selection coefficients from data that includes correlations (counts)
    between pairs of mutant sites. This is not possible with data that includes
    single codon/amino acid counts only.
    
    Required arguments:
        - name: String that will be associated with saved files, serving as an
            identifier for the data set
        - n_replicates: Number of experimental replicates
        - corr_cutoff_pct: Cutoff on the maximum allowed drop in correlation
            between replicates, used to determine the optimal regularization
            strength

    Optional arguments:
        - freq_dir (default: '.'): Directory where variant frequency data has been 
            saved; assumes data saved in the format of compute_variant_frequencies
        - output_dir (default: '.'): Directory where selection coefficients will be 
            saved
        - with_epistasis (default: False): Flag to infer epistatic interactions; 
            this requires correlation data for sets of three and four codons/amino 
            acids; not currently functional (implemented separately in C++ code in 
            the epistasis_inference directory)
        - plot_gamma (default: True): Plot correlation between replicates as a 
            function of the regularization strength gamma
    '''
    ################################################################################

    # Get codon frequencies and auxiliary information
    freq_types = ['single', 'double']
    if with_epistasis:
        # freq_types = freq_types + ['triple', 'quadruple']
        pass
        
    replicates = [i+1 for i in range(n_replicates)]
    aa_freqs = get_aa_freqs(freq_dir, name, freq_types, replicates)

    # Get frequency change and covariance, used to compute selection coefficients, and map to indices
    dx, icov, p2i = compute_dx_covariance(aa_freqs)
    
    # Compute optimal regularization value
    ## Get correlations for each value of gamma
    max_reads = get_max_reads(freq_dir, name, replicates)
    gamma_values = np.logspace(np.log10(1/max_reads), 4, num=20)
    corrs = []
    for g in gamma_values:
        s = np.zeros_like(dx)
        for r_idx in range(n_replicates):
            s[r_idx] = np.inner(np.linalg.inv(icov[r_idx] + g*np.eye(len(icov[r_idx]))), dx[r_idx])
        corrs.append(np.mean([st.pearsonr(s[i], s[j]).statistic for i in range(n_replicates) for j in range(i+1, n_replicates)]))

    ## (Optional) plot the results
    if plot_gamma:
        plot_regularization(corrs, gamma_values)

    ## Select best regularization value
    gamma_opt = get_best_regularization(corrs, gamma_values, corr_cutoff_pct)
    print('Found best regularization strength gamma = %.1e, R = %.2f' % (gamma_opt, corrs[list(gamma_values).index(gamma_opt)]))
    
    ## Compute selection coefficients at optimal gamma
    s = np.zeros_like(dx)
    for r_idx in range(n_replicates):
        s[r_idx] = np.inner(np.linalg.inv(icov[r_idx] + gamma_opt*np.eye(len(icov[r_idx]))), dx[r_idx])
        
    s_joint = np.inner(np.linalg.inv(np.sum(icov, axis=0) + gamma_opt*np.eye(len(icov[0]))), np.sum(dx, axis=0))

    # Convert selection coefficients to a data frame and save to file
    df_ref = aa_freqs['single'][0][aa_freqs['single'][0]['WT_indicator']==True]
    ref_aa = {}
    for df_iter, row in df_ref.iterrows():
        ref_aa[row['site']] = row['aa']
        
    sel_cols = ['site', 'amino_acid', 'WT_indicator'] + ['rep_%d' % r for r in replicates] + ['joint']
    sel_data = []
    for pair, loc in p2i.items():
        sel_data.append([pair[0], pair[1], pair[1]==ref_aa[pair[0]]] + [s[r][loc] for r in range(n_replicates)] + [s_joint[loc]])

    path = get_selection_file(output_dir, name, file_ext='.csv.gz')
    df_temp = pd.DataFrame(data=sel_data, columns=sel_cols)
    df_temp.to_csv(path, index=False, compression='gzip')


def infer_independent(name, n_replicates, corr_cutoff_pct, freq_dir='.', output_dir='.', plot_gamma=True):
    '''
    Infer selection coefficients from data that consists of single amino acid
    frequencies only.
    
    Required arguments:
        - name: String that will be associated with saved files, serving as an
            identifier for the data set
        - n_replicates: Number of experimental replicates
        - corr_cutoff_pct: Cutoff on the maximum allowed drop in correlation
            between replicates, used to determine the optimal regularization
            strength

    Optional arguments:
        - freq_dir (default: '.'): Directory where variant frequency data has been
            saved; assumes data saved in the format of compute_variant_frequencies
        - output_dir (default: '.'): Directory where selection coefficients will be
            saved
        - plot_gamma (default: True): Plot correlation between replicates as a 
            function of the regularization strength gamma
    '''
    ################################################################################

    # Get amino acid frequencies and auxiliary information
    freq_types = ['single']
    replicates = [i+1 for i in range(n_replicates)]
    aa_freqs = get_aa_freqs(freq_dir, name, freq_types, replicates)

    # Get frequency change and covariance, used to compute selection coefficients, and map to indices
    dx, icov, aa2i = compute_dx_covariance_independent(aa_freqs)
    
    # Compute optimal regularization value
    ## Get correlations for each value of gamma
    max_reads = get_max_reads(freq_dir, name, replicates)
    gamma_values = np.logspace(np.log10(1/max_reads), 4, num=20)

    ## Get correlations for each value of gamma
    L = len(dx[0])
    corrs = []
    for g in gamma_values:
        s = np.zeros_like(dx)
        for r_idx in range(n_replicates):
            for seq_i in range(L):
                s[r_idx][seq_i] = np.inner(np.linalg.inv(icov[r_idx][seq_i] + g*np.eye(len(icov[r_idx][seq_i]))), dx[r_idx][seq_i])
            
        corrs.append(np.mean([st.pearsonr(s[i].flatten(), s[j].flatten()).statistic for i in range(n_replicates) for j in range(i+1, n_replicates)]))

    ## (Optional) plot the results
    if plot_gamma:
        plot_regularization(corrs, gamma_values)

    ## Select best regularization value
    gamma_opt = get_best_regularization(corrs, gamma_values, corr_cutoff_pct)
    print('Found best regularization strength gamma = %.1e, R = %.2f' % (gamma_opt, corrs[list(gamma_values).index(gamma_opt)]))
    
    ## Compute selection coefficients at optimal gamma
    s = np.zeros_like(dx)
    for r_idx in range(n_replicates):
        for seq_i in range(L):
            s[r_idx][seq_i] = np.inner(np.linalg.inv(icov[r_idx][seq_i] + gamma_opt*np.eye(len(icov[r_idx][seq_i]))), dx[r_idx][seq_i])
    
    s_joint = np.zeros_like(dx[0])
    for seq_i in range(L):
        s_joint[seq_i] = np.inner(np.linalg.inv(np.sum([icov[r_idx][seq_i] for r_idx in range(n_replicates)], axis=0) + gamma_opt*np.eye(len(icov[0][seq_i]))), np.sum([dx[r_idx][seq_i] for r_idx in range(n_replicates)], axis=0))

    # Convert selection coefficients to a data frame and save to file
    df_ref = aa_freqs['single'][0][aa_freqs['single'][0]['WT_indicator']==True]
    ref_aa = {}
    for df_iter, row in df_ref.iterrows():
        ref_aa[row['site']] = row['aa']

    sites = list(np.unique(aa_freqs['single'][0]['site']))
        
    sel_cols = ['site', 'amino_acid', 'WT_indicator'] + ['rep_%d' % r for r in replicates] + ['joint']
    sel_data = []
    for seq_i in range(L):
        for aa, loc in aa2i.items():
            sel_data.append([sites[seq_i], aa, aa==ref_aa[sites[seq_i]]] + [s[r][seq_i][loc] for r in range(n_replicates)] + [s_joint[seq_i][loc]])

    path = get_selection_file(output_dir, name, file_ext='.csv.gz')
    df_temp = pd.DataFrame(data=sel_data, columns=sel_cols)
    df_temp.to_csv(path, index=False, compression='gzip')
