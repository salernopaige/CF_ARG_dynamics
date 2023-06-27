
import subprocess
import csv
import argparse
import re
import os
import sys
import tqdm
import pickle
import multiprocessing
import pandas as pd
from pathlib import Path

def read_total_reads(infile_fastq):
    cmd = ["zcat", infile_fastq]
    cmd2 = ["wc", "-l"]

    ps = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    output = subprocess.check_output(cmd2, stdin=ps.stdout)
    ps.wait()

    return int(output.decode('ascii')) / 4

def frac_sam_aligned(infile_sam):
    print("Currently on {}".format(infile_sam))
    unique_reads = set()
    unique_reads_hit = set()
    with open(infile_sam) as f:
        for line in csv.reader(f, delimiter = "\t"):
            if line[0][0] == "@": #Header
                continue
            unique_reads.add(line[0]) 
            if line[2] == "*": #Unmapped
                continue
            if line[5] == "*": #No CIGAR string
                continue
            #if int(line[4]) <= 20: #Bad quality
            #    continue
            read_len = sum([int(x) for x in re.split("[A-Z=]{1}",line[5]) if x != ""])
            m = re.search("([0-9]+)M",line[5])
            if m:
                hits = int(m.group(1))
            else:
                hits = 0
            if hits / read_len < 0.9:
                continue
            unique_reads_hit.add(line[0])
    return len(unique_reads_hit) / len(unique_reads) * 100


def load_micop_output(infile):
    ab_dict = {}
    df = pd.read_csv(infile, sep = "\t")
    levels = ["k","p","c","o","f","g","s", "s"]
    ab_dict = {}
    for (i,dat) in df.iterrows():
        id_split = dat["TAXPATH"].split("|")
        tax_split = dat["TAXPATHSN"].split("|")
        assert len(id_split) == len(tax_split)
        tax_str = ""
        for j in range(len(id_split)):
            tax_str = tax_str + "%s_%s_%s" %(levels[j],id_split[j],tax_split[j]) + "|"
        tax_str = tax_str.rstrip("|")
        ab_dict[tax_str] = dat["PERCENTAGE"]
    return ab_dict


def load_micop_full(indir = "/dartfs/rc/lab/R/RossB/DartCF_infant_meta/Results/"):
    results_all = {}
    for indir_results in Path(indir).glob("*"):
        if not indir_results.is_dir():
            continue
        infile = indir_results / ("{}_micop_fungi_abundance.tsv".format(indir_results.name))
        infile_raw = indir_results / ("{}_micop_fungi_abundance_rawcounts.tsv".format(indir_results.name))
        if (not infile.exists()) | (not infile_raw.exists()):
            print("Missing file for {}".format(indir_results.name))
            continue
        results_all[indir_results.name] = load_micop_output(infile)
    return pd.DataFrame(results_all).fillna(0)

def count_raw(infile_raw):
    df = pd.read_csv(infile_raw, sep = "\t")
    if df.shape[0] == 0:
        print("Problem with {}".format(infile_raw.name))
        return 0
    assert df.iloc[0,3] == "Eukaryota"
    return df.iloc[0,4]

def check_sams(indir):
    results_all = {}
    for indir_results in Path(indir).glob("*"):
        if not indir_results.is_dir():
            continue
        infile_sam = indir_results / ("{}_micop_fungi.sam".format(indir_results.name))
        if (not infile_sam.exists()):
            continue
        if os.stat(infile_sam).st_size == 0:
            print("Problem with {}".format(infile_sam))



def mutli_proc_calcs(col, indir, outdir):
    print("Currently on {}".format(col))
    infile_fastq = indir / "{}/{}_R1.qc.humanDecontaminated.fastq.gz".format(col,col)
    infile_raw = outdir / "{}/{}_micop_fungi_abundance_rawcounts.tsv".format(col,col)
    return count_raw(infile_raw) / read_total_reads(infile_fastq) * 100

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Combine the micop results from different samples.')
    parser.add_argument("--indir",type=str)
    parser.add_argument("--outdir",type=str)
    parser.add_argument("--outfile",type=str)
    parser.add_argument("--infile_pickle",type=str)
    parser.add_argument("--n_cpu",type=int)
    args = parser.parse_args()
    indir = Path(args.indir)
    outdir = Path(args.outdir)
    outfile = indir / args.outfile
    infile_pickle = indir / args.infile_pickle

    check_sams(outdir)

    df_full = load_micop_full(outdir)
    print(df_full)
    if not Path(infile_pickle).exists():
        pool = multiprocessing.Pool(processes = args.n_cpu)
        all_frac = pool.starmap(mutli_proc_calcs, [(col, indir, outdir) for col in df_full.columns])
        pickle.dump( all_frac, open( infile_pickle, "wb" ) )
    else:
        all_frac = pickle.load( open( infile_pickle, "rb" ) )
    print(all_frac)
    df_frac = pd.DataFrame(all_frac).T
    df_frac.columns = df_full.columns
    df_frac.index = ["PercentageReads"]
    pd.concat([df_frac,df_full]).to_csv(outfile, sep = "\t")
