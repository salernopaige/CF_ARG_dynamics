
import re
import csv
from tqdm import tqdm
import pandas as pd
from Bio import SeqIO
from pathlib import Path
from optparse import OptionParser

def screen_for_secondary(bitwise_flag):
    bitwise_flag_bin = list(reversed(list('{0:012b}'.format(bitwise_flag))))
    if bitwise_flag_bin[8] == "1":
        # This is a secondary alignment
        return True
    return False


def parse_sam_file(infile, READ_LEN = 100.0, MATCH_CUTOFF = 0.0, MISMATCH_CUTOFF = 0.00, QUAL_CUTOFF = 20, keep_multi_reads=False):
    counts = {}
    count_skipped = 0
    if re.search("gz$", infile.name):
        open_fn = gzip.open
    else:
        open_fn = open
    with open_fn(infile) as f:
        reader = csv.reader(f, delimiter = "\t", quoting=csv.QUOTE_NONE)
        for line in reader:
            # Skip through the header
            if line[0][0] == "@":
                continue
            seq_name = line[0]
            bitwise_flag = int(line[1])
            hit_seq = line[2]
            map_qual = int(line[4])
            cigar_string = line[5]

            # Ignore this if we have no hits
            if hit_seq == "*":
                continue
            if cigar_string == "*":
                count_skipped += 1
                continue

            m = re.search("([0-9]*)M",cigar_string)
            assert m, "Can't find matches in the cigar_string"
            n_match = float(m.group(1))

            # Ignore those which don't hit enough baseas
            # TODO: Fix this, READ_LEN should not be a paramter
            if n_match / READ_LEN < MATCH_CUTOFF:
                count_skipped += 1
                continue
            n_mismatch = float(re.search("XM:i:([0-9]*)","\t".join(line)).group(1))
            if n_mismatch / n_match > MISMATCH_CUTOFF:
                count_skipped += 1
                continue

            if not keep_multi_reads:
                # When there are a series of identical hits, this will be low
                if map_qual < QUAL_CUTOFF:
                    count_skipped += 1
                    continue
                # Ignore this is this is a secondary alignment of a read that is otherwise counted
                if screen_for_secondary(bitwise_flag):
                    count_skipped += 1
                    continue

            # Now count these up
            if line[2] not in counts:
                counts[line[2]] = 0

            counts[line[2]] += 1
    return counts, count_skipped


def create_full_count_table(infile_list, infile_seqs, infile_total_reads = None, 
                            match_cutoff=0.50, mismatch_fraction = 0.05, qual_cutoff = 20,
                            keep_multi_reads=False):
    #Need to include all missing sequences
    all_seqs = set()
    with open(infile_seqs) as f:
        for seq in SeqIO.parse(f, "fasta"):
            all_seqs.add(seq.id)
    full_results = {}

    for infile in tqdm(infile_list, desc="processing sam files"):
        infile = Path(infile)
        if not infile.exists():
            print("Can't find file %s" %(infile))
            continue

        counts, counts_skipped = parse_sam_file(infile, MATCH_CUTOFF = match_cutoff, MISMATCH_CUTOFF = mismatch_fraction, QUAL_CUTOFF = qual_cutoff, keep_multi_reads=keep_multi_reads)
        counts["total"] = sum(counts.values())
        counts["skipped"] = counts_skipped

        # Index by the directory name
        full_results[infile.parents[0].name] = counts
        
        assert all([ x in all_seqs for x in counts if x not in ["total","skipped"]]), "I'm seeing sequences in a sam that are not in your sequence file"

    for d in full_results:
        for key in all_seqs:
            if key not in full_results[d]:
                full_results[d][key] = 0

    df = pd.DataFrame(full_results).T
    df["Tag"] = df.index
    cols_order = ["Tag"] + [x for x in df.columns if x not in ["total","skipped","Tag"]] + ["total","skipped"]
    df = df[cols_order]
    if infile_total_reads is not None:
        df_reads = pd.DataFrame.from_csv(infile_total_reads)
        df_merge = df_reads.merge(df, how = "right", right_index = True, left_on = "name")
        return df_merge
    else:
        return df

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-f", "--infile_fasta")
    parser.add_option("-o", "--outfile")
    (options, args) = parser.parse_args()

    assert all([x[-3:] == "sam" for x in args]), "You should only be passing .sam files in the args"

    df = create_full_count_table(args, options.infile_fasta)
    df.to_csv(options.outfile, sep = ",", index=False)
