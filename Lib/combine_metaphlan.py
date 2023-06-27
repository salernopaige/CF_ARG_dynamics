import pandas as pd
from pathlib import Path
import logging
import re
import tqdm
from optparse import OptionParser


def load_single_metaphlan(infile, col_oi="relative_abundance",mode="rel_ab_w_read_stats", tax_level="s"):
    header=3
    if mode == "rel_ab_w_read_stats":
        header=4
    df_meta = pd.read_csv(infile, sep = "\t", header=header)

    if df_meta.index[0] == 0:
        cols = list(df_meta.columns)
        cols[0] = cols[0].replace("#","")
        df_meta.columns = cols
    else:
        # Ocasionally the metaphlan file does not have "additional_species" and it loads the clade_name as the index
        # We need to correct that
        assert mode!="rel_ab_w_read_stats"
        assert df_meta.index[0] == "k__Bacteria"
        df_meta.insert(0, "clade_name", df_meta.index)
        df_meta.columns = ["clade_name", "NCBI_tax_id", "relative_abundance", "additional_species"]

    expr_dict = {}
    for (i, dat) in df_meta.iterrows():
        if (col_oi == "estimated_number_of_reads_from_the_clade") & (dat["clade_name"] == "k__Bacteria"):
            expr_dict[dat["clade_name"]] = dat[col_oi]
        elif re.search("{}__[A-Za-z_]+$".format(tax_level), dat["clade_name"]):
            expr_dict[dat["clade_name"]] = dat[col_oi]
    return expr_dict



if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    parser = OptionParser()
    parser.add_option("-m", "--mode", default="tax_level")
    parser.add_option("-i", "--infiles", default=[], action="append")  #Not used for mode=tax_level. Use args instead
    parser.add_option("-t", "--tax_level", default="s")
    parser.add_option("-p", "--metaphlan_mode", default="rel_ab_w_read_stats")
    parser.add_option("-c", "--col_oi", default="relative_abundance")
    parser.add_option("-o", "--outfile")
    (options, args) = parser.parse_args()
    if options.mode == "tax_level":
        master_dict = {}
        for infile in args:
            infile = Path(infile)
            try:
                master_dict[infile.name] = load_single_metaphlan(infile, options.col_oi, options.metaphlan_mode, options.tax_level)
            except:
                logging.info("problem with {}".format(infile.name))
        df_metaphlan = pd.DataFrame(master_dict).fillna(0)
        df_metaphlan.insert(0, "clade", df_metaphlan.index)
        df_metaphlan.to_excel(options.outfile, index=False)
    if options.mode == "full_bind":
        infiles_use = options.infiles
        df_full = pd.read_excel(infiles_use[0], engine='openpyxl')
        for infile in infiles_use[1:]:
            df = pd.read_excel(infile, engine='openpyxl')
            assert all(df_full.columns == df.columns)
            df_full = df_full.append(df)
        df_full.to_excel(options.outfile, index=False, engine='openpyxl')
