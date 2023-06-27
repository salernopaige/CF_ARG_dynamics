#!/usr/bin/env python2.7.9

def run_shortbred(marker_db, infile, outfile, temp_file):
    shortbred_quantify.py --markers marker_db --wgs infile --results outfile --tmp temp_file

run_shortbred(snakemake.config["markers"], snakemake.input, snakemake.output[0], snakemake.output[1])
