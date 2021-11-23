#!/usr/bin/env python

import logging
import argparse
import subprocess
import sys
import os
import glob
import pathlib
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
 
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file", required=True, default=None,
                        help="Fastq file. Default: None.")
    parser.add_argument("--sample_id", required=True, default=None,
                        help="Sample ID. Default: None.")    
    parser.add_argument("--ref", required=True, default=None,
                        help="Reference FASTA. Default: None.")    
    parser.add_argument("--loglevel", required=False, default="DEBUG",
                        help="Set logging level to INFO, WARNING or DEBUG (default).")
    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        logging.error("Input fastq file not found: %s" % args.input_file)
        sys.exit(1)
    if not os.path.exists(args.ref):
        logging.error("Input Ref file not found: %s" % args.ref)
        sys.exit(1)
    return args

def parse_ref(ref):
    logging.info("Parsing Ref file: %s" % ref )

    ref_Dict = {}
    for records in SeqIO.parse(ref, "fasta"):
        ref_Dict[records.seq] = records.id
    return ref_Dict    

def count_reads( input_fastq ):
    logging.info("Counting reads of %s" % input_fastq )
    count_Dict = {}
    count = 0
    for records in SeqIO.parse(input_fastq, "fastq"):
        if count_Dict.get(records.seq):
            count_Dict[records.seq] += 1
            continue 
        else:
            count_Dict[records.seq] = 1
    return count_Dict

def output_count_table( sample_id , ref_Dict , counts_Dict ):
    output_csv = "%s.countMatrix.csv" % ( sample_id )
    logging.info("Writing output to %s" % output_csv )
    output_file = open( output_csv, "a")
    output_file.write("featureID,%s\n" % sample_id )
    
    summary_count = [0,0,0]
    inserts = {}
    for seq, count in counts_Dict.items():
        revComp = Seq(seq).reverse_complement()
        # count of ref_Pos_strand, ref_Neg_strand, not_ref
        if seq in ref_Dict.keys():
            summary_count[0] += count
            continue
        elif revComp in ref_Dict.keys() and revComp in counts_Dict.keys():
            counts_Dict[revComp] += count
            summary_count[1] += count
            #print("%s\t%s\tIN_REF" % (seq , str(count) ) )
            # del counts_Dict[seq]
            continue
        elif revComp in ref_Dict.keys() and revComp not in counts_Dict.keys():
            inserts[revComp] = count
            summary_count[1] += count
            # print("%s\t%s\tREVCOMP_ONLY" % (seq , str(count) ) )
            # del counts_Dict[seq]
            continue
        else:
            summary_count[2] += count
            # print("%s\t%s\tNOT_FOUND" % (seq , str(count) ) )
            # del counts_Dict[seq]
        
    for seq, count in inserts.items():
        counts_Dict[seq] = count
    
    print(*summary_count , sep = "\t")
    
    for seq, seq_id in ref_Dict.items():
        out_string = ''
        if seq in counts_Dict.keys():
            out_string = "%s,%s" % ( seq_id , str(counts_Dict[seq]) )
            del counts_Dict[seq]
        else:
            out_string = "%s,0" % ( seq_id )
        output_file.write("%s\n" % out_string )
    
    return output_file
    
def done():
    logging.info("DONE")


if __name__ == "__main__":
    args = parse_args()
    
    ref_Dict = parse_ref(args.ref)
    counts_Dict = count_reads( args.input_file )
    out_file = output_count_table( args.sample_id , ref_Dict , counts_Dict )
    done()
