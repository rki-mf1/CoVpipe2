#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Stephan Fuchs (Robert Koch Institute, MF-1, fuchss@rki.de)

VERSION = "0.0.9"
import os
import argparse
import subprocess
import re
import sys
import gzip

def parse_args(CMD=None):
    parser = argparse.ArgumentParser(prog="rename_in_gff3.py", description="changes genotype in VCFs", )
    parser.add_argument('vcf', metavar="FILE", help="vcf file", type=str)
    parser.add_argument('-o', help="output file (will be overwritten!)", type=str, required=True)
    parser.add_argument('--gz', help="bgzip compressed input", action="store_true")
    return parser.parse_args(CMD)

def get_cmd_from_snake(snakemake):
    cmd = ["-o", snakemake.output[0],
          "--gz", snakemake.input[0]]
    cmd = [str(arg) for arg in cmd]
    return cmd

# open file handles considering compression state
def get_filehandle(in_fname, gz):
    if not gz:
        inhandle = open(in_fname, "r")
    else:
        inhandle = gzip.open(in_fname, "rt")
    return inhandle
                
def process(in_fname, out_fname, gz=False):
    out_gz = out_fname.endswith(".gz")
    intermediate = re.sub("\.gz$","", out_fname)


    with get_filehandle(in_fname, gz) as inhandle: 
        with open(intermediate, "w") as outhandle: 
            for l, line in enumerate(inhandle):
       
                #skip empty or comment lines
                if len(line.strip()) == 0 or line.startswith("#"):
                    outhandle.write(line)
                    continue
                
                #checking line for mixed variants
                fields = line.split("\t")
                ref = fields[3]
                if len(ref) == 1: 
                    outhandle.write(line)
                    continue
                fields[4] = ",".join([x + "?"*(len(ref)-len(x)) for x in fields[4].split(",")])
                outhandle.write("\t".join(fields))

        if out_gz:
            bgzip_outname(intermediate, out_fname)
        
def bgzip_outname(_file, outfile=None):
    if outfile is not None:
        with open(outfile, "wb") as out_fh:
            with subprocess.Popen(["bgzip", _file, "-c"], 
                                   stdout=out_fh) as bg_proc:
                ret = bg_proc.wait() 
        os.remove(_file)  
    else:
        with subprocess.Popen(["bgzip", _file], 
                               stderr=subprocess.PIPE) as bg_proc:
            out, err = bg_proc.communicate()
            ret = bg_proc.wait() 
                
def main(CMD=None):
        args = parse_args(CMD)
        process(args.vcf, args.o, args.gz)

if __name__ == "__main__":
        CMD = None
        if "snakemake" in globals(): 
            CMD = get_cmd_from_snake(snakemake)
        main(CMD=CMD)
