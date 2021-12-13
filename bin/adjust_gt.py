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
    parser.add_argument('--ao', metavar="STR", help="tag for read count supporting the respective variant (default: AO)", type=str, default="AO")
    parser.add_argument('--dp', metavar="STR", help="tag for total read count at the repsective position (default: DP)", type=str, default="DP")
    parser.add_argument('--gt', metavar="STR", help="tag for genotype (default: GT)", type=str, default="GT")
    parser.add_argument('-o', help="output file (will be overwritten!)", type=str, required=True)
    parser.add_argument('--gz', help="bgzip compressed input", action="store_true")
    parser.add_argument('--vf', metavar="FLOAT", help="minimal variant fraction to set a homogeneous genotype (default: 0.9)", type=float, default=0.9)
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    return parser.parse_args(CMD)

def get_cmd_from_snake(snakemake):
    cmd = ["-o", snakemake.output[0],
           "--vf", snakemake.params.frac,
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
                
def process(in_fname, out_fname, min_vf, ao_tag="ao", dp_tag="AO", gt_tag="GT", gz=False):
    #sanity checks
    if min_vf <= 0.5:
        sys.exit("error: min_vf has to be greater than 0.5")
    out_gz = out_fname.endswith(".gz")
    intermediate = re.sub("\.gz$","", out_fname)

    #regex generation
    ao_pattern = re.compile(r"(?:^|\t|;)" + re.escape(ao_tag) + "=([0-9,]+)(?:$|\t|;)")
    dp_pattern = re.compile(r"(?:^|\t|;)" + re.escape(dp_tag) + "=([0-9]+)(?:$|\t|;)")
    gt_pattern = re.compile(r"(?:^|\t|;)" + re.escape(gt_tag) + "=([0-9]+/[0-9]+)(?:$|\t|;)")

    with get_filehandle(in_fname, gz) as inhandle: 
        with open(intermediate, "w") as outhandle: 
            for l, line in enumerate(inhandle):
       
                #skip empty or comment lines
                if len(line.strip()) == 0 or line.startswith("#"):
                    outhandle.write(line)
                    continue
                
                #checking line for mixed variants
                fields = line.split("\t")
                if "," not in fields[4]:
                    outhandle.write(line)
                    continue
              
                #find ao and dp
                ao = ao_pattern.findall(fields[7])
                dp = dp_pattern.findall(fields[7])
                
                if len(ao) > 1:
                    sys.exit("error: multiple occurences of " + ao_tag + " tag in line " + str(l+1))
                if len(dp) > 1:
                    sys.exit("error: multiple occurences of " + dp_tag + " tag in line " + str(l+1))
                
                #calc fractions and check threshold
                fracs = [int(x)/int(dp[0]) for x in ao[0].split(",")]
                m = max(fracs)
                
                if m < min_vf:
                    outhandle.write(line)
                    continue
                
                #generate new GT
                gt = str(fracs.index(m))
                gt = gt + "/" + gt
                
                #find GT position
                gt_pos = fields[8].split(":").index(gt_tag)
                
                #replacing GT info (considering line eventual breaks at end)
                cols = fields[9].split(":")
                cols[gt_pos] = gt
                if gt_pos == len(cols)-1:
                    cols[gt_pos] += "\n"
                fields[9] = ":".join(cols)
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
        process(args.vcf, args.o, args.vf, args.ao, args.dp, args.gt, args.gz)

if __name__ == "__main__":
        CMD = None
        if "snakemake" in globals(): 
            CMD = get_cmd_from_snake(snakemake)
        main(CMD=CMD)
