#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Marie Lataretu (Robert Koch Institute, MF-1, lataretum@rki.de)
#adapted from: Stephan Fuchs (Robert Koch Institute, MF-1, fuchss@rki.de)

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
    parser.add_argument('--tt', metavar="STR", help="tag for type (default: TYPE)", type=str, default="TYPE")
    parser.add_argument('-o', help="output file (will be overwritten!)", type=str, required=True)
    parser.add_argument('--gz', help="bgzip compressed input", action="store_true")
    parser.add_argument('--vf', metavar="FLOAT", help="minimal variant fraction to set include an indel into the consensus (default: 0.9)", type=float, default=0.9)
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    return parser.parse_args(CMD)

# open file handles considering compression state
def get_filehandle(in_fname, gz):
    if not gz:
        inhandle = open(in_fname, "r")
    else:
        inhandle = gzip.open(in_fname, "rt")
    return inhandle
                
def process(in_fname, out_fname, min_vf, ao_tag="AO", dp_tag="DP", type_tag="TYPE", gz=False):
    #sanity checks
    if min_vf <= 0:
        sys.exit("error: min_vf has to be greater than 0")
    out_gz = out_fname.endswith(".gz")
    intermediate = re.sub("\.gz$","", out_fname)

    #regex generation
    ao_pattern = re.compile(r"(?:^|\t|;)" + re.escape(ao_tag) + "=([0-9,]+)(?:$|\t|;)")
    dp_pattern = re.compile(r"(?:^|\t|;)" + re.escape(dp_tag) + "=([0-9]+)(?:$|\t|;)")

    with get_filehandle(in_fname, gz) as inhandle: 
        with open(intermediate, "w") as outhandle: 
            for l, line in enumerate(inhandle):
                #skip empty or comment lines
                if len(line.strip()) == 0 or line.startswith("#"):
                    outhandle.write(line)
                    continue

                fields = line.split("\t")
                
                #find TYPE
                info = { i.split('=')[0]: i.split('=')[1] for i in fields[7].split(";") }

                if not ( 'ins' in info[type_tag] or 'del' in info[type_tag] ):
                    # not a indel
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
                
                if m > min_vf:
                    outhandle.write(line)
                    continue
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
    process(args.vcf, args.o, args.vf, args.ao, args.dp, args.tt, args.gz)

if __name__ == "__main__":
    main()
