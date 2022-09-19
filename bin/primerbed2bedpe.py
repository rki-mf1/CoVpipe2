#!/usr/bin/env python3

import argparse
import pandas as pd
import re

version="0.0.1"

# Author: https://gitlab.com/rekm in the project https://gitlab.com/RKIBioinformaticsPipelines/ncov_minipipe
# small adaptions from https://gitlab.com/MarieLataretu

def main():
    parser = get_argparser()
    args = parser.parse_args()
    reformat_bedpe(args)


def get_argparser():
    parser = argparse.ArgumentParser("bed2bedpe")
    parser.description = (
            "Tool to quickly transform a primer bed file to a bedpe\n"
            "It uses the naming sheme of primers to asign them.\n "
            "amp1(_LEFT)\n"
            "amp2(_RIGHT)\n"
            "{_LEFT|_RIGHT} are the default identifiers the tool looks for. "
            "When multiple pairs are possible for an amplicon,"
            " the cartesian product is used.\n"
            "Multiple paired end records with the same id will be outputted\n"
            " and the score field will be incremented for each duplicate")
    parser.add_argument("primer_bed", help="bed file with primer positions")
    parser.add_argument("--forward_identifier", default="_LEFT",
                        help="identifier used to identify left amplicon primer."
                             " default: _LEFT")
    parser.add_argument("--reverse_identifier", default="_RIGHT",
                        help="identifier used to identify right amplicon primer."
                             " default: _RIGHT")
    parser.add_argument("-o","--output",default="primer.bedpe", 
            help="Output file name.  default: primer.bedpe")
    return parser

def reformat_bedpe(args):
    bed_file = args.primer_bed
    outfile = args.output
    for_id = args.forward_identifier
    rev_id = args.reverse_identifier

    with open(bed_file, 'r') as bed:
        for line in bed:
            assert len(line.split("\t")) >= 4, "ERROR: Check BED file format. Found less than four tab seperated entries in line."

    bed =  pd.read_csv(bed_file, sep="\t",  usecols=[0, 1, 2, 3],
                       names=["chrom","start","end","name"])
    pattern = re.compile("(.+)(%s|%s)(.*)" % (for_id,rev_id))
    func = lambda x: split_pattern(x, pattern, for_id, rev_id)
    parse_bed_ids = pd.concat(list(bed.iloc[:,3].apply(func))) 
    parse_bed_ids = pd.merge(bed, parse_bed_ids,how="left", left_on="name", right_on="key")
    amplicons = parse_bed_ids.groupby("id")
    assert amplicons.ngroups > 0, "ERROR: No primer paires to write. Check name foramt in BED file."
    with open(outfile,"w") as o_fh:
        for amplicon_id, keys in amplicons:
            assert keys.shape[0] >= 2, f'ERROR: Cloud not find at least two primers for amplicon {amplicon_id}.'
            if keys.shape[0] > 2:
                print(f'LOG: Found more than two primers for amplicon {amplicon_id}')
            (left, l_data), (right, r_data) = keys.groupby("reverse") 
            for (idx,(li,ri)) in enumerate([(x,y) for x in range(len(l_data)) for y in range(len(r_data))]):
                out = []
                l_temp = l_data.iloc[li,:]
                r_temp = r_data.iloc[ri,:]
                out.extend(list(l_temp[["chrom","start","end"]]))
                out.extend(list(r_temp[["chrom","start","end"]]))
                out.append(amplicon_id)
                out.append(idx)
                # out.extend(list(l_temp["strand"])) # not respected by bamclipper
                # out.extend(list(r_temp["strand"])) # not respected by bamclipper
                print(*out, sep="\t", file=o_fh)
    assert sum(1 for line in open(outfile)) >= sum(1 for line in open(bed_file))/2, "ERROR: BED to BEDPE conversion might be corrupted. Wrote less primer pairs than expected."


def split_pattern(val, pattern, for_id, rev_id):
    ret = {"key":val,
           "id": None,
           "reverse": False,
           "orient_id": None,
           "add": None}
    match = pattern.match(val)
    if match is not None:
        ret["id"] = match.groups()[0]
        ret["reverse"] = match.groups()[1] in rev_id
        ret["orient_id"] = match.groups()[1]
        ret["add"] = match.groups()[2]
    return pd.DataFrame(ret, index=[ret["key"]])


if __name__ == "__main__":
    main()
