#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author: Stephan Fuchs (Robert Koch Institute, MF-1, fuchss@rki.de)

VERSION = "0.0.9"
import os
import argparse
import subprocess
import re
import sys
import gzip
import warnings


def parse_args(CMD=None):
    parser = argparse.ArgumentParser(
        prog="rename_in_gff3.py",
        description="changes genotype in VCFs",
    )
    parser.add_argument("vcf", metavar="FILE", help="vcf file", type=str)
    parser.add_argument(
        "--ao",
        metavar="STR",
        help="tag for read count supporting the respective variant (default: AO)",
        type=str,
        default="AO",
    )
    parser.add_argument(
        "--ro",
        metavar="STR",
        help="tag for read count supporting REF (default: RO)",
        type=str,
        default="RO",
    )
    parser.add_argument(
        "--dp",
        metavar="STR",
        help="tag for total read count at the repsective position (default: DP)",
        type=str,
        default="DP",
    )
    parser.add_argument(
        "--gt",
        metavar="STR",
        help="tag for genotype (default: GT)",
        type=str,
        default="GT",
    )
    parser.add_argument(
        "-o", help="output file (will be overwritten!)", type=str, required=True
    )
    parser.add_argument("--gz", help="bgzip compressed input", action="store_true")
    parser.add_argument(
        "--vf",
        metavar="FLOAT",
        help="minimal alternative variant fraction to set a homozygous genotype (default: 0.9)",
        type=float,
        default=0.9,
    )
    parser.add_argument(
        "--rf",
        metavar="FLOAT",
        help="maximal reference fraction to set a homozygous genotype (default: 0.1)",
        type=float,
        default=0,
    )
    parser.add_argument("--version", action="version", version="%(prog)s " + VERSION)
    return parser.parse_args(CMD)


# open file handles considering compression state
def get_filehandle(in_fname, gz):
    if not gz:
        inhandle = open(in_fname, "r")
    else:
        inhandle = gzip.open(in_fname, "rt")
    return inhandle


def process(
    in_fname,
    out_fname,
    min_vf,
    max_rf,
    ao_tag="AO",
    ro_tag="RO",
    dp_tag="DP",
    gt_tag="GT",
    gz=False,
):
    # sanity checks
    if min_vf > 0 and min_vf <= 0.5:
        warnings.warn(
            f"[WARNING] Minimal variant fraction to set a homozygous genotype (--vf) is below 0.5 ({min_vf}). Assuming you know what you are doing."
        )
    if max_rf > 0 and max_rf >= 0.5:
        warnings.warn(
            f"[WARNING] Maximal reference fraction to set a homozygous genotype (--rf) is below 0.5 ({max_rf}). Assuming you know what you are doing."
        )
    out_gz = out_fname.endswith(".gz")
    intermediate = re.sub("\.gz$", "", out_fname)

    # regex generation
    ao_pattern = re.compile(r"(?:^|\t|;)" + re.escape(ao_tag) + "=([0-9,]+)(?:$|\t|;)")
    ro_pattern = re.compile(r"(?:^|\t|;)" + re.escape(ro_tag) + "=([0-9]+)(?:$|\t|;)")
    dp_pattern = re.compile(r"(?:^|\t|;)" + re.escape(dp_tag) + "=([0-9]+)(?:$|\t|;)")

    with get_filehandle(in_fname, gz) as inhandle:
        with open(intermediate, "w") as outhandle:
            for line_index, line in enumerate(inhandle):
                # skip empty or comment lines
                if len(line.strip()) == 0 or line.startswith("#"):
                    outhandle.write(line)
                    continue

                fields = line.split("\t")

                # find GT position
                gt_pos = fields[8].split(":").index(gt_tag)
                # split FORMAT values
                cols = fields[9].split(":")

                # check line for homo/heterozygous
                gt1, gt2 = cols[gt_pos].split("/")
                if gt1 == gt2:
                    # homozygous
                    outhandle.write(line)
                    continue

                # find values
                ao = ao_pattern.findall(fields[7])
                ro = ro_pattern.findall(fields[7])
                dp = dp_pattern.findall(fields[7])

                if len(ao) > 1:
                    sys.exit(
                        "error: multiple occurrences of "
                        + ao_tag
                        + " tag in line "
                        + str(line_index + 1)
                    )
                if len(ro) > 1:
                    sys.exit(
                        "error: multiple occurrences of "
                        + ro_tag
                        + " tag in line "
                        + str(line_index + 1)
                    )
                if len(dp) > 1:
                    sys.exit(
                        "error: multiple occurrences of "
                        + dp_tag
                        + " tag in line "
                        + str(line_index + 1)
                    )

                # calc fractions and check threshold
                alt_fracs = [int(x) / int(dp[0]) for x in ao[0].split(",")]
                max_alt_fraq = max(alt_fracs)

                ref_frac = int(ro[0]) / int(dp[0])

                # old GT, possibly overwritten
                gt = cols[gt_pos]
                if max_rf > 0:
                    # gt adjustment for RO activated
                    if ref_frac <= max_rf and "0" in cols[gt_pos]:
                        # RO low -> remove from GT
                        # generate new GT
                        if "," in fields[4]:
                            gt = "1/2"
                        else:
                            gt = "1/1"
                if min_vf > 0:
                    # gt adjustment for AO activated
                    if max_alt_fraq >= min_vf:
                        # adjust GT to homozygous call of highest AO
                        # generate new GT
                        gt = str(alt_fracs.index(max_alt_fraq) + 1)  # REF == 0 -> ++1
                        gt = gt + "/" + gt
                # replacing GT info (considering line eventual breaks at end)
                cols[gt_pos] = gt
                if gt_pos == len(cols) - 1:
                    cols[gt_pos] += "\n"
                fields[9] = ":".join(cols)
                outhandle.write("\t".join(fields))
        if out_gz:
            bgzip_outname(intermediate, out_fname)


def bgzip_outname(_file, outfile=None):
    if outfile is not None:
        with open(outfile, "wb") as out_fh:
            with subprocess.Popen(["bgzip", _file, "-c"], stdout=out_fh) as bg_proc:
                ret = bg_proc.wait()
        os.remove(_file)
    else:
        with subprocess.Popen(["bgzip", _file], stderr=subprocess.PIPE) as bg_proc:
            out, err = bg_proc.communicate()
            ret = bg_proc.wait()


def main(CMD=None):
    args = parse_args(CMD)
    process(
        args.vcf, args.o, args.vf, args.rf, args.ao, args.ro, args.dp, args.gt, args.gz
    )


if __name__ == "__main__":
    main()
