#!/usr/bin/env python3

import pandas as pd
import os
import re
import sys

DEFAULT_FILENAME_PATTERNS = [
    ("illumina_some_institute_1",{
        "regex":r"(?P<date>\d+)_(?P<lab_id>\d{2}-\d{4,5}(-[^_]+)?)_(?P<sample_id>.+)"
                r"_(?P<snum>S\d+)_(?P<lane>L\d{3})_(?P<read>R[12])"
                r"_(?P<running>\d{3})",
        "ambig": ["lane","running"]
        }),
    ("illumina1",{
        "regex":r"(?P<sample_id>.+)_(?P<snum>S[\d]+)_(?P<lane>L[\d]{3})_"
                r"(?P<read>R[12])_(?P<running>\d{3})",
        "ambig":["lane","running"] }),
    ("illumina2",{
        "regex":r"(?P<sample_id>.+)_(?P<snum>S[\d]+)_"
                r"(?P<read>R[12])_(?P<running>\d{3})",
        "ambig":["running"] }),
    ("illumina3",{
        "regex":r"(?P<sample_id>.+)_(?P<snum>S[0-9]+)_(?P<lane>L[0-9]{3})_"
                r"(?P<read>R[12])",
        "ambig":["lane"] }),
    ("illumina4",{
        "regex":r"(?P<sample_id>.+)_S(?P<snum>[0-9]+)_"
                r"(?P<read>R[12])",
        "ambig":[] }),
    ("illumina_fallback",{
        "regex":r"(?P<sample_id>.+)_"
                r"(?P<read>R[12])(?P<residual>_.+)?",
        "ambig":[] }),
    ("SRA",{
         "regex":r"(?P<sample_id>SR.+)_(?P<read>[12])",
         "ambig":[] }),

    ("fallback1",{
         "regex":r"(?P<sample_id>.+)_(?P<read>[A-Za-z0-9]+)",
         "ambig":[] }),
    ("fallback2",{
         "regex":r"(?P<sample_id>.+)\.(?P<read>[A-Za-z0-9]+)",
         "ambig":[],
         "sep":"."})
    ]

PAIRED_READ_REF = {
        "1":1,
        "2":2,
        "R1":1,
        "R2":2,
        "FORWARD":1,
        "REVERSE":2,
        "FWD":1,
        "REV":2,
        "F":1,
        "R":2,
        "P":1,
        "M":2,
        "PLUS":1,
        "MINUS":2,
        "SENSE":1,
        "ANTI":2
        }

def check_and_add_fastq(_files, res, pdir=None, sample_names=None, alt_ids=None):
    single_mode = len(_files) == 1 and sample_names is not None
    patterns=dict(DEFAULT_FILENAME_PATTERNS)
    patterns_in_order = [key for key, _ in DEFAULT_FILENAME_PATTERNS]
    end=r'.(fq|fnq|fastq)(.gz)?'
    for _file in (fl for fl in _files if re.match(f".*{end}", fl)):
        if pdir is None:
            path = os.path.dirname(os.path.realpath(_file))
        else:
            path = os.path.realpath(pdir)
        _file = os.path.basename(_file)
        for pnum, pattern in enumerate(patterns_in_order):
            regex_string = ('{patternregex}{end}'.format(
                                patternregex=patterns[pattern]["regex"],
                                end=end))
            match = re.match(regex_string, _file)
            if match is None:
                if pnum-1 == len(patterns):
                    print(f"FastQ does not meet any known spec: file: {_file}")
                continue
            sep = default_if_not("sep", patterns[pattern],"_")
            ambig_keys=default_if_not("ambig",patterns[pattern], [])
            match_dict = match.groupdict()
            ambig = sep.join(match_dict[a]
                             for a in sorted(ambig_keys))
            nonambig = sep.join(match_dict[na]
                                for na in sorted(match_dict.keys())
                                if na not in ambig_keys + ["read"] 
                                and match_dict[na] is not None)
            sub_sample_id = sep.join(
                    [val if val is not None else ""  
                     for (_id, val) in list(match.groupdict().items())
                     if _id not in ["read"]])
            if pattern.startswith("illumina_some_institute"):
                matchdict = match.groupdict()
                sub_sample_id = "{sid}_{lid}".format(
                        sid=matchdict["sample_id"],
                        lid=matchdict["lab_id"])


            read = "R1"
            read_id = 1
            try:
                read = match_dict["read"]
            except KeyError:
                pass
            try:
                read_id = PAIRED_READ_REF[read.upper()]
            except KeyError:
                print("Warning: Read name \"",read,"\" not known")
                print(f"         Using Regex-pattern: {pattern} {regex_string} File: {_file}")


            key=(path,nonambig,ambig,read)
            if key in res["sample_data"].index:
                break
            new_entry=pd.DataFrame(dict(zip(res["sample_data"].columns,
                    ([a]
                     for a
                     in  [path, nonambig, ambig, read, read_id,
                          "PLACE_HOLDER",
                          sample_names[0] if single_mode
                                         and sample_names is not None
                                         else sub_sample_id,
                          "PLACE_HOLDER",
                          alt_ids[0] if single_mode
                                     and alt_ids is not None
                                     else "PLACE_HOLDER",
                          _file, pattern,
                          '{regexpattern}{end}'.format(
                              regexpattern=patterns[pattern]["regex"],end=end),
                          "new"]))))
            new_entry.set_index(["path","unambig_id","ambig_id","read"],inplace=True, drop=False)
            res["sample_data"]= res["sample_data"].append(new_entry)
            # print(read, nonambig, ambig,  pattern, _file, path, sep='\t')
            break
    #eprint(res["sample_data"][["sub_sample_id", "alt_sub_sample_id"]])

def default_if_not(key, _dict, default):
    try:
        return _dict[key]
    except KeyError:
        return default

def resolve_sample_id_conflicts(sample_data):
    generate_alternative_ids(sample_data)
    prelim_groups = sample_data.groupby(["sub_sample_id"]).groups
    for sub_id in prelim_groups:
        paths_subgroup = sample_data[
                sample_data.sub_sample_id == sub_id].groupby(level="path").groups
        num_subids = len(paths_subgroup)
        if  num_subids == 1:
            continue
        print("Warning: Found multiple samples ({num_subids}) with sample id {sub_id}! DUPLICATE_[Number]_ will be prepended to duplicates")
        for (_id, path) in enumerate(paths_subgroup):
            if _id == 0: # Skip first occurance
                continue
            sample_data.loc[((sample_data.sub_sample_id == sub_id) &
                             (sample_data.path == path)) ,
                             "sub_sample_id"] = "Duplicate_{_id}_{samp}".format(
                                     _id=_id, samp=sub_id)

def generate_alternative_ids(sample_data):
    if not ((sample_data["alt_sub_sample_id"]=="PLACE_HOLDER").any()):
        # All Samples known and will just be passed on
        return
    old_samples = sample_data.loc[
            sample_data.alt_sub_sample_id != "PLACE_HOLDER"].copy()
    old_sample_groups = old_samples.groupby(
            level=["path", "unambig_id"]).groups
    num_old_samples = len(old_sample_groups)
    new_samples = sample_data.loc[
            sample_data.alt_sub_sample_id == "PLACE_HOLDER"].copy()
    new_sample_groups = new_samples.groupby(
            level=["path","unambig_id"]).groups
    for key in (key
                for key in new_sample_groups if key in old_sample_groups):
        new_sample_data = sample_data.loc[new_sample_groups[key]].copy()
        old_sample_data = sample_data.loc[old_sample_groups[key]].copy()
        old_sample_sub_groups = old_sample_data.groupby(
                level=["ambig_id"]).groups
        new_sample_sub_groups = new_sample_data.groupby(
                level=["ambig_id"]).groups
        for sub_key in (s_key
                    for s_key in new_sample_sub_groups
                    if s_key in old_sample_sub_groups):
            old_alt_id = old_sample_data.loc[
                    old_sample_sub_groups[sub_key]].alt_sub_sample_id[0]
            sample_data.loc[
                    new_sample_sub_groups[sub_key],
                    "alt_sub_sample_id"] = old_alt_id

        for sub_key in (skey
                    for skey in new_sample_sub_groups
                    if skey not in old_sample_sub_groups):
            raise RuntimeError(
                    "Use of given sample names and conflict resolution through"
                    " longest common prefix is not implemented yet!!!")

    new_samples = sample_data.loc[
            sample_data.alt_sub_sample_id == "PLACE_HOLDER"].copy()
    new_sample_groups = new_samples.groupby(
            level=["path","unambig_id"]).groups
    for (_id,(path, unambig)) in enumerate(new_sample_groups):
        main_name = "Sample_{nid}".format(nid=_id + 1 + num_old_samples)
        sample_data.loc[
                new_sample_groups[(path, unambig)],
                                  "alt_sample_id" ] = main_name
        sub_groups = sample_data[
                sample_data.alt_sample_id == main_name].groupby(
                        level=["ambig_id"]).groups
        if len(sub_groups)==1:
            # Easy Case when Sample not split between lanes or
            # in multiple files with differing running number
            sample_data.loc[sample_data.alt_sample_id == main_name,
                            "alt_sub_sample_id"] = main_name
        else:
            #import pdb; pdb.set_trace()
            for (_id, ambig) in enumerate(sub_groups):
                sample_data.loc[((sample_data.alt_sample_id == main_name) &
                               (sample_data.ambig_id == ambig)),
                              "alt_sub_sample_id"] = "{sid}.{nid}".format(
                                      sid=main_name, nid=_id+1)
    if (sample_data["alt_sub_sample_id"]=="PLACE_HOLDER").any():
        raise RuntimeError("Some unique ids could not be assigned")

def generate_sample_config(sample_data):
    samples = dict(
          (key, dict(("read{read}".format(read=read_id), os.path.join(
                            sample_data[(
                                (sample_data.sub_sample_id==key) &
                                (sample_data.read_id == read_id))].path[0],
                            sample_data[(
                                (sample_data.sub_sample_id==key) &
                                (sample_data.read_id == read_id))].file[0]
                            ))
                     for read_id in sample_data[
                         sample_data.sub_sample_id==key].read_id ))
          for key in sample_data.sub_sample_id)
    for key in samples:
        samples[key]["alt_id"] = sample_data[
                sample_data.sub_sample_id==key].alt_sub_sample_id[0]
    return samples


def main(fastq_dir):
    assert os.path.isdir(fastq_dir), f"{fastq_dir} does not exist"
    fastq_all = [f for f in os.listdir(fastq_dir) if re.search(r'.(fq|fnq|fastq)(.gz)?', f)]

    res={"sample_data":pd.DataFrame(
                    columns=["path", "unambig_id", "ambig_id","read", "read_id",
                            "sample_id", "sub_sample_id",
                            "alt_sample_id", "alt_sub_sample_id",
                            "file","match","regex", "state"],
                    ),
            "delta_files":[]}
    res['sample_data'].set_index(["path","unambig_id","ambig_id","read"], inplace=True, drop=False)

    check_and_add_fastq(fastq_all, res, fastq_dir)
    samp_data = res["sample_data"]
    resolve_sample_id_conflicts(samp_data)
    foo = generate_sample_config(samp_data)

    print(f"sample,fastq_1,fastq_2")
    for sample in foo:
        print(f"{sample},{foo[sample]['read1']},{foo[sample]['read2']}")

if __name__ == "__main__":
    fastq_dir = sys.argv[1]
    main(fastq_dir)
