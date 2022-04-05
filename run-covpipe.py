from re import sub
import yaml
import subprocess


yaml_path = '/home/lataretum/scratch/projects/covpipenext/params.yml'

with open(yaml_path, "r") as stream:
    try:
        yaml = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

params_string = []
for param in yaml:
    if yaml[param]:
        if yaml[param] != True:
            params_string.append(f'--{param}')
            params_string.append(f'{yaml[param]}')
        else:
            params_string.append(f'--{param}')
    else:
        params_string.append(f'--{param}')
        params_string.append(f'false')

# print(["nextflow", "run", "-profile", "mamba,slurm", "-resume"] + params_string)

subprocess.run(["nextflow", "run", "CoVpipeNext.nf", "-profile", "mamba,slurm", "-resume"] + params_string)

report_out_dir = yaml['output'] + '/' + yaml['runinfo_dir'] 
subprocess.run(['cp', '.nextflow.log', f'{report_out_dir}/nextflow.log'])