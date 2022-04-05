WORKDIR='/home/lataretum/scratch/projects/covpipenext/'
DATE='2022-42-42'
ID='test-id-bla'
FASTQ='/home/lataretum/scratch/projects/covpipenext/sample_sheet.csv'

# conda create -n nextflow -c conda-forge -c anaconda -c bioconda nextflow=21.04.0 dos2unix xlrd pandas csvkit
eval "$(conda shell.bash hook)"
conda activate nf

# should be the dir where the git repo lives
CURRENT_DIR=$PWD
DIR=${WORKDIR}/${DATE}/${ID}

mkdir -p $DIR && cd $DIR

# check if sample sheet is empty
if [[ ! -s $FASTQ ]]; then
    echo "ERROR! Sample sheet $FASTQ seems to be empty, please check!"
    exit 1
fi

# run pipeline 

# first pull the latest revisions and select the newest per default
# or select a specific release (see full list via `nextflow info RKIBioinformaticsPipelines/covpipenext`)
nextflow pull RKIBioinformaticsPipelines/covpipenext -hub gitlab
REVISION=$(nextflow info  RKIBioinformaticsPipelines/covpipenext -hub gitlab | sed 's/ [*]//' | sed 's/ //g' | sed 's/\[t\]//g' | awk 'BEGIN{FS=" "};{if($0 ~ /^ *0/ || $0 ~ /^ *1/ ){print $0}}' | sort -Vr | head -1)
nextflow run RKIBioinformaticsPipelines/covpipenext -hub gitlab -r $REVISION -profile slurm,singularity -w $DIR/work \
    --output $DIR/results --list --fastq $FASTQ --run_id $ID \
    -params-file '/home/lataretum/scratch/projects/covpipenext/params.yml'

# set permissions
chmod -R g+rwX ${WORKDIR}/${DATE}/ --silent
chmod 775 $CACHEDIR/* --silent
