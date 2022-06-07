process sc2rf {
    label 'sc2rf'

    input:
    tuple val(name), path(fasta)

    output:
    tuple val(name), path("${name}_sc2rf.csv"), emit: csv

    script:
    """
    git clone https://github.com/lenaschimmel/sc2rf.git
    cd sc2rf
    python3 sc2rf.py --csvfile ../${name}_sc2rf.csv --parents 1-35 --breakpoints 0-10 \
                        --max-intermission-count 3 --max-intermission-length 1 \
                        --unique 1 --max-ambiguous 10000 --max-name-length 55 \
                        ../${fasta} > /dev/null 2>&1
    cd ..

    # add empty csv line, if no output
    if [[ \$(wc -l <${name}_sc2rf.csv) -eq 1 ]] 
    then
        echo ${name},,,, >> ${name}_sc2rf.csv
    fi
    """
    stub:
    """
    touch ${name}_sc2rf.csv
    """
}
