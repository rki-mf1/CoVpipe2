process sc2rf {
    label 'sc2rf'

    input:
    tuple val(name), path(fasta)

    output:
    tuple val(name), path("${name}_sc2rf.csv"), emit: csv

    script:
    """
    git clone https://github.com/ktmeaton/ncov-recombinant.git
    cd ncov-recombinant/sc2rf
    # --clades can't be the last argument before the fasta input, else script fails without an error, see https://github.com/lenaschimmel/sc2rf/issues/35
    python3 sc2rf.py --csvfile ../../${name}_sc2rf.csv --parents 1-1000 --breakpoints 1-2 \
                        --ansi \
                        --max-intermission-count 3 --max-intermission-length 1 \
                        --clades 'all' \
                        --unique 1 --max-ambiguous 10000 --max-name-length 55 \
                         ../../${fasta}
    cd ../..

    if [[ \$(wc -l <${name}_sc2rf.csv) -eq 0 ]] 
    then
        # add empty csv line, if no output at all (e.g. script teminates with parameter problem, but exits with 0)
        echo ${name},,,,,, >> ${name}_sc2rf.csv
    elif [[ \$(wc -l <${name}_sc2rf.csv) -eq 1 ]] 
    then
        # only header as output
        # remove header because of weired newline character after header that breaks the R report
        # overwrite with empty csv line
        echo ${name},,,,,, > ${name}_sc2rf.csv
    else
        # remove header because of weired newline character after header that breaks the R report
        sed -i 1d ${name}_sc2rf.csv
    fi
    """
    stub:
    """
    touch ${name}_sc2rf.csv
    """
}
