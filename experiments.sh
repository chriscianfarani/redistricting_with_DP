#!/bin/bash
if [ $# -ne 1 ]; then
    echo "Usage: $0 <experiment #>"
    exit 1
fi

set -e

EXP=${1}

# CHANGE THIS PATH TO THE JULIA BINARY ON YOUR SYSTEM
JULIA="${HOME}/.julia/bin/julia"

JULIA_ENV="."

SRC_DIR="./julia"

if [ ${EXP} -eq 1 ]; then 
    echo "Running experiment ${EXP}:"
    echo "Generate ensembles to reproduce Figure 1 for Lousiana State Senate"

    dist_type="SS"
    state="LA"
    num_steps=100000

    exp_name="LA_SS_offsets"

    for thresh in 0.001 0.005 0.01 0.05 0.1 0.2 0.3
    do
        for offset in 0.0 0.0005 0.001 0.002 
        do 
            eps=$(echo "${thresh} - ${offset}" | bc)
            if [ $(echo "${eps} <= 0" | bc) -eq 1 ]; then
                continue
            fi
            $JULIA --project=$JULIA_ENV $SRC_DIR/run_recom.jl ${state} ${num_steps} all dp bg ${dist_type} black ${exp_name} --epsilon ${eps} --burn_time 100000 --thresh_name --nowrite --save_every 10 
        done
    done
elif [ ${EXP} -eq 2 ]; then
    echo "Running experiment ${EXP}:"
    echo "Generate ensemble for Louisiana State Senate testing population balance (OPOV test)"

    state="LA"
    dist_type="SS"
    level="bg"
    num_steps=1000000
    exp_name="LA_SS_pop_bal"

    $JULIA --project=$JULIA_ENV $SRC_DIR/run_recom.jl ${state} ${num_steps} all dp ${level} ${dist_type} black ${exp_name} --epsilon 0.05 --burn_time 100000 --nowrite --save_every 10
elif [ ${EXP} -eq 3 ]; then
    echo "Running experiment ${EXP}:"
    echo "Generate ensemble for Georgia State House testing majority-minority districts (Gingles 1 test)"

    state="GA"
    dist_type="SH"
    level="bg"
    group="black"
    num_steps=1000000
    exp_name="GA_SH_gingles"

    $JULIA --project=$JULIA_ENV $SRC_DIR/run_recom.jl ${state} ${num_steps} burst dp ${level} ${dist_type} ${group} ${exp_name} --epsilon 0.05 --burn_time 100000 --vap --nowrite --save_every 10
elif [ ${EXP} -eq 4 ]; then
    echo "Running experiment ${EXP}:"
    echo "Compute critical offset and no-offset discrepancy rate for Georgia State House"

    state="GA"
    dist_type="SH"
    level="bg"
    num_steps=100000

    $JULIA --project=$JULIA_ENV $SRC_DIR/run_recom_critical_offset.jl ${state} ${num_steps} all dp ${level} ${dist_type} black --burn_time 100000 --nowrite
fi