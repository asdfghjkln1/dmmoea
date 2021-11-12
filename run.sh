#!/bin/bash
set +e
# PATH=$PATH:/bin/usr:Rscript
export $WD_PATH=$(pwd)
export PARAMS_NSGA="parameters_nsga.txt"
export PARAMS_DNSGA="parameters_dnsga.txt"
export PARAMS_DMNSGA="parameters_dmnsga.txt"
export SCENARIO="scenario.txt"
export SCENARIO_PARALLEL="scenario_parallel.txt"
export N_RUNS=30

Rscript dmmoea.R "nsga2" $WD_PATH $PARAMS_NSGA $SCENARIO_PARALLEL || true
Rscript dmmoea.R "dnsga2" $WD_PATH $PARAMS_DNSGA $SCENARIO_PARALLEL || true
Rscript dmmoea.R "dmnsga2" $WD_PATH $PARAMS_DMNSGA $SCENARIO || true
Rscript dmmoea_post_tuning.R $WD_PATH $N_RUNS || true

#read -rn1