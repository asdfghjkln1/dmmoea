#!/bin/bash
set +e
# PATH=$PATH:/bin/usr:Rscript
export WD_PATH=$(pwd)
export PARAMS_NSGA="parameters_nsga.txt"
export PARAMS_DNSGA="parameters_dnsga.txt"
export PARAMS_DMNSGA="parameters_dmnsga.txt"
export SCENARIO="scenario.txt"
export SCENARIO_PARALLEL="scenario_parallel.txt"
export N_RUNS=30
export OBJ_FUN="XieBeni"

Rscript dmmoea.R "nsga2" $WD_PATH $OBJ_FUN $PARAMS_NSGA $SCENARIO_PARALLEL > nsga_output.txt 2>&1
Rscript dmmoea.R "dnsga2" $WD_PATH $OBJ_FUN $PARAMS_DNSGA $SCENARIO_PARALLEL > dnsga_output.txt 2>&1
Rscript dmmoea.R "dmnsga2" $WD_PATH $OBJ_FUN $PARAMS_DMNSGA $SCENARIO > dmnsga_output.txt 2>&1
Rscript dmmoea_post_tuning.R $WD_PATH $N_RUNS $OBJ_FUN > post_output.txt 2>&1

#read -rn1
