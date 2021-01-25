#!/bin/bash
GOMC=/Projects/msoroush/proposal/NIH/wsu-uiuc/codes/GOMC/bin/GOMC_GPU_GCMC
cudaFile=/Common/linux/cuda-10.2
LOG="-out.log"
CPU=8
SYSTEM=(010k 100k)
CBMC=(cbmc-001-001 cbmc-005-005 cbmc-010-010 cbmc-020-010 cbmc-020-020 cbmc-050-010 cbmc-050-020 cbmc-050-030 cbmc-100-010 cbmc-100-020 cbmc-100-030)

BASE_DIR=`pwd`

for i in "${!SYSTEM[@]}"
do
    sys=${SYSTEM[i]}
    echo "Working on: ${sys}"
    echo ""

    SYS_DIR=${BASE_DIR}/${sys}
    if [ ! -d ${SYS_DIR} ]; then
        echo "There is no ${SYS_DIR} directory!"
        exit 1
    fi
    cd ${SYS_DIR}

    for j in "${!CBMC[@]}"
    do
        cbmc=${CBMC[j]}
        echo "Running CBMC trials: ${cbmc}"

        CBMC_DIR=${SYS_DIR}/${cbmc}
        if [ ! -d ${CBMC_DIR} ]; then
            echo "There is no ${CBMC_DIR} directory!"
            exit 1
        fi
        cd ${CBMC_DIR}

        FILE="${sys}-${cbmc}"
        LOG_FILE="${FILE}${LOG}"
        ${cudaFile}/bin/nsys profile -t cuda,nvtx -o ${FILE} --stats=true -f true -w true ${GOMC} +p${CPU} in.conf |& tee ${LOG_FILE}

    done

done






echo ""
echo "Finished"
