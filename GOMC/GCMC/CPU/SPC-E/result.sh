#!/bin/bash
GOMC=/Projects/msoroush/proposal/NIH/wsu-uiuc/codes/GOMC/bin/GOMC_GPU_GCMC
cudaFile=/Common/linux/cuda-10.2
LOG="-out.log"
CPU=4
SYSTEM=(010k 050k 100k)
CBMC=(cbmc-001-001 cbmc-005-005 cbmc-010-010 cbmc-020-010 cbmc-020-020 cbmc-050-010 cbmc-050-020 cbmc-050-030 cbmc-100-010 cbmc-100-020 cbmc-100-030)

BASE_DIR=`pwd`

for i in "${!SYSTEM[@]}"
do
    sys=${SYSTEM[i]}
    echo "Working on: ${sys}"
    echo ""

    SYS_DIR=${BASE_DIR}/${sys}

    OUT_FILE_ACC=${SYS_DIR}/"${sys}-transfer-acceptance.dat"
    OUT_FILE_TIME=${SYS_DIR}/"${sys}-transfer-profile.dat"
    echo -e "#CBMC"   "\t %Acceptance_AVG" "\t %Acceptance_insert" "\t %Acceptance_delete" > ${OUT_FILE_ACC}
    echo -e "#CBMC"   "\t %Average-CPU (nano sec)" "\t Average-CPU (nano sec)" > ${OUT_FILE_TIME}
    echo -e "#CBMC"   "\t %Short-Range" "\t Long-range" >> ${OUT_FILE_TIME}

    if [ ! -d ${SYS_DIR} ]; then
        echo "There is no ${SYS_DIR} directory!"
        continue
    fi
    cd ${SYS_DIR}

    for j in "${!CBMC[@]}"
    do
        cbmc=${CBMC[j]}
        echo "Running CBMC trials: ${cbmc}"

        CBMC_DIR=${SYS_DIR}/${cbmc}
        if [ ! -d ${CBMC_DIR} ]; then
            echo "There is no ${CBMC_DIR} directory!"
            continue
        fi
        cd ${CBMC_DIR}

        FILE="${sys}-${cbmc}"
        LOG_FILE="${FILE}${LOG}"
        if [ ! -f ${LOG_FILE} ]; then
            echo "There is no ${LOG_FILE} file!"
            continue
        fi

        box_0=`grep -w "Accepted Mol-Transfer" ${LOG_FILE} | awk '{print $5}'`
        box_1=`grep -w "Accepted Mol-Transfer" ${LOG_FILE} | awk '{print $6}'`
        acceptance=$( echo "0.5 * (${box_0} + ${box_1})" | bc )
        echo -e "${cbmc} \t ${acceptance} \t ${box_0} \t ${box_1}" >> ${OUT_FILE_ACC}

        cpu=`grep -w "transform_swap_move" ${LOG_FILE} | awk '{print $4}'`
        gpu=`grep -w "ewald_molecule_swap_recip_energy" ${LOG_FILE} | awk '{print $4}'`
        echo -e "${cbmc} \t ${cpu} \t ${gpu}" >> ${OUT_FILE_TIME}
    done

done


echo ""
echo "Finished"
