#!/bin/bash

EVENT=$1
SAMPLE=$2

RBAMS=$(find bams -name '*.bam' | sort | fgrep $SAMPLE)
DBAMS=$(fgrep $SAMPLE /ifs/dmprequest/12-245/bam_fs_locs.txt  | cut -d, -f2 | sort | awk '{print "/ifs"$1}')

echo "Snap R $SAMPLE"
echo ./igv_plotter \
    --igv-jar-path IGV_2.3.97/igv.jar \
    -o isnap_R_${SAMPLE}___ \
    --collapse \
    $RBAMS \
    $EVENT

echo "Snap D"
echo ./igv_plotter \
    --igv-jar-path IGV_2.3.97/igv.jar \
    -o isnap_D_${SAMPLE}___ \
    --collapse \
    $DBAMS \
    $EVENT

echo $DBAMS

convert \
    isnap_D_${SAMPLE}___s1__${EVENT/:/_}.png \
    isnap_R_${SAMPLE}___s1__${EVENT/:/_}.png \
    +append isnap_${SAMPLE}___s1__${EVENT/:/_}.png

echo "Done $SAMPLE"
# rm isnap_D_${SAMPLE}___s1__${EVENT/:/_}.png \
#     isnap_R_${SAMPLE}___s1__${EVENT/:/_}.png \