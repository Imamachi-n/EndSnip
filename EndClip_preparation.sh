#!/bin/bash
echo "-- Prepare map file from annotation gtf file..."
python3 A_make_map_from_gtf.py /mnt/hgfs/Projects/EndClip/EndClip/data/gencode.v24.basic.annotation.gtf /mnt/hgfs/Projects/EndClip/EndClip/data/gencode.v24.basic.annotation.map

echo "-- Prepare 3UTR database file from annotation gtf file..."
#python3 EndClip_prep.py data/EndClip_prep_configure_file.txt
