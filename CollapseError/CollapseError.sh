mosdepth -t 64 -b 50000 -n -f $1 mapping_50k $2

gzip -d 50k.regions.bed.gz

python CalculateCollapseRatio.py YZ.ngs_50k.regions.bed


