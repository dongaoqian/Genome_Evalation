mosdepth -t 64 -b 50000 -n -f $1 50k $2

gzip -d 50k.regions.bed.gz

python CalculateCollapseRatio.py 50k.regions.bed


