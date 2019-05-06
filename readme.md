# ascat to binMatrix
to align cnv segemnts to each genomic window, for heatmap.

## getting start
1. get heatmap matrix:

`
python ascat_to_binMatrix.py filename.ascat.txt
`

output matrix will be files named with "filename.ascat.txt.gainloss.binMat.txt" and "filename.ascat.txt.relcn.binMat.txt".

also, file "chromosome.colorband.txt" and 'chromosome.nameprobe.txt'

2. modify r script to draw heatmap, using ComplexHeatmap package.
## file format
> *.ascat cnv file foramt, 6 columns:
`
sample  chr  startpos  endpos  nMajor nMinor
`

