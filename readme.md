# ascat to binMatrix
to align cnv segemnts to each genomic window, for heatmap.

## file format
1. *.ascat cnv file foramt, 6 columns:
`
sample  chr  startpos  endpos  nMajor nMinor
`

## getting start
`
python ascat_to_binMatrix.py filename.ascat.txt
`
output will be files named with "filename.ascat.txt.gainloss.binMat.txt" and "filename.ascat.txt.relcn.binMat.txt".