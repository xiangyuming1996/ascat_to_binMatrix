# Heatmap of CNV ( or CNA, copy number alterations )
# with ComplexHeatmap package

#library('gplots')
library('ComplexHeatmap')
#library('circlize')
#setwd('/Users/ysh/tmp/heat_map_test_2')
#setwd('D:/xiongys/3_project/20180804_??ѵ??ʿ_OncoScan/heat_map_test_2')

mat = read.csv2('all.renamed.ASCAT.segments.txt.gainloss.binMat.txt',
                sep='\t', header=TRUE, row.names=1)

# for black and white color band of chromosome at top of heatmap.
chromosome.anno_df = data.frame(
    chromosome = t(read.csv2('chromosome.colorband.txt',
    header = FALSE, sep = '\t')))

# to indicate chr Number text.
chromosome.name.loci = read.csv2('chromosome.nameprobe.txt', sep='\t', header= FALSE)
chromosome.name.loci = as.numeric(chromosome.name.loci)

chromosome.name.link = column_anno_link(at = chromosome.name.loci,
    labels = as.character(1:22), side = 'top', labels_gp = gpar(fontsize=7),
    link_width = unit(1, 'mm'))

chromosome.ha = HeatmapAnnotation(df = chromosome.anno_df,
    col = list(chromosome = c("odd"="black","even"="gray80")),
    show_legend = FALSE, gap = unit(2,'mm'), link = chromosome.name.link)

# Sample annotation at right side, with color bar.
sample.anno_df = data.frame(group = c(rep('N_tumor', 24),rep('N_normal', 24),
    rep('R_tumor', 24), rep('R_normal', 24)))
sample.ha = rowAnnotation(df = sample.anno_df,
    col = list(group = c('N_tumor' = 'darkorchid1', 'N_normal' = 'yellow',
                         'R_tumor' = 'chocolate', 'R_normal' = 'chartreuse')),
    gp = gpar(col = 'black'),  # for the color of grid border.
    show_legend = TRUE)

# show stack bar at right side of heatmap
#bar.matrix = as.matrix(read.csv('test', header=TRUE, sep = '\t'))
#sample.stackbar.ha = rowAnnotation(stacked_bar = row_anno_barplot(bar.matrix,
#    border = FALSE, gp = gpar(fill = c('red', 'blue'))),
#    width = unit(1, 'cm'))

col = c('gain' = 'red', 'loss' = 'blue')

pdf('heatmap.pdf', width = 10, height = 5)
Heatmap(mat,
        name = "CNA",
        col = col,   # color
        na_col = 'gray95',
        cluster_rows = FALSE,
        show_row_name=FALSE,
        #row_dend_side = 'left',
        #show_row_dend = FALSE,
        #row_names_side = 'right',
#        split = c(rep('a',24), rep('b', 24), rep('c', 24), rep('d', 24)),
        row_title_gp=gpar(col = c('white')),
        
        
        cluster_columns=FALSE,
        show_column_names = FALSE,
        top_annotation = chromosome.ha,
        top_annotation_height = unit(2, 'mm')
        #heatmap_legend_param = list(title = 'CNA', at=c(2,1,-2,-1,0), labels=c('somatic gain', 'germline gain', 'somatic loss', 'germline loss', 'copy neutral'))
        )
# +
#sample.ha +
#sample.stackbar.ha

dev.off()
