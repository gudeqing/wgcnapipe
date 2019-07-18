# coding=utf-8
import subprocess
import argparse
import pandas as pd

"""
1. Correlate the module eigengenes with the trait
2. Correlate each gene with the trait
"""
# get arguments
parser = argparse.ArgumentParser(description="wgcna step3: relate module/gene to traits")
parser.add_argument('-datExpr', type=str, metavar="exp_matrix_file", required=True,
                    help="expression matrix file, for gene/traits relationship calculation")
parser.add_argument('-MEs', type=str, metavar="module eigengenes", required=True,
                    help="module eigengenes, for module/traits relationship calculation")
parser.add_argument('-traits', type=str, required=True, metavar="phenotype_data",
                    help="sample name in row, traits data in column; "
                         "或者为样本分组信息文件,第一行是header，第一列为样本id, 第二列为样本分组, "
                         "样本信息必须和输入的表达矩阵完全匹配,不多不少")
parser.add_argument('-corType', type=str, metavar="correlation_type", default="pearson",
                    help="correlation type, 'pearson' or 'spearman', 'kendall'")
parser.add_argument('-nThreads', type=int, default=16, )
args = parser.parse_args()


# read expr
r_cmds = \
"""
# module and traits relation
library('WGCNA')
enableWGCNAThreads()
datME = read.table('{eigengenes}', header=T, row.names=1)
datME = as.data.frame(t(datME))
MEs = orderMEs(datME)
traits = read.table("{traits}", header=T, row.names=1)
if (dim(traits)[2] == 1 & class(traits[1,1])=="factor"){bracket1}
    tmp = model.matrix(~0+ traits[,1])
    colnames(tmp) = levels(traits[,1])
    traits = tmp
{bracket2}
traits = as.data.frame(traits)
correlation = signif(cor(MEs, traits, use="p", method="{cor_type}", nThreads={threads}), 3)
pvalues = signif(corPvalueStudent(correlation, nSamples = dim(traits)[1]), 3)

pdf(file='Module-Trait.pdf', width = 12, height = 9)
textMatrix = paste(signif(correlation, 2), "\n(", signif(pvalues, 1), ")", sep = "")
dim(textMatrix) = dim(correlation)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = correlation, xLabels = names(traits), yLabels = names(MEs), 
    ySymbols = names(MEs), colorLabels = FALSE, colors = greenWhiteRed(50), 
    textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, 
    zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()

write.table(correlation, 'module_trait.correlation.xls', col.names=NA, quote=F, sep='\\t', row.names=T)
write.table(pvalues, 'module_trait.correlation_pvalues.xls', col.names=NA, quote=F, sep='\\t', row.names=T)
# gene and traits relation
exp = read.table('{exp_matrix}', header=T, row.names=1)
correlation2 = signif(cor(t(exp), traits, use="p", method="{cor_type}"), 3)
write.table(correlation2, 'gene_trait.correlation.xls', col.names=NA, quote=F, sep='\\t', row.names=T)
""".format(
    eigengenes=args.MEs,
    traits=args.traits,
    cor_type=args.corType,
    threads=args.nThreads,
    exp_matrix=args.datExpr,
    bracket1="{",
    bracket2="}"
)

with open('wgcna_relate_analysis.r', 'w') as f:
    f.write(r_cmds)
subprocess.check_call("Rscript wgcna_relate_analysis.r", shell=True)
