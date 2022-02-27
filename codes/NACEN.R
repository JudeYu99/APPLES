#! /usr/bin/env Rscript

# Rscript NACEN.R -d /opt/anaconda3/bin/mkdssp -p protein.pdb

library(argparse)

# 创建参数解析对象
parser <- ArgumentParser()

parser$add_argument("-d", "--dssp_file", type = "character", default = "/opt/anaconda3/bin/mkdssp", help = "The directory of DSSP executive file.")
parser$add_argument("-p", "--pdb_path", type = "character", help = "The directory of protein complex pdb files.")

args <- parser$parse_args()

library(NACEN)

dssp_file <- args$dssp_file
pdb_path <- args$pdb_path

Net <- NACENConstructor(pdb_path, WeightType = "none", exefile = dssp_file)
print(" *** Network constructed! *** ")

Net_analyzed <- NACENAnalyzer(Net$AM, Net$NodeWeights)
print(" *** Network analyzed! *** ")

NetP <- data.frame(Net_analyzed$NetP)
output <- data.frame(cbind(as.character(NetP$ID), as.character(NetP$B), as.character(NetP$chain), as.character(NetP$Resid), as.character(NetP$Res)))
colnames(output) <- c("ID", "Betweenness", "Chain", "ResID", "ResName")
write.table(output, "../APPLES_OUT/NACEN_output.csv", quote = F, row.names = F, sep = "\t")
