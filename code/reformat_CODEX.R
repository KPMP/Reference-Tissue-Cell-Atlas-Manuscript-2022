#Reformat CODEX data for input to VorteX

library(data.table)
a <- fread("CODEX_76_reg3.csv", sep=",")
library(plyr)
library(dplyr)
a <- as.data.frame(a)[, !duplicated(colnames(a))]
a <- filter(a, a$QC2_reg3 == 1)
a <- select(a, -matches("DAPI"))
a <- select(a, -matches("Empty"))
a <- select(a, -matches("Blank"))
a_expr <- a[,9:39]
a_expr_scaled <- scale(a_expr)
a_expr_scaled <- cbind(a$cell_id, a_expr_scaled)
write.table(a_expr_scaled, "a_input_with_size.csv",  sep=",", quote=FALSE, row.names=FALSE)