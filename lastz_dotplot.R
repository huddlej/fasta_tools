
args <- commandArgs(trailingOnly=TRUE)

input.file <- args[1]
output.file <- args[2]

data <- read.table(input.file, header=TRUE)
pdf(output.file, width=6, height=6)
plot(data, type="l")
dev.off()
