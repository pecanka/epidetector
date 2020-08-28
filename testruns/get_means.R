source("d:/Dropbox/Projects/R/general_functions.R")

files = list.files(pattern=".epi$")
f = file_sort(files)[1]

cat("Analyzing file ",f," ...\n")

d = read.table(f, header=TRUE, na="*")

cc = c("AS4","DS4","CS")
cc = cc[cc %in% colnames(d)]

m = apply(d[,cc], 2, mean, na.rm=TRUE)
v = apply(d[,cc], 2, var, na.rm=TRUE)

print(rbind(m,v))

