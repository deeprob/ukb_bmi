
args = commandArgs(trailingOnly=TRUE)
# check to see that one argument is given
if (length(args)!=2) {
  stop("contingency table filename and output filename must be given", call.=FALSE)
}

contingency_file = args[1]
out_file = args[2]

# get the table and store as a df
tabla = read.csv(contingency_file, header=TRUE, row.names=1)

# run fisher's exact
test <- fisher.test(tabla) # alternative="greater"
oddsratio = unname(test$estimate)
pvalue = test$p.value
conf_int_lower = test$conf.int[1]
conf_int_upper = test$conf.int[2]

results = cbind(
    oddsratio, 
    pvalue, 
    conf_int_lower,
    conf_int_upper
    )

colnames(results) = c("oddsratio", "pvalue", "conf_int_lower", "conf_int_upper")

# save to file .. 
write.table(results, file=out_file, sep=",", row.names=FALSE, col.names=TRUE)
