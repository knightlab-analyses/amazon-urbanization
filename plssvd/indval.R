library(labdsv)

set.seed(42)

metabolites <- read.csv('metabolites.csv', row.names=1, header=TRUE)
microbes <- read.csv('microbes.csv', row.names=1, header=TRUE)
mapping <- read.csv('mapping.csv')

metabolite.iv <- indval(metabolites, as.numeric(mapping$category), numitr=1000)
microbe.iv <- indval(microbes, as.numeric(mapping$category), numitr=1000)

gr <- microbe.iv$maxcls
iv <- microbe.iv$indcls
pv <- microbe.iv$pval

microbe.summary <- data.frame(group=gr, indval=iv, pval=pv)
write.table(microbe.summary, 'microbe.indicator.values.csv', sep=',')

gr <- metabolite.iv$maxcls
iv <- metabolite.iv$indcls
pv <- metabolite.iv$pval

metabolite.summary <- data.frame(group=gr, indval=iv, pval=pv)
write.table(metabolite.summary, 'metabolite.indicator.values.csv', sep=',')

