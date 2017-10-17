herv.expr.overlap <- findOverlaps(herv.ranges, expr.ranges)

herv.meth.overlap <- findOverlaps(herv.ranges, meth.ranges)


expr.herv.ids <- expr.ranges$ids[subjectHits(herv.expr.overlap)]
expr.herv.data <- expr.data[expr.herv.ids,]

expr.mean <- apply(expr.data, 1, mean)
expr.sd <- apply(expr.data, 1, sd)

expr.herv.mean <- apply(expr.herv.data, 1, mean)
expr.herv.sd <- apply(expr.herv.data, 1, sd)

pdf('/media/data/Masterarbeit/Plots/expr.var.pdf', width = 11, height = 11)
layout(matrix(c(1:4), 2, 2, byrow = F))
plot(expr.mean, expr.sd/expr.mean * 100, pch = 3, xlab = 'mean expression', ylab = 'coefficent of variation', main = 'Expression per probe')
hist(expr.sd/expr.mean * 100, breaks = seq(0, 35, 0.25), xlab = 'coefficent of variation', main = '')

plot(expr.herv.mean, expr.herv.sd/expr.herv.mean * 100, pch = 3, xlab = 'mean expression', ylab = 'coefficent of variation', main = 'Expression per probe overlapping with herv elements (191)')
hist(expr.herv.sd/expr.herv.mean * 100, breaks = seq(1, 10, 0.25), xlab = 'coefficent of variation', main = '')
dev.off()

meth.herv.ids <- intersect(names(meth.ranges)[subjectHits(herv.meth.overlap)], rownames(meth.data))
meth.herv.data <- meth.data[meth.herv.ids,]

meth.mean <- apply(meth.data, 1, mean)
meth.sd <- apply(meth.data, 1, sd)

meth.herv.mean <- apply(meth.herv.data, 1, mean, na.rm = T)
meth.herv.sd <- apply(meth.herv.data, 1, sd, na.rm = T)

pdf('/media/data/Masterarbeit/Plots/meth.var.pdf', width = 11, height = 11)
layout(matrix(c(1:6), 2, 3, byrow = F))
plot(meth.mean, meth.var/meth.mean * 100, pch = 3, xlab = 'mean methylation', ylab = 'coefficent of variation', main = 'methylation per probe')
hist(meth.var/meth.mean * 100, breaks = seq(0, 35, 0.25), xlab = 'coefficent of variation', main = '')

plot(meth.herv.mean, meth.herv.sd/meth.herv.mean * 100, pch = 3, xlab = 'mean methylation', ylab = 'coefficent of variation', main = 'methylation per probe overlapping with herv elements (1595)')
hist(meth.herv.sd/meth.herv.mean * 100, pch = 3, xlab = 'coefficient of variation', main = '')
plot(meth.herv.mean, meth.herv.sd, pch = 3, xlab = 'mean methylation', ylab = 'standard deviation', main = 'methylation per probe overlapping with herv elements (1595)')
dev.off()