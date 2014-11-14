# Goodness of fit analysis of all sites of all proteins pooled together by
# performing global fits of empirical and predicted rates
# this is better than averaging protein-specific goodness-of-fit measures 
# because it automatically deals with the issue that different proteins have different number of sites.
# (averaging should use a weighting scheme to deal with this).

# read data
if(T) d = read.csv("foldx_rate_profiles.csv")

# R2 values for single-variable and two-variable models.
fit = lm(nr4s~rate.ddg+rate.lms,data=d)
sfit = summary(fit)
r2.total = sfit$r.squared

fit.ddg = lm(nr4s~rate.ddg,data=d)
sfit.ddg = summary(fit.ddg)
r2.ddg = sfit.ddg$r.squared

fit.lms = lm(nr4s~rate.lms,data=d)
sfit.lms = summary(fit.lms)
r2.lms = sfit.lms$r.squared

# Variance partitioning (absolute contributions to explained variance)
unique.ddg = r2.total - r2.lms
unique.lms = r2.total - r2.ddg
common.ddgORlms=r2.total-unique.ddg-unique.lms

table = data.frame("fit"=c("nr4s~rate.ddg","nr4s~rate.lms","nr4s~rate.ddg+rate.lms","","",""),
          "contribution"=c("","","total","common","unique ddG", "unique MLmS"),
          "R2"=c(r2.ddg,r2.lms,r2.total,common.ddgORlms,unique.ddg,unique.lms))

# Variance partitioning (relative to total explained variance of empirical rates)
table$nR2 = table$R2/r2.total
table$nR2[1:2]=c("","")
table$R2=table$R2

print(table)
write.csv(table,"global-variance-partitioning-table.csv")

# correlations and partial correlations
require(ppcor)
c.test.ddg=cor.test(d$nr4s,d$rate.ddg)
pc.test.ddg=pcor.test(d$nr4s,d$rate.ddg,d$rate.lms)
spc.test.ddg=spcor.test(d$nr4s,d$rate.ddg,d$rate.lms)
c.test.lms=cor.test(d$nr4s,d$rate.lms)
pc.test.lms=pcor.test(d$nr4s,d$rate.lms,d$rate.ddg)
spc.test.lms=spcor.test(d$nr4s,d$rate.lms,d$rate.ddg)

global.correlation.table = data.frame("measure"=c("correlation","partial","semi-partial"),
                                      "ddg.r"=c(c.test.ddg$estimate,pc.test.ddg$estimate,spc.test.ddg$estimate),
                                      "lms.r"=c(c.test.lms$estimate,pc.test.lms$estimate,spc.test.lms$estimate),
                                      "ddg.pvalue"=c(c.test.ddg$p.value,pc.test.ddg$p.value,spc.test.ddg$p.value),
                                      "lms.pvalue"=c(c.test.lms$p.value,pc.test.lms$p.value,spc.test.lms$p.value))
print(global.correlation.table)
write.csv(global.correlation.table,"global-fit-table.csv")

