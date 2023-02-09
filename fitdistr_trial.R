library("fitdistrplus")
library("VGAM")
k=data("groundbeef")
str(groundbeef)
plotdist(groundbeef$serving,histo=TRUE, demp = TRUE)
plotdist(data_ss$sm, histo = TRUE, demp=TRUE)
descdist(data_ss$sm, boot=1000)

fb=fitdist(data_ss$sm, "kumar")
summary(fb)

fg=fitdist(da$sm, "gamma")
fln=fitdist(da$sm,"lnorm")
par(mffrow=c(2,2))
plot.legend=c("beta","lognormal","gamma")
denscomp(list(fb,fln,fg), legendtext = plot.legend)
qqcomp(list(fb,fln,fg), legendtext = plot.legend)
cdfcomp(list(fb,fln,fg),legendtext = plot.legend)
ppcomp(list(fb,fln,fg),legendtext = plot.legend)
gofstat(list(fb,fg,fln),fitnames=c("beta","gamma","lnorm"))



## kumaraswamy distribution
shape1 <- exp(1); shape2 <- exp(2)
kdata <- data.frame(y = rkumar(n = 1000, shape1, shape2))
fit <- vglm(y ~ 1, kumar, data = kdata, trace = TRUE)
c(with(kdata, mean(y)), head(fitted(fit), 1))
coef(fit, matrix = TRUE)
Coef(fit)
summary(fit)

# distribution 
