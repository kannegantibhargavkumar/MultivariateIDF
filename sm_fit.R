### 
#distribution to soil moisture
library(readxl)
library(extraDistr)
## soil moisture at warunji
data=read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/data_input/sm_warunji.xlsx")
da=na.omit(data)
hist(as.matrix(da[,2]))
plot(density(as.matrix(da[,2])))


## soil moisture at c bridge
data_s=read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Uk_work/Data/Data_at_cbridge/soil_moisture_cbridge.xlsx")
data_ss=na.omit(data_s)
plot(density(as.matrix(data_ss[,2])))



rf_data=read_excel("F:/Research Work/Synchronization work/Multivariate IDF/Krishna basin work/data_input/rf_warunji.xlsx")
rf=rf_data[rf_data[,4]!=0,]
hist(as.matrix(rf[,4]))


x <- rkumar(1e5, 5, 16)
hist(x, 100, freq = FALSE)
curve(dkumar(x, 5, 16), 0, 1, col = "red", add = TRUE)
hist(pkumar(x, 5, 16))
plot(ecdf(x))
curve(pkumar(x, 5, 16), 0, 1, col = "red", lwd = 2, add = TRUE)
