### script for contour plot along with scatter
library(plotly)
library(stringr)
library(reshape2)



data.loess <- loess(qsec ~ wt * hp, data = mtcars)

# Create a sequence of incrementally increasing (by 0.3 units) values for both wt and hp
xgrid <-  seq(min(mtcars$wt), max(mtcars$wt), 0.3)
ygrid <-  seq(min(mtcars$hp), max(mtcars$hp), 0.3)

# Generate a dataframe with every possible combination of wt and hp
data.fit <-  expand.grid(wt = xgrid, hp = ygrid)
# Feed the dataframe into the loess model and receive a matrix output with estimates of
# acceleration for each combination of wt and hp
mtrx3d <-  predict(data.loess, newdata = data.fit)
# Abbreviated display of final matrix
mtrx3d[1:4, 1:4]
# Transform data to long form
mtrx.melt <- melt(mtrx3d, id.vars = c('wt', 'hp'), measure.vars = 'qsec')
names(mtrx.melt) <- c('wt', 'hp', 'qsec')
# Return data to numeric form
mtrx.melt$wt <- as.numeric(str_sub(mtrx.melt$wt, str_locate(mtrx.melt$wt, '=')[1,1] + 1))
mtrx.melt$hp <-  as.numeric(str_sub(mtrx.melt$hp, str_locate(mtrx.melt$hp, '=')[1,1] + 1))

fig <- plot_ly(mtrx.melt, x = ~wt, y = ~hp, z = ~qsec, type = "contour",
               width = 600, height = 500)

fig
### 