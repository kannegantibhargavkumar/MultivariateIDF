## ternary diagram for return periods
sf_rp=uni_fit[["sf"]][["RP"]]
rf_rp=uni_fit[["rf"]][["RP"]]
surge_rp=uni_fit[["surge"]][["RP"]]
sm_rp=uni_fit[["sm"]][["RP"]]
rp=cbind(sf_rp,rf_rp,surge_rp,sm_rp)
rp=data.frame(rp)
colnames(rp)=c("sf","rf","surge","sm")

# normalize the series
normalize=function(data){
  data_n=(data-min(data))/(max(data)-min(data))
}
# normalizing the data
rp_n=matrix(data=NA, nrow=nrow(rp), ncol=(ncol(rp))) # frist column is year
for(i in 1:(ncol(rp))){
  
  rp_n[,i]=normalize(as.matrix(rp[,(i)]))
}
colnames(rp_n)={c('sf','rf','surge','sm')}

# checking whether dependence is changed bcoz of normalization
s_rp=rowSums(rp_n[,2:4])
rp_in=rp_n[,2:4] # SF is dependent variable
rp_ternary=matrix(data=NA, nrow=nrow(rp_in), ncol=(ncol(rp_in)))
for (j in 1:nrow(rp_in)){
  s_temp=s_rp[j]
  rp_ternary[j,]=as.matrix(rp_in[j,])/s_temp
  
}

## ternary diagram
## input data for ternany plot should be dataframe
rp_ternary=as.data.frame(rp_ternary)
colnames(data_ternary)=c('rf','surge','sm')

### ternary plot
fig_rp=plot_ly(data=rp_ternary, a=rp_ternary[,1],b=rp_ternary[,2],c=rp_ternary[,3],
              type="scatterternary",
              mode="markers",marker=list(size=15,color=rp$sf, colorscale="Blues", showscale=TRUE, reversescale=TRUE),alpha = 1.0) 
fig_rp=fig_rp%>%  layout(title="",showlegend=F,
                       xaxis=list(title="", showgrid=F, zeroline=F, showticklabels=F),
                       yaxis=list(title="", showgrid=F, zeroline=F, showticklabels=F),
                       sum=1,
                       ternary=list(
                         aaxis=list(title="Rainfall", tickformat=".0", tickfont=list(size=15)),
                         baxis = list(title = "Surge", tickformat = ".0", tickfont = list(size = 15)),
                         caxis = list(title = "Soil Moisture", tickformat = ".0", tickfont = list(size = 15))))
#annotations = list(
# list(xref = "paper", yref = "paper", align = "center",
#      x = 0.1, y = 1, text = "", ax = 0, ay = 0,
#      font = list(family = "serif", size = 20, color = "white"),
#      bgcolor = "#b3b3b3", bordercolor = "black", borderwidth = 2)))

fig_rp
