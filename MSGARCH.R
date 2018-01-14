rm(list = ls())

require('rugarch')
require('RHmm')
require('reshape2')
require('ggplot2')

data = EuStockMarkets[1:500, 4]

for(x in 1 : 5)
{
y= data
Fit = HMMFit(data, nStates=2)
Path = viterbi(Fit, data)
fb = forwardBackward(Fit, data)

plot(Path$states, type='s', main='Implied States', xlab='', ylab='State')

matplot(fb$Gamma, type='l', main='Probabilities', ylab='Probability')

legend(x='topright', c('State1','State2'),  fill=1:2, bty='n')

diff = Path$states[-1] - Path$states[-length(Path$states)]
state2 = which(Path$states %in% unique(Path$states)[1])
state1 = which(Path$states %in% unique(Path$states)[2])
vals = unique(diff) 
state21 = length(which(diff %in% vals[2]))
state12 = length(which(diff %in% vals[3]))
state11 = (length(state1) - length(state12)) / length(state1)
state22 = (length(state2) - length(state21)) / length(state2)
state21 = state21 / length(state2)
state12 = state12 / length(state1)

PMat = matrix(c(state11,state12,state21,state22),ncol=2)

#Determining the most optimal p and q parameters of GARCH for each regime

p1 = 0
q1 = 0
p2 = 0
q2 = 0
for(p in 1 : 2)
{
  for(q in 1 : 2)
  {
    garch.spec.norm = ugarchspec(variance.model=list(garchOrder=c(p,q)),mean.model= list(armaOrder=c(0,0)),distribution.model = "norm")

    garch.fit.norm1 = ugarchfit(spec=garch.spec.norm,data=data[state2],solver.control=list(trace = 1))

    garch.fit.norm2 = ugarchfit(spec=garch.spec.norm,data=data[state1],solver.control=list(trace = 1))
   {
      
    if(!exists('max1'))
    {
      max1 = garch.fit.norm1@fit$LLH
      p1 = p
      q1 = q
    }
    else  
    {
      if(garch.fit.norm1@fit$LLH > max1)
      {
        max1 = garch.fit.norm1@fit$LLH
        p1 = p
        q1 = q
      }
    }
  }
  {
    if(!exists('max2'))
    {
      max2 = garch.fit.norm2@fit$LLH
      p2 = p
      q2 = q
    }
    
    else
    {
      if(garch.fit.norm2@fit$LLH > max1)
      {
        max2 = garch.fit.norm2@fit$LLH
        p2 = p
        q2 = q
      }
    }
  }
  
  }
}

#Applying the most optimal p nd q parameters for GARCH models in both the regimes and building the GARCH models. 
# Building the GARCH model for regime 1 using the specification optimat p and q, i.e.,  p1 and q1 derived above.
garch.spec.norm = ugarchspec(variance.model=list(garchOrder=c(p1,q1)),mean.model= list(armaOrder=c(0,0)),distribution.model = "norm")
garch.fit.norm1 = ugarchfit(spec=garch.spec.norm,data=data[state2],solver.control=list(trace = 1))

# Building the GARCH model for regime 2 using the specification optimat p and q, i.e.,  p2 and q2 derived above.
garch.spec.norm2 = ugarchspec(variance.model=list(garchOrder=c(p2,q2)),mean.model= list(armaOrder=c(0,0)),distribution.model = "norm")
garch.fit.norm2 = ugarchfit(spec=garch.spec.norm2,data=data[state1],solver.control=list(trace = 1))

statecurrent = Path$states[length(Path$states)]

coef1 = garch.fit.norm1@fit$coef
coef2 = garch.fit.norm2@fit$coef

m11 = c()
m12 =c()
m21 = c()
m22 =c()
for( p in 1 : p1)
{
  m11[p] =(garch.fit.norm1@fit$residuals[length(garch.fit.norm1@fit$residuals) - (p - 1)]) ^ 2
}
for( p in 1 : q1)
{
  m12[p] =(garch.fit.norm1@fit$var[length(garch.fit.norm1@fit$var) - (p - 1)])
}
for( p in 1 : p2)
{
  m21[p] =(garch.fit.norm2@fit$residuals[length(garch.fit.norm2@fit$residuals) - (p - 1)]) ^ 2
}
for( p in 1 : q2)
{
  m22[p] =(garch.fit.norm2@fit$var[length(garch.fit.norm2@fit$var) - (p - 1)])
}

m1 = c(1,1,m11,m12)
m2 = c(1,1,m21,m22)

futureval = ((PMat[statecurrent,1] * sum(coef1 * m1)) + (PMat[statecurrent,2] * sum(coef2 * m2))) ^ 0.5 + mean(data)

data = c(data,futureval)

}

##########################################################################

RMSE = mean(((data[500:505] - EuStockMarkets[500:505,4])^2))^0.5
RMSE
xaxis = c(1:505,501:505)
yaxis = c(data,EuStockMarkets[501:505,4])
col = as.character(c(rep(1,500),rep(2,5),rep(1,5)))

df = data.frame(X = xaxis, Y = yaxis, Z = col)

plot = ggplot(df,aes(x = X, y = Y, col = Z,shape = Z)) + geom_point(size= 3)
plot = plot + scale_shape_manual(values = c("1" = 1, "2" = 2),labels = c("1" = "Actual","2" = "Predicted"))
plot = plot + scale_color_manual(values = c("1" = "red", "2" = "blue"),labels = c("1" = "Actual","2" = "Predicted"))
plot = plot + annotate("text", label = paste("RMSE:",round(RMSE,2)), x = 200, y = 2800) + theme(legend.position = "bottom",legend.title = element_blank())
plot = plot + xlab("Days from Jan 2009") + ylab("FTSE value") 
plot 



residueactual = EuStockMarkets[501:505] - EuStockMarkets[500:504]
residuefcst = data[501:505] - data[500:504]

RMSE = (mean((residueactual - residuefcst)^2)) ^ 0.5
xaxis = c(1:504,500:504)
yaxis = c((data[-1] - data[-length(data)]),(EuStockMarkets[501:505,4] - EuStockMarkets[500:504,4]))
col = as.character(c(rep(1,504),rep(2,5)))
df = data.frame(X = xaxis, Y = yaxis, Z = col)

plot = ggplot(df,aes(x = X, y = Y, col = Z,shape = Z)) + geom_point(size =3)
plot = plot + scale_shape_manual(values = c("1" = 1, "2" = 2),labels = c("1" = "Actual","2" = "Predicted"))
plot = plot + scale_color_manual(values = c("1" = "red", "2" = "blue"),labels = c("1" = "Actual","2" = "Predicted"))
plot = plot + annotate("text", label = paste("RMSE:",round(RMSE,2)), x = 200, y = 100) + theme(legend.position = "bottom",legend.title = element_blank())
plot = plot + xlab("Days from Jan 2009") + ylab("Residual values") 
plot
