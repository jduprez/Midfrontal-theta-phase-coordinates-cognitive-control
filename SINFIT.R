## SINFIT

## This script fits sine waves to power-RT correlation according to phase bins data obtained from the first part of the data analysis (using the midfrontal_theta_phase.m script)
## It uses non linear least square fitting (nls, {minpack} package) to estimate amplitude, frequency, and phase parameters of the fitted sine wave
## These parameters are then entered in a lm model combining a sine and a cosine function to evaluate the goodness of the model


# Load and reshape data--------------------------------------------------------
# Loads and reshapes data in 5 dataframes (1 by condition)

data = read.table("rtcorrphase_R_13-2_300-1200.txt", header = TRUE)

data1 = droplevels(data[which(data$condition == 1),])
data2 = droplevels(data[which(data$condition == 2),])
data3 = droplevels(data[which(data$condition == 3),])
data4 = droplevels(data[which(data$condition == 4),])
data5 = droplevels(data[which(data$condition == 5),])

data1 = data.frame(data1[,-1], row.names=data1[,1])
data1 = data1[,-21]
data1 = as.matrix(data1)

data2 = data.frame(data2[,-1], row.names=data2[,1])
data2 = data2[,-21]
data2 = as.matrix(data2)

data3 = data.frame(data3[,-1], row.names=data3[,1])
data3 = data3[,-21]
data3 = as.matrix(data3)

data4 = data.frame(data4[,-1], row.names=data4[,1])
data4 = data4[,-21]
data4 = as.matrix(data4)

data5 = data.frame(data5[,-1], row.names=data5[,1])
data5 = data5[,-21]
data5 = as.matrix(data5)

# Now we reshape all the dataframes as one 3D array

dat = list(data1, data2, data3, data4, data5)
dat = simplify2array(dat) # Gives a  n (subject) * 20 (values) * 5 (conditions) array



# Initialize output arrays ----------------------------------------------

# For the nlsLM coefficents

coeff <- array(rep(NA, max(data[,1])*12*5), dim=c(max(data[,1]),12,5)) # Creates a n (subjects)*12(coefficients)*5(conditions) array
# 16 Values :"amp","amp_std","amp_t", "amp_pval", "freq", "freq_std", 
# "freq_t", "freq_pval", "phase", "phase_std", "phase_t", "phase_pval"

# For lm fit coefficients
lmcoeff = array(rep(NA, max(data[,1])*10*5), dim=c(max(data[,1]),10,5)) # Creates a n (subjects)*10(coefficients)*5(conditions) array
# 10 Values : 'xs','xs_std','xs_t','xs_pval','xc','xc_std','xc_t','xc_pval',
# 'modelF','modelpval'

# For Rsquared and adjusted Rsquared of the lm fit

lmRsq = array(rep(NA,max(data[,1]*2*5)), dim=c(max(data[,1]),2,5))# Creates a n (subjects)*2(coefficients)*5(conditions) array
# 2 Values : 'Rsquared','adjusted Rsquared'


# For fitted values
fittedval = array(rep(NA, max(data[,1])*20*5), dim=c(max(data[,1]),20,5)) # n (subjects)*20(fitted values)*5(conditions)


# Create phase vector (which is here a simple 20 elements vector from 0 to 1)
phasebin = seq(0,1,length.out=20)



# Model and Fit -----------------------------------------------------------

require(minpack.lm) # needed for nls
require(OpenMx) # needed for rvectorize
require(broom) # needed for glance(model)$p.value

# Loop through conditions

for (condi in 1:5){
  
  # Loop through subjects
  
  for (subno in 1:28){
    
    tempdat = dat[subno,,condi] # temporary vector containing data from subno and condi
    
    fitsine<-nlsLM(tempdat~A*sin(2*pi*W*phasebin+phi), start=list(A = -0.2, W = 1, phi =-pi), upper =c(1,2,pi)) # nonlinear least squares model. finds estimated sine parameters coefficients. Need to specify a starting point for each parameter
    
    coeff[subno, ,condi] = t(rvectorize(summary(fitsine)$coefficients[1:3,1:4])) # store coefficients for group level analysis
    
    xs<-(summary(fitsine)$coefficients[1,1])*sin(2*pi*(summary(fitsine)$coefficients[2,1])*phasebin+(summary(fitsine)$coefficients[3,1])) # sin part for the linear model (A*sin*[2*pi*freq]*phasebin+phase)
    xc<-(summary(fitsine)$coefficients[1,1])*cos(2*pi*(summary(fitsine)$coefficients[2,1])*phasebin+(summary(fitsine)$coefficients[3,1])) # cos part for the linear model (A*cos*[2*pi*freq]*phasebin+phase)
    
    fitsine2<-lm(tempdat ~ xs + xc) # linear model to test the overall fit of the parameters to the data
    
    lmcoeff[subno,,condi] = c(c(t(rvectorize(summary(fitsine2)$coefficients[2:3,1:4]))),c(glance(fitsine2)$statistic),c(glance(fitsine2)$p.value)) # store coefficients of the lm fit
    
    lmRsq[subno,,condi] = c(summary(fitsine2)$r.squared, summary(fitsine2)$adj.r.squared) # stores Rsquared and adjusted Rsquared of the lm model
    
    fittedval[subno,,condi] = fitsine2$fitted.values # store fitted values for plotting
    
  } # end subject loop
  
} # end condition loop

# Plots --------------------------------------------------------------------


# Plot for one condition (change x in dat[subno,,x] to plot condition 1 or 2 or...or 5)

# Add text using the following Corner_text function:
Corner_text <- function(text, location="topright"){
  par(font=2)
  legend(location,legend=text, bty ="n", pch=NA,cex=1.1) 
}

  
col_plot= c("#0072BD",  "#ED1220","#D95319", "#7E2F8E", "#77AC30")# same colors as the matlab plots


phase= seq(-pi,pi,length.out=20)

par(mfrow=c(3,3),
        oma = c(5,4,0,0) + 0.05,
        mar = c(0,0,0.5,0.5) + 0.05)

lablist=c(expression(-pi), expression(-pi/2), 0, expression(pi/2), expression(pi))

for (subno in 1:1){ # Plot 9 sub max at a time
  if ( (subno == 26) | (subno ==27) | (subno ==28)) {
  plot(dat[subno,,5]~phase, xlab="Phase (rad)", yaxt="n", pch =16, col = "dimgrey", ylab='', xaxt = "n")
  axis(1, seq(-pi,pi,length.out=5), labels=lablist)
  par(new=TRUE)
  plot(fittedval[subno,,5], type = "l", lwd = 2, col = col_plot[5], ann=FALSE, xaxt="n", yaxt="n")
  #Corner_text(paste(subno))
  }else{
  plot(dat[subno,,5]~phase, xaxt="n", yaxt="n", pch =16, col = "dimgrey", ylab='', xlab='')
  par(new=TRUE)
  plot(fittedval[subno,,5], type = "l", lwd = 2, col = col_plot[5], ann=FALSE, xaxt="n", yaxt="n", ylab='', xlab='')  
  }
  box(lwd=0.1, col ="grey")
}

# Use 400 247 for image size !!

# Significant amplitude estimated coefficients - subject level

conditions = c('no_conflict','little_conflict','lot_of_conflict','mixed_correct','error') # Creates a character vector with the conditions name

# Loop through conditions

for (condi in 1:5){

print(conditions[condi]) # prints the condition
print("NLS amplitude non-significant")
print(which(coeff[,4,condi]>0.05, arr.ind = FALSE, useNames = TRUE)) # returns the index of NS coefficients (p>0.05). 4th coefficients : amplitude p-value
print("LM sine part non-significant")
print(which(lmcoeff[,10,condi]>0.05, arr.ind = FALSE, useNames = TRUE)) # returns the index of NS lm model (p>0.05). 10th coefficients : overall model p-value
  
} # end of loop




# Group level statistics --------------------------------------------------

# all pairwise non parametric comparisons of Rsquared between conditions

for (condi in 1:5){
  
  if (condi <=5){
  
    for (condi2 in 1:5)
      
      if (condi2 != condi){
  pairtest = wilcox.test(lmRsq[,1,condi],lmRsq[,1,condi2], paired=TRUE)
  print(c(condi, 'vs', condi2))
  print(pairtest$p.value)
  
      } else {}
  
  } else {}
  
}

# Plot Rsquared as a function of condition

require(plyr) # for the join function (joins dataframe)
require(plotrix) # for the std.error function
require(ggplot2) # for plots


rsqlm_mean1 = data.frame(n=c(1:max(data[,1])), rsq=lmRsq[,1,1], condition = rep(1,max(data[,1])))
rsqlm_mean2 = data.frame(n=c(1:max(data[,1])), rsq=lmRsq[,1,2], condition = rep(2,max(data[,1])))
rsqlm_mean3 = data.frame(n=c(1:max(data[,1])), rsq=lmRsq[,1,3], condition = rep(3,max(data[,1])))
rsqlm_mean4 = data.frame(n=c(1:max(data[,1])), rsq=lmRsq[,1,4], condition = rep(4,max(data[,1])))
rsqlm_mean5 = data.frame(n=c(1:max(data[,1])), rsq=lmRsq[,1,5], condition = rep(5,max(data[,1])))

rsqlm_mean = join(rsqlm_mean1, rsqlm_mean2, type='full')
rsqlm_mean = join(rsqlm_mean, rsqlm_mean3, type='full')
rsqlm_mean = join(rsqlm_mean, rsqlm_mean4, type='full')
rsqlm_mean = join(rsqlm_mean, rsqlm_mean5, type='full')


mean_rsq = with(rsqlm_mean, aggregate(rsqlm_mean$rsq, list(rsqlm_mean$condition), mean))
stderr_rsq = with(rsqlm_mean, aggregate(rsqlm_mean$rsq, list(rsqlm_mean$condition), std.error))

gg_rsq = cbind(mean_rsq,stderr_rsq)
gg_rsq = gg_rsq[,-3]
colnames(gg_rsq) = c('condition','avg_rsq','std')

# Plot amplitude according to conditions 

rsq_plot = qplot(x = condition, y = avg_rsq, data = gg_rsq, geom = "point")
rsq_plot + geom_point(size=3)+geom_errorbar(aes(ymin=avg_rsq-std, ymax=avg_rsq+std), width=0.1)

# Linear mixed model on amplitude coefficients 

require(nlme)
require(plyr) # for the join function (joins dataframe)
require(plotrix) # for the std.error function
require(ggplot2) # for plots


amp1 = data.frame(n=c(1:max(data[,1])), amp=abs(coeff[,1,1]), condition = rep(1,max(data[,1])))
amp2 = data.frame(n=c(1:max(data[,1])), amp=abs(coeff[,1,2]), condition = rep(2,max(data[,1])))
amp3 = data.frame(n=c(1:max(data[,1])), amp=abs(coeff[,1,3]), condition = rep(3,max(data[,1])))
amp4 = data.frame(n=c(1:max(data[,1])), amp=abs(coeff[,1,4]), condition = rep(4,max(data[,1])))
amp5 = data.frame(n=c(1:max(data[,1])), amp=abs(coeff[,1,5]), condition = rep(5,max(data[,1])))

ampgp = join(amp1, amp2, type='full')
ampgp = join(ampgp, amp3, type='full')
ampgp = join(ampgp, amp4, type='full')
ampgp = join(ampgp, amp5, type='full')

mean_amp = with(ampgp, aggregate(ampgp$amp, list(ampgp$condition), mean))
stderr_amp = with(ampgp, aggregate(ampgp$amp, list(ampgp$condition), std.error))

gg_amp = cbind(mean_amp,stderr_amp)
gg_amp = gg_amp[,-3]
colnames(gg_amp) = c('condition','avg_amp','std')


# all pairwise non parametric comparisons of amplitude between conditions

for (condi in 1:5){
  
  if (condi <=5){
    
    for (condi2 in 1:5)
      
      if (condi2 != condi){
        pairtest = wilcox.test(abs(coeff[,1,condi]),abs(coeff[,1,condi2]), paired=TRUE)
        print(c(condi, 'vs', condi2))
        print(pairtest$p.value)
        
      } else {}
    
  } else {}
  
}


# Plot amplitude according to conditions 

amp_plot = qplot(x = condition, y = avg_amp, data = gg_amp, geom = "point")
amp_plot + geom_point(size=3)+geom_errorbar(aes(ymin=avg_amp-std, ymax=avg_amp+std), width=0.1)



# Determine and export max and min points phases for CFC analyses -------

# Initialize 3D array to export

maxminph = array(rep(NA, max(data[,1])*2*5), dim=c(max(data[,1]),2,5))

# Export empirical phase of max correlation

for (condi in 1:5){
  
  # Loop through subjects
  
  for (subno in 1:28){
    
    maxminph[subno,1,condi] = which.min(dat[subno,,condi]) # gets the index of the min fitted value
    maxminph[subno,2,condi] = which.max(dat[subno,,condi]) # gets the index of the max fitted value
    
  }
  
  file = paste("emp_maxminph13_300-1200_",condi,".txt", sep="") # filename
  write.table(maxminph[,,condi], file, append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE) # export one txt file per condition
}


## END

