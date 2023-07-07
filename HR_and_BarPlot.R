##All code for Figure 1

##Calculating Cox Proportional Hazard Ratios (HR)  
#INPUT FILE should include the survival time (hrs), status (1=dead, 0=alive), and the experimental groups
  #Note- by default the first group listed will be treated as the control group 
  #example of input file: 
  #70.5	1	no_preconditioning
  #100.5	1	no_preconditioning
  #100.5	1	no_preconditioning
  # ...
  #76.5	1	preconditioning 
  #92.5	1	preconditioning 
  #94.5	1	preconditioning 
  # ...

#Code to obtain HR
library(survival)
mydata<- read.delim("INPUT FILE NAME", header=FALSE)
names(mydata)=c("time", "status", "x")
mydata$x<-as.character(mydata$x)
mydata$x <- factor(mydata$x, levels=unique(mydata$x))
cox<-coxph(Surv(time, status)~x, data=mydata)
summary(cox)
#In the output, the HR is the exp(coef) value and the p-value is Pr(>|z|) value 

#Plot the ln(HR) of all DGRP tested as a boxplot 
#Note the natural log was taken to normalize the HR data 
#Example of INPUT FILE2, ln(HR) listed in ascending orde: 
  #DGRP  ln(HR)
  #26  0.417920549
  #38  0.301880792
  # ...

require(ggplot2)
data<- read.delim("INPUT FILE2", header=TRUE)
data$DGRP.line <-as.character(data$DGRP.line)
data$DGRP.line <- factor(data$DGRP.line, levels=unique(data$DGRP.line))
ggplot(data, aes(x = DGRP.line, y = ln.HR.)) + geom_bar(stat = "identity")


RStudio 2022.02.2+485 "Prairie Trillium" Release (8acbd38b0d4ca3c86c570cf4112a8180c48cc6fb, 2022-04-19) for Windows
Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) QtWebEngine/5.12.8 Chrome/69.0.3497.128 Safari/537.36

survival package version 3.3.1
ggplot2 package version 3.3.6
