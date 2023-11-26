#require(compareGroups)
#require(survival)
# load REGICOR data
if (!require("foreach")) install.packages("foreach")
library(foreach)
if (!require("doParallel")) install.packages("doParallel")
library(doParallel)
#require randtoolbox for random number generations
if (!require("randtoolbox")) install.packages("randtoolbox")
library(randtoolbox)
if (!require("Rcpp")) install.packages("Rcpp")
library(Rcpp)
if (!require("Rfast")) install.packages("Rfast")
library(Rfast)

if (!require("matrixStats")) install.packages("matrixStats")
library(matrixStats)
if (!require("parallel")) install.packages("parallel")
library(parallel)

require(tidyverse)
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)




num_cores <- detectCores() - 1  # Use all available cores except one

set.seed(1)
locationdiff1 <- function(data, indices1,indices2) {
  group1_resampled <- data[indices1]
  group2_resampled <- data[indices2]
  tHL(x=group1_resampled) - tHL(x=group2_resampled)  # Example: difference in means
}

bootstrap_parallel_location <- function(group1, group2, R) {
  combined_data <- c(group1, group2)
  locationdiff1 <- function(data, indices1,indices2) {
    group1_resampled <- data[indices1]
    group2_resampled <- data[indices2]
    tHL(x=group1_resampled) - tHL(x=group2_resampled)  # Example: difference in means
  }
  
  mediansorted<-function(sortedx,lengthx){
    if (lengthx%%2==0){
      median1<-(sortedx[lengthx/2]+sortedx[(lengthx/2)+1])/2
    }
    else{median1<-sortedx[(lengthx+1)/2]}
    names(median1)<-NULL
    median1
  }
  tHL<-function (x,sorted=FALSE,bootsize=100000,boot=TRUE,releaseall=FALSE){
    if(sorted){
      sortedx<-x
    }else{
      sortedx<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    }
    lengthx<-length(sortedx)
    set.seed(1)
    if(lengthx>40 || boot){
      subtract<-t(replicate(bootsize, sort(sample(sortedx, size = 2))))
      getlm<-function(vector){ 
        ((vector[1]+vector[2]))/2
      }
      dp2HLx<-apply(subtract,MARGIN=1,FUN=getlm)
      dp2HLx<-Rfast::Sort(x=dp2HLx,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
      
      bootstrappedsample2<-c()
      
    }else{
      combinationsample2<-t(as.data.frame(Rfast::comb_n(sortedx,2)))
      getlm<-function(vector){ 
        ((vector[1]+vector[2]))/2
      }
      dp2HLx<-apply(subtract,MARGIN=1,FUN=getlm)
      dp2HLx<-Rfast::Sort(x=dp2HLx,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
      
      combinationsample2<-c()
      
    }
    
    #HL1<-mean(x=sortedx)
    mean1<-(mean(dp2HLx))
    tHL1<-(mean(dp2HLx,0.05))
    tHL2<-(mean(dp2HLx,0.1))
    tHL3<-(mean(dp2HLx,0.15))
    tHL4<-(mean(dp2HLx,0.2))
    tHL5<-(mean(dp2HLx,0.25))
    tHL6<-(mean(dp2HLx,0.3))
    tHL7<-(mean(dp2HLx,0.35))
    tHL8<-(mean(dp2HLx,0.4))
    tHL9<-(mean(dp2HLx,0.45))
    mHL1<-(mediansorted(sortedx=dp2HLx,lengthx=length(dp2HLx)))
    if(releaseall){
      sdall<-c(sdvar=sd(dp2HLx))
      return(c(sdall,msd1))
    }
    return(c(tHL0=mean1,tHL05=tHL1,tHL10=tHL2,tHL15=tHL3,tHL20=tHL4,tHL25=tHL5,tHL30=tHL6,tHL35=tHL7,tHL40=tHL8,tHL45=tHL9,msd50=mHL1))
  }
  
  indices_list <- replicate(R, list(indices1 = sample(1:length(group1), replace = TRUE),
                                    indices2 = length(group1) + sample(1:length(group2), replace = TRUE)), simplify = FALSE)
  numCores <- detectCores()-4 # Detect the number of available cores
  cl <- makeCluster(numCores) # Create a cluster with the number of cores
  registerDoParallel(cl) # Register the parallel backend
  
  results <- foreach(indices = indices_list, .packages = c("stats")) %dopar% {
    locationdiff1(data = combined_data, indices1 = indices$indices1, indices2 = indices$indices2)
  }
  
  stopCluster(cl)
  registerDoSEQ()
  
  #results <- parLapply(cl = makeCluster(num_cores), X = indices_list, fun = function(indices) {
  #  
  #  locationdiff1(data = combined_data, indices1 = indices$indices1, indices2 = indices$indices2)
  #})
  
  
  return(results)
}

scalediff1 <- function(data, indices1,indices2) {
  group1_resampled <- data[indices1]
  group2_resampled <- data[indices2]
  tsd(group1_resampled) - tsd(group2_resampled)  # Example: difference in means
}
msd<-function (x,sorted=FALSE,bootsize=100000,boot=TRUE,releaseall=FALSE){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  lengthx<-length(sortedx)
  set.seed(1)
  if(lengthx>40 || boot){
    subtract<-t(replicate(bootsize, sort(sample(sortedx, size = 2))))
    getlm<-function(vector){ 
      ((vector[1]-vector[2])^2)/2
    }
    dp2varx<-apply(subtract,MARGIN=1,FUN=getlm)
    dp2varx<-Rfast::Sort(x=dp2varx,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    
    bootstrappedsample2<-c()
    
  }else{
    combinationsample2<-t(as.data.frame(Rfast::comb_n(sortedx,2)))
    getlm<-function(vector){ 
      ((vector[1]-vector[2])^2)/2
    }
    dp2varx<-apply(subtract,MARGIN=1,FUN=getlm)
    dp2varx<-Rfast::Sort(x=dp2varx,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    
    combinationsample2<-c()
    
  }
  
  #HL1<-mean(x=sortedx)
  
  msd1<-sqrt(mediansorted(sortedx=dp2varx,lengthx=length(dp2varx)))
  
  if(releaseall){
    sdall<-c(sdvar=sd(dp2varx))
    return(c(sdall,msd1))
  }
  return(c(msd1))
}
mediansorted<-function(sortedx,lengthx){
  if (lengthx%%2==0){
    median1<-(sortedx[lengthx/2]+sortedx[(lengthx/2)+1])/2
  }
  else{median1<-sortedx[(lengthx+1)/2]}
  names(median1)<-NULL
  median1
}
bootstrap_parallel_scale <- function(group1, group2, R) {
  combined_data <- c(group1, group2)
  
  
  indices_list <- replicate(R, list(indices1 = sample(1:length(group1), replace = TRUE),
                                    indices2 = length(group1) + sample(1:length(group2), replace = TRUE)), simplify = FALSE)
  
  numCores <- detectCores()-4 # Detect the number of available cores
  cl <- makeCluster(numCores) # Create a cluster with the number of cores
  registerDoParallel(cl) # Register the parallel backend
  
  
  scalediff1 <- function(data, indices1,indices2) {
    group1_resampled <- data[indices1]
    group2_resampled <- data[indices2]
    tsd(group1_resampled) - tsd(group2_resampled)  # Example: difference in means
  }
  
  mediansorted<-function(sortedx,lengthx){
    if (lengthx%%2==0){
      median1<-(sortedx[lengthx/2]+sortedx[(lengthx/2)+1])/2
    }
    else{median1<-sortedx[(lengthx+1)/2]}
    names(median1)<-NULL
    median1
  }
  
  tsd<-function (x,sorted=FALSE,bootsize=100000,boot=TRUE,releaseall=FALSE){
    if(sorted){
      sortedx<-x
    }else{
      sortedx<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    }
    lengthx<-length(sortedx)
    set.seed(1)
    if(lengthx>40 || boot){
      subtract<-t(replicate(bootsize, sort(sample(sortedx, size = 2))))
      getlm<-function(vector){ 
        ((vector[1]-vector[2])^2)/2
      }
      dp2varx<-apply(subtract,MARGIN=1,FUN=getlm)
      dp2varx<-Rfast::Sort(x=dp2varx,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
      
      bootstrappedsample2<-c()
      
    }else{
      combinationsample2<-t(as.data.frame(Rfast::comb_n(sortedx,2)))
      getlm<-function(vector){ 
        ((vector[1]-vector[2])^2)/2
      }
      dp2varx<-apply(subtract,MARGIN=1,FUN=getlm)
      dp2varx<-Rfast::Sort(x=dp2varx,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
      
      combinationsample2<-c()
      
    }
    
    #HL1<-mean(x=sortedx)
    sd1<-sqrt(mean(dp2varx))
    tsd1<-sqrt(mean(dp2varx,0.05))
    tsd2<-sqrt(mean(dp2varx,0.1))
    tsd3<-sqrt(mean(dp2varx,0.15))
    tsd4<-sqrt(mean(dp2varx,0.2))
    tsd5<-sqrt(mean(dp2varx,0.25))
    tsd6<-sqrt(mean(dp2varx,0.3))
    tsd7<-sqrt(mean(dp2varx,0.35))
    tsd8<-sqrt(mean(dp2varx,0.4))
    tsd9<-sqrt(mean(dp2varx,0.45))
    msd1<-sqrt(mediansorted(sortedx=dp2varx,lengthx=length(dp2varx)))
    if(releaseall){
      sdall<-c(sdvar=sd(dp2varx))
      return(c(sdall,msd1))
    }
    return(c(tsd0=sd1,tsd05=tsd1,tsd10=tsd2,tsd15=tsd3,tsd20=tsd4,tsd25=tsd5,tsd30=tsd6,tsd35=tsd7,tsd40=tsd8,tsd45=tsd9,msd50=msd1))
  }
  
  
  results <- foreach(indices = indices_list, .packages = c("stats")) %dopar% {
    
    
    scalediff1(data = combined_data, indices1 = indices$indices1, indices2 = indices$indices2)
  }
  
  stopCluster(cl)
  registerDoSEQ()
  
  return(results)
}
MSD_process_compare_location <- function(group1, group2, R) {
  bootstrap_results <- bootstrap_parallel_location(group1, group2, R)
  
  bootstrap <- as.matrix(as.data.frame(bootstrap_results))
  
  # Calculate the observed statistic
  observed_statistic <- locationdiff1(data = c(group1, group2), indices1 = (1:length(c(group1))),indices2 = length(c(group1))+(1:length(c(group2))))
  
  # Calculate the bias
  bias <- rowMeans2(bootstrap) - observed_statistic
  
  # Calculate the standard error
  std_error <- rowSds(bootstrap)
  
  # Calculate the confidence interval
  ci <- apply(bootstrap,1,quantile,probs = c(0.05, 0.95))
  
  all<-rbind(observed_statistic,bias,std_error,ci)
  
  all
}
MSD_process_compare_scale <- function(group1, group2, R) {
  bootstrap_results <- bootstrap_parallel_scale(group1, group2, R)
  
  bootstrap <- as.matrix(as.data.frame(bootstrap_results))
  
  # Calculate the observed statistic
  observed_statistic <- scalediff1(data = c(group1, group2), indices1 = (1:length(c(group1))),indices2 = length(c(group1))+(1:length(c(group2))))
  
  # Calculate the bias
  bias <- rowMeans2(bootstrap) - observed_statistic
  
  # Calculate the standard error
  std_error <- rowSds(bootstrap)
  
  # Calculate the confidence interval
  ci <- apply(bootstrap,1,quantile,probs = c(0.05, 0.95))
  
  all<-rbind(observed_statistic,bias,std_error,ci)
  
  all
}


tsd<-function (x,sorted=FALSE,bootsize=100000,boot=TRUE,releaseall=FALSE){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  lengthx<-length(sortedx)
  set.seed(1)
  if(lengthx>40 || boot){
    subtract<-t(replicate(bootsize, sort(sample(sortedx, size = 2))))
    getlm<-function(vector){ 
      ((vector[1]-vector[2])^2)/2
    }
    dp2varx<-apply(subtract,MARGIN=1,FUN=getlm)
    dp2varx<-Rfast::Sort(x=dp2varx,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    
    bootstrappedsample2<-c()
    
  }else{
    combinationsample2<-t(as.data.frame(Rfast::comb_n(sortedx,2)))
    getlm<-function(vector){ 
      ((vector[1]-vector[2])^2)/2
    }
    dp2varx<-apply(subtract,MARGIN=1,FUN=getlm)
    dp2varx<-Rfast::Sort(x=dp2varx,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    
    combinationsample2<-c()
    
  }
  
  #HL1<-mean(x=sortedx)
  sd1<-sqrt(mean(dp2varx))
  tsd1<-sqrt(mean(dp2varx,0.05))
  tsd2<-sqrt(mean(dp2varx,0.1))
  tsd3<-sqrt(mean(dp2varx,0.15))
  tsd4<-sqrt(mean(dp2varx,0.2))
  tsd5<-sqrt(mean(dp2varx,0.25))
  tsd6<-sqrt(mean(dp2varx,0.3))
  tsd7<-sqrt(mean(dp2varx,0.35))
  tsd8<-sqrt(mean(dp2varx,0.4))
  tsd9<-sqrt(mean(dp2varx,0.45))
  msd1<-sqrt(mediansorted(sortedx=dp2varx,lengthx=length(dp2varx)))
  if(releaseall){
    sdall<-c(sdvar=sd(dp2varx))
    return(c(sdall,msd1))
  }
  return(c(tsd0=sd1,tsd05=tsd1,tsd10=tsd2,tsd15=tsd3,tsd20=tsd4,tsd25=tsd5,tsd30=tsd6,tsd35=tsd7,tsd40=tsd8,tsd45=tsd9,msd50=msd1))
}
tmeanbatch<-function (x,sorted=FALSE){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  lengthx<-length(sortedx)
  set.seed(1)
  return(c(tm0=mean(sortedx),tm05=mean(sortedx,0.05),tm10=mean(sortedx,0.1),tm15=mean(sortedx,0.15),tm20=mean(sortedx,0.2),tm25=mean(sortedx,0.25),tm30=mean(sortedx,0.3),tm35=mean(sortedx,0.35),tm40=mean(sortedx,0.4),tm45=mean(sortedx,0.45),tm50=mean(sortedx,0.5)))
}


tHL<-function (x,sorted=FALSE,bootsize=100000,boot=TRUE,releaseall=FALSE){
  if(sorted){
    sortedx<-x
  }else{
    sortedx<-Rfast::Sort(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
  }
  lengthx<-length(sortedx)
  set.seed(1)
  if(lengthx>40 || boot){
    subtract<-t(replicate(bootsize, sort(sample(sortedx, size = 2))))
    getlm<-function(vector){ 
      ((vector[1]+vector[2]))/2
    }
    dp2HLx<-apply(subtract,MARGIN=1,FUN=getlm)
    dp2HLx<-Rfast::Sort(x=dp2HLx,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    
    bootstrappedsample2<-c()
    
  }else{
    combinationsample2<-t(as.data.frame(Rfast::comb_n(sortedx,2)))
    getlm<-function(vector){ 
      ((vector[1]+vector[2]))/2
    }
    dp2HLx<-apply(subtract,MARGIN=1,FUN=getlm)
    dp2HLx<-Rfast::Sort(x=dp2HLx,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL)
    
    combinationsample2<-c()
    
  }
  
  #HL1<-mean(x=sortedx)
  mean1<-(mean(dp2HLx))
  tHL1<-(mean(dp2HLx,0.05))
  tHL2<-(mean(dp2HLx,0.1))
  tHL3<-(mean(dp2HLx,0.15))
  tHL4<-(mean(dp2HLx,0.2))
  tHL5<-(mean(dp2HLx,0.25))
  tHL6<-(mean(dp2HLx,0.3))
  tHL7<-(mean(dp2HLx,0.35))
  tHL8<-(mean(dp2HLx,0.4))
  tHL9<-(mean(dp2HLx,0.45))
  mHL1<-(mediansorted(sortedx=dp2HLx,lengthx=length(dp2HLx)))
  if(releaseall){
    sdall<-c(sdvar=sd(dp2HLx))
    return(c(sdall,msd1))
  }
  return(c(tHL0=mean1,tHL05=tHL1,tHL10=tHL2,tHL15=tHL3,tHL20=tHL4,tHL25=tHL5,tHL30=tHL6,tHL35=tHL7,tHL40=tHL8,tHL45=tHL9,msd50=mHL1))
}





require(tidyverse)
library(dplyr)
HKBU_yuanzhu_COVID_19<- read.csv(("Yang_COVID_19_Plasma_MSD.csv"))


HKBU_yuanzhu_COVID_19 <- HKBU_yuanzhu_COVID_19 %>%
  pivot_longer(cols = starts_with("X"), names_to = "X_number", values_to = "X") %>%
  dplyr::select(-X_number) %>%
  arrange(Sample, Group)
# 
# #There are two main groups, control and biosolarization, 
# #each group has three subgroups, control, nitrogen stress, and salt stress
# #control
# 
# #colnames(Chen_mice)[colnames(Chen_mice)=="Factors"]<-"group"
colnames(HKBU_yuanzhu_COVID_19)[colnames(HKBU_yuanzhu_COVID_19)=="X"]<-"value"
# 
HKBU_yuanzhu_COVID_19$Group <- as.character(HKBU_yuanzhu_COVID_19$Group)
HKBU_yuanzhu_COVID_19$Group <- as.factor(HKBU_yuanzhu_COVID_19$Group)






require(multcomp)
HKBU_yuanzhu_COVID_19$Group <- as.character(HKBU_yuanzhu_COVID_19$Group)
HKBU_yuanzhu_COVID_19$Group <- as.factor(HKBU_yuanzhu_COVID_19$Group)

HKBU_yuanzhu_COVID_19_H <- HKBU_yuanzhu_COVID_19 %>%
  filter(Group == "Group:H")
HKBU_yuanzhu_COVID_19_CA <- HKBU_yuanzhu_COVID_19 %>%
  filter(Group == "Group:CA")
HKBU_yuanzhu_COVID_19_CO <- HKBU_yuanzhu_COVID_19 %>%
  filter(Group == "Group:CO")
# Chen_mice_wt <- HKBU_yuanzhu_COVID_19 %>%
#   filter(Group == "wt")

#Trimmed_Mean1<-trimpb(x=HKBU_yuanzhu_COVID_19_H$value,y=HKBU_yuanzhu_COVID_19_CA$value, tr=0.2,alpha=0.05,nboot=2000,WIN=FALSE,plotit=TRUE,win=0.1,pop=1)


# Function to calculate the difference or statistic of interest
group1 <-HKBU_yuanzhu_COVID_19_H$value
group2 <-HKBU_yuanzhu_COVID_19_CA$value


# Set the number of bootstrap iterations
R <- 1000

all<-MSD_process_compare_location(group1, group2, R)


write.csv(all,paste("location_COVID_H_CA.csv", sep = ","), row.names = FALSE)


group1 <-HKBU_yuanzhu_COVID_19_H$value
group2 <-HKBU_yuanzhu_COVID_19_CO$value


# Set the number of bootstrap iterations
R <- 1000

all<-MSD_process_compare_location(group1, group2, R)


write.csv(all,paste("location_COVID_H_CO.csv", sep = ","), row.names = FALSE)


HL_H<-mean(x=HKBU_yuanzhu_COVID_19_H$value)

HL_CA<-mean(x=HKBU_yuanzhu_COVID_19_CA$value)

HL_CO<-mean(x=HKBU_yuanzhu_COVID_19_CO$value)

#HL_wt<-mean(x=Chen_mice_wt$value)


HL_H
HL_CA
HL_CO

#no significant changes among health and recovery patients with antibodies (CA) or rapid-faded antibodies (CO)

#reasonable considering the key role of vitamin A in metabolism.

#however, the trend is more severe in wile-type, since the TRACK group already has vitamin A deficiency even under normal diet.


group1 <-HKBU_yuanzhu_COVID_19_H$value
group2 <-HKBU_yuanzhu_COVID_19_CA$value


# Set the number of bootstrap iterations

R <- 1000

all<-MSD_process_compare_scale(group1, group2, R)


write.csv(all,paste("scale_COVID_H_CA.csv", sep = ","), row.names = FALSE)


#control vs. salt stress
group1 <-HKBU_yuanzhu_COVID_19_H$value
group2 <-HKBU_yuanzhu_COVID_19_CO$value

R <- 1000

all<-MSD_process_compare_scale(group1, group2, R)


write.csv(all,paste("scale_COVID_H_CO.csv", sep = ","), row.names = FALSE)



msd2_H<-sd(x=HKBU_yuanzhu_COVID_19_H$value)

msd2_CA<-sd(x=HKBU_yuanzhu_COVID_19_CA$value)

msd2_CO<-sd(x=HKBU_yuanzhu_COVID_19_CO$value)

msd2_H
msd2_CA
msd2_CO
#no significant changes among health and recovery patients with antibodies (CA) or rapid-faded antibodies (CO)

grouped_statistics <- HKBU_yuanzhu_COVID_19 %>%
  group_by(Group) %>%
  summarise(HL = mean(value), msd = sd(value))
print(grouped_statistics)

write.csv(grouped_statistics,paste("All_COVID-19_1.csv", sep = ","), row.names = FALSE)

library(dplyr)

grouped_statistics <- HKBU_yuanzhu_COVID_19 %>%
  group_by(Group) 


group_tHL <- aggregate(grouped_statistics$value, list(grouped_statistics$Group), tHL)

print(group_tHL)

group_tsd <- aggregate(grouped_statistics$value, list(grouped_statistics$Group), tsd)

# Print the result
print(group_tsd)

ranked_group_tsd <- apply(group_tsd, 2, rank)
print(ranked_group_tsd)

write.csv(rbind(group_tHL,group_tsd),paste("All_COVID-19_comparison.csv", sep = ","), row.names = FALSE)



#%>%
#  summarise(HL = tsd(x=value,sorted=FALSE,bootsize=100000,boot=TRUE,releaseall=FALSE), msd = sd(value))
#print(grouped_statistics)
#tsd(value,sorted=FALSE,bootsize=100000,boot=TRUE,releaseall=FALSE)








HKBU_yuanzhu_COVID_19_MSD<- read.csv(("Yang_COVID_19_Plasma_MSD_Class.csv"))
HKBU_yuanzhu_COVID_19_MSD_list <- split(HKBU_yuanzhu_COVID_19_MSD, HKBU_yuanzhu_COVID_19_MSD$X)
for (i in seq_along(HKBU_yuanzhu_COVID_19_MSD_list)) {
  first_row <- HKBU_yuanzhu_COVID_19_MSD[1, ]
  if (!identical(HKBU_yuanzhu_COVID_19_MSD_list[[i]][1, ], first_row)) {
    HKBU_yuanzhu_COVID_19_MSD_list[[i]] <- rbind(first_row, HKBU_yuanzhu_COVID_19_MSD_list[[i]])
  }
}
grouped_statisticsall<-c()
ciall<-c()
for (i in seq_along(HKBU_yuanzhu_COVID_19_MSD_list)) {
  HKBU_yuanzhu_COVID_19_MSD_list_1q <- HKBU_yuanzhu_COVID_19_MSD_list[[i]]
  HKBU_yuanzhu_COVID_19_MSD_list_1 <- t(HKBU_yuanzhu_COVID_19_MSD_list_1q)
  if(nrow(HKBU_yuanzhu_COVID_19_MSD_list_1)<=2){
    next
  }
  if(ncol(HKBU_yuanzhu_COVID_19_MSD_list_1)<=2){
    next
  }
  HKBU_yuanzhu_COVID_19_MSD_list_1<-HKBU_yuanzhu_COVID_19_MSD_list_1[3:nrow(HKBU_yuanzhu_COVID_19_MSD_list_1),]
  colnames(HKBU_yuanzhu_COVID_19_MSD_list_1)<-c("Group",paste("X", 1:(ncol(HKBU_yuanzhu_COVID_19_MSD_list_1)-1), sep = ""))
  require(multcomp)
  HKBU_yuanzhu_COVID_19_MSD_list_1<-as.data.frame(HKBU_yuanzhu_COVID_19_MSD_list_1)
  require(tidyverse)
  library(dplyr)
  HKBU_yuanzhu_COVID_19_MSD_list_10 <- HKBU_yuanzhu_COVID_19_MSD_list_1 %>%
    pivot_longer(cols = starts_with("X"), names_to = "X_number", values_to = "X") %>%
    dplyr::select(-X_number) %>%
    arrange(Group)
  
  #There are two main groups, solarization and biosolarization, 
  #each group has three subgroups, control, nitrogen stress, and salt stress
  #control
  
  
  HKBU_yuanzhu_COVID_19_MSD_list_10<-as.data.frame(HKBU_yuanzhu_COVID_19_MSD_list_10)
  HKBU_yuanzhu_COVID_19_MSD_list_10$Group <- as.character(HKBU_yuanzhu_COVID_19_MSD_list_10$Group)
  HKBU_yuanzhu_COVID_19_MSD_list_10$Group <- as.factor(HKBU_yuanzhu_COVID_19_MSD_list_10$Group)
  colnames(HKBU_yuanzhu_COVID_19_MSD_list_10)[colnames(HKBU_yuanzhu_COVID_19_MSD_list_10)=="X"]<-"value"
  
  
  if (!is.numeric(HKBU_yuanzhu_COVID_19_MSD_list_10$value)) {
    HKBU_yuanzhu_COVID_19_MSD_list_10$value <- as.numeric(HKBU_yuanzhu_COVID_19_MSD_list_10$value)
  }
  #ff5<-((c(HKBU_yuanzhu_COVID_19_MSD_list_1q[2,1],"observed_statistic","bias","std_error","ci.5%","ci.95%")))
  #rep(HKBU_yuanzhu_COVID_19_MSD_list_1q[2,1],22)
  ff5<-rep(HKBU_yuanzhu_COVID_19_MSD_list_1q[2,1],24)
  ciall<-rbind(ciall,ff5)

  HKBU_yuanzhu_COVID_19_MSD_list_10_H <- HKBU_yuanzhu_COVID_19_MSD_list_10 %>%
    filter(Group == "Group:H")
  HKBU_yuanzhu_COVID_19_MSD_list_10_CA <- HKBU_yuanzhu_COVID_19_MSD_list_10 %>%
    filter(Group == "Group:CA")
  HKBU_yuanzhu_COVID_19_MSD_list_10_CO <- HKBU_yuanzhu_COVID_19_MSD_list_10 %>%
    filter(Group == "Group:CO")
  # Chen_mice_wt <- HKBU_yuanzhu_COVID_19 %>%
  #   filter(Group == "wt")
  
  #Trimmed_Mean1<-trimpb(x=HKBU_yuanzhu_COVID_19_H$value,y=HKBU_yuanzhu_COVID_19_CA$value, tr=0.2,alpha=0.05,nboot=2000,WIN=FALSE,plotit=TRUE,win=0.1,pop=1)
  
  
  group1 <-HKBU_yuanzhu_COVID_19_MSD_list_10_H$value
  group2 <-HKBU_yuanzhu_COVID_19_MSD_list_10_CA$value
  
  R <- 1000
  
  all5<-cbind(rep("H-CA",5),MSD_process_compare_location(group1, group2, R),rep("H-CA",5),MSD_process_compare_scale(group1, group2, R))
  
  
  #HCA1<-c("H-CA",MSD_process_compare_location(group1, group2, R),"H-CA",MSD_process_compare_scale(group1, group2, R))
  
  ciall<-rbind(ciall,all5)
  
  group1 <-HKBU_yuanzhu_COVID_19_MSD_list_10_H$value
  group2 <-HKBU_yuanzhu_COVID_19_MSD_list_10_CO$value
  
  R <- 1000
  
  all5<-cbind(rep("H-CO",5),MSD_process_compare_location(group1, group2, R),rep("H-CO",5),MSD_process_compare_scale(group1, group2, R))
  
  #HCA1<-c("H-CA",MSD_process_compare_location(group1, group2, R),"H-CA",MSD_process_compare_scale(group1, group2, R))
  
  ciall<-rbind(ciall,all5)
  
  group1 <-HKBU_yuanzhu_COVID_19_MSD_list_10_CA$value
  group2 <-HKBU_yuanzhu_COVID_19_MSD_list_10_CO$value
  
  R <- 1000
  
  all5<-cbind(rep("CA-CO",5),MSD_process_compare_location(group1, group2, R),rep("CA-CO",5),MSD_process_compare_scale(group1, group2, R))
  
  #HCA1<-c("H-CA",MSD_process_compare_location(group1, group2, R),"H-CA",MSD_process_compare_scale(group1, group2, R))
  
  ciall<-rbind(ciall,all5)
  
  #no significant changes among health and recovery patients with antibodies (CA) or rapid-faded antibodies (CO)
  
  grouped_statistics <- HKBU_yuanzhu_COVID_19_MSD_list_10 %>%
    group_by(Group) 
  
  
  group_tHL <- aggregate(grouped_statistics$value, list(grouped_statistics$Group), tHL)
  
  print(group_tHL)
  
  group_tsd <- aggregate(grouped_statistics$value, list(grouped_statistics$Group), tsd)
  
  # Print the result
  print(group_tsd)
  
  ranked_group_tsd <- apply(group_tsd, 2, rank)
  print(ranked_group_tsd)
  grouped_statistics1<-rbind(group_tHL,group_tsd)

  ff4<-t(t(rep(HKBU_yuanzhu_COVID_19_MSD_list_1q[2,1],nrow(grouped_statistics1))))
  
  grouped_statistics1<-cbind(ff4,grouped_statistics1)

  grouped_statisticsall<-rbind(grouped_statisticsall,grouped_statistics1)
  
  
}
write.csv(grouped_statisticsall,paste("grouped_statisticsall.csv", sep = ","), row.names = FALSE)
write.csv(ciall,paste("ciall.csv", sep = ","), row.names = FALSE)
    

