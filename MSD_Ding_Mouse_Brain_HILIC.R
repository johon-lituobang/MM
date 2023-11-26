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
Junding_Mice_HILIC_positive_MSD<- read.csv(("Ding_Mouse_Brain_HILIC_positive_MSD.csv"))


Junding_Mice_HILIC_positive_MSD <- Junding_Mice_HILIC_positive_MSD %>%
  pivot_longer(cols = starts_with("X"), names_to = "X_number", values_to = "X") %>%
  dplyr::select(-X_number) %>%
  arrange(Sample, Group,Region,Age,Gender,AG)
colnames(Junding_Mice_HILIC_positive_MSD)[colnames(Junding_Mice_HILIC_positive_MSD)=="X"]<-"value"
# 
Junding_Mice_HILIC_positive_MSD$AG <- as.character(Junding_Mice_HILIC_positive_MSD$AG)
Junding_Mice_HILIC_positive_MSD$AG <- as.factor(Junding_Mice_HILIC_positive_MSD$AG)






require(multcomp)
Junding_Mice_HILIC_positive_MSD$AG <- as.character(Junding_Mice_HILIC_positive_MSD$AG)
Junding_Mice_HILIC_positive_MSD$AG <- as.factor(Junding_Mice_HILIC_positive_MSD$AG)
Junding_Mice_HILIC_positive_MSD_3 <- Junding_Mice_HILIC_positive_MSD %>%
  filter(AG == "Age:3weeksGender:Female")

Junding_Mice_HILIC_positive_MSD_59 <- Junding_Mice_HILIC_positive_MSD %>%
  filter(AG == "Age:59weeksGender:Female")
Junding_Mice_HILIC_positive_MSD_3_M <- Junding_Mice_HILIC_positive_MSD %>%
  filter(AG == "Age:3weeksGender:Male")

Junding_Mice_HILIC_positive_MSD_59_M <- Junding_Mice_HILIC_positive_MSD %>%
  filter(AG == "Age:59weeksGender:Male")



# Function to calculate the difference or statistic of interest

group1 <-Junding_Mice_HILIC_positive_MSD_3$value
group2 <-Junding_Mice_HILIC_positive_MSD_59$value



R <- 5

all<-MSD_process_compare_location(group1, group2, R)

write.csv(all,paste("location_mouse_HILIC_positive_3_59_female.csv", sep = ","), row.names = FALSE)




group1 <-Junding_Mice_HILIC_positive_MSD_3_M$value
group2 <-Junding_Mice_HILIC_positive_MSD_59_M$value


R <- 5

all<-MSD_process_compare_location(group1, group2, R)

write.csv(all,paste("location_mouse_HILIC_positive_3_59_male.csv", sep = ","), row.names = FALSE)


group1 <-Junding_Mice_HILIC_positive_MSD_3$value
group2 <-Junding_Mice_HILIC_positive_MSD_3_M$value

R <- 5

all<-MSD_process_compare_location(group1, group2, R)

write.csv(all,paste("location_mouse_HILIC_positive_3.csv", sep = ","), row.names = FALSE)


group1 <-Junding_Mice_HILIC_positive_MSD_59$value
group2 <-Junding_Mice_HILIC_positive_MSD_59_M$value

R <- 5

all<-MSD_process_compare_location(group1, group2, R)

write.csv(all,paste("location_mouse_HILIC_positive_59.csv", sep = ","), row.names = FALSE)


#no significant changes among health and recovery patients with antibodies (CA) or rapid-faded antibodies (CO)

#reasonable considering the key role of vitamin A in metabolism.

#however, the trend is more severe in wile-type, since the TRACK group already has vitamin A deficiency even under normal diet.

# Install and load the 'boot' package
if (!require("boot")) install.packages("boot")
library(boot)

# Function to calculate the difference or statistic of interest
group1 <-Junding_Mice_HILIC_positive_MSD_3$value
group2 <-Junding_Mice_HILIC_positive_MSD_59$value

R <- 5

all<-MSD_process_compare_scale(group1, group2, R)

write.csv(all,paste("scale_mouse_HILIC_positive_3_59_female.csv", sep = ","), row.names = FALSE)


group1 <-Junding_Mice_HILIC_positive_MSD_3_M$value
group2 <-Junding_Mice_HILIC_positive_MSD_59_M$value

R <- 5

all<-MSD_process_compare_scale(group1, group2, R)

write.csv(all,paste("scale_mouse_HILIC_positive_3_59_male.csv", sep = ","), row.names = FALSE)


group1 <-Junding_Mice_HILIC_positive_MSD_3$value
group2 <-Junding_Mice_HILIC_positive_MSD_3_M$value

R <- 5

all<-MSD_process_compare_scale(group1, group2, R)

write.csv(all,paste("scale_mouse_HILIC_positive_3.csv", sep = ","), row.names = FALSE)

#no significant changes among health and recovery patients with antibodies (CA) or rapid-faded antibodies (CO)


group1 <-Junding_Mice_HILIC_positive_MSD_59$value
group2 <-Junding_Mice_HILIC_positive_MSD_59_M$value

R <- 5

all<-MSD_process_compare_scale(group1, group2, R)

write.csv(all,paste("scale_mouse_HILIC_positive_59.csv", sep = ","), row.names = FALSE)


#no significant changes among health and recovery patients with antibodies (CA) or rapid-faded antibodies (CO)

grouped_statistics <- Junding_Mice_HILIC_positive_MSD %>%
  group_by(AG) %>%
  summarise(HL = wilcox.test(x=value,conf.int = TRUE)$estimate, msd = msd(value))
print(grouped_statistics)

write.csv(grouped_statistics,paste("Junding_Mice_HILIC_positive_MSD_statistics.csv", sep = ","), row.names = FALSE)


library(dplyr)

grouped_statistics <- Junding_Mice_HILIC_positive_MSD %>%
  group_by(Group) 


group_tHL <- aggregate(grouped_statistics$value, list(grouped_statistics$AG), tHL)

print(group_tHL)
ranked_group_tHL <- apply(group_tHL, 2, rank)
print(ranked_group_tHL)
group_tsd <- aggregate(grouped_statistics$value, list(grouped_statistics$AG), tsd)

# Print the result
print(group_tsd)

ranked_group_tsd <- apply(group_tsd, 2, rank)
print(ranked_group_tsd)

write.csv(rbind(group_tHL,group_tsd),paste("Junding_Mice_HILIC_positive_comparison.csv", sep = ","), row.names = FALSE)


write.csv(rbind(ranked_group_tHL,ranked_group_tsd),paste("Junding_Mice_HILIC_positive_comparison_rank.csv", sep = ","), row.names = FALSE)










require(tidyverse)
library(dplyr)
Junding_Mice_HILIC_negative_MSD<- read.csv(("Ding_Mouse_Brain_HILIC_negative_MSD.csv"))


Junding_Mice_HILIC_negative_MSD <- Junding_Mice_HILIC_negative_MSD %>%
  pivot_longer(cols = starts_with("X"), names_to = "X_number", values_to = "X") %>%
  dplyr::select(-X_number) %>%
  arrange(Sample, Group,Region,Age,Gender,AG)
colnames(Junding_Mice_HILIC_negative_MSD)[colnames(Junding_Mice_HILIC_negative_MSD)=="X"]<-"value"
# 
Junding_Mice_HILIC_negative_MSD$AG <- as.character(Junding_Mice_HILIC_negative_MSD$AG)
Junding_Mice_HILIC_negative_MSD$AG <- as.factor(Junding_Mice_HILIC_negative_MSD$AG)






require(multcomp)
Junding_Mice_HILIC_negative_MSD$AG <- as.character(Junding_Mice_HILIC_negative_MSD$AG)
Junding_Mice_HILIC_negative_MSD$AG <- as.factor(Junding_Mice_HILIC_negative_MSD$AG)
Junding_Mice_HILIC_negative_MSD_3 <- Junding_Mice_HILIC_negative_MSD %>%
  filter(AG == "Age:3weeksGender:Female")

Junding_Mice_HILIC_negative_MSD_59 <- Junding_Mice_HILIC_negative_MSD %>%
  filter(AG == "Age:59weeksGender:Female")
Junding_Mice_HILIC_negative_MSD_3_M <- Junding_Mice_HILIC_negative_MSD %>%
  filter(AG == "Age:3weeksGender:Male")

Junding_Mice_HILIC_negative_MSD_59_M <- Junding_Mice_HILIC_negative_MSD %>%
  filter(AG == "Age:59weeksGender:Male")



# Function to calculate the difference or statistic of interest

group1 <-Junding_Mice_HILIC_negative_MSD_3$value
group2 <-Junding_Mice_HILIC_negative_MSD_59$value



R <- 5

all<-MSD_process_compare_location(group1, group2, R)

write.csv(all,paste("location_mouse_HILIC_negative_3_59_female.csv", sep = ","), row.names = FALSE)




group1 <-Junding_Mice_HILIC_negative_MSD_3_M$value
group2 <-Junding_Mice_HILIC_negative_MSD_59_M$value


R <- 5

all<-MSD_process_compare_location(group1, group2, R)

write.csv(all,paste("location_mouse_HILIC_negative_3_59_male.csv", sep = ","), row.names = FALSE)


group1 <-Junding_Mice_HILIC_negative_MSD_3$value
group2 <-Junding_Mice_HILIC_negative_MSD_3_M$value

R <- 5

all<-MSD_process_compare_location(group1, group2, R)

write.csv(all,paste("location_mouse_HILIC_negative_3.csv", sep = ","), row.names = FALSE)


group1 <-Junding_Mice_HILIC_negative_MSD_59$value
group2 <-Junding_Mice_HILIC_negative_MSD_59_M$value

R <- 5

all<-MSD_process_compare_location(group1, group2, R)

write.csv(all,paste("location_mouse_HILIC_negative_59.csv", sep = ","), row.names = FALSE)


#no significant changes among health and recovery patients with antibodies (CA) or rapid-faded antibodies (CO)

#reasonable considering the key role of vitamin A in metabolism.

#however, the trend is more severe in wile-type, since the TRACK group already has vitamin A deficiency even under normal diet.

# Install and load the 'boot' package
if (!require("boot")) install.packages("boot")
library(boot)

# Function to calculate the difference or statistic of interest
group1 <-Junding_Mice_HILIC_negative_MSD_3$value
group2 <-Junding_Mice_HILIC_negative_MSD_59$value

R <- 5

all<-MSD_process_compare_scale(group1, group2, R)

write.csv(all,paste("scale_mouse_HILIC_negative_3_59_female.csv", sep = ","), row.names = FALSE)


group1 <-Junding_Mice_HILIC_negative_MSD_3_M$value
group2 <-Junding_Mice_HILIC_negative_MSD_59_M$value

R <- 5

all<-MSD_process_compare_scale(group1, group2, R)

write.csv(all,paste("scale_mouse_HILIC_negative_3_59_male.csv", sep = ","), row.names = FALSE)


group1 <-Junding_Mice_HILIC_negative_MSD_3$value
group2 <-Junding_Mice_HILIC_negative_MSD_3_M$value

R <- 5

all<-MSD_process_compare_scale(group1, group2, R)

write.csv(all,paste("scale_mouse_HILIC_negative_3.csv", sep = ","), row.names = FALSE)

#no significant changes among health and recovery patients with antibodies (CA) or rapid-faded antibodies (CO)


group1 <-Junding_Mice_HILIC_negative_MSD_59$value
group2 <-Junding_Mice_HILIC_negative_MSD_59_M$value

R <- 5

all<-MSD_process_compare_scale(group1, group2, R)

write.csv(all,paste("scale_mouse_HILIC_negative_59.csv", sep = ","), row.names = FALSE)


#no significant changes among health and recovery patients with antibodies (CA) or rapid-faded antibodies (CO)

grouped_statistics <- Junding_Mice_HILIC_negative_MSD %>%
  group_by(AG) %>%
  summarise(HL = wilcox.test(x=value,conf.int = TRUE)$estimate, msd = msd(value))
print(grouped_statistics)

write.csv(grouped_statistics,paste("Junding_Mice_HILIC_negative_MSD_statistics.csv", sep = ","), row.names = FALSE)







library(dplyr)

grouped_statistics <- Junding_Mice_HILIC_negative_MSD %>%
  group_by(AG) 


group_tHL <- aggregate(grouped_statistics$value, list(grouped_statistics$AG), tHL)

print(group_tHL)
ranked_group_tHL <- apply(group_tHL, 2, rank)
print(ranked_group_tHL)
group_tsd <- aggregate(grouped_statistics$value, list(grouped_statistics$AG), tsd)

# Print the result
print(group_tsd)

ranked_group_tsd <- apply(group_tsd, 2, rank)
print(ranked_group_tsd)

write.csv(rbind(group_tHL,group_tsd),paste("Junding_Mice_HILIC_negative_comparison.csv", sep = ","), row.names = FALSE)


write.csv(rbind(ranked_group_tHL,ranked_group_tsd),paste("Junding_Mice_HILIC_negative_comparison_rank.csv", sep = ","), row.names = FALSE)

