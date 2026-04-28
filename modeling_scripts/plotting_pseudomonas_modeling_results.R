library(tidyr)
library(data.table)
library(plyr)
library(dplyr)
library(ggplot2)
library(readxl)
library(pals)
library(ggthemes)
library(gggenes)
library(ggpubr)
library(gridExtra)
library(stringr)
library(stringdist)
library(ggbeeswarm)
library(vegan)
library(Hmisc)


use=fread("240323_ABXUSEDATAFORMODELING.csv")[,-1]
op_pdf=fread("240614-op_metadata-for-modeling.csv")
mcp=fread("240614-crassphage-for-modeling.csv")
crass=mcp[,c("sample","minmax_crassphage")]%>%unique()
data=fread("240616-reformatted-mag-abundance-data-pseudomonas-modeling.csv")[,-1]
inf_flow=op_pdf%>%subset(Stage=="Influent"&variable=="Flow (MGD)")%>%mutate(x=(ordered_week/4)-0.25,IF=mean)
sxt_use=use%>%subset(k==4)%>%subset(grepl("Folate Antag",drug_class))%>%mutate(x=(ordered_week/4)-0.25,use=lagged_value)
use_data=sxt_use[,c("x","use")]
model_data=data%>%#subset(genomeid==gnid)%>%
    merge(.,use_data,by="x")%>%
    merge(.,inf_flow[,c("Date","IF")]%>%unique(),by="Date")%>%merge(.,crass,by="sample")
#model_data%>%write.csv("240828-FormattedPseudomonasModelingData.csv")
model_data=fread("240828-FormattedPseudomonasModelingData.csv")[,-1]

spkey=fread("species_name_key_pseudomonas.csv")
cos_func = function(month, amplitude, phase, omega,intercept) {
  amplitude * cos(omega *(month - phase)) + intercept
}
getAdjustedXcorrs=function(gnid){
  #gnid="gs35"
  period=12
  omega = 2*pi/period
  working_model=model3%>%subset(genomeid==gnid)
  working_data=model_data%>%subset(genomeid==gnid)
  
  phase=(working_model%>%subset(raw_or_adjusted=="adj")%>%subset(term=="phase"))$estimate  
  amplitude=(working_model%>%subset(raw_or_adjusted=="adj")%>%subset(term=="amplitude"))$estimate
  alow=(working_model%>%subset(raw_or_adjusted=="adj")%>%subset(term=="amplitude"))$ci.lower
  ahigh=(working_model%>%subset(raw_or_adjusted=="adj")%>%subset(term=="amplitude"))$ci.upper
  intercept=(working_model%>%subset(raw_or_adjusted=="raw")%>%subset(term=="intercept"))$estimate
  BSXT=(working_model%>%subset(raw_or_adjusted=="raw")%>%subset(term=="B_SXT"))$estimate
  stat_function(fun = cos_func, args = list(a = amplitude, phase = phase, omega = omega, intercept = 0), size = 0.7)   

working_data$model_predicted_data=amplitude*cos(omega*(working_data$x-phase))+(BSXT*working_data$use)+intercept
working_data$model_predicted_data_no_use=amplitude*cos(omega*(working_data$x-phase))+intercept
working_data$seasonality_val=(amplitude*cos(omega*(working_data$x-phase))+intercept)
working_data$bsxt_contribution=(BSXT*working_data$use)+intercept

  ci = data.frame(month=seq(0,12,0.01)) %>%
    mutate(lower_ci = purrr::map_dbl(month, ~cos_func(., alow, phase, omega, 0))) %>%
    mutate(upper_ci = purrr::map_dbl(month, ~cos_func(., ahigh, phase, omega, 0))) 
  
  
res=working_data%>%
  ggplot(aes(x)) + 
  geom_point(aes(x = x, y = y-mean(y)),shape=21,size=2)+ 
  #geom_point(aes(x=x,y=model_predicted_data_no_use-mean(model_predicted_data_no_use)),shape=21,fill="black",size=2)+
  geom_point(aes(x=x,y=model_predicted_data-mean(model_predicted_data)),shape=21,fill="black",size=2)+
  stat_function(fun = cos_func, args = list(a = amplitude, phase = phase, omega = omega, intercept = 0), size = 0.7)+
  geom_ribbon(data = ci, aes(x = month, ymin = lower_ci, ymax = upper_ci), alpha = 0.3)+cleanup+
  labs(x="",y="Seasonal deviate",x="Month",title=gnid)+
  #scale_x_continuous(breaks=c(1,3,5,7,9,11),labels=c("Sept 2020","Nov 2020","Jan 2021","Mar 2021","May 2021","July 2021"))+ 
  scale_x_continuous(breaks=c(1,3,5,7,9,11),labels=c("09-20","11-20","01-21","03-21","05-21","07-21"))+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme(axis_title.y=element_blank())

return(res)
}

list_of_plots=lapply(unique(m3sig$genomeid),getAdjustedXcorrs)

gp=ggarrange(plotlist = list_of_plots)


