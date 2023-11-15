## model 2 - model ARG abundance given crAssphage, mexK, flow, and day of week
model2_summary=function(prd,data,drug_class,drug_class_use){
  
  period=prd
  omega = 2*pi/period
  
  data$y=data$sum
  start_intercept=mean(data$y)
  start_slope=0
  
  model = nls(y ~ amplitude * cos(omega * (month_ordered - phase)) +  # sinusoidal component
               # (use_scaled*use_slope)              + # linear component for usage 
                (flow_mgd_scaled*flow_slope)              + # linear component for flow. 
                (crassphagespergbp_scaled*crass_slope) +  # linear component for crAssphage
                (DOW*day_slope)+
                (mexK_scaled*mexK_slope)+intercept, # linear component for mexK
              
              start = list(amplitude = 0.1,
                           phase = 0,
                           day_slope=0,
                           use_slope=0,
                           mexK_slope=0,
                           flow_slope=0,
                           crass_slope = start_slope,
                           intercept = start_intercept),data=data)
  
  coeffs=summary(model)$coefficients%>%as_tibble(rownames = "term") %>% magrittr::set_colnames(c("term","estimate","std.error","statistic","p.value")) 
  cints=confint.default(model) %>% magrittr::set_colnames(c("ci.lower", "ci.upper")) %>% as_tibble(rownames = "term")
  aic=AIC(model)  
  
  out=data.frame(coeffs,cints)[,c("term","estimate","std.error","statistic","p.value","ci.lower","ci.upper")]
  
  out$AIC=aic
  out$period=period
  out.adj=convert_a_phases_func(out)
  out$p.adj=p.adjust(out$p.value)
  out$raw_or_adjusted="raw"
  out.adj$raw_or_adjusted="adj"
  
  SSE=sum((data$y - fitted(model))^2)
  SST=sum((data$y - mean(data$y))^2)
  R2=1-SSE/SST
  
  fo=rbind.fill(out,out.adj)
  fo$R2=R2
  fo$period=period
  fo$drug_class=drug_class
  return(fo)
}


## model 2.1 includes a linear component for antibiotic usage - otherwise - same as 2
model2.1_summary=function(prd,data,drug_class,drug_class_use){
  
  period=prd
  omega = 2*pi/period
  
  data$y=data$sum
  start_intercept=mean(data$y)
  start_slope=0
  
  model = nls(y ~ amplitude * cos(omega * (month_ordered - phase)) +  # sinusoidal component
                (use_scaled*use_slope)              + # linear component for usage 
                (flow_mgd_scaled*flow_slope)              + # linear component for flow. 
                (crassphagespergbp_scaled*crass_slope) +  # linear component for crAssphage
                (DOW*day_slope)+
                (mexK_scaled*mexK_slope)+intercept, # linear component for mexK
              
              start = list(amplitude = 0.1,
                           phase = 0,
                           day_slope=0,
                           use_slope=0,
                           mexK_slope=0,
                           flow_slope=0,
                           crass_slope = start_slope,
                           intercept = start_intercept),data=data)
  
  coeffs=summary(model)$coefficients%>%as_tibble(rownames = "term") %>% magrittr::set_colnames(c("term","estimate","std.error","statistic","p.value")) 
  cints=confint.default(model) %>% magrittr::set_colnames(c("ci.lower", "ci.upper")) %>% as_tibble(rownames = "term")
  aic=AIC(model)  
  
  out=data.frame(coeffs,cints)[,c("term","estimate","std.error","statistic","p.value","ci.lower","ci.upper")]
  
  out$AIC=aic
  out$period=period
  out.adj=convert_a_phases_func(out)
  out$p.adj=p.adjust(out$p.value)
  out$raw_or_adjusted="raw"
  out.adj$raw_or_adjusted="adj"
  
  SSE=sum((data$y - fitted(model))^2)
  SST=sum((data$y - mean(data$y))^2)
  R2=1-SSE/SST
  
  fo=rbind.fill(out,out.adj)
  fo$R2=R2
  fo$period=period
  fo$drug_class=drug_class
  return(fo)
}







