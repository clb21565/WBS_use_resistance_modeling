
# data start in august 2020 and goes to august 2021; make datafarme in the appropriate order for modeling.  

month_ordered=c(1,2,3,4,5,6,7,8,9,10,11,12,13)
month=c(8,9,10,11,12,1,2,3,4,5,6,7,8)
Year=c(2020,2020,2020,2020,2020,2021,2021,2021,2021,2021,2021,2021,2021)
order=data.frame(month_ordered,month,Year)

#function for calculating minmax scaling of data 
minmax_scale=function(values){
  normalized_values=(values - (min(values))) / (max(values) - min(values))
  return(normalized_values)
}

add_minmax=function(df){
  df$use = minmax_scale(df$total_class_use)
  df$flow_mgd_scaled_max=minmax_scale(df$flow_mgd)
  df$crassphagespergbp_scaled_max=minmax_scale(df$crassphagespergbp)
  df$mexK_scaled_max=minmax_scale(df$mexK)
  df$sum_scaled=minmax_scale(df$sum)
  df$`Abx Usage`=minmax_scale(df$total_class_use)
  df$`Influent flow`=minmax_scale(df$flow_mgd)
  df$`crAssphage`=minmax_scale(df$crassphagespergbp)
  df$`mexK rel abund`=minmax_scale(df$mexK)
  return(df)
}
