# Applies Jacobian LE chaos detection method to time series

library(rEDM)
library(dplyr)
library(tidyr)
library(purrr)
#to run in parallel
#library(furrr)

source("scripts/chaos_localstability_methods.R")

ts_d=read.csv("./data/alldata_filtered.csv")
ts_d=ts_d %>% group_by(Site, Name, Level) %>% nest(.key="data_rescale") %>% 
  mutate(data_rescale=map(data_rescale, as.data.frame))
#for updating
#dnestfilt=ts_d

ts_results=select(ts_d, Site, Name, Level)

# for testing
# ts_d=ts_d[1:2,]
# ts_results=select(ts_d, Site, Name, Level)
# ts_d=ts_d %>% filter(Site=="Port Erin Bay")
# ts_results=select(ts_d, Site, Name, Level)
# ts_results$modelresults3=map(ts_d$data_rescale, smap_model_options, y="AbundanceScale", model=3, seaspred=F) #fd-ut
# ts_results$modelresults4=map(ts_d$data_rescale, smap_model_options, y="AbundanceScale", model=4, seaspred=F) #gr-ut
# ts_results$modelresults5=map(ts_d$data_rescale, smap_model_options, y="AbundanceScale", model=5, seaspred=F) #gr-log

#original model
# ts_results$modelresults5=future_map(ts_d$data_rescale, smap_model_options, y="AbundanceScale", model=5, Efix=4, taufix=1) #gr-log
# ts_results$modelresultsbest=ts_results$modelresults5
# ts_d$bestR2=map_dbl(ts_results$modelresults5, ~.x$modelstats$R2abund) #R2 for abundance
# ts_d$bestR2m=map_dbl(ts_results$modelresults5, ~.x$modelstats$R2model) #R2 for growth rate

#### Jacobian LE Method ####

#fit models
#models 1 and 2 are dropped because the results of 1&3 and 2&5 are identical
ts_results$modelresults3=map(ts_d$data_rescale, smap_model_options, y="AbundanceScale", model=3, seaspred=F) #fd-ut
ts_results$modelresults4=map(ts_d$data_rescale, smap_model_options, y="AbundanceScale", model=4, seaspred=F) #gr-ut
ts_results$modelresults5=map(ts_d$data_rescale, smap_model_options, y="AbundanceScale", model=5, seaspred=F) #gr-log

#pull out R2 values for each model
ts_d$R3m=map_dbl(ts_results$modelresults3, ~.x$modelstats$R2model)
ts_d$R4m=map_dbl(ts_results$modelresults4, ~.x$modelstats$R2model)
ts_d$R5m=map_dbl(ts_results$modelresults5, ~.x$modelstats$R2model)
ts_d$R3a=map_dbl(ts_results$modelresults3, ~.x$modelstats$R2abund)
ts_d$R4a=map_dbl(ts_results$modelresults4, ~.x$modelstats$R2abund)
ts_d$R5a=map_dbl(ts_results$modelresults5, ~.x$modelstats$R2abund)

#get best model of the 3 forms
ts_d$bestmodel=select(ts_d,R3a,R4a,R5a) %>% apply(1,which.max)
ts_d$bestR2=select(ts_d,R3a,R4a,R5a) %>% apply(1,max) #R2 for abundance
ts_d$bestR2m=select(ts_d,R3m,R4m,R5m) %>% apply(1,max) #R2 for growth rate
ts_results$modelresultsbest=cbind(select(ts_results, modelresults3, modelresults4, modelresults5),ts_d$bestmodel) %>% apply(1, function(x) {m=as.numeric(x["ts_d$bestmodel"]); x[m][[1]]})
#get E, tau, theta values from best model
ts_d$E=map_dbl(ts_results$modelresultsbest, ~.x$modelstats$E)
ts_d$tau=map_dbl(ts_results$modelresultsbest, ~.x$modelstats$tau)
ts_d$theta=map_dbl(ts_results$modelresultsbest, ~.x$modelstats$theta)
#get best model form
ts_d$modelform=map_chr(ts_results$modelresultsbest, ~.x$form)

#get jacobian matrices
ts_results$jacobians=map(ts_results$modelresultsbest, getJacobians)

#get local LE (eigenvalue)
ts_d$LocalLE = map2(ts_results$modelresultsbest, ts_results$jacobians, LElocal)

#get annual LE
ts_d$AnnualLE = pmap(list(ts_d$data_rescale, ts_results$modelresultsbest, ts_results$jacobians), LEannual)

#get global LE
ts_results$GlobalLE=map2(ts_results$modelresultsbest, ts_results$jacobians, LEshift)
ts_d$JLEmin=map_dbl(ts_results$GlobalLE, ~.x$minci) #LE lower confidence bound (*this is the LE estimate to use!*)
ts_d$JLEmean=map_dbl(ts_results$GlobalLE, ~.x$minmean) #LE mean (*this is the LE estimate to use!*)
ts_d$JLEsign=ifelse(ts_d$JLEmin>0.01, "chaotic", "not chaotic")
# length(which(ts_d$JLE>0.01))/length(which(!is.na(ts_d$JLE)))

#get power spectrum for local LE
ts_d$LocalLEpower = map(ts_d$LocalLE, getspectrum, ts="lle")

#get power for period 12 (annual)
ts_d$LocalLEpower12 = map_dbl(ts_d$LocalLE, getspectrum12, ts="lle")

#predictability of time series (abundance, growth rate, both, neither)
predthreshold=0.2
ts_d$predictable_ag=ifelse(ts_d$bestR2>predthreshold & ts_d$bestR2m>predthreshold, "ag",
                             ifelse(ts_d$bestR2>predthreshold & ts_d$bestR2m<=predthreshold, "a",
                                    ifelse(ts_d$bestR2<=predthreshold & ts_d$bestR2m>predthreshold, "g", "none")))

#proportion positive lle
ts_d$pplle=map_dbl(ts_d$LocalLE, function(d) {
  ts=d$lle
  length(which(ts>0.01))/length(which(!is.na(ts)))
})

#number of months with at least 4 lle estimates
ts_d$momin=map2_dbl(ts_d$data_rescale, ts_d$LocalLE, function(d1,d2) {
  m=d1$Month
  ts=!is.na(d2$lle)
  tab=table(m,ts)
  length(which(tab[,2]>4))
})

#month with max mean lle
ts_d$llemax_mo=map2_dbl(ts_d$data_rescale, ts_d$LocalLE, function(d1,d2) {
  df=data.frame(Month=d1$Month,lle=d2$lle)
  tab=aggregate(lle~Month, data=df, FUN=mean, na.rm=T)
  tab$Month[which.max(tab$lle)]
})

#lle range, monthly median
ts_d$llerange=map2_dbl(ts_d$data_rescale, ts_d$LocalLE, function(d1,d2) {
  df=data.frame(Month=d1$Month,lle=d2$lle)
  tab=aggregate(lle~Month, data=df, FUN=median, na.rm=T)
  max(tab$lle)-min(tab$lle)
})

#lle range, monthly median, scaled
ts_d$llerange_rel=map2_dbl(ts_d$data_rescale, ts_d$LocalLE, function(d1,d2) {
  df=data.frame(Month=d1$Month,lle=d2$lle)
  df$lle=(df$lle-mean(df$lle,na.rm=T))/sd(df$lle,na.rm=T)
  tab=aggregate(lle~Month, data=df, FUN=median, na.rm=T)
  max(tab$lle)-min(tab$lle)
})

#residual
ts_d$residuals=map(ts_results$modelresultsbest, function(x) {
  res=x$resultsdf$Obs_abund-x$resultsdf$Pred_abund
  data.frame(res=res,logabsres=log(abs(res)))
})

#VER
ts_d$instability=map(ts_results$jacobians, function(x) {
  res=apply((x[1,,,drop=F])^2,3,sum)
})

#update sites used
# load("./output/ts_results_updated.Rdata")
# ts_d=left_join(select(dnestfilt, Site, Name, Level,data_rescale), select(ts_d,-data_rescale))
# ts_results=left_join(select(dnestfilt, Site, Name, Level), ts_results)

#### Exclude Wadden Sea extra species ####
ws=filter(ts_d, Site=="Wadden Sea")
exclude=c("Acartia_sp_c","Acartia_sp_nau","C_hamatus_c","C_hamatus_nau",
          "Oithona_sp_c","Oithona_sp_nau","T_longicornis_c","T_longicornis_naup")
ts_d2=filter(ts_d, !(Name %in% exclude))
ts_results2=filter(ts_results, !(Name %in% exclude))

#### Export Results ####

#save results

#updated results
save(ts_d, ts_d2, ts_results, ts_results2, file = "./output/ts_results_updated3.Rdata")
#load("./output/ts_results_updated2.Rdata")

#global LE and other site-level values
exportres1=ts_d2 %>% select(Site, Name, Level, R2abund=bestR2, R2gr=bestR2m, predictable_ag, modelform, E, tau, theta, JLEmin, JLEmean, JLEsign, 
                            LocalLEpower12, pplle, momin, llemax_mo, llerange, llerange_rel)
#annual LE
exportres2=ts_d2 %>% select(Site, Name, Level, AnnualLE) %>% unnest()
#local LE
exportres3=ts_d2 %>% select(Site, Name, Level, data_rescale, LocalLE, residuals, instability) %>% unnest(data_rescale, LocalLE, residuals, instability) %>% 
  select(Site, Name, Level, Year, Month, lle, det, ver=instability, res, logabsres)
#local LE power spectrum
exportres4=ts_d2 %>% select(Site, Name, Level, LocalLEpower) %>% unnest()

write.csv(exportres1, "./output/results_global_updated.csv", row.names = F)
write.csv(exportres2, "./output/results_annual_updated.csv", row.names = F)
write.csv(exportres3, "./output/results_local_updated.csv", row.names = F)
write.csv(exportres4, "./output/results_localpower_updated.csv", row.names = F)

#rerun
# load("./output/ts_results_updated2.Rdata")
# rerunsites=c("Lake Dora (FL)","Lake Eustis (FL)","Lake Griffin (FL)","Lake Harris (FL)")
# rerunindex=which(ts_d$Site %in% rerunsites)
# ts_d_original=ts_d
# ts_results_original=ts_results
# 
# ts_d1=read.csv("./data/alldata_filtered.csv")
# ts_d1=ts_d1 %>% group_by(Site, Name, Level) %>% nest(.key="data_rescale") %>% 
#   mutate(data_rescale=map(data_rescale, as.data.frame))
# 
# ts_d$data_rescale=ts_d1$data_rescale
# 
# #fit models
# #models 1 and 2 are dropped because the results of 1&3 and 2&5 are identical
# ts_results$modelresults3[rerunindex]=map(ts_d$data_rescale[rerunindex], smap_model_options, y="AbundanceScale", model=3, seaspred=F) #fd-ut
# ts_results$modelresults4[rerunindex]=map(ts_d$data_rescale[rerunindex], smap_model_options, y="AbundanceScale", model=4, seaspred=F) #gr-ut
# ts_results$modelresults5[rerunindex]=map(ts_d$data_rescale[rerunindex], smap_model_options, y="AbundanceScale", model=5, seaspred=F) #gr-log
# #rerun rest
# 
# ts_d$JLEmean[rerunindex]
# ts_d_original$JLEmean[rerunindex]
