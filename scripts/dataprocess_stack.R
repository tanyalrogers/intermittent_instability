#stack together raw data files in datainterm/
#some processing, filtering
#files are output to data/ or data_public/

filenames <- list.files("datainterm/","*.comm.csv",full.names = T)
data <- lapply(filenames, read.csv, stringsAsFactors = F)
alldata <- do.call(rbind, data)
#write.csv(alldata,"data/alldata.csv", row.names = F)

library(dplyr)
library(tidyr)
library(purrr)

dnest=alldata %>% group_by(Site,Name,Level) %>% nest()

#trim off and remove long periods with no dynamics
tsprocess=function(data) {
  AbundancePos=ifelse(data$Abundance<=0,NA,data$Abundance)
  for(i in 1:(length(AbundancePos)-12)) {
    if(all(is.na(AbundancePos[i:(i+12)]))) {
      data$Abundance[i:(i+12)]<-NA
    }
  }  
  #data$AbundanceScale=data$AbundancePos/sd(data$AbundancePos, na.rm = T)
  #data$AbundanceLogScale=log(data$AbundanceScale)
  data=data[(first(which(!is.na(data$Abundance)))):(last(which(!is.na(data$Abundance)))),]

  return(as.data.frame(data))
}
dnest$data2=map(dnest$data,tsprocess)

#figure out what to do with zeros
# dnest$mina=map_dbl(dnest$data2, ~min(.x$Abundance[which(.x$Abundance>0)], na.rm=T))
# dnest$minb=map_dbl(dnest$data2, ~min(.x$Abundance, na.rm=T))
# dnest$minc=round(dnest$mina-dnest$minb,2)
# dnest$maxa=map_dbl(dnest$data2, ~max(.x$Abundance, na.rm=T))
# dnest=dnest %>% arrange(desc(minc))

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
rescale=function(data, diag=F) {
  mind=min(data$Abundance, na.rm=T)
  minnz=min(data$Abundance[data$Abundance>0], na.rm=T)
  if(mind>0) {
    data$AbundanceScale=data$Abundance/sd(data$Abundance, na.rm=T)
    case=1 #no zeros, nothing added
  } else {
    if(all(is.wholenumber(data$Abundance), na.rm=T) & minnz<100) {
      data$AbundanceScale=(data$Abundance+1)/sd(data$Abundance, na.rm=T)
      case=2 #zeros, integers, plus 1, unless minnz is >100
    } else {
      data$AbundanceScale=(data$Abundance+minnz)/sd(data$Abundance, na.rm=T)
      case=3 #zeros, non-integers, plus min non-zero value
    }
  }
  data$AbundanceLogScale=log(data$AbundanceScale)
  if(diag) return(case) 
  else return(as.data.frame(data))
}
dnest$data_rescale=map(dnest$data2, rescale)
dnest$data_rescale_case=map_dbl(dnest$data2, rescale, diag=T)


#filtering metrics
dnest$tslength=map_int(dnest$data_rescale,nrow)
dnest$ndatapoints=map_int(dnest$data_rescale, function(d) {
  length(which(!is.na(d$Abundance)))
  })
dnest$longestrun=map_int(dnest$data_rescale, function(d) {
  len=rle(!is.na(d$Abundance))
  max(len$lengths[len$values==TRUE])
})
dnest$propmissing=map_dbl(dnest$data_rescale, function(d) {
  ts=d$Abundance #[(first(which(!is.na(d$AbundancePos)))):(last(which(!is.na(d$AbundancePos))))]
  length(which(is.na(ts)))/length(ts)
})
#min run 10 allows testing of E/tau: 1/1,2/1,1/2
minrun=10
dnest$nruns10=map_int(dnest$data_rescale, function(d) {
  len=rle(!is.na(d$Abundance))
  length(which(len$lengths[len$values==TRUE]>=minrun))
})
#effective sample size for E=2
dnest$effectivenE2=map_dbl(dnest$data_rescale, function(d) {
  len=rle(!is.na(d$Abundance))
  len2=len$lengths[len$values==TRUE]
  if(any(len2>2)) {sum(len2[len2>2]-2)}
  else {0}
})

#filter
dnestfilt=dnest %>% 
  filter(longestrun>=24) %>% 
  filter(effectivenE2>=40) 
  #filter(!(longestrun<12 & runs_gtmin<2)) %>% 
  #arrange(effectivenE2)

alldata2=unnest(dnestfilt, data_rescale) %>% select(Site,Name,Level,Year:AbundanceLogScale)
write.csv(alldata2,"data/alldata_filtered.csv", row.names = F)

seriesmeta=select(dnestfilt,Site,Name,Level,data_rescale_case:effectivenE2)
write.csv(seriesmeta,"data/seriesmetadata.csv", row.names = F)

#####################

# env data (temperature, photoperiod, productivity))
filenames <- list.files("datainterm/","*.abiotic.csv",full.names = T)
data2 <- lapply(filenames, read.csv, stringsAsFactors = F)
tempdata <- lapply(data2,FUN = function(d){d[,c("Site","Year","Month","Temp")]})
tempdata2 <- do.call(rbind, tempdata)
tempdata3 <- filter(tempdata2, Site %in% unique(dnestfilt$Site))
write.csv(tempdata3,"data/temperaturedata.csv", row.names = F)

# site level data (monthly and total)
library(geosphere)
library(lubridate)

meta=read.csv("data/sitemetadata.csv", stringsAsFactors = F)
dates=data.frame(Month=1:12,Day=yday(dmy(paste(15,1:12,2021,sep = "-"))))
sites=unique(meta$Site)

monthlyphper=expand.grid(Site=sites, Month=1:12) %>% 
  left_join(dates) %>% 
  left_join(select(meta,Site,Lat))
monthlyphper$Photoperiod=daylength(monthlyphper$Lat,monthlyphper$Day)

monthlytemps=aggregate(Temp~Site*Month, tempdata2, mean, na.rm=T) %>% arrange(Site,Month)
monthlydata=monthlytemps %>% right_join(select(monthlyphper,Site,Month,Photoperiod)) %>% 
  arrange(Site, Month)

annualdata=monthlydata %>% group_by(Site) %>% 
  summarise(Temp_mean=mean(Temp, na.rm = T),Temp_range=max(Temp)-min(Temp),
            Photoperiod_range=max(Photoperiod)-min(Photoperiod)) %>% 
  ungroup() %>% as.data.frame()
meta=select(meta,1:8) %>% left_join(annualdata)

chlmean=alldata %>% filter(Name=="Chla") %>% group_by(Site) %>% 
  summarise(Chla_mean=mean(Abundance, na.rm = T))
meta=meta %>% left_join(chlmean)
meta$Chla_mean[meta$Chla_status=="cell density"]=NA

#write.csv(monthlydata,"data/monthly_temp_photoper_data.csv", row.names = F)
write.csv(meta, "data/sitemetadata.csv", row.names = F)

#Remove restricted data
#################
keep=meta$Site[meta$Redistributable=="yes"]

alldata_public=filter(alldata2,Site %in% keep)
write.csv(alldata_public,"data_public/alldata_filtered.csv", row.names = F)

tempdata_public=filter(tempdata3,Site %in% keep)
write.csv(tempdata_public,"data_public/temperaturedata.csv", row.names = F)

write.csv(seriesmeta,"data_public/seriesmetadata.csv", row.names = F)
write.csv(meta, "data_public/sitemetadata.csv", row.names = F)

#Diagnostics
#################

table(dnest$Site,dnest$Level)
table(dnestfilt$Site,dnestfilt$Level)

table(dnestfilt$Level)
length(unique(dnestfilt$Site))

test=dnest$data2[[2]]$AbundancePos
len=rle(!is.na(test))

freq=table(dnestfilt$longestrun)
cumfreq0 = c(0, cumsum(freq)) 
x=as.numeric(names(cumfreq0));x[1]=0
plot(x,cumfreq0/max(cumfreq0), type = "o",
     xlab="max run length", ylab="cumulative frequency")
abline(v=c(10,20,30,40),lty=2, col="gray")

########

plotts=function(dnest,table, index, y) {
  ts=dnest[index,table][[1]][[1]]
  plot(ts[,y], type="o", ylab=y)
}

plotts(dnestfilt,"data_rescale",249,"AbundanceLogScale")
