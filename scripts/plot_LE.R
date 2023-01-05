# Plots results

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(purrr)
library(ggmosaic)
library(lme4)

source("scripts/ggplot_themes.R")

#load results
rset="_updated"
global=read.csv(paste0("output/results_global",rset,".csv"), stringsAsFactors = F)
local=read.csv(paste0("output/results_local",rset,".csv"), stringsAsFactors = F)
localpower=read.csv(paste0("./output/results_localpower",rset,".csv"), stringsAsFactors = F)
meta=read.csv("data/sitemetadata.csv", stringsAsFactors = F)
seriesmeta=read.csv("data/seriesmetadata.csv", stringsAsFactors = F)
seriests=read.csv("data/alldata_filtered.csv", stringsAsFactors = F)
global=left_join(global,meta) %>% left_join(seriesmeta)
local=left_join(local,seriests)

#set site order
global=arrange(global, desc(abs(Lat)))
siteorder=unique(global$Site)
global$Site=factor(global$Site, levels = siteorder)
localpower$Site=factor(localpower$Site, levels = siteorder)
local$Site=factor(local$Site, levels = siteorder)
meta$Site=factor(meta$Site, levels = siteorder)

#compute CV for each series
global=local %>% group_by(Site,Name,Level) %>% 
  summarise(AbundSD=sd(Abundance, na.rm=T),
            AbundCV=AbundSD/mean(Abundance, na.rm=T)) %>% right_join(global)

#recode functional groups and set order
global$Level=factor(global$Level, levels = c("SP","FG","Total"))
global$Level2=recode(global$Level, SP="Species", FG="Functional Group", Total="Trophic Level")
global$Level3=recode(global$Level, Total="TL")

# data frames for plotting ####

#series with year round sampling (for seasonality)
#momin=number of months with at least 4 lle estimates, must be at least 10 months
globalseas=global %>% filter(momin>9) %>% 
  mutate(lnChla=log(Chla_mean), absLat=abs(Lat))

#sites with all 3 levels and at least 6 total series
sitelev=names(which(apply(table(globalseas$Site,globalseas$Level), 1, function(x) {length(which(x>0))==3})==T))
sitelev2=names(which(apply(table(globalseas$Site,globalseas$Level), 1, function(x) {sum(x)>6})==T))
globallevels=filter(globalseas, Site %in% intersect(sitelev, sitelev2))

#site level averages (all global quantities)
globalseassummary=globalseas %>% group_by(Site,Type,Level2,Level3) %>% 
  summarise_if(is.numeric,mean)

#site level medians  (all global quantities)
globalsummarymed=global %>% group_by(Site,Type,Level2,Level3) %>% 
  summarise_if(is.numeric,median)

#prop chaotic series
globalseassummary2=globalseas %>% group_by(Site,Type,Level2) %>% 
  summarise(pchaotic=sum(JLEsign=="chaotic")/n())
globalseassummary = left_join(globalseassummary,globalseassummary2)

##################

#Count of series by level
table(globalseas$Level)

#Effect of taxonomic resolution on stability ####

#Overall chaos classification by level
as.data.frame(table(Level=globalseas$Level,JLEsign=globalseas$JLEsign)) %>% spread(JLEsign,Freq) %>% 
  mutate(propchaotic=round(chaotic/(chaotic+`not chaotic`),3),
         total=(chaotic+`not chaotic`))
#pp eigenvalues (average)
aggregate(pplle~Level, data=globalseas, FUN=mean)

#Chaos classification and pp eigenvalues
length(which(globalseas$pplle>0 & globalseas$JLEsign=="not chaotic"))/length(which(globalseas$JLEsign=="not chaotic"))
length(which(globalseas$pplle>0.25 & globalseas$JLEsign=="not chaotic"))/length(which(globalseas$JLEsign=="not chaotic"))
min(globalseas$pplle[globalseas$JLEsign=="chaotic"])
max(globalseas$pplle[globalseas$JLEsign=="chaotic"])

#Dominant period of power spectrum, chaos class, pplle
localpowerdom=localpower %>% group_by(Site, Name, Level) %>% 
  filter(Power==max(Power)) %>% ungroup() %>% 
  mutate(Seasonal=ifelse(Period>=11.5 & Period<=12.5,"yes","no"),
         Level=factor(Level, levels = c("SP","FG","Total")))
globalseas2=left_join(globalseas,localpowerdom)
globalseas2$ppllepos=ifelse(globalseas2$pplle>0,"yes","no")
globalseas2$JLEsign2=ifelse(globalseas2$JLEsign=="chaotic","chaotic",
                           ifelse(globalseas2$ppllepos=="yes",
                                  "not chaotic with\nlocal instability","not chaotic,\nalways stable"))
globalseas2=globalseas2 %>% mutate(Seasonal=factor(Seasonal,levels = c("yes","no")))
#overall prop seasonal
as.data.frame(table(Level=globalseas2$Level,Seasonal=globalseas2$Seasonal)) %>% spread(Seasonal,Freq) %>% 
  mutate(propseasonal=round(yes/(yes+no),3),
         total=(yes+no))
#local instability given non-chaotic
as.data.frame(table(Level=globalseas2$Level,JLEsign2=globalseas2$JLEsign2)) %>% spread(JLEsign2,Freq) %>% 
  mutate(total=(chaotic+`not chaotic with\nlocal instability`+`not chaotic,\nalways stable`),
         propchaotic=round(chaotic/total,3),
         propinst=round(`not chaotic with\nlocal instability`/(`not chaotic with\nlocal instability`+`not chaotic,\nalways stable`),3))

#mosaic plot (all series)
ggplot(globalseas2) +
  facet_grid(.~Level2, scales = "free_x") +
  geom_mosaic(aes(x=product(Seasonal, JLEsign2), fill = Seasonal),show.legend = F) +  
  labs(y="Seasonal period dominant", x="Stability") +
  classic + removefacetbackground +
  theme(panel.border = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  xlabvert + scale_fill_brewer(palette = "Set2") +
  scale_y_productlist(expand = expand_scale(mult = 0)) +
  scale_x_productlist(expand = expand_scale(mult = 0))
ggsave("figures/stability_mosaic.png", width = 6, height = 3, units = "in")
ggsave("figures/Figure_1.pdf", width = 6, height = 3, units = "in")

#mosaic plot version (individual sites)
serpersite=as.data.frame(table(Site=droplevels(globallevels$Site)))
ggplot(globallevels) +
  facet_wrap(Site~., ncol=4, scales = "free_x") +
  geom_mosaic(aes(x=product(JLEsign, Level3), fill = JLEsign),show.legend = F) +
  geom_text(data=serpersite, aes(label=Freq,x=0.1,y=0.9),size=3) +
  labs(y="Stability", x="Taxonomic resolution") +
  classic + removefacetbackground +
  theme(panel.border = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
  xlabvert + scale_fill_brewer(palette = "Paired", direction = -1) +
  scale_y_productlist(expand = expand_scale(mult = c(0, 0))) +
  scale_x_productlist(expand = expand_scale(mult = c(0, 0))) +
  theme(panel.spacing.x = unit(10,"pt"),axis.text.x = element_text(size=7))
ggsave("figures/classification_aggeffects_mosaic.png", width = 8, height = 4, units = "in")

# Effect of taxonomic resolution on CV, R2, E, and LE magnitude ####

#site level medians (uses global, not global seas)
CV=ggplot(globalsummarymed, aes(x=Level2, y=AbundCV, group=Site, color=Type)) +
  geom_line() +
  geom_point(alpha=0.5) + 
  labs(y="CV for abundance", x="Taxonomic resolution", color=NULL) +
  classic + scale_color_brewer(palette = "Set1", direction = -1) +
  theme(legend.position = c(1,1), legend.justification = c(1,1), legend.background = element_blank())
R2=ggplot(globalsummarymed, aes(x=Level2, y=R2abund, group=Site, color=Type)) +
  geom_line(show.legend = F) +
  geom_point(alpha=0.5, show.legend = F) + 
  labs(y=expression(R^2~'for'~abundance), x="Taxonomic resolution", color=NULL) +
  classic + scale_color_brewer(palette = "Set1", direction = -1)
plot_grid(CV,R2,labels = "auto", align = "hv", ncol = 2)
ggsave("figures/CV_R2.png", width = 7, height = 3, units = "in")
ggsave("figures/Figure_5.pdf", width = 7, height = 3, units = "in")

#embedding dimension
ggplot(globalsummarymed, aes(x=Level2, y=E, group=Site, color=Type)) +
  geom_line() +
  geom_point(alpha=0.5) + 
  geom_hline(yintercept=0, lty=2, color="gray") +
  labs(y="E", x="Taxonomic resolution") +
  classic + scale_color_brewer(palette = "Set1", direction = -1)
#LE value
ggplot(globalsummarymed, aes(x=Level2, y=JLEmean, group=Site, color=Type)) +
  geom_line() +
  geom_point(alpha=0.5) + 
  geom_hline(yintercept=0, lty=2, color="gray") +
  labs(y="LE", x="Taxonomic resolution") +
  classic + scale_color_brewer(palette = "Set1", direction = -1)

#faceted by site, showing all points
exclude=c("Blelham Tarn (UK)","Esthwaite Water (UK)") #only 1 level
#R2
ggplot(filter(global,!(Site %in% exclude)), aes(x=Level3, y=R2abund)) +
  facet_wrap(Site~., ncol=4) +
  geom_hline(yintercept=0, lty=2, color="gray") +
  geom_jitter(alpha=0.6,width = 0.05, height = 0, size=1, color="gray") +
  geom_point(data=filter(globalsummarymed,!(Site %in% exclude)), size=1) +
  geom_line(data=filter(globalsummarymed,!(Site %in% exclude)), aes(group=Site)) +
  labs(y=expression(R^2~'for'~abundance), x="Taxonomic resolution") +
  classic + removefacetbackground
#CV
ggplot(filter(global,!(Site %in% exclude)), aes(x=Level3, y=AbundCV)) +
  facet_wrap(Site~., ncol=4) +
  geom_hline(yintercept=0, lty=2, color="gray") +
  geom_jitter(alpha=0.6,width = 0.05, height = 0, size=1, color="gray") +
  geom_point(data=filter(globalsummarymed,!(Site %in% exclude)), size=1) +
  geom_line(data=filter(globalsummarymed,!(Site %in% exclude)), aes(group=Site)) +
  labs(y="CV for abundance", x="Taxonomic resolution") +
  classic + removefacetbackground

#sites with all 3 levels
ggplot(globallevels, aes(x=Level3, y=R2abund)) +
  facet_wrap(Site~., ncol=3) +
  geom_hline(yintercept=0, lty=2, color="gray") +
  geom_jitter(aes(color=JLEsign),alpha=0.5,width = 0.1, height = 0) + 
  geom_line(data=filter(globalsummarymed, Site %in% intersect(sitelev, sitelev2)), aes(group=Site)) +
  labs(y="R2 for abundance", x="Taxonomic resolution") +
  classic + removefacetbackground +
  scale_color_brewer(palette = "Paired", direction = -1)
ggplot(globallevels, aes(x=Level3, y=E)) +
  facet_wrap(Site~., ncol=3) +
  geom_hline(yintercept=0, lty=2, color="gray") +
  geom_jitter(aes(color=JLEsign),alpha=0.5,width = 0.1, height = 0) + 
  geom_line(data=filter(globalsummarymed, Site %in% intersect(sitelev, sitelev2)), aes(group=Site)) +
  labs(y="R2 for abundance", x="Taxonomic resolution") +
  classic + removefacetbackground +
  scale_color_brewer(palette = "Paired", direction = -1)
ggplot(globallevels, aes(x=Level3, y=JLEmean)) +
  facet_wrap(Site~., ncol=3) +
  geom_hline(yintercept=0, lty=2, color="gray") +
  geom_jitter(aes(color=JLEsign),alpha=0.5,width = 0.1, height = 0) + 
  geom_line(data=filter(globalsummarymed, Site %in% intersect(sitelev, sitelev2)), aes(group=Site)) +
  labs(y="R2 for abundance", x="Taxonomic resolution") +
  classic + removefacetbackground +
  scale_color_brewer(palette = "Paired", direction = -1)

#Distribution of E, tau, theta by level of aggregation ####

ggplot(global, aes(x=factor(E),fill=JLEsign)) +
  facet_wrap(~Level2, scales = "free_y") +
  geom_bar(position = "stack") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  classic + removefacetbackground + labs(x="E")
ggplot(global, aes(x=factor(tau), fill=JLEsign)) +
  facet_wrap(~Level2, scales = "free_y") +
  geom_bar(position = "stack") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  classic + removefacetbackground + labs(x="tau")
ggplot(global, aes(x=factor(theta),fill=JLEsign)) +
  facet_wrap(~Level2, scales = "free_y") +
  geom_bar(position = "stack") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  classic + removefacetbackground + labs(x="theta")

#LE vs. pp eigenvalues vs. seasonality of eigenvalues ####

# select individual time series to circle and plot
# test1=globalseas %>% filter(Level=="Total")
# test2=globalseas %>% filter(Level=="FG")
# test3=globalseas %>% filter(Level=="SP")
# plotly::plot_ly(test1, x= ~pplle, y= ~llerange, name = ~Name, text=~Site)
# plotly::plot_ly(test2, x= ~pplle, y= ~llerange, name = ~JLEsign, text=~paste(Site,Name))
# plotly::plot_ly(test3, x= ~pplle, y= ~llerange, name = ~JLEsign, text=~paste(Site,Name))
# plotly::plot_ly(test3, x= ~JLEmean, y= ~llerange, name = ~JLEsign, text=~paste(Site,Name))
globalsub=globalseas %>% filter((Site=="Lake Geneva" & Name=="Cyclops prealpinus") | 
                                  (Site=="Lake Muggelsee" & Name=="Leptodora kindtii") |
                                  (Site=="Loch Leven" & Name=="Eudiaptomus.gracilis") |
                                  (Site=="Port Erin Bay" & Name=="Paracalanus.parvus") |
                                  (Site=="Port Erin Bay" & Name=="Biddulphia.mobiliensis") |
                                  (Site=="Narragansett Bay" & Name=="Acartia.tonsa")) %>% arrange(desc(llerange)) %>% ungroup()
localsub=globalsub %>% left_join(local, by=c("Site", "Name")) %>% 
  mutate(Name=factor(Name, levels = unique(Name)))
localpowersub=globalsub %>% left_join(localpower, by=c("Site","Name") )%>% 
  mutate(Name=factor(Name, levels = unique(Name))) %>% 
  select(Site, Name, Frequency, Power) 
temp=as.data.frame(unique.data.frame(select(localpowersub,Site,Name)))
temp$Frequency=1/12
temp$Power=globalsub$LocalLEpower12
localpowersub=rbind(localpowersub, temp)
globalsub=globalsub %>% mutate(Name=factor(Name, levels = unique(Name)), Number=1:nrow(globalsub))
localsub=localsub %>% group_by(Name) %>% mutate(Year_rel=(Year-min(Year))/(max(Year)-min(Year)))

#circling individual timeseries
#LE vs seaslle
gvgv12a=ggplot(globalseas, aes(x=JLEmean,y=llerange,color=pplle)) +
  facet_grid(.~Level2, scales = "free_y") +
  annotate("rect", xmin=0.01,xmax=Inf,ymin=-Inf,ymax=Inf,fill="gray90") +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_point(aes(color=pplle, shape=Type), size=2, alpha=0.7) + 
  geom_point(data=globalsub, pch=21, size=4, color="firebrick") +
  geom_text(data=globalsub, aes(label=Number), nudge_x = 0.15, pch=21, size=4, color="firebrick") +
  classic + removefacetbackground + scale_color_viridis_c() + #scale_color_distiller(palette = "Spectral", direction = 1) + legalpha +
  labs(x="LE", y="Seasonality (annual)\nof local eigenvalue", color="Local eigenvalue\n(proportion\npositive)") +
  guides(shape="none")
#pplle vs seaslle
gvgv12b=ggplot(globalseas, aes(y=llerange,x=pplle)) +
  facet_grid(.~Level2, scales = "free_y") +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_vline(xintercept = c(0,1), lty=2, color="gray") +
  geom_point(aes(color=JLEsign, shape=Type), size=2, alpha=0.7) + 
  geom_point(data=globalsub, pch=21, size=4, color="firebrick") +
  geom_text(data=globalsub, aes(label=Number), nudge_x = 0.08, pch=21, size=4, color="firebrick") +
  classic + removefacetbackground + scale_color_brewer(palette = "Paired", direction = -1) + legalpha +
  labs(color="LE\nclassification", y="Seasonality (annual)\nof local eigenvalue", x="Local eigenvalue (proportion positive)", shape="System")
plot_grid(gvgv12a,gvgv12b,nrow = 2,labels = "auto",align = "hv")
ggsave("figures/GLEvLLEalt_llerange_update.png", width = 8, height = 5, units = "in")
ggsave("figures/Figure_2.pdf", width = 8, height = 5, units = "in")

# Plots of individual time series ####

#Log Abundance
abplot=ggplot(localsub, aes(x=Month, y=AbundanceLogScale)) +
  facet_grid(Name~., scales="free_y")+
  geom_line(aes(group=Year), alpha=0.2, show.legend = F,lwd=0.5, color="firebrick") +
  geom_point(size=1, alpha=0.2, show.legend = F, color="firebrick") +
  labs(y="Log Scaled Abundance",x="Month") + 
  scale_x_continuous(breaks=1:12) + 
  classic + removefacets
#LE and monthly eigenvalue
leplot=ggplot(localsub, aes(x=Month, y=lle)) +
  facet_grid(Name~.)+
  annotate("rect", ymin=0.01,ymax=Inf,xmin=-Inf,xmax=Inf,fill="gray90",lwd=0.25) +
  geom_line(aes(group=Year), alpha=0.2, show.legend = F,lwd=0.5, color="firebrick") +
  geom_point(size=1, alpha=0.2, show.legend = F, color="firebrick") +
  #geom_boxplot(data=filter(zooann, !is.na(AnnualLE)),aes(x=0,y=AnnualLE), show.legend = F, outlier.size = 0.75, lwd=0.25) +
  geom_point(data=globalsub, aes(x=0,y=JLEmean), show.legend = F, size=3, pch=18, color="firebrick") +
  labs(y="Eigenvalue",x="   Month") + 
  scale_x_continuous(breaks=0:12,labels = c("LE",1:12)) + classic + removefacets
#Eigenvalue power spectrum
specplot=ggplot(localpowersub, aes(x=Frequency, y=Power)) +
  facet_grid(Name~.) +
  geom_vline(xintercept = 1/12, color="gray80", lwd=1, lty=1) +
  geom_vline(xintercept = 1/6, color="gray80", lwd=1, lty=1) +
  geom_line(lwd=0.5, color="firebrick", show.legend = F) + 
  geom_text(data=globalsub, aes(x=0.5,y=1,label=Number)) +
  classic + removefacets
plot_grid(abplot,leplot,specplot,nrow = 1,labels = "auto",align = "hv", rel_widths = c(0.9,1,0.9))
ggsave("figures/LEpwrsample.png", width = 10, height = 7.5, units = "in")
ggsave("figures/Figure_3.pdf", width = 10, height = 7.5, units = "in")

#color lines by year
ggplot(localsub, aes(x=Month, y=AbundanceLogScale, color=Year_rel)) +
  facet_grid(Name~., scales="free_y")+
  geom_line(aes(group=Year), alpha=0.2, show.legend = T,lwd=0.5) +
  geom_point(size=1, alpha=0.2, show.legend = F) +
  labs(y="Log Scaled Abundance",x="Month", color="Year") + 
  scale_x_continuous(breaks=1:12) + 
  classic + removefacets + 
  theme(legend.position = "top", legend.key.height = unit(4, "mm"), legend.box.margin = margin(-5,0,-10,0), legend.title = element_text(vjust = 1)) + 
  scale_color_viridis_c(option = "A",begin = 0.1, end=0.9, breaks=c(0,1),labels=c("min","max"),direction = -1)
ggplot(localsub, aes(x=Month, y=lle, color=Year_rel)) +
  facet_grid(Name~.)+
  annotate("rect", ymin=0.01,ymax=Inf,xmin=-Inf,xmax=Inf,fill="gray90",lwd=0.25) +
  geom_line(aes(group=Year), alpha=0.2, show.legend = T,lwd=0.5) +
  geom_point(size=1, alpha=0.2, show.legend = F) +
  geom_point(data=globalsub, aes(x=0,y=JLEmean), show.legend = F, size=3, pch=18, color="firebrick") +
  labs(y="Eigenvalue",x="   Month", color="Year") + 
  scale_x_continuous(breaks=0:12,labels = c("LE",1:12)) + 
  classic + removefacets + 
  theme(legend.position = "top", legend.key.height = unit(4, "mm"), legend.box.margin = margin(-5,0,-10,0), legend.title = element_text(vjust = 1)) + 
  scale_color_viridis_c(option = "A",begin = 0.1, end=0.9, breaks=c(0,1),labels=c("min","max"),direction = -1)

# pp eig vs. eig seasonality by name or functional group ####
ggplot(filter(globalseas,Level=="Total"), aes(y=llerange,x=pplle)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_vline(xintercept = c(0,1), lty=2, color="gray") +
  geom_point(aes(color=Name, shape=Type), size=2, alpha=0.9) + 
  classic + removefacetbackground + scale_fill_brewer(palette = "Set1") + legalpha +
  geom_point(data=filter(globalseas,Level=="Total",JLEsign=="chaotic"), pch=21, size=4, color="black") +
  labs(color="Group", y="Seasonality (annual)\nof local eigenvalue", x="Local eigenvalue (proportion positive)", shape="System")
ggplot(filter(globalseas,Level=="FG",Name %in% c("large.clad","small.clad","copepod","pred","rotifer")), aes(y=llerange,x=pplle)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_vline(xintercept = c(0,1), lty=2, color="gray") +
  geom_point(aes(color=Name, shape=Type), size=2, alpha=0.9) + 
  classic + removefacetbackground + scale_color_brewer(palette = "Set1") + legalpha +
  geom_point(data=filter(globalseas,Level=="FG",Name %in% c("large.clad","small.clad","copepod","pred","rotifer"),JLEsign=="chaotic"), pch=21, size=4, color="black") +
  labs(color="Group", y="Seasonality (annual)\nof local eigenvalue", x="Local eigenvalue (proportion positive)", shape="System")

#chaos by functional group ####
ggplot(filter(global,Level=="FG",Name %in% c("large.clad","small.clad","copepod","pred","rotifer")), aes(x=Name,fill=JLEsign)) +
  #facet_wrap(~Level2, scales = "free_y") +
  geom_bar(position = "stack") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  classic + removefacetbackground + labs(x="Functional Group")

#Season with peak instability ####
globalseas$peak=case_when(globalseas$llemax_mo %in% c(12,1,2) ~ "winter",
                          globalseas$llemax_mo %in% c(3,4,5) ~ "spring",
                          globalseas$llemax_mo %in% c(6,7,8) ~ "summer",
                          globalseas$llemax_mo %in% c(9,10,11) ~ "fall")
globalseas$peak=factor(globalseas$peak, levels=c("winter","spring","summer","fall"))

ggplot(filter(globalseas, llerange>0.25, pplle>0), aes(x=peak, fill=Site)) +
  facet_grid(Type~Level2) +
  geom_bar(position = "stack",show.legend = F) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  classic + removefacetbackground  + scale_fill_viridis_d() +
  labs(x="Season with peak instability", y="Count")
ggsave("figures/LEpeakseason.png", width = 6, height = 4, units = "in")

#individual months
ggplot(filter(globalseas, llerange>0.25, pplle>0), aes(x=factor(llemax_mo), fill=Site)) +
  facet_grid(Type~Level) +
  geom_bar(position = "stack", show.legend = F) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  classic + labs(x="Month")

#total phyto and zoo
ggplot(filter(globalseas, Level=="Total", llerange>0.25, pplle>0), aes(x=peak)) +
  facet_grid(Type~Name) +
  geom_bar(position = "stack",show.legend = F) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  classic + removefacetbackground  + scale_fill_viridis_d() +
  labs(x="Season with peak instability", y="Count")

#Distribution of species level LEs by site ####
ggplot(filter(global,Level=="SP"), aes(x=JLEmean)) +
  facet_wrap(.~Site, scales = "free") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  geom_histogram(bins=10,boundary=0) + classic + removefacetbackground

# CCF max for each series (eigenvalue vs abundance or gr) ####
localccf=full_join(local,select(globalseas, Site, Name, Level, Level2, pplle, llerange)) %>% 
  filter(llerange>0.25, !is.na(det), pplle>0)
ccfmax=function(x,y) {
  z=ccf(x,y,na.action = na.pass,plot=F,lag.max = 6)
  return(z$lag[which.max(z$acf)])
}
localccf=localccf %>% group_by(Site,Name,Level2) %>% 
  mutate(gr=c(NA,diff(AbundanceLogScale))) %>% 
  summarise(ccfmax=ccfmax(lle,AbundanceLogScale),
            ccfmaxgr=ccfmax(lle,gr),
            ccfabgr=ccfmax(AbundanceLogScale,gr)) %>% 
  left_join(global)

ggplot(localccf, aes(x=ccfmax,fill=JLEsign)) +
  facet_grid(Type~Level2, scales = "free_y") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  geom_histogram(bins=13) + 
  geom_vline(xintercept = 0) + 
  classic + removefacetbackground +
  labs(x="Lag with max correlation between local eigenvalue and log abundance\n(neg: local eigenvalue precedes abundance, pos: local eigenvalue lags abundance)",
       fill="Classification", y="Count") +
  scale_fill_brewer(palette = "Paired", direction = -1)
#ggsave("figures/CCFplots1.png", width = 6, height = 3, units = "in")

ggplot(localccf, aes(x=ccfmaxgr,fill=JLEsign)) +
  facet_grid(Type~Level2, scales = "free_y") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  geom_histogram(bins=13) + 
  geom_vline(xintercept = 0) +
  classic + removefacetbackground + 
  labs(x="Lag with max correlation between local eigenvalue and growth rate\n(neg: local eigenvalue precedes growth rate, pos: local eigenvalue lags growth rate)",
       fill="Classification", y="Count") +
  scale_fill_brewer(palette = "Paired", direction = -1)
#ggsave("figures/CCFplots2.png", width = 6, height = 3, units = "in")

ggplot(localccf, aes(x=ccfabgr,fill=JLEsign)) +
  facet_grid(Type~Level2, scales = "free_y") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  geom_histogram(bins=13) + 
  geom_vline(xintercept = 0) +
  classic + removefacetbackground + 
  labs(x="Lag with max correlation log abundance and growth rate\n(neg: abundance precedes growth rate, pos: abundance lags growth rate)",
       fill="Classification", y="Count") +
  scale_fill_brewer(palette = "Paired", direction = -1)

#not separating lake and marine
ccfgr=ggplot(localccf, aes(x=ccfmaxgr,fill=JLEsign)) +
  facet_wrap(.~Level2, scales = "free_y") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  geom_histogram(bins=13) + 
  geom_vline(xintercept = 0) +
  classic + removefacetbackground + 
  labs(x="Lag with max correlation between local eigenvalue and growth rate\n(neg: local eigenvalue precedes growth rate, pos: local eigenvalue lags growth rate)",
       fill="Classification", y="Count") +
  scale_fill_brewer(palette = "Paired", direction = -1)
ccfab=ggplot(localccf, aes(x=ccfmax,fill=JLEsign)) +
  facet_wrap(.~Level2, scales = "free_y") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  geom_histogram(bins=13) + 
  geom_vline(xintercept = 0) +
  classic + removefacetbackground + 
  labs(x="Lag with max correlation between local eigenvalue and log abundance\n(neg: local eigenvalue precedes abundance, pos: local eigenvalue lags abundance)",
       fill="Classification", y="Count") +
  scale_fill_brewer(palette = "Paired", direction = -1)

plot_grid(ccfgr,ccfab,labels = "auto", nrow=2)
ggsave("figures/CCFplots.png", width = 8, height = 5, units = "in")

#Lake-level stability vs covariates  ####

crosssiteplot=function(data,x,y,pointsize=3) {
  ylab=case_when(y=="llerange" ~ "Seasonality (annual)\nof local eigenvalue",
                 y=="llerange_rel" ~ "Relative seasonality\n(annual) of local eigenvalue",
                 y=="JLEmean" ~ "LE")
  xlab=case_when(x=="Temp_mean" ~ "Mean temperature",
                 x=="Temp_range" ~ "Temperature range",
                 x=="absLat" ~ "Latitude (abs)",
                 x=="lnChla" ~ "ln Chla")
  ggplot(data=data, aes_string(color="Level2", shape="Type", x=x, y=y)) +
    facet_grid(.~Level2) +
    geom_smooth(aes(fill=Level2, group=Level2),alpha=0.1,method="lm",se = T,show.legend = F) +
    geom_point(size=pointsize, alpha=0.5, show.legend = F) + 
    scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2") +
    classic + removefacetbackground +
    labs(x=xlab,y=ylab,shape="System",color="Taxonomic resolution")
}

crosssiteplot2=function(x,y,pointsize=3) {
  ylab=case_when(y=="llerange" ~ "Seasonality (annual)\nof local eigenvalue",
                 y=="llerange_rel" ~ "Relative seasonality\n(annual) of local eigenvalue",
                 y=="JLEmean" ~ "LE")
  xlab=case_when(x=="Temp_mean" ~ "Mean temperature",
                 x=="Temp_range" ~ "Temperature range",
                 x=="absLat" ~ "Latitude (abs)",
                 x=="lnChla" ~ "ln Chla")
  ggplot(data=globalseassummary, aes_string(color="Level2", shape="Type", x=x, y=y)) +
    facet_grid(.~Level2) +
    geom_smooth(aes(fill=Level2, group=Level2),alpha=0.2,method="lm",se = T,show.legend = F) +
    geom_point(data=globalseas, size=pointsize-1, alpha=0.2, color="gray50", show.legend = F) +
    geom_point(size=pointsize, alpha=0.9, show.legend = F) + 
    scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2") +
    classic + removefacetbackground +
    labs(x=xlab,y=ylab,shape="System",color="Taxonomic resolution")
}

#Individual series

#Absolute seasonality

cs1=crosssiteplot(globalseas,x="Temp_mean", y="llerange")
cs2=crosssiteplot(globalseas,x="Temp_range", y="llerange")
cs3=crosssiteplot(globalseas,x="absLat", y="llerange")
cs4=crosssiteplot(globalseas,x="lnChla", y="llerange")
plot_grid(cs1,cs2,cs3,cs4,nrow = 2, labels = "auto")
#ggsave("figures/crosssite_seas_indiv.png", width = 9, height = 6, units = "in")

#Relative seasonality
cs1=crosssiteplot(globalseas,x="Temp_mean", y="llerange_rel")
cs2=crosssiteplot(globalseas,x="Temp_range", y="llerange_rel")
cs3=crosssiteplot(globalseas,x="absLat", y="llerange_rel")
cs4=crosssiteplot(globalseas,x="lnChla", y="llerange_rel")
plot_grid(cs1,cs2,cs3,cs4,nrow = 2, labels = "auto")
#ggsave("figures/crosssite_relseas_indiv.png", width = 9, height = 6, units = "in")

#GLE
cs1=crosssiteplot(globalseas,x="Temp_mean", y="JLEmean")
cs2=crosssiteplot(globalseas,x="Temp_range", y="JLEmean")
cs3=crosssiteplot(globalseas,x="absLat", y="JLEmean")
cs4=crosssiteplot(globalseas,x="lnChla", y="JLEmean")
plot_grid(cs1,cs2,cs3,cs4,nrow = 2, labels = "auto")
#ggsave("figures/crosssite_gle_indiv.png", width = 9, height = 6, units = "in")

#Site means

#Absolute seasonality
cs1=crosssiteplot(globalseassummary,x="Temp_mean", y="llerange")
cs2=crosssiteplot(globalseassummary,x="Temp_range", y="llerange")
cs3=crosssiteplot(globalseassummary,x="absLat", y="llerange")
cs4=crosssiteplot(globalseassummary,x="lnChla", y="llerange")
plot_grid(cs1,cs2,cs3,cs4,nrow = 2, labels = "auto")
#ggsave("figures/crosssite_seas_mean.png", width = 9, height = 6, units = "in")

#Relative seasonality
cs1=crosssiteplot(globalseassummary,x="Temp_mean", y="llerange_rel")
cs2=crosssiteplot(globalseassummary,x="Temp_range", y="llerange_rel")
cs3=crosssiteplot(globalseassummary,x="absLat", y="llerange_rel")
cs4=crosssiteplot(globalseassummary,x="lnChla", y="llerange_rel")
plot_grid(cs1,cs2,cs3,cs4,nrow = 2, labels = "auto")
#ggsave("figures/crosssite_relseas_mean.png", width = 9, height = 6, units = "in")

#GLE
cs1=crosssiteplot(globalseassummary,x="Temp_mean", y="JLEmean")
cs2=crosssiteplot(globalseassummary,x="Temp_range", y="JLEmean")
cs3=crosssiteplot(globalseassummary,x="absLat", y="JLEmean")
cs4=crosssiteplot(globalseassummary,x="lnChla", y="JLEmean")
plot_grid(cs1,cs2,cs3,cs4,nrow = 2, labels = "auto")
#ggsave("figures/crosssite_gle_mean.png", width = 9, height = 6, units = "in")

#All response variables, individual series
# cs11=crosssiteplot(globalseas,x="Temp_mean", y="llerange", pointsize = 2)
# cs21=crosssiteplot(globalseas,x="Temp_range", y="llerange", pointsize = 2)
# cs12=crosssiteplot(globalseas,x="Temp_mean", y="llerange_rel", pointsize = 2)
# cs22=crosssiteplot(globalseas,x="Temp_range", y="llerange_rel", pointsize = 2)
# cs13=crosssiteplot(globalseas,x="Temp_mean", y="JLEmean", pointsize = 2)
# cs23=crosssiteplot(globalseas,x="Temp_range", y="JLEmean", pointsize = 2)
# plot_grid(cs11,cs21,cs12,cs22,cs13,cs23, nrow = 3, labels = "auto")
# ggsave("figures/crosssite_indiv.png", width = 9, height = 8, units = "in")

#All response variables, individual series
cs11=crosssiteplot2(x="Temp_mean", y="llerange", pointsize = 2)
cs21=crosssiteplot2(x="Temp_range", y="llerange", pointsize = 2)
cs12=crosssiteplot2(x="Temp_mean", y="llerange_rel", pointsize = 2)
cs22=crosssiteplot2(x="Temp_range", y="llerange_rel", pointsize = 2)
cs13=crosssiteplot2(x="Temp_mean", y="JLEmean", pointsize = 2)
cs23=crosssiteplot2(x="Temp_range", y="JLEmean", pointsize = 2)
plot_grid(cs11,cs21,cs12,cs22,cs13,cs23, nrow = 3, labels = "auto")
ggsave("figures/crosssite_indiv.png", width = 9, height = 8, units = "in")
ggsave("figures/Figure_6.pdf", width = 9, height = 8, units = "in")


# Lake-level stability vs covariates (regressions) ####

#lat, mean temp, and chla are very correlated
#temp range is not correlated with the other predictors
cor(select(ungroup(globalseassummary),Temp_range,Temp_mean,absLat,Photoperiod_range,lnChla),use="p")

#LE
m1s=lm(JLEmean~Temp_mean+Temp_range, data=filter(globalseassummary,Level2=="Species"))
m1f=lm(JLEmean~Temp_mean+Temp_range, data=filter(globalseassummary,Level2=="Functional Group"))
m1t=lm(JLEmean~Temp_mean+Temp_range, data=filter(globalseassummary,Level2=="Trophic Level"))
anova(m1s);anova(m1f);anova(m1t)
#mean temp is sig for FG and marginal for TL. range marginal for TL

#seas
m1s=lm(llerange~Temp_mean+Temp_range, data=filter(globalseassummary,Level2=="Species"))
m1f=lm(llerange~Temp_mean+Temp_range, data=filter(globalseassummary,Level2=="Functional Group"))
m1t=lm(llerange~Temp_mean+Temp_range, data=filter(globalseassummary,Level2=="Trophic Level"))
anova(m1s);anova(m1f);anova(m1t)
#nothing is sig.

#rel seas
m1s=lm(llerange_rel~Temp_mean+Temp_range, data=filter(globalseassummary,Level2=="Species"))
m1f=lm(llerange_rel~Temp_mean+Temp_range, data=filter(globalseassummary,Level2=="Functional Group"))
m1t=lm(llerange_rel~Temp_mean+Temp_range, data=filter(globalseassummary,Level2=="Trophic Level"))
anova(m1s);anova(m1f);anova(m1t)
#mean temp is sig for TL and marginal for FG, range marginal for TL

#pplle
m1s=lm(pplle~Temp_mean+Temp_range, data=filter(globalseassummary,Level2=="Species"))
m1f=lm(pplle~Temp_mean+Temp_range, data=filter(globalseassummary,Level2=="Functional Group"))
m1t=lm(pplle~Temp_mean+Temp_range, data=filter(globalseassummary,Level2=="Trophic Level"))
anova(m1s);anova(m1f);anova(m1t)
#nothing is sig

#with all data points and random effects
#LE
m1s=lmer(JLEmean~Temp_mean+Temp_range+(1|Site), data=filter(globalseas,Level2=="Species"))
m1f=lmer(JLEmean~Temp_mean+Temp_range+(1|Site), data=filter(globalseas,Level2=="Functional Group"))
m1t=lmer(JLEmean~Temp_mean+Temp_range+(1|Site), data=filter(globalseas,Level2=="Trophic Level"))
#lmerTest::rand(m1s)
car::Anova(m1s);car::Anova(m1f);car::Anova(m1t)
#mean temp is sig for FG and TL

#seas
m1s=lmer(llerange~Temp_mean+Temp_range+(1|Site), data=filter(globalseas,Level2=="Species"))
m1f=lmer(llerange~Temp_mean+Temp_range+(1|Site), data=filter(globalseas,Level2=="Functional Group"))
m1t=lmer(llerange~Temp_mean+Temp_range+(1|Site), data=filter(globalseas,Level2=="Trophic Level"))
#lmerTest::rand(m1s)
car::Anova(m1s);car::Anova(m1f);car::Anova(m1t)
#nothing is sig.

#rel seas
m1s=lmer(llerange_rel~Temp_mean+Temp_range+(1|Site), data=filter(globalseas,Level2=="Species"))
m1f=lmer(llerange_rel~Temp_mean+Temp_range+(1|Site), data=filter(globalseas,Level2=="Functional Group"))
m1t=lmer(llerange_rel~Temp_mean+Temp_range+(1|Site), data=filter(globalseas,Level2=="Trophic Level"))
#lmerTest::rand(m1s)
car::Anova(m1s);car::Anova(m1f);car::Anova(m1t)
#mean temp is sig for FG and TL

# Predictability: Eigenvalue or VER vs. residual ####

local3=full_join(local,select(globalseas, Site, Name, Level, Level2, pplle, llerange)) %>% 
  filter(!is.na(lle),!is.na(logabsres), pplle>0, llerange>0.25, det>-20)

#eigenvalue
ggplot(local3,aes(x=lle, y=logabsres)) +
  facet_grid(.~Level2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept =0) +
  geom_point(aes(color=Site),alpha=0.05, show.legend = F) +
  geom_smooth(show.legend = F, se = F, method = "lm") +
  classic + labs(x="LLE", y="Log abs residual") + removefacetbackground

m1=lme4::lmer(logabsres~lle*Level2+(1|Site),data=local3)
car::Anova(m1)
summary(m1)

#det
ggplot(local3,aes(x=det, y=logabsres)) +
  facet_grid(.~Level2, scales = "free_x") +
  geom_hline(yintercept = 0) + geom_vline(xintercept =0) +
  geom_point(aes(color=Site),alpha=0.05, show.legend = F) +
  geom_smooth(show.legend = F, se = F, method = "lm", lwd=1, color="navy") +
  classic + labs(x="Log abs Det", y="Log abs residual") + removefacetbackground
#ggsave("figures/pred_det_alldata_update.png",width = 6,height = 3)

m1=lme4::lmer(logabsres~det*Level2+(1|Site),data=local3)
car::Anova(m1)
summary(m1)

#VER
ggplot(local3,aes(x=log(ver), y=logabsres)) +
  facet_grid(.~Level2) +
  geom_hline(yintercept = 0) + geom_vline(xintercept =0) +
  geom_point(aes(color=Site),alpha=0.05, show.legend = F) +
  geom_smooth(show.legend = F, se = F, method = "lm", lwd=1, color="navy") +
  classic + labs(x="Log VER", y="Log abs residual") + removefacetbackground
ggsave("figures/pred_coefs_alldata_update.png",width = 6,height = 3)
ggsave("figures/Figure_4.pdf",width = 6,height = 3)

m1=lme4::lmer(logabsres~log(ver)*Level2+(1|Site),data=local3)
car::Anova(m1)
summary(m1)

#distribution of correlations between det and logabsres
#for each series
local4=local3 %>% group_by(Site,Name,Level2) %>% 
  summarise(pcor=cor(det,logabsres,use = "p")) %>% 
  left_join(global)
ggplot(local4, aes(x=pcor,fill=JLEsign)) +
  facet_wrap(.~Level2, scales = "free_y") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  geom_histogram(bins=6,boundary=0) + 
  geom_vline(xintercept = 0) +
  classic + removefacetbackground +
  labs(x="Correlation between log abs det and log abs residual")
#for each site
local5=local3 %>% group_by(Site,Level2) %>% 
  summarise(pcor=cor(det,logabsres,use = "p"))
ggplot(local5, aes(x=pcor)) +
  facet_wrap(.~Level2, scales = "free_y") +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  geom_histogram(bins=6,boundary=0) + 
  geom_vline(xintercept = 0) +
  classic + removefacetbackground +
  labs(x="Correlation between log abs det and log abs residual")

#individual LLE vs temperature ####
# temperature=read.csv("data/temperaturedata.csv")
# local4=left_join(local3,temperature)
# 
# ggplot(local4,aes(y=lle, x=Temp)) +
#   facet_grid(.~Level2, scales = "free_x") +
#   geom_hline(yintercept = 0) + #geom_vline(xintercept =0) +
#   geom_point(aes(color=Site),alpha=0.05, show.legend = F) +
#   geom_smooth(show.legend = F, se = F, method = "lm") +
#   classic + labs(x="Temperature", y="LLE") + removefacetbackground
# 
# ggplot(local4,aes(y=log(ver), x=Temp)) +
#   facet_grid(.~Level2, scales = "free_x") +
#   geom_hline(yintercept = 0) + #geom_vline(xintercept =0) +
#   geom_point(aes(color=Site),alpha=0.05, show.legend = F) +
#   geom_smooth(show.legend = F, se = F, method = "lm") +
#   classic + labs(x="Temperature", y="Det") + removefacetbackground

# plot an indiv time series (for diagnostics) #####
plotts=function(Sitef, Namef, y) {
  ts=filter(local,Site==Sitef, Name==Namef)
  #plot(ts[,"Month"], ts[,y], type="o", ylab=y, main=paste(Sitef, Namef))
  plot(ts[,y], type="o", ylab=y, main=paste(Sitef, Namef))
  print(ts[,y])
}

plotts(Sitef="Lake Zurich", Namef="Daphn.hyali","AbundanceScale")
plotts(Sitef="Port Erin Bay", Namef="Chaetoceras.debile","Abundance")
plotts(Sitef="Lake Geneva", Namef="Bythotrephes longimanus","lle")
plotts(Sitef="Lake Geneva", Namef="Cyclops prealpinus","det")
plotts(Sitef="Lake Mendota (WI)", Namef="ACANTHOCYCLOPS SP","AbundanceLogScale")

# metadata table for paper ####

metat=select(meta, Site, Country, Type, Lat, Lon, Chla_status)
ssyears=local %>% group_by(Site) %>% 
  summarise(minyear=min(Year), maxyear=max(Year))
ndatarange=global %>% group_by(Site) %>% 
  summarise(minn=min(ndatapoints), maxn=max(ndatapoints),
            #mint=min(tslength), maxt=max(tslength),
            species=paste(Name[Level=="SP"], collapse = ", "),
            functional=paste(Name[Level=="FG"], collapse = ", "),
            total=paste(Name[Level=="Total"], collapse = ", "))
metat=left_join(metat, ssyears) %>% left_join(ndatarange) %>% arrange(Type, desc(abs(Lat)))
write.csv(metat, "archive/metatable.csv", row.names = F)
