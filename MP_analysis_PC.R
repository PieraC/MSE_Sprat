rm(list=ls())

library(FLCore)
library(ggplotFL)
library(Cairo)
library(reshape2)
library(data.table)
library(dplyr)
### set up parallel computing environment
library(doParallel)
cl <- makeCluster(28)
registerDoParallel(cl)
getDoParWorkers()
getDoParName()


### ------------------------------------------------------------------------ ###
### list files
### ------------------------------------------------------------------------ ###
### plot all scenarios from external HDD
files <- list.files("C:/Users/PC09/Documents/MSE_Sprat/output/Piera", pattern=".RData")
#files=files[2]
stkTmp <- list()
for(i in files){
  stkTmp[[i]] <- mget(load(paste0("./output/Piera/", i))) 
}

# create list
OMs <- list()
# add 1 over 2
OMs[["OM1.1o2.noUC"]] <- stkTmp[[1]][[1]][[1]]$pstk
OMs[["OM1.1o2.UC"]] <- stkTmp[[1]][[2]][[1]]$pstk
OMs[["OM2.1o2.noUC"]] <- stkTmp[[1]][[1]][[2]]$pstk
OMs[["OM2.1o2.UC"]] <- stkTmp[[1]][[2]][[2]]$pstk
# add harvest strategies with different %
OMs[["OM1.HR.01.noUC"]] <- stkTmp[[2]][[1]][[1]]$pstk
OMs[["OM1.HR.02.noUC"]] <- stkTmp[[2]][[2]][[1]]$pstk
OMs[["OM1.HR.03.noUC"]] <- stkTmp[[2]][[3]][[1]]$pstk
OMs[["OM1.HR.04.noUC"]] <- stkTmp[[2]][[4]][[1]]$pstk
OMs[["OM2.HR.01.noUC"]] <- stkTmp[[2]][[1]][[2]]$pstk
OMs[["OM2.HR.02.noUC"]] <- stkTmp[[2]][[2]][[2]]$pstk
OMs[["OM2.HR.03.noUC"]] <- stkTmp[[2]][[3]][[2]]$pstk
OMs[["OM2.HR.04.noUC"]] <- stkTmp[[2]][[4]][[2]]$pstk
# add harvest strategies with different %
OMs[["OM1.HR.015.noUC"]] <- stkTmp[[3]][[1]][[1]]$pstk
OMs[["OM1.HR.016.noUC"]] <- stkTmp[[3]][[2]][[1]]$pstk
OMs[["OM1.HR.017.noUC"]] <- stkTmp[[3]][[3]][[1]]$pstk
OMs[["OM1.HR.018.noUC"]] <- stkTmp[[3]][[4]][[1]]$pstk
OMs[["OM1.HR.019.noUC"]] <- stkTmp[[3]][[5]][[1]]$pstk
OMs[["OM2.HR.015.noUC"]] <- stkTmp[[3]][[1]][[2]]$pstk
OMs[["OM2.HR.016.noUC"]] <- stkTmp[[3]][[2]][[2]]$pstk
OMs[["OM2.HR.017.noUC"]] <- stkTmp[[3]][[3]][[2]]$pstk
OMs[["OM2.HR.018.noUC"]] <- stkTmp[[3]][[4]][[2]]$pstk
OMs[["OM2.HR.019.noUC"]] <- stkTmp[[3]][[5]][[2]]$pstk



### ------------------------------------------------------------------------ ###
### statistics ####
### ------------------------------------------------------------------------ ###
B_lim <- 400

# #stats <- foreach(scenario = STK_HR_1o2_NoUC, .packages = "FLCore", .export = "B_lim") %dopar% {
stats <- list()

# Statistics from 33 to 48 when more stable
for(i in 1:length(OMs)){
    #stats <- foreach(stk_tmp = res, .packages = "FLCore", .export = "B_lim") %dopar% {
  #lapply(seq_along(res), function(i){
  ### load stk
  stk_tmp <- OMs[[i]]#get(load(paste0("output/Piera/", scenario, ".RData")))
  res_temp <- list()
  n_iter <- dim(stk_tmp)[6]
  
  ### proportion where SSB was below B_lim
  res_temp$ssb_below_blim_total <-
    sum(ssb(stk_tmp)[, ac(33:48)] < B_lim) / (n_iter*16)
  ### proportion of iterations where SSB dropped below B_lim
  res_temp$ssb_below_blim_iter <-
    sum(yearSums(ssb(stk_tmp)[, ac(33:48)] < B_lim) > 0) / n_iter
  
  ### stock collapse = ssb < 1
  res_temp$collapse_total <-
    sum(ssb(stk_tmp)[, ac(33:48)] < 1) / (n_iter*24)
  ### proportion of iterations with collapse
  res_temp$collapse_iter <-
    sum(yearSums(ssb(stk_tmp)[, ac(33:48)] < 1) > 0) / n_iter
  
  ### how frequently is max F reached?
  res_temp$fmaxed_total <- sum(fbar(stk_tmp)[, ac(33:48)] == 2) / (n_iter*16)
  ### in how many iterations did this happen?
  res_temp$fmaxed_iter <- sum(yearSums(fbar(stk_tmp)[, ac(33:48)] == 2) > 0)/
    n_iter
  
  ### yield
  res_temp$yield <- mean(catch(stk_tmp[,ac(33:48)]))
  res_temp$rel_yield <- mean(catch(stk_tmp[,ac(33:48)])) / mean(catch(stk_tmp[,ac(1:24)]))
  ### scenario definition
  #res_temp$scenario <- allscenarios[19:27]
  res_temp <- as.data.frame(res_temp)
  stats[[i]] <- res_temp
  #return()
  
}#)

STATS <- do.call(rbind, stats)
STATS$scenarios <- names(OMs)


### save
saveRDS(STATS, file = "./output/Piera/AllOMs_STATS_w15-19.RDS")
write.csv(STATS, file = "./output/Piera/AllOMs_STATS_w15-19.csv")


### ------------------------------------------------------------------------ ###
### PLOTS
### ------------------------------------------------------------------------ ###
# ********************************************************
# plot of SSB relative to Bmsy per each operating model
OMs.SSB <- lapply(OMs, function(x){
  as.data.frame(ssb(x))
})

OMs.SSB <- data.table(do.call(rbind, OMs.SSB))
OMs.SSB$scenario <- rep(names(OMs), each=50*500)
OMs.SSB$OM <- substr(OMs.SSB$scenario,1,3)

myOMs <- c("OM1.1o2.noUC", "OM1.1o2.UC", "OM2.1o2.noUC", "OM2.1o2.UC",
           "OM1.HR.01.noUC", "OM1.HR.02.noUC", "OM1.HR.03.noUC", "OM1.HR.04.noUC",
           "OM2.HR.01.noUC", "OM2.HR.02.noUC", "OM2.HR.03.noUC", "OM2.HR.04.noUC")
           
P.SIM.SSB <- ggplot(OMs.SSB[year>32 & year<50 & scenario %in% myOMs,], aes(x=factor(year), y=data, fill=scenario)) + geom_boxplot() +
  theme_light(35) + ylab("SSB") + xlab("Year") +
  theme(legend.position = "bottom", legend.title=element_blank()) +
  geom_hline(yintercept=B_lim, col="red") +
  theme(legend.text=element_text(size=20), legend.spacing.x = unit(1.0, 'cm'),legend.spacing.y = unit(2, 'cm')) +
  #geom_hline(yintercept=OMs.catch[year<25 & substr(scenario,1,3)=="OM2",mean(data)], col="blue") +
  facet_grid(OM~.)

# Plot
png("./output/Piera/plots/SIM_SSB.png", width=6000, height=4000, res=250)
P.SIM.SSB
dev.off()

# ********************************************************
# Catch
# plot of Catch relative to Bmsy per each operating model
OMs.catch <- lapply(OMs, function(x){
  as.data.frame(catch(x))
})

OMs.catch <- data.table(do.call(rbind, OMs.catch))
OMs.catch$scenario <- rep(names(OMs), each=50*500)
OMs.catch$OM <- substr(OMs.catch$scenario,1,3)
#OMs.SSB <- rbind(as.data.frame(ssb(stk.tmp.om1)), as.data.frame(ssb(stk.tmp.om2)))
#OMs.catch$OMs <- c(rep("OM1", 22*250),rep("OM2", 22*250))

P.SIM.CATCH <- ggplot(OMs.catch[year>32 & year<50 & scenario %in% myOMs,], aes(x=factor(year), y=data, fill=scenario)) + geom_boxplot() +
  theme_light(38) + ylab("Catch") + xlab("Year") +
  theme(legend.position = "bottom", legend.title=element_blank()) +
  geom_hline(yintercept=OMs.catch[year<25,mean(data)], col="red") +
  theme(legend.text=element_text(size=20), legend.spacing.x = unit(1.0, 'cm'),legend.spacing.y = unit(2, 'cm')) +
  #geom_hline(yintercept=OMs.catch[year<25 & substr(scenario,1,3)=="OM2",mean(data)], col="blue") +
  facet_grid(OM~.)


png("./output/Piera/plots/SIM_CATCH.png", width=6000, height=4000, res=250)
P.SIM.CATCH
dev.off()

###
# Calculate SD of Yield to check how much it varies
YieldSD <- OMs.catch[,sd(data), by=.(year,scenario)]
P.YieldSD <-ggplot(YieldSD[year>32 & year < 49 & scenario %in% myOMs,], aes(year, V1, col=scenario, shape=scenario)) + 
  geom_line(size=1.2) + geom_point() +
  theme_light(25) + theme(legend.position = "bottom") + 
  scale_y_continuous(breaks=seq(0,300,50)) +
  ylab("Yield SD")

png("./output/Piera/plots/SIM_YieldSD.png", width=3000, height=2000, res=250)
P.YieldSD
dev.off()


# ********************************************************
# Yield vs Risk
STATS$OMs <- substr(STATS$scenarios, 1,3)
STATS <- data.table(STATS)
STATS$scenarios_f <- factor(STATS$scenarios)
P.YvsR <- ggplot(STATS[scenarios %in% myOMs,], aes(ssb_below_blim_total, yield, col=scenarios)) + 
  geom_point(size=2, stroke=2) +
  #scale_shape_manual(values=1:nlevels(STATS$scenarios_f)) +
  facet_grid(.~OMs) + theme_light(25) + xlab("Risk") + ylab("Yield") +
  theme(legend.text=element_text(size=20), legend.title=element_blank(),legend.spacing.x = unit(1.0, 'cm'),legend.spacing.y = unit(2, 'cm')) +
  theme(legend.position="bottom")

png("./output/Piera/plots/SIM_YvsR.png", width=4000, height=2000, res=250)
P.YvsR
dev.off()

### ------------------------------------------------------------------------ ###
### Plot true VS perceived stock
survey_time=5/12
# For 1o2 and 20% HR Management strategy
STK_MS <- FLStocks(OM1_1o2 = stkTmp[[1]][[1]][[1]]$pstk, OM2_1o2 = stkTmp[[1]][[1]][[2]]$pstk, OM1_HR20=stkTmp[[2]][[2]][[1]]$pstk, OM2_HR20=stkTmp[[2]][[2]][[2]]$pstk,
                   OM1_HR40=stkTmp[[2]][[4]][[1]]$pstk, OM2_HR40=stkTmp[[2]][[4]][[2]]$pstk)
#IDX_MS <- FLStocks(OM_1o2 = stkTmp[[1]][[1]][[1]]$pstk, OM_1o2 = stkTmp[[1]][[1]][[2]]$pstk, OM1_HR20=stkTmp[[1]][[1]][[2]]$pstk, OM2_HR20=stkTmp[[2]][[2]][[2]]$pstk)

TrueBio <- list()
Perceived_1o2 <- list()
Perceived_HR20 <- list()
Perceived_HR40 <- list()
for(i in 1:length(STK_MS)) TrueBio[[i]] <- data.table(as.data.frame(quantSums(stock.n(STK_MS[[i]])*stock.wt(STK_MS[[i]])*exp(-(harvest(STK_MS[[i]])*survey_time + m(STK_MS[[i]])*survey_time)))))
for(i in 1:2) Perceived_1o2[[i]] <- data.table(as.data.frame(stkTmp[[1]][[1]][[i]]$idx.bio))
for(i in 1:2) Perceived_HR20[[i]] <- data.table(as.data.frame(quantSums(stkTmp[[2]][[2]][[i]]$idx@index*stock.wt(stkTmp[[2]][[2]][[i]]$pstk))))
for(i in 1:2) Perceived_HR40[[i]] <- data.table(as.data.frame(quantSums(stkTmp[[2]][[4]][[i]]$idx@index*stock.wt(stkTmp[[2]][[4]][[i]]$pstk))))

OM <- c("OM1", "OM2")
TrueBio <- data.table(do.call(rbind.data.frame, mapply(cbind, TrueBio, "OM"=OM, SIMPLIFY=F)))
Perceived_1o2 <- data.table(do.call(rbind.data.frame, mapply(cbind, Perceived_1o2, "OM"=OM, SIMPLIFY=F)))
Perceived_HR20 <- data.table(do.call(rbind.data.frame, mapply(cbind, Perceived_HR20, "OM"=OM, SIMPLIFY=F)))
Perceived_HR40 <- data.table(do.call(rbind.data.frame, mapply(cbind, Perceived_HR40, "OM"=OM, SIMPLIFY=F)))

TrueBio[,c("type","HR"):=list("True",c(rep("HR1o2",50000),rep("HR20",50000),rep("HR40",50000)))]
Perceived_1o2[,c("type","HR"):=list("Perceived", rep("HR1o2",50000))]
Perceived_HR20[,c("type","HR"):=list("Perceived", rep("HR20",50000))]
Perceived_HR40[,c("type","HR"):=list("Perceived", rep("HR40",50000))]

names(Perceived_1o2)[1] <- "age"
STK_PERC <- rbind(TrueBio, Perceived_1o2, Perceived_HR20, Perceived_HR40)


P_STK_PERC <- ggplot(STK_PERC[year<50 & HR!="HR40",], aes(year, data, color=type, fill=type)) + 
  stat_summary(fun.y="mean", geom="line", size=1.2) + ylab("Total Biomass at survey time") +
  stat_summary(fun.data="mean_se", geom="ribbon", alpha=.3, colour=NA) +
  theme_bw(25) + theme(legend.title = element_blank(), legend.position="top") +
  facet_grid(HR~OM)

png("./output/Piera/plots/TruePerc_Biomass.png", width=4000, height=3500, res=300)
P_STK_PERC 
dev.off()



### ------------------------------------------------------------------------ ###
### PLOT STKS
### ------------------------------------------------------------------------ ###
my_OMs <- OMs[myOMs]
OM1_1o2_STKs <- FLStocks(OM1_HCR1o2.UC=my_OMs[["OM1.1o2.UC"]], OM1_HCR1o2.noUC=my_OMs[["OM1.1o2.noUC"]])

OM2_1o2_STKs <- FLStocks(OM2_HCR1o2.UC=my_OMs[["OM2.1o2.UC"]], OM2_HCR1o2.noUC=my_OMs[["OM2.1o2.noUC"]])

OM1_HR_STKs <- FLStocks(OM1.HR.10=my_OMs[["OM1.HR.01.noUC"]], OM1.HR.20=my_OMs[["OM1.HR.02.noUC"]],
                   OM1.HR.30=my_OMs[["OM1.HR.03.noUC"]], OM1.HR.40=my_OMs[["OM1.HR.04.noUC"]])

OM2_HR_STKs <- FLStocks(OM2.HR.10=my_OMs[["OM2.HR.01.noUC"]], OM2.HR.20=my_OMs[["OM2.HR.02.noUC"]],
                   OM2.HR.30=my_OMs[["OM2.HR.03.noUC"]], OM2.HR.40=my_OMs[["OM2.HR.04.noUC"]])


OM1_1o2 <- plot(OM1_1o2_STKs) + ggtitle("OM1 - 1 over 2") + theme_light(25) + xlab("Year") +
  theme(legend.position="bottom", legend.title=element_blank(),plot.title = element_text(hjust = 0.5))
OM2_1o2 <- plot(OM2_1o2_STKs) + ggtitle("OM2 - 1 over 2") + theme_light(25) + 
  theme(legend.position="bottom", legend.title=element_blank(),plot.title = element_text(hjust = 0.5))
OM1_HR <- plot(OM1_HR_STKs) + ggtitle("OM1 - Fixed Harvest Rate") + theme_light(25) + 
  theme(legend.position="bottom", legend.title=element_blank(),plot.title = element_text(hjust = 0.5))
OM2_HR <- plot(OM2_HR_STKs) + ggtitle("OM2 - Fixed Harvest Rate") + theme_light(25) + 
  theme(legend.position="bottom", legend.title=element_blank(),plot.title = element_text(hjust = 0.5))

png("./output/Piera/plots/STK_OM1_1o2.png", width=4000, height=3500, res=300)
OM1_1o2 
dev.off()

png("./output/Piera/plots/STK_OM1_HR.png", width=4000, height=3500, res=300)
OM1_HR
dev.off()

png("./output/Piera/plots/STK_OM2_1o2.png", width=4000, height=3500, res=300)
OM2_1o2 
dev.off()

png("./output/Piera/plots/STK_OM2_HR.png", width=4000, height=3500, res=300)
OM2_HR
dev.off()



### ------------------------------------------------------------------------ ###
### PLOT CATCHABILITIES
### ------------------------------------------------------------------------ ###
# Plot catchability
idxQ <- data.table(rbind(as.data.frame(stkTmp[[1]][[1]][[1]]$observations[[1]]$idx@index.q),
                         as.data.frame(stkTmp[[1]][[1]][[2]]$observations[[1]]$idx@index.q),
                         as.data.frame(stkTmp[[1]][[1]][[3]]$observations[[1]]$idx@index.q),
                         as.data.frame(stkTmp[[1]][[1]][[4]]$observations[[1]]$idx@index.q),
                         as.data.frame(stkTmp[[1]][[1]][[5]]$observations[[1]]$idx@index.q),
                         as.data.frame(stkTmp[[1]][[1]][[6]]$observations[[1]]$idx@index.q)))
# OM names              
idxQ$OM <- c(rep("OM1",500*6*50),rep("OM2",500*6*50),rep("OM3",500*8*50),
             rep("OM4",500*8*50),rep("OM5",500*6*50),rep("OM6",500*6*50))

# summarise
myQ_sum <- idxQ %>%
  group_by(age,year,OM) %>%
  summarize(mean_size = mean(data), stDev=sd(data))

# Calculate confidence intervals
myQ_sum <- data.table(myQ_sum)
myQ_sum[, se:=1.96*stDev/sqrt(500)]
myQ_sum[,c("lCI","uCI","OMs"):=list(mean_size-se,mean_size+se,substr(OM,1,3))]
# Plot
P_Q <- ggplot(myQ_sum[year %in% seq(1,50,10) & OM %in% c("OM1","OM3","OM5"),], 
              aes(age, mean_size, col=factor(year), group=year)) + 
  geom_line() + 
  scale_x_continuous(breaks=seq(0,7,1)) + ylab("q") +
  geom_ribbon(aes(ymax=uCI,ymin=lCI,fill=factor(year)),alpha=.25, colour=NA) +
  facet_wrap(~OMs) + 
  guides(colour=guide_legend(title="Year"),fill=guide_legend(title="Year")) +
  theme_bw(22) +
  theme(legend.position = "top", legend.direction = "horizontal")
# save plot
png("./Input/Piera/OMs_Q.png", width=3000, height=1500, res=300)
P_Q
dev.off()

#



# Plot catchability
idxQ <- as.data.frame(OM_list[[1]]$observations$idx@.Data[[1]]@index.q)
q_modeled <- as.data.frame(q_modeled)
q_modeled$ages <- 0:5
names(q_modeled)[1] <- "q"

Qplot <- ggplot(q_modeled, aes(ages, q)) + geom_line(col="red", size=1.2) +
  theme_classic(20)
 
png("./output/Piera/plots/Qplot.png", width=1500, height=1500, res=200)
Qplot
dev.off()

plot_lst <- function(scenario, name = ""){
  stk_lst <- lapply(scenario, function(x){
    get(load(paste0("output/Piera/", x, ".RData")))
  })
  # stk_lst <- res[as.character(select$scenario)]
  names(stk_lst) <- scenario
  ### median over all iterations
  stk_lst <- lapply(stk_lst, function(x){
    qapply(x, iterMedians)
  })
  stk_lst <- FLStocks(stk_lst)
  p <- plot(stk_lst, col=names) + theme_bw() +
    xlab("year")
  #return(p)
  ggsave(filename = paste0("output/Piera/plots/", name, ".png"),
         width = 18, height = 15, units = "cm", dpi = 300, type = "cairo-png",
         plot = p)
}

plot_lst(allscenarios, name="median")

