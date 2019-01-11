# Use code from Simon Fisher ()

### ------------------------------------------------------------------------ ###
### create biological stocks from life-history parameters ####
### as used for WKLIFE VII 2017
### ------------------------------------------------------------------------ ###
### author: Simon Fischer (Cefas), simon.fischer@cefas.co.uk
### FRAMEWORK based on the a4a standard MSE developed at JRC
### by Ernesto Jardim, Iago Mosqueira, Finlay Scott, et al.
### additional contributors:
### Karin Olsson
### ------------------------------------------------------------------------ ###
### created 08/2017
### last modifications:
### 2017-12 Simon Fischer
### 2018-12 Piera Carpi
### ------------------------------------------------------------------------ ###
# source("http://flr-project.org/R/instFLR.R")

rm(list=ls())
### load packages
library(FLife)
library(FLBRP)
library(FLasher)
library(foreach)
library(doParallel)
library(data.table)
### set up cluster for parallel computing
#cl <- makeCluster(28)
#registerDoParallel(cl)

### load additional functions
### select R scripts from functions folder
load_files <- list.files("functions/")
load_files <- load_files[grepl(pattern = "*.R$", x = load_files)]
### source the scripts
invisible(lapply(paste0("functions/", load_files), source))


### ------------------------------------------------------------------------ ###
### create STK from life history charactheristics ####
### ------------------------------------------------------------------------ ###
# 1 OPERATING MODELS
biol.params.OM1 <- data.frame(linf=16, k=0.6, t0=-0.8, a = 0.0000048, b=3.19, s = 0.5, v=1000) #, sr=1, s=0.5)
biol.params.OM2 <- data.frame(linf=16, k=0.4, t0=-0.8, a = 0.0000048, b=3.19, s = 0.5, v=1000) #, sr=1, s=0.5)
biol.params.OM3 <- data.frame(linf=13, k=0.6, t0=-0.8, a = 0.0000048, b=3.19, s = 0.5, v=1000) #, sr=1, s=0.5)

#i <- biol.params.OM1
# brps <- foreach(i = split(stocks_lh, 1:nrow(stocks_lh)),
#                 .errorhandling = "pass", 
#                 .packages = c("FLife", "FLasher", "FLBRP")) %dopar% {
### create brp
brps.OM1 <- getBrp(i=biol.params.OM1)
brps.OM2 <- getBrp(i=biol.params.OM2)
brps.OM3 <- getBrp(i=biol.params.OM3)

### get reference points
ref_pts <- refpts(brps.OM1)
#names(ref_pts) <- stocks_lh[2,"stock"]

brps <- list(brps.OM1, brps.OM2, brps.OM3)
names(brps) <- seq_along(brps)
### save stock list
#write.csv(file = "input/stock_list_full2.csv", x = stocks_lh)


### ------------------------------------------------------------------------ ###
### create FLStock & project forward  
### ------------------------------------------------------------------------ ###
### number of iterations and projection years
its <- 500
ny_sim <-25

# empty list for OMs
myOM <- list()
myOMs <- list()
# options for sd in rectruitment
recVar <- c(0.3, 0.5)

#OMs <- foreach(i = brps, .errorhandling = "pass", .packages = c("FLife", "FLasher", "FLBRP")) %dopar% {
for(y in 1:length(brps)){
  set.seed(12*y)
  for(rr in recVar){
    #i <- brps.OM1
    # convert FLBRP into FLStock
    i = brps[[y]]
    # create FLStock
    stk <- as(i, "FLStock")[, 2]
    # f and m before spawning = 75%
    harvest.spwn(stk) <- 0.75 
    m.spwn(stk) <- 0.75
    ### name first year "1"
    stk <- qapply(stk, function(x) {dimnames(x)$year <- "1"; return(x)})
    ### extend object to year ny_sim
    stk <- fwdWindow(stk, i, end = ny_sim)
    ### propagate with requested number of iterations
    stk <- propagate(stk, its)
    ### create FLSR object
    #stk_sr <- FLSR(params=params(i), model = model(i))
    stk_sr <- FLSR(params=params(i), model = "segreg")
    ### create residuals for (historical) projection
    #set.seed(123)
    residuals <- rlnoise(its, rec(stk) %=% 0, sd = rr, b = 0)
    ### project forward
    years_target <- ((dims(stk)$minyear+1):ny_sim)
    ## Calculate f corresponding to Patterson E
    Ftarget <- 0.4*quantMeans(m(stk[ac(1:3),1,,,,1]))/0.6
    # set up control
    my_control <- fwdControl(year = years_target, value = an(Ftarget), quant = "f")   
    # Apply control to simulated stock
    STK <- fwd(stk, sr = stk_sr, control=my_control, residuals = residuals[, ac(years_target)])
    ### project last 25 years with 2 scenarios: one-way & roller-coaster
    years_target <- 2:ny_sim
    ### one-way: linear increase to Ftarget, which corresponds to E=0.4
    stk_proj <- oneWayTrip(STK, sr = stk_sr, brp = brps, fmax = Ftarget, 
                           years = years_target,
                           residuals = residuals[, ac(years_target)],
                           f0 = 0.1) #refpts(brps)["msy", "harvest"]*0.5
    myOM[[ac(rr)]] <- list(sr=stk_sr, sr.res=residuals, stk=stk_proj, brp=i)

  }
  myOMs[[y]] <-  myOM  
}

# Plot my OMs
STK_OMs <- FLStocks(OM1_03 = myOMs[[1]][["0.3"]]$stk, OM1_05 = myOMs[[1]][["0.5"]]$stk, 
                    OM2_03 = myOMs[[2]][["0.3"]]$stk, OM2_05 = myOMs[[2]][["0.5"]]$stk,
                    OM3_03 = myOMs[[3]][["0.3"]]$stk, OM3_05 = myOMs[[3]][["0.5"]]$stk)
P.OMs <- plot(STK_OMs) + theme_light(30) 
# Print plot of operating models
png("./Input/Piera/OMs_STK.png", width=4000, height=3300, res=250)
P.OMs
dev.off()

# *************************************
# Plot other biological characteristics
# MATURITY
STK_mat <- data.table(getQuant(STK_OMs[c(1,3)], mat))
P.MAT <- ggplot(STK_mat[iter==1 & year==1,], aes(age, data, col=OM)) + geom_line(size=1.2) +
  theme_classic(25) + ylab("% mat") +
  scale_x_continuous(breaks=seq(0,7,1)) +
  scale_color_manual(labels = c("OM1 & OM3", "OM2"), values = c("salmon", "deepskyblue4")) +
  theme(legend.position = c(0.7,0.3), legend.title = element_blank())

STK_m <- data.table(getQuant(STK_OMs[c(1,3,5)], m))
P.M <- ggplot(STK_m[iter==1 & year==1,], aes(age, data, col=OM)) + geom_line(size=1.2) +
  theme_classic(25) + ylab("M") + scale_y_continuous(breaks=seq(0.5,1.8,0.2)) +
  scale_x_continuous(breaks=seq(0,7,1)) +
  scale_color_manual(labels = c("OM1", "OM2", "OM3"), values=c("salmon", "deepskyblue4", "seagreen3")) +
  theme(legend.position = c(0.7,0.7), legend.title = element_blank())

STK_wt <- data.table(getQuant(STK_OMs[c(1,3,5)], catch.wt))
P.WT <- ggplot(STK_wt[iter==1 & year==1,], aes(age, data, col=OM)) + geom_line(size=1.2) +
  theme_classic(25) + ylab("kg") + scale_y_continuous(breaks=seq(0.004,0.032,0.003)) +
  scale_x_continuous(breaks=seq(0,7,1)) +
  scale_color_manual(labels = c("OM1", "OM2", "OM3"), values=c("salmon", "deepskyblue4", "seagreen3")) +
  theme(legend.position = c(0.7,0.2), legend.title = element_blank())

STK_f <- data.table(getQuant(STK_OMs[c(1,3,5)], harvest))
STK_sel <- STK_f[iter==1 & year==1,]
refAge <- STK_f[iter==1 & year==1 & age==2,data]
refAge <- c(rep(refAge[1], 6), rep(refAge[2],8),rep(refAge[3],6))
STK_sel[,sel:=data/refAge]#,rep(0.04036131,6)]
P.SEL <- ggplot(STK_sel, aes(age, sel, col=OM)) + geom_line(size=1.2) +
  theme_classic(25) + ylab("sel") + scale_y_continuous(breaks=seq(0.2,1,0.1)) +
  scale_x_continuous(breaks=seq(0,7,1)) +
  scale_color_manual(labels = c("OM1", "OM2", "OM3"), values=c("salmon", "deepskyblue4", "seagreen3")) +
  theme(legend.position = c(0.7,0.2), legend.title = element_blank())

# Print all biological assumptions
png("./Input/Piera/OMs_BIOL_%3d.png", width=2000, height=1500, res=300)
P.MAT
P.M
P.WT
P.SEL
dev.off()

### ------------------------------------------------------------------------ ###
### Save OMs and life history parameters used
### ------------------------------------------------------------------------ ###

### extract life-history parameters
lhpar.OM1 <- attr(brps.OM1, "lhpar")
lhpar.OM2 <- attr(brps.OM2, "lhpar")
lhpar.OM3 <- attr(brps.OM3, "lhpar")

attr(brps, "lhpar") <- NULL

OMs <- list()                 
### return list
OMs[["OM1_0.3"]] <- list(lhpar = lhpar.OM1, sr=myOMs[[1]][["0.3"]]$sr, sr.res=myOMs[[1]][["0.3"]]$sr.res, brp=myOMs[[1]][["0.3"]]$brp, stk=myOMs[[1]][["0.3"]]$stk)
OMs[["OM1_0.5"]] <- list(lhpar = lhpar.OM1, sr=myOMs[[1]][["0.5"]]$sr, sr.res=myOMs[[1]][["0.5"]]$sr.res, brp=myOMs[[1]][["0.5"]]$brp, stk=myOMs[[1]][["0.5"]]$stk)
OMs[["OM2_0.3"]] <- list(lhpar = lhpar.OM2, sr=myOMs[[2]][["0.3"]]$sr, sr.res=myOMs[[2]][["0.3"]]$sr.res, brp=myOMs[[2]][["0.3"]]$brp, stk=myOMs[[2]][["0.3"]]$stk)
OMs[["OM2_0.5"]] <- list(lhpar = lhpar.OM2, sr=myOMs[[2]][["0.5"]]$sr, sr.res=myOMs[[2]][["0.5"]]$sr.res, brp=myOMs[[2]][["0.5"]]$brp, stk=myOMs[[2]][["0.5"]]$stk)
OMs[["OM3_0.3"]] <- list(lhpar = lhpar.OM3, sr=myOMs[[3]][["0.3"]]$sr, sr.res=myOMs[[3]][["0.3"]]$sr.res, brp=myOMs[[3]][["0.3"]]$brp, stk=myOMs[[3]][["0.3"]]$stk)
OMs[["OM3_0.5"]] <- list(lhpar = lhpar.OM3, sr=myOMs[[3]][["0.5"]]$sr, sr.res=myOMs[[3]][["0.5"]]$sr.res, brp=myOMs[[3]][["0.5"]]$brp, stk=myOMs[[3]][["0.5"]]$stk)

### save list
saveRDS(OMs, "Input/Piera/OMs.rds")
#               

# END
#stopCluster(cl)
myrec <- data.table(rbind(as.data.frame(rec(myOMs[[1]][["0.3"]]$stk[,1:25])), 
                          as.data.frame(rec(myOMs[[1]][["0.5"]]$stk[,1:25])),
                          as.data.frame(rec(myOMs[[2]][["0.3"]]$stk[,1:25])),
                          as.data.frame(rec(myOMs[[2]][["0.5"]]$stk[,1:25])),
                          as.data.frame(rec(myOMs[[3]][["0.3"]]$stk[,1:25])),
                          as.data.frame(rec(myOMs[[3]][["0.5"]]$stk[,1:25]))))
myrec$OM <- c(rep("OM1_03",500*25),rep("OM1_05",500*25),rep("OM2_03",500*25),
              rep("OM2_05",500*25),rep("OM3_03",500*25),rep("OM3_05",500*25)) 
myrec$OMid <- substr(myrec$OM, 1,3)
myrec$recvarId <- substr(myrec$OM, 5,6)
P_RECit1 <- ggplot(myrec[iter==15,], aes(year, data, col=recvarId)) + geom_line(size=1.2) +
  facet_grid(.~OMid) + 
  theme_bw(22) + 
  ylab("Rec") + xlab("Year") + 
  theme(legend.position = c(0.1, 0.8), legend.title = element_blank())

png("./Input/Piera/OMs_REC_1it.png", width=3000, height=1200, res=300)
P_RECit1
dev.off()
