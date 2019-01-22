
rm(list=ls())

### ------------------------------------------------------------------------ ###
### libraries ####
### ------------------------------------------------------------------------ ###
required_pckgs <- c("FLash", "FLAssess", "ggplotFL", "FLBRP", "data.table")
### save as object in order to avoid output to screen
. <- lapply(required_pckgs, function(x){
  suppressMessages(library(x, character.only = TRUE))
})
library(doParallel)

cl <- makeCluster(28)
registerDoParallel(cl)

### ------------------------------------------------------------------------ ###
### functions ####
### ------------------------------------------------------------------------ ###
### select R scripts from functions folder
load_files <- list.files("functions/")
load_files <- load_files[grepl(pattern = "*.R$", x = load_files)]
### source the scripts
invisible(lapply(paste0("functions/", load_files), source))


# list all OMs
myOMs <- list.files("Input/Piera/stocks/", pattern=".RData")
# read in
OM_list <- list()
for(i in 1:length(myOMs)){
  OM_list[[i]] <- mget(load(paste0("Input/Piera/stocks/", myOMs[i])))
}


### START WITH SIMULATIONS 

########################################################
# ======================================================
# FIXED PERCENTAGE OF INDEX BIOMASS
# ======================================================
# percentages to be tested for HS
perc <- 0.2
survey_time <- 5/12
# Uncertainty cap yes or no?
uncertainty_Cap = FALSE

HR.02.noUC.20idx <- foreach(x = OM_list, .packages = required_pckgs) %dopar% {
  
  pstk <- x$stk
  idx <- x$observations$idx@.Data[[1]]
  stk <- x$stk
  sr.om <- x$sr.om
  sr.om.res <- x$sr.om.res
  
  TAC <- FLQuant(NA, dimnames=list(TAC="all", year=c(x$dy,x$vy), iter=1:x$it))
  TAC[,ac(x$dy)] <- catch(stk)[,ac(x$dy)]
  TAC[,ac(x$iy)] <- TAC[,ac(x$dy)] #assume same TAC in the first intermediate year
  ctrl <- getCtrl(c(TAC[,ac(x$iy)]), "catch", x$iy, x$it)
  
  for(i in x$vy[-length(x$vy)]){
    ## i <- vy[-length(vy)][1]
    print(i)
    gc()
    ay <- an(i)
    cat(i, "> ")
    vy0 <- 1:(ay-x$y0) # data years (positions vector) - one less than current year
    sqy <- ac((ay-1):(ay-x$nsqy)) # years for status quo computations
    
    stk0 <- pstk[,vy0]
    #catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
    ## note that vy0 is changing below so index is being updated
    # for (index_counter in 1:length(idx)){
    idx0 <- idx[,vy0]
    #index(idx[[index_counter]])[,i] <- stock.n(stk)[,i]*index.q(idx[[index_counter]])[,i] #+ 1
    # calculation at time of the survey
    index(idx)[,i] <- index.q(idx)[,i]*stock.n(pstk)[,i]*exp(-(harvest(pstk)[,i]*survey_time+m(pstk)[,1]*survey_time)) #+ 1
    # Calculate index tot bio
    idx.bio <- quantSums(index(idx)[,i]*stock.wt(pstk)[,i])*1.2
    #}
    ##
    
    # Calculate TAC from index
    advice <- idx.bio*perc
    #   rmodel <- sr
    # fmod <- ~ s(age, k=4) + s(year, k=16)
    # fit <- sca(stk0, FLIndices(idx0), fmodel=fmod, srmodel=rmodel)
    # stk0 <- stk0 + fit
    # fwd control
    fsq0 <- yearMeans(fbar(pstk)[,sqy])
    csq0 <- yearMeans(landings(pstk)[,sqy])
    dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max")) #, ay + 1
    arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
    ## ftrg.vec <- rep(ftrg, it) ## original
    refpt <- data.frame(ssb = 1, harvest = 1, catch=1)
    # Estimate catch vector
    catch.vec <- an(advice)
    #Bescape <- blim
    # if(uncertainty_Cap == TRUE){
    #   new_catchVec <- hcr.check.catch(catch(stk0)[, ac(an(i) - 1)], catch.vec=catch.vec)}else {
        new_catchVec=catch.vec
     # }
    # OM proj
    TAC[,ac(ay+1)] <- catch.vec
    
    # apply the TAC to the operating model stock
    ctrl <- getCtrl(c(TAC[,ac(ay+1)]), "catch", ay+1, it)
    # ##ctrl@trgtArray[,"val",] <- c(TAC[,ac(ay+1)]) #+ BB[,ac(ay)])
    pstk <- fwd(pstk, ctrl=ctrl, sr=sr.om, sr.residuals = sr.om.res, sr.residuals.mult = TRUE) #, exp(sr.res[,ac(ay:(ay+1))]) sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
    #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
  }
  
  x$idx <- idx
  x$pstk <- pstk
  return(x)

}



trueBio <- quantSums(stock.n(HR.04.noUC[[2]]$pstk)*stock.wt(HR.04.noUC[[2]]$pstk)*exp(-(harvest(HR.04.noUC[[2]]$pstk)*survey_time + m(HR.04.noUC[[2]]$pstk)*survey_time)))
perceived <- quantSums(HR.04.noUC[[2]]$idx@index*stock.wt(HR.04.noUC[[2]]$pstk))
STK_perception <- FLQuants(perceived=perceived, true=trueBio)
STK_perc_DF <- data.table(as.data.frame(STK_perception))

# Plot with standard error
ggplot(STK_perc_DF, aes(year, data, color=qname, fill=qname)) + 
  stat_summary(fun.y="mean", geom="line", size=1.2) + 
  stat_summary(fun.data="mean_se", geom="ribbon", alpha=.3, colour=NA)
survey_time=5/12


plot(HR.02.noUC[[1]]$pstk)

save(HR.01.noUC, HR.02.noUC, HR.03.noUC, HR.04.noUC, file = paste0("Output/Piera/HCR_HR_NoUC.RData"))
save(HR.02.noUC.20idx, file = paste0("Output/Piera/HCR_HR20_noUC_20idx.RData"))

# Create FLSTocks for OM1 and OM2 
STKs_OM1 <- FLStocks(OM1.01 = HR.01.noUC[[1]]$pstk, OM1.02=HR.02.noUC[[1]]$pstk,OM1.03=HR.03.noUC[[1]]$pstk,OM1.04=HR.04.noUC[[1]]$pstk)
STKs_OM2 <- FLStocks(OM2.01=HR.01.noUC[[2]]$pstk,OM2.02=HR.02.noUC[[2]]$pstk,OM2.03=HR.03.noUC[[2]]$pstk,OM2.04=HR.04.noUC[[2]]$pstk)
# stk_20perc_abrupt <- stk
# stk_20perc_uncCap <- stk
png("./output/Piera/plots/OM1_HCR_HR.png", width=4000, height=3500, res=300)
plot(STKs_OM1) + theme_light(25)
dev.off()
png("./output/Piera/plots/OM2_HCR_HR.png", width=4000, height=3500, res=300)
plot(STKs_OM2) + theme_light(25)
dev.off()



########################################################
# ======================================================
# FIXED PERCENTAGE OF INDEX BIOMASS
# SENSITIVITY TO CHECK EFFECT OF OVERESTIMATING BIOMASS WITH IDX
# ======================================================
# percentages to be tested for HS
perc <- 0.15
survey_time <- 5/12
# Uncertainty cap yes or no?
uncertainty_Cap = FALSE

HR.015.noUC <- foreach(x = OM_list, .packages = required_pckgs) %dopar% {
  
  pstk <- x$stk
  idx <- x$observations$idx@.Data[[1]]
  stk <- x$stk
  sr.om <- x$sr.om
  sr.om.res <- x$sr.om.res
  
  TAC <- FLQuant(NA, dimnames=list(TAC="all", year=c(x$dy,x$vy), iter=1:x$it))
  TAC[,ac(x$dy)] <- catch(stk)[,ac(x$dy)]
  TAC[,ac(x$iy)] <- TAC[,ac(x$dy)] #assume same TAC in the first intermediate year
  ctrl <- getCtrl(c(TAC[,ac(x$iy)]), "catch", x$iy, x$it)
  
  for(i in x$vy[-length(x$vy)]){
    ## i <- vy[-length(vy)][1]
    print(i)
    gc()
    ay <- an(i)
    cat(i, "> ")
    vy0 <- 1:(ay-x$y0) # data years (positions vector) - one less than current year
    sqy <- ac((ay-1):(ay-x$nsqy)) # years for status quo computations
    
    stk0 <- pstk[,vy0]
    #catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
    ## note that vy0 is changing below so index is being updated
    # for (index_counter in 1:length(idx)){
    idx0 <- idx[,vy0]
    #index(idx[[index_counter]])[,i] <- stock.n(stk)[,i]*index.q(idx[[index_counter]])[,i] #+ 1
    # calculation at time of the survey
    index(idx)[,i] <- index.q(idx)[,i]*stock.n(pstk)[,i]*exp(-(harvest(pstk)[,i]*survey_time+m(pstk)[,1]*survey_time)) #+ 1
    # Calculate index tot bio
    idx.bio <- quantSums(index(idx)[,i]*stock.wt(pstk)[,i])
    #}
    ##
    
    # Calculate TAC from index
    advice <- idx.bio*perc
    #   rmodel <- sr
    # fmod <- ~ s(age, k=4) + s(year, k=16)
    # fit <- sca(stk0, FLIndices(idx0), fmodel=fmod, srmodel=rmodel)
    # stk0 <- stk0 + fit
    # fwd control
    fsq0 <- yearMeans(fbar(pstk)[,sqy])
    csq0 <- yearMeans(landings(pstk)[,sqy])
    dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max")) #, ay + 1
    arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
    ## ftrg.vec <- rep(ftrg, it) ## original
    refpt <- data.frame(ssb = 1, harvest = 1, catch=1)
    # Estimate catch vector
    catch.vec <- an(advice)
    #Bescape <- blim
    if(uncertainty_Cap == TRUE){
      new_catchVec <- hcr.check.catch(catch(stk0)[, ac(an(i) - 1)], catch.vec=catch.vec)}else {
        new_catchVec=catch.vec
      }
    #TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
    # OM proj
    TAC[,ac(ay+1)] <- catch.vec
    
    # apply the TAC to the operating model stock
    ctrl <- getCtrl(c(TAC[,ac(ay+1)]), "catch", ay+1, it)
    # ##ctrl@trgtArray[,"val",] <- c(TAC[,ac(ay+1)]) #+ BB[,ac(ay)])
    pstk <- fwd(pstk, ctrl=ctrl, sr=sr.om, sr.residuals = sr.om.res, sr.residuals.mult = TRUE) #, exp(sr.res[,ac(ay:(ay+1))]) sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
    #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
  }
  
  x$idx <- idx
  x$pstk <- pstk
  return(x)
  
}


save(HR.015.noUC, HR.016.noUC, HR.017.noUC, HR.018.noUC, HR.019.noUC, file = paste0("Output/Piera/HCR_HR_noUC_015-019.RData"))

# Create FLSTocks for OM1 and OM2 
STKs_OM1 <- FLStocks(OM1.01 = HR.01.noUC.30percIDX[[1]]$pstk, OM1.02=HR.02.noUC.30percIDX[[1]]$pstk,OM1.03=HR.03.noUC.30percIDX[[1]]$pstk,OM1.04=HR.04.noUC.30percIDX[[1]]$pstk)
STKs_OM2 <- FLStocks(OM2.01=HR.01.noUC[[2]]$pstk,OM2.02=HR.02.noUC[[2]]$pstk,OM2.03=HR.03.noUC[[2]]$pstk,OM2.04=HR.04.noUC[[2]]$pstk)
# stk_20perc_abrupt <- stk
# stk_20perc_uncCap <- stk
png("./output/Piera/plots/OM1_HCR_HR.png", width=4000, height=3500, res=300)
plot(STKs_OM1) + theme_light(25)
dev.off()
png("./output/Piera/plots/OM2_HCR_HR.png", width=4000, height=3500, res=300)
plot(STKs_OM2) + theme_light(25)
dev.off()




########################################################
# ======================================================
# RULE 1 OVER 2 
# ======================================================
########################################################
# Uncertainty cap?
uncertainty_Cap = FALSE 

# Start simulations
HCR.1o2.noUC <- foreach(x = OM_list, .packages = required_pckgs) %dopar% {
  
  pstk <- x$stk
  idx <- x$observations$idx@.Data[[1]]
  stk <- x$stk
  sr.om <- x$sr.om
  sr.om.res <- x$sr.om.res
  
  TAC <- FLQuant(NA, dimnames=list(TAC="all", year=c(x$dy,x$vy), iter=1:x$it))
  TAC[,ac(x$dy)] <- catch(stk)[,ac(x$dy)]
  TAC[,ac(x$iy)] <- TAC[,ac(x$dy)] #assume same TAC in the first intermediate year
  ctrl <- getCtrl(c(TAC[,ac(x$iy)]), "catch", x$iy, x$it)
  
  # create idx.bio
  idx.bio <- FLQuant(dim=c(1,dim(observations$idx@.Data[[1]])[-1]))
  idx.bio[,1:24] <- quantSums(index(idx)[,1:24]*stock.wt(stk)[,1:24])
  scaler <- list()
  
  # Start with loop
  for(i in x$vy[-length(x$vy)]){
    ## i <- vy[-length(vy)][1]
    print(i)
    gc()
    ay <- an(i)
    cat(i, "> ")
    vy0 <- 1:(ay-x$y0) # data years (positions vector) - one less than current year
    sqy <- ac((ay-1):(ay-x$nsqy)) # years for status quo computations
    
    stk0 <- pstk[,vy0]
    #catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
    ## note that vy0 is changing below so index is being updated
    # for (index_counter in 1:length(idx)){
    idx0 <- idx[,vy0]
    #index(idx[[index_counter]])[,i] <- stock.n(stk)[,i]*index.q(idx[[index_counter]])[,i] #+ 1
    # calculation at time of the survey
    index(idx)[,i] <- index.q(idx)[,i]*stock.n(pstk)[,i]*exp(-(harvest(pstk)[,i]*survey_time+m(pstk)[,1]*survey_time)) #+ 1
    # Calculate index tot bio
    idx.bio[,i] <- quantSums(index(idx)[,i]*stock.wt(pstk)[,i])
    #}
    # Calculate 1 over 2 rule
    scaler[[i]] <- an(idx.bio[,i]/yearMeans(idx.bio[,ac(c(ay-(1:2)))]))
    if(uncertainty_Cap == TRUE){
      new_scaler <- hcr.check.scaler(scaler=scaler[[i]])} else {
        new_scaler=scaler[[i]]
      }
    # Calculate catch advice
    catch.vec <- an(catch(pstk[,ac(ay-1)]))*new_scaler
    # fwd control
    fsq0 <- yearMeans(fbar(pstk)[,sqy])
    csq0 <- yearMeans(landings(pstk)[,sqy])
    dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max")) #, ay + 1
    arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
    ## ftrg.vec <- rep(ftrg, it) ## original
    refpt <- data.frame(ssb = 1, harvest = 1, catch=1)
    # Populate TAC object with target Catch
    TAC[,ac(ay+1)] <- catch.vec
    # apply the TAC to the operating model stock
    ctrl <- getCtrl(c(TAC[,ac(ay+1)]), "catch", ay+1, it)
    # Project forward applying the HCR
    pstk <- fwd(pstk, ctrl=ctrl, sr=sr.om, sr.residuals = sr.om.res, sr.residuals.mult = TRUE) #, exp(sr.res[,ac(ay:(ay+1))]) sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
    #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
  }
  x$idx.bio <- idx.bio
  x$pstk <- pstk
  return(x)
  
}

plot(stk)



########################################################
# ======================================================
# RULE 1 OVER 2 - with UNCERTAINTY CAP
# ======================================================
########################################################
uncertainty_Cap = TRUE
#HR.1o2.UC <- list()

HCR.1o2.UC <- foreach(x = OM_list, .packages = required_pckgs) %dopar% {
  
  pstk <- x$stk
  idx <- x$observations$idx@.Data[[1]]
  stk <- x$stk
  sr.om <- x$sr.om
  sr.om.res <- x$sr.om.res
  
  TAC <- FLQuant(NA, dimnames=list(TAC="all", year=c(x$dy,x$vy), iter=1:x$it))
  TAC[,ac(x$dy)] <- catch(stk)[,ac(x$dy)]
  TAC[,ac(x$iy)] <- TAC[,ac(x$dy)] #assume same TAC in the first intermediate year
  ctrl <- getCtrl(c(TAC[,ac(x$iy)]), "catch", x$iy, x$it)
  
  # create idx.bio
  idx.bio <- FLQuant(dim=c(1,dim(observations$idx@.Data[[1]])[-1]))
  idx.bio[,1:24] <- quantSums(index(idx)[,1:24]*stock.wt(stk)[,1:24])
  scaler <- list()
  
  # Start with loop
  for(i in x$vy[-length(x$vy)]){
    ## i <- vy[-length(vy)][1]
    print(i)
    gc()
    ay <- an(i)
    cat(i, "> ")
    vy0 <- 1:(ay-x$y0) # data years (positions vector) - one less than current year
    sqy <- ac((ay-1):(ay-x$nsqy)) # years for status quo computations
    
    stk0 <- pstk[,vy0]
    #catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
    ## note that vy0 is changing below so index is being updated
    # for (index_counter in 1:length(idx)){
    idx0 <- idx[,vy0]
    #index(idx[[index_counter]])[,i] <- stock.n(stk)[,i]*index.q(idx[[index_counter]])[,i] #+ 1
    # calculation at time of the survey
    index(idx)[,i] <- index.q(idx)[,i]*stock.n(pstk)[,i]*exp(-(harvest(pstk)[,i]*survey_time+m(pstk)[,1]*survey_time)) #+ 1
    # Calculate index tot bio
    idx.bio[,i] <- quantSums(index(idx)[,i]*stock.wt(pstk)[,i])
    #}
    # Calculate 1 over 2 rule
    scaler[[i]] <- an(idx.bio[,i]/yearMeans(idx.bio[,ac(c(ay-(1:2)))]))
    if(uncertainty_Cap == TRUE){
      new_scaler <- hcr.check.scaler(scaler=scaler[[i]])} else {
        new_scaler=scaler[[i]]
      }
    # Calculate catch advice
    catch.vec <- an(catch(pstk[,ac(ay-1)]))*new_scaler
    # fwd control
    fsq0 <- yearMeans(fbar(pstk)[,sqy])
    csq0 <- yearMeans(landings(pstk)[,sqy])
    dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max")) #, ay + 1
    arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
    ## ftrg.vec <- rep(ftrg, it) ## original
    refpt <- data.frame(ssb = 1, harvest = 1, catch=1)
    # Populate TAC object with target Catch
    TAC[,ac(ay+1)] <- catch.vec
    # apply the TAC to the operating model stock
    ctrl <- getCtrl(c(TAC[,ac(ay+1)]), "catch", ay+1, it)
    # Project forward applying the HCR
    pstk <- fwd(pstk, ctrl=ctrl, sr=sr.om, sr.residuals = sr.om.res, sr.residuals.mult = TRUE) #, exp(sr.res[,ac(ay:(ay+1))]) sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
    #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
  }
  
  x$pstk <- pstk
  x$idx.bio <- idx.bio
  return(x)
  
}

plot(stk)

save(HCR.1o2.noUC, HCR.1o2.UC, file = paste0("Output/Piera/HCR_1o2.RData"))


# Create FLSTocks for OM1 and OM2 
STKs_OM1_1o2 <- FLStocks(OM1.1o2.noUC = HCR.1o2.noUC[[1]]$pstk, OM1.1o2.UC=HCR.1o2.UC[[1]]$pstk)
STKs_OM2_1o2 <- FLStocks(OM2.1o2.noUC = HCR.1o2.noUC[[2]]$pstk, OM2.1o2.UC=HCR.1o2.UC[[2]]$pstk)
# stk_20perc_abrupt <- stk
# stk_20perc_uncCap <- stk
png("./output/Piera/plots/OM1_HCR_1o2.png", width=4000, height=3500, res=300)
plot(STKs_OM1_1o2) + theme_light(25)
dev.off()
png("./output/Piera/plots/OM2_HCR_1o2.png", width=4000, height=3500, res=300)
plot(STKs_OM2_1o2) + theme_light(25)
dev.off()



idxBio <- data.table(as.data.frame(HCR.1o2.noUC[[1]]$idx.bio))
totBio <- data.table(as.data.frame(quantSums(HCR.1o2.noUC[[1]]$pstk@mat*HCR.1o2.noUC[[1]]$pstk@stock.wt*HCR.1o2.noUC[[1]]$pstk@stock.n*exp(-(HCR.1o2.noUC[[1]]$pstk@harvest*survey_time+HCR.1o2.noUC[[1]]$pstk@m*survey_time)))))
ssb <- data.table(as.data.frame(ssb(HCR.1o2.noUC[[1]]$pstk)))
#totBio <- data.table(as.data.frame(quantSums((HCR.1o2.noUC[[2]]$pstk@stock.n)*(HCR.1o2.noUC[[2]]$pstk@stock.wt))))
idxBio[,type:="perceived"]
totBio[,type:="true"]
ssb[,type:="true"]
Bio <- rbind(idxBio[,-1], totBio[,-1])
BioSum <- Bio %>%
  group_by(year,type) %>%
  summarize(mean_size = mean(data), stDev=sd(data))
BioSum <- data.table(BioSum)
BioSum[, se:=1.96*stDev/sqrt(500)]
BioSum[,c("lCI","uCI"):=list(mean_size-se,mean_size+se)]
# Plot
ggplot(BioSum[year!=50,], aes(year, mean_size, col=type, group=type)) + 
  geom_line() + 
  scale_x_continuous(breaks=seq(1,50,5)) + ylab("t") +
  geom_ribbon(aes(ymax=uCI,ymin=lCI,fill=type),alpha=.25) +
  theme_bw(22) +
  theme(legend.position = "top", legend.title = element_blank(), legend.direction = "horizontal")

ggplot(Bio)



########################################################
# ======================================================
# HARVEST RATE WITH OBS ERROR DIRECTLY ON IDX.BIO
# ======================================================
########################################################
# percentages to be tested for HS
OM_list_obsErr <- OM_list[c(2,4)]
perc <- 0.2
survey_time <- 5/12
# Uncertainty cap yes or no?
uncertainty_Cap = FALSE

HR.02.ObsErr <- foreach(x = OM_list_obsErr, .packages = required_pckgs) %dopar% {
  
  pstk <- x$stk
  idx <- x$observations$idx@.Data[[1]]
  stk <- x$stk
  sr.om <- x$sr.om
  sr.om.res <- x$sr.om.res
  
  TAC <- FLQuant(NA, dimnames=list(TAC="all", year=c(x$dy,x$vy), iter=1:x$it))
  TAC[,ac(x$dy)] <- catch(stk)[,ac(x$dy)]
  TAC[,ac(x$iy)] <- TAC[,ac(x$dy)] #assume same TAC in the first intermediate year
  ctrl <- getCtrl(c(TAC[,ac(x$iy)]), "catch", x$iy, x$it)
  
  for(i in x$vy[-length(x$vy)]){
    ## i <- vy[-length(vy)][1]
    print(i)
    gc()
    ay <- an(i)
    cat(i, "> ")
    vy0 <- 1:(ay-x$y0) # data years (positions vector) - one less than current year
    sqy <- ac((ay-1):(ay-x$nsqy)) # years for status quo computations
    
    stk0 <- pstk[,vy0]
    #catch.n(stk0) <- catch.n(stk0) + 1 # avoid zeros
    ## note that vy0 is changing below so index is being updated
    # for (index_counter in 1:length(idx)){
    idx0 <- idx[,vy0]
    #index(idx[[index_counter]])[,i] <- stock.n(stk)[,i]*index.q(idx[[index_counter]])[,i] #+ 1
    # calculation at time of the survey
    index(idx)[,i] <- index.q(idx)[,i]*stock.n(pstk)[,i]*exp(-(harvest(pstk)[,1]*survey_time+m(pstk)[,1]*survey_time)) #+ 1
    # Calculate index tot bio
    idx.bio <- quantSums(index(idx)[,i]*stock.wt(pstk)[,i])*rlnorm(n = it, sdlog = 0.5)
    #}
    ##
    
    # Calculate TAC from index
    advice <- idx.bio*perc
    #   rmodel <- sr
    # fmod <- ~ s(age, k=4) + s(year, k=16)
    # fit <- sca(stk0, FLIndices(idx0), fmodel=fmod, srmodel=rmodel)
    # stk0 <- stk0 + fit
    # fwd control
    fsq0 <- yearMeans(fbar(pstk)[,sqy])
    csq0 <- yearMeans(landings(pstk)[,sqy])
    dnms <- list(iter=1:it, year=c(ay, ay + 1), c("min", "val", "max")) #, ay + 1
    arr0 <- array(NA, dimnames=dnms, dim=unlist(lapply(dnms, length)))
    ## ftrg.vec <- rep(ftrg, it) ## original
    refpt <- data.frame(ssb = 1, harvest = 1, catch=1)
    # Estimate catch vector
    catch.vec <- an(advice)
    #Bescape <- blim
    if(uncertainty_Cap == TRUE){
      new_catchVec <- hcr.check.catch(catch(stk0)[, ac(an(i) - 1)], catch.vec=catch.vec)}else {
        new_catchVec=catch.vec
      }
    # arr0[,,"val"] <- c(rep(NA, it), new_catchVec)
    # #arr0[,,"min"] <- c(fsq0, rep(NA, it))
    # arr0 <- aperm(arr0, c(2,3,1))
    # ctrl <- fwdControl(data.frame(year=c(ay, ay+1), quantity=c('f', 'catch'), val=NA)) # , ay + 1
    # ctrl@trgtArray <- arr0
    ## 
    #stkTmp <- stf(stk0, 2)
    #stkTmp <- fwd(stkTmp, ctrl=ctrl, sr=sr.om, sr.residuals = sr.om.res, sr.residuals.mult = TRUE) #, exp(sr.res[,ac(ay:(ay+1))])sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE
    ## USING F
    #TAC[,ac(ay+1)] <- catch(stkTmp)[,ac(ay+1)]
    # OM proj
    TAC[,ac(ay+1)] <- catch.vec
    
    # apply the TAC to the operating model stock
    ctrl <- getCtrl(c(TAC[,ac(ay+1)]), "catch", ay+1, it)
    
    # 
    # ctrl@target <- ctrl@target[2,]
    # ## original was catch
    # ##ctrl@target[,"quantity"] <- "catch"
    # ctrl@trgtArray <- ctrl@trgtArray[2,,,drop=FALSE]
    # ## original was catch
    # ##ctrl@trgtArray[,"val",] <- c(TAC[,ac(ay+1)]) #+ BB[,ac(ay)])
    pstk <- fwd(pstk, ctrl=ctrl, sr=sr.om, sr.residuals = sr.om.res, sr.residuals.mult = TRUE) #, exp(sr.res[,ac(ay:(ay+1))]) sr.residuals = exp(sr.res[,ac(ay:(ay+1))]), sr.residuals.mult = TRUE) #, sr.residuals = exp(sr.res[,ac(ay+1)]), sr.residuals.mult = TRUE
    #BB[,ac(ay+1)] <- iterMedians(TAC[,ac(ay+1)]) - catch(pstk)[,ac(ay+1)]
  }
  
  idx_new <- idx
  x$pstk <- pstk
  x$idx_new <- idx_new
  return(x)
  
}

save(HR.02.ObsErr, file = paste0("Output/Piera/HCR_HR02_obsErr.RData"))

HR.sim <- mget(load("Output/Piera/HCR_HR_NoUC.RData"))
# Create FLSTocks for OM1 and OM2 
STKs_OM1_HS_02_obsErr <- FLStocks(OM1.HR.02.obsErr = HR.02.ObsErr[[1]]$pstk, OM1.HR.02=HR.sim[[2]][[1]]$pstk)
STKs_OM2_HS_02_obsErr <- FLStocks(OM2.HR.02.obsErr = HR.02.ObsErr[[2]]$pstk, OM2.HR.02=HR.sim[[2]][[2]]$pstk)
# stk_20perc_abrupt <- stk
# stk_20perc_uncCap <- stk
png("./output/Piera/plots/OM1_HR_02_obsErr.png", width=4000, height=3500, res=300)
plot(STKs_OM1_HS_02_obsErr) + theme_light(25)
dev.off()
png("./output/Piera/plots/OM2_HR_02_obsErr.png", width=4000, height=3500, res=300)
plot(STKs_OM2_HS_02_obsErr) + theme_light(25)
dev.off()

# Check idx 
IDX <- window(HR.02.ObsErr[[2]]$idx_new, end=49)
IDX.BIO <- quantSums(index(IDX)*catch.wt(window(HR.02.ObsErr[[2]]$pstk, end=49)))
SSB <- ssb(window(HR.02.ObsErr[[2]]$pstk, end=49))

Cor <- data.table(cbind(as.data.frame(IDX.BIO),as.data.frame(SSB)$data))
names(Cor)[8] <- "SSB"

plot(Cor[,cor(data,SSB), by=iter])
