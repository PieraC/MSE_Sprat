### ------------------------------------------------------------------------ ###
### set-up OMs for MSE simulation ####
### set additional parameters for simulation
### From Simon Fisher git hub
### ------------------------------------------------------------------------ ###
### the R session needs to be restarted as FLash and FLasher environments
### interfere
rm(list=ls())

required_pckgs <- c("FLash", "FLAssess", "ggplotFL", "FLBRP", "data.table")
### save as object in order to avoid output to screen
. <- lapply(required_pckgs, function(x){
  suppressMessages(library(x, character.only = TRUE))
})
library(doParallel)

cl <- makeCluster(28)
registerDoParallel(cl)


### load additional functions
### select R scripts from functions folder
load_files <- list.files("functions/")
load_files <- load_files[grepl(pattern = "*.R$", x = load_files)]
### source the scripts
invisible(lapply(paste0("functions/", load_files), source))

set.seed(0)

#houseKeeping <- OMs
### ------------------------------------------------------------------------ ###
### load OMs
### ------------------------------------------------------------------------ ###

### load list with OMs
OM_list <- readRDS("Input/Piera/OMs.rds")
# housekeeping OM_names
OM_names <- names(OM_list)
# > OM_names
# [1] "OM1_0.3" "OM1_0.5" "OM2_0.3" "OM2_0.5" "OM3_0.3" "OM3_0.5"

#for(i in 1:length(OM_list)) names(OM_list[[i]]) <- c("lhpar", "sr", "stk")
### names to numbers of OMs
names(OM_list) <- seq_along(OM_list)


### ------------------------------------------------------------------------ ###
### specify dimensions of simulation ####
### ------------------------------------------------------------------------ ###
it <- dims(OM_list[[1]]$stk)$iter # iterations
fy <- dims(OM_list[[1]]$stk)$maxyear + 25 # final year
y0 <- range(OM_list[[1]]$stk)[["minyear"]] # initial data year
dy <- range(OM_list[[1]]$stk)[["maxyear"]] # final data year
iy <- dims(OM_list[[1]]$stk)$maxyear # initial year of projection (also intermediate)
ny <- fy - iy + 1 # number of years to project from intial year
nsqy <- 3 # number of years to compute status quo metrics
vy <- ac(iy:fy) # vector of years to be projected


### ------------------------------------------------------------------------ ###
### "loop" through all stocks ####
### ------------------------------------------------------------------------ ###

OM_list <- foreach(x = seq_along(OM_list), .export = ls(), .packages = required_pckgs) %dopar% {
  
  set.seed(0)
  ### ---------------------------------------------------------------------- ###
  ### create residuals for SRR
  ### ---------------------------------------------------------------------- ###
  
  ### create residuals
  # different values for sd
  if(substr(OM_names[x],5,7)=="0.3"){
    srbh.res <- rlnoise(n = it, FLQuant(0, dimnames = list(year = vy)), sd = 0.3, b = 0)
  } else {
    srbh.res <- rlnoise(n = it, FLQuant(0, dimnames = list(year = vy)), sd = 0.5, b = 0)
  } 
  # add residuals to OM model and assign name
  OM_list[[x]]$srbh.res <- srbh.res
  names(OM_list[[x]])[names(OM_list[[x]]) == "sr"] <- "srbh"
  
  ### ---------------------------------------------------------------------- ###
  ### add life-history parameters to stock ####
  ### ---------------------------------------------------------------------- ###
  ### extract lhpar
  lhpar <- OM_list[[x]][["lhpar"]]
  ### add some more parameters
  ### natural mortality and max age
  add_pars <- FLPar(M = mean(m(OM_list[[x]][["stk"]][, 1])),
                    max_age = range(OM_list[[x]][["stk"]])[["max"]])
  ### add
  lhpar <- rbind2(lhpar, add_pars)
  ### change some names
  dimnames(lhpar)$params[dimnames(lhpar)$params %in% c("linf", "k")] <- 
    c("L_inf", "K")
  ### propagate
  lhpar <- propagate(lhpar, dims(OM_list[[1]][["stk"]])$iter)
  ### save as attribute in stk
  attr(OM_list[[x]][["stk"]], "lhpar") <- lhpar
  
  ### ---------------------------------------------------------------------- ###
  ### get reference points from FLBRP
  ### ---------------------------------------------------------------------- ###
  
  #attr(OM_list[[x]]$stk, "refpts") <- refpts(OM_list[[x]]$brp)
  #attr(OM_list[["stk_one_way"]], "refpts") <- refpts(OM_list[["brp"]])
  
  ### ---------------------------------------------------------------------- ###
  ### set up operating model: extend ####
  ### ---------------------------------------------------------------------- ###
  OM_list[[x]]$stk <- stf(OM_list[[x]]$stk, fy-dy, nsqy, nsqy)
  OM_list[[x]]
}

names(OM_list) <- OM_names

# Housekeeping just in case...
housekeeping <- OM_list
### ------------------------------------------------------------------------ ###
### create initial observations and add uncertainty ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### "loop" trough stocks
### ------------------------------------------------------------------------ ###

OM_list <- foreach(x = OM_list, .export = setdiff(ls(), "OM_list"),
                   .packages = required_pckgs) %dopar% {
  
  ### set seed 
  set.seed(1)
  
  ### ---------------------------------------------------------------------- ###
  ### create index, based on stock object
  ### ---------------------------------------------------------------------- ###
  survey_time <- 5/12
  ### get ages
  ages <- an(dimnames(stock.n(x$stk))[["age"]])
  ### define model for selectivity: logistic function
  ### inflection point of curve = 10% of max age
  q_model <- FLModelSim(model = ~ max_q/(1+exp(-steepness*(age - age50))), 
                       params = FLPar(max_q = 1.5, steepness = 1.1, 
                                      age50 = 0.1))#max(ages)/max(ages+1)))#max(ages)/50))
  ### model selectivity
  q_modeled <- predict(q_model, age = ages)
  
  ### create index template
  idx <- FLIndex(index = stock.n(x$stk))
  
  # CHECK OUT COHORT IN PELTIC
  ### insert catchability (observation error added later on)
  index.q(idx) <- c(q_modeled)
  
  ### calculate historical index values
  index(idx) <- index.q(idx) * stock.n(x$stk) * exp(-(harvest(x$stk)*survey_time + m(x$stk)*survey_time))# * stock.wt(x$stk)
  
  # plot(FLQuants(idx = quantSums(index(idx)),
  #               tsb = quantSums(stock.n(x$stk)*stock.wt(x$stk)),
  #               ssb = ssb(x$stk)))
  
  ### save in observations object
  x$observations <- list(idx = FLIndices(idx = idx))
  
  ### ---------------------------------------------------------------------- ###
  ### estimate some required reference points
  ### ---------------------------------------------------------------------- ###
  ### I_lim = lowest observed index value
  I_lim <- apply(quantSums(index(x$observations[[1]][["idx"]])), 6, min, 
                na.rm = TRUE)
  ### I_trigger not defined yet
  I_trigger <- NA
  
  ### add noise
  #I_lim <- I_lim * rlnorm(n = length(I_lim), sdlog = 0.05)
  #I_trigger <- I_trigger * rlnorm(n = length(I_trigger), sdlog = 0.05)
  
  ### add reference points to stk
  attr(x = x$stk, which = "refpts") <- FLPar(I_lim = I_lim,
                                            I_trigger = I_trigger)

  return(x)
  
}

### ---------------------------------------------------------------------- ###
### calculate F0.1 with YPR
### ---------------------------------------------------------------------- ###
OM_list <- foreach(x = OM_list, .packages = required_pckgs) %dopar% {

  ### needed for catch rule 3.2.1 from WKMSYCat34, option c for factor f

  ### run YPR for range of F's
  F_list <- seq(0, 3, 0.01)

  ### inverse von Bertalanffy growth function
  inv_vB <- function(L, L_inf, K, t0){
    res <- -log(1 - (L / L_inf)) / K + t0
    return(res)
  }

  ### length at first capture
  L_c <- 8# From self sampling data yearMeans(calc_Lc(attr(x$stk, "catch_len")))
  # ### convert into age  ### calculate t_c from L_c
  t_c <- inv_vB(L = L_c,
                L_inf = attr(x$stk, "lhpar")["L_inf"],
                K = attr(x$stk, "lhpar")["K"], t0 = attr(x$stk, "lhpar")["t0"])
  t_c <- round(t_c)

  ### calculate F0.1 for 1 iteration
  F0.1 <- lapply((1:dims(x$stk)$iter)[1], function(y){

    ### conduct YPR for list of F's
    res_list <- lapply(F_list, function(f){

      YPR(a =	attr(x$stk, "lhpar")["a", y], b =	attr(x$stk, "lhpar")["b", y],
          Linf =	attr(x$stk, "lhpar")["L_inf", y],
          K =	attr(x$stk, "lhpar")["K", y], t0 =	attr(x$stk, "lhpar")["t0", y],
          tc =	c(FLCore::iter(t_c, y)),
          M =	attr(x$stk, "lhpar")["M", y], F =	f,
          max_age = attr(x$stk, "lhpar")["max_age", y])

    })
    ### format
    res_list <- do.call(rbind, res_list)

    ### search for F0.1
    ### get slope (step size is equal between all elements)
    slope <- diff(res_list)

    ### find first F where slope is < 0.1 * initial slope
    res <- F_list[which(slope < 0.1 * slope[1])[1]]

    return(res)

  })
  F0.1 <- unlist(F0.1)

  ### save in refpts
  attr(x$stk, "refpts") <- rbind2(attr(x$stk, "refpts"),
                                  FLPar(F0.1YPR = F0.1, iter = dims(x$stk@lhpar)$iter))

  return(x)

}


# ### ------------------------------------------------------------------------ ###
### observation error
### ------------------------------------------------------------------------ ###

OM_list <- foreach(x = OM_list, .packages = required_pckgs) %dopar% {
  
  ### set seed 
  set.seed(1)
  # survey time
  survey_time = 5/12
  ### ---------------------------------------------------------------------- ###
  ### add index uncertainty
  ### ---------------------------------------------------------------------- ###
  
  ### get index
  idx <- x$observations$idx$idx

  ### add uncertainty to catchability
  ### log-normal noise, cv = 0.5
  index.q(idx) <- index.q(idx) * rlnorm(n = length(index.q(idx)), sdlog = 0.5)
  
  ### update index values
  index(idx) <- index.q(idx) * stock.n(x$stk) * exp(-(harvest(x$stk)*survey_time+m(x$stk)*survey_time)) #+ 1# * stock.wt(x$stk)
  
  ### save in observations object
  x$observations <- list(idx = FLIndices(idx = idx))
  
  ### ---------------------------------------------------------------------- ###
  ### update index reference points
  ### ---------------------------------------------------------------------- ###
  ### uncertainty already implemented in index values
  ### the new minimum values already include this uncertainty
  
  attr(x$stk, "refpts")["I_lim"] <- apply(quantSums(index(idx)), 6, min, 
                                          na.rm = TRUE)
  
  ### ---------------------------------------------------------------------- ###
  ### add uncertainty to life-history parameters and reference points
  ### ---------------------------------------------------------------------- ###
  ### reference points
  # attr(x$stk, "refpts")["LFeFmsy"] <- attr(x$stk, "refpts")["LFeFmsy"] * 
  #   rlnorm(n = length(attr(x$stk, "refpts")["LFeFmsy"]), sdlog = 0.1)
  # attr(x$stk, "refpts")["LFeM"] <- attr(x$stk, "refpts")["LFeM"] * 
  #   rlnorm(n = length(attr(x$stk, "refpts")["LFeM"]), sdlog = 0.1)
  # attr(x$stk, "refpts")["F_MSY"] <- attr(x$stk, "refpts")["F_MSY"] * 
  #   rlnorm(n = length(attr(x$stk, "refpts")["F_MSY"]), sdlog = 0.1)
  attr(x$stk, "refpts")["F0.1YPR"] <- attr(x$stk, "refpts")["F0.1YPR"] *
    rlnorm(n = length(attr(x$stk, "refpts")["F0.1YPR"]), sdlog = 0.1)
  # attr(x$stk, "refpts")["I_F_proxy"] <- attr(x$stk, "refpts")["I_F_proxy"] * 
  #   rlnorm(n = length(attr(x$stk, "refpts")["I_F_proxy"]), sdlog = 0.1)
  # 
  ### lhist
  attr(x$stk, "lhpar")["L_inf"] <- attr(x$stk, "lhpar")["L_inf"] * 
    rlnorm(n = length(attr(x$stk, "lhpar")["L_inf"]), sdlog = 0.1)
  attr(x$stk, "lhpar")["K"] <- attr(x$stk, "lhpar")["K"] * 
    rlnorm(n = length(attr(x$stk, "lhpar")["K"]), sdlog = 0.1)
  attr(x$stk, "lhpar")["t0"] <- attr(x$stk, "lhpar")["t0"] * 
    rlnorm(n = length(attr(x$stk, "lhpar")["t0"]), sdlog = 0.1)
  attr(x$stk, "lhpar")["a"] <- attr(x$stk, "lhpar")["a"] *
    rlnorm(n = length(attr(x$stk, "lhpar")["a"]), sdlog = 0.1)
  attr(x$stk, "lhpar")["b"] <- attr(x$stk, "lhpar")["b"] *
    rlnorm(n = length(attr(x$stk, "lhpar")["b"]), sdlog = 0.1)
  attr(x$stk, "lhpar")["M"] <- attr(x$stk, "lhpar")["M"] * 
    rlnorm(n = length(attr(x$stk, "lhpar")["M"]), sdlog = 0.5)
  # 
  return(x)
  
}

### ------------------------------------------------------------------------ ###
### save
### ------------------------------------------------------------------------ ###

. <- foreach(x = seq_along(OM_list)) %dopar% {
  
  stk <- OM_list[[x]]$stk
  observations <- OM_list[[x]]$observations
  sr.om <- OM_list[[x]]$srbh
  sr.om.res <- OM_list[[x]]$srbh.res
  
  save(stk, observations, sr.om, sr.om.res, it, fy, y0, dy, iy, ny, nsqy, vy, #
       file = paste0("Input/Piera/stocks/",  
                     x, "_new.RData"))
  
}


stopCluster(cl)

