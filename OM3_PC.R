### ------------------------------------------------------------------------ ###
### Check OMs with plots ####
### set additional parameters for simulation
### 
### ------------------------------------------------------------------------ ###
### the R session needs to be restarted as FLash and FLasher environments
### interfere
library(FLife)
library(FLBRP)
library(FLasher)
library(foreach)
library(data.table)
library(ggplot2)

### load additional functions
### select R scripts from functions folder
load_files <- list.files("functions/")
load_files <- load_files[grepl(pattern = "*.R$", x = load_files)]
### source the scripts
invisible(lapply(paste0("functions/", load_files), source))

# load OMs
OMs <- list.files("Input/Piera/stocks/", pattern=".RData")

OM_list <- list()
for(i in OMs) OM_list[[i]] <- mget(load(paste0("Input/Piera/stocks/", i))) 


# ================================================
# CORRELATION STOCK NUMBERS vs IDX NUMBERS IN OMs
# ================================================
# Check correlation between stock.n and idx.n for main ages (0, 1, 2)
STK.n <- data.table()
IDX.n <- data.table()

for(i in 1:length(OM_list)){
  STK.n <- rbind(STK.n, data.table(as.data.frame(stock.n(window(OM_list[[i]]$stk, end=25)))))
  IDX.n <- rbind(IDX.n, data.table(as.data.frame(window(index(OM_list[[i]]$observations$idx[[1]]),end=25))))
}
# add column for OMs
# my OMs 
OM_names <- c("OM1","OM2","OM3","OM4","OM5","OM6")
STK.n[,OM:=c(rep(OM_names[1:2],each=75000),rep(OM_names[3:4],each=100000), rep(OM_names[5:6],each=75000))]
# merge into single dataframe and rename cols
num <- cbind(STK.n, IDX.n[,data])
names(num)[c(7,9)] <- c("STK", "IDX")
# select random iterations on which check stats
myIt <- seq(1, 500, 10)
# plot for age 1 STK vs IDX 
ggplot(num[age==1 & iter %in% myIt,], aes(STK, IDX, group=iter, col=iter)) + geom_point()

# Estimate correlation
myCor_age0 <- num[age==0,cor(STK,IDX),by=.(iter,OM)]
myCor_age1 <- num[age==1,cor(STK,IDX),by=.(iter,OM)]
myCor_age2 <- num[age==2,cor(STK,IDX),by=.(iter,OM)]
# assign id and bind together
myCor_age0$age <- "age0"
myCor_age1$age <- "age1"
myCor_age2$age <- "age2"
myCor <- rbind(myCor_age0,myCor_age1,myCor_age2)
names(myCor)[3] <- "cor"
# calculate summary statistics
myCor_sum <- summarySE(myCor, measurevar="cor", groupvars=c("OM","age"))

# Add column with recVar specifications for plotting
myCor_sum$OM_recVar <- substr(myCor_sum$OM, 5,6)
# Plot with CI
# Use 95% confidence intervals

P.COR <- ggplot(myCor_sum, aes(x=age, y=cor, colour=OM, group=OM)) + 
  geom_errorbar(aes(ymin=cor-ci, ymax=cor+ci), size=1.2,width=.2, position=position_dodge(width=0.5)) +
  geom_point(size=3, position=position_dodge(width=0.5)) + theme_bw(25) + 
  #theme(legend.position = ) +
  facet_grid(.~OM_recVar)

# save plot
png("./Input/Piera/OMs_CORR_IdxStk.png", width=3000, height=1000, res=250)
P.COR
dev.off()



# ================================================
# PLOT IDX BIO AND GET CV
# ================================================
IDX.bio <- list()
IDX.CV <- list()
for(ii in 1:length(OM_list)){
  IDX.bio[[ii]] <- quantSums(window(index(OM_list[[ii]]$observations$idx[[1]]), end=25) * stock.wt(window(OM_list[[ii]]$stk, end=25)))
  IDX.CV[[ii]] <- sqrt(iterVars(IDX.bio[[ii]]))/iterMeans(IDX.bio[[ii]])
  }


IDX.CV.df <- do.call(rbind.data.frame, lapply(IDX.CV, as.data.frame))
IDX.CV.df$OM <- rep(OM_names, each=25)

png("./Input/Piera/OMs_IDX_CV.png", width=2000, height=1500, res=200)
ggplot(IDX.CV.df, aes(year, data, col=OM, group=OM)) + geom_line(size=1.2) +
  theme_bw(25) + theme(legend.position=c(0.3,0.85), legend.direction="horizontal", legend.title=element_blank()) +
  ylab("CV")
dev.off()
png("./Input/Piera/OMs_IDX_%3d.png", width=3000, height=2500, res=250)
plot(IDX.bio) + theme_light(25) + xlab("Year") + ylab("t")
plot(IDX.CV.df) + theme_light(25) + xlab("Year") + ylab("CV") 
dev.off()


# ================================================
# PLOT CATCHABILITY
# ================================================
IDX.q <- as.data.frame(index.q(OM_list[[2]]$observations$idx[[1]][,25,,,,1]))
P.Q <- ggplot(IDX.q, aes(age,data)) + geom_line(col="blue", size=1.2) +
  theme_classic(25) + scale_y_continuous(breaks=seq(0.2, 1.4, 0.2)) +
  ylab("q")

png("./Input/Piera/OMs_IDX_Q.png", width=2500, height=2000, res=250)
P.Q
dev.off()


# ================================================
# PLOT THE TWO OPERATING MODELS 
# ================================================
OM1 <- OM_list[["1_new.RData"]]$stk
OM2 <- OM_list[["2_new.RData"]]$stk
OM3 <- OM_list[["3_new.RData"]]$stk
OM4 <- OM_list[["4_new.RData"]]$stk
OM5 <- OM_list[["5_new.RData"]]$stk
OM6 <- OM_list[["6_new.RData"]]$stk

# set my it
myIt <- seq(1, 500, 100)
OM1_rec <- as.data.frame(rec(OM1[,1:25,,,,myIt]))
OM2_rec <- as.data.frame(rec(OM2[,1:25,,,,(myIt+1)]))

OM1_ssb <- as.data.frame(ssb(OM1[,1:25,,,,myIt]))
OM2_ssb <- as.data.frame(ssb(OM2[,1:25,,,,(myIt+1)]))

P.OM1rec <- ggplot(OM1_rec, aes(year, data, col=iter)) + geom_line(size=1.2) +
  theme_light(25) + ylab("Rec") + ggtitle("OM1 rec")
P.OM2rec <- ggplot(OM2_rec, aes(year, data, col=iter)) + geom_line(size=1.2) +
  theme_light(25) + ylab("Rec")+ ggtitle("OM2 rec")
P.OM1ssb <- ggplot(OM1_ssb, aes(year, data, col=iter)) + geom_line(size=1.2) +
  theme_light(25) + ylab("SSB")+ ggtitle("OM1 SSB")
P.OM2ssb <- ggplot(OM2_ssb, aes(year, data, col=iter)) + geom_line(size=1.2) +
  theme_light(25) + ylab("SSB")+ ggtitle("OM2 SSB")

png("./Input/Piera/OMs_iterExample_%3d.png", width=2500, height=1800, res=250)
P.OM1rec
P.OM2rec
P.OM1ssb
P.OM2ssb
dev.off()

# Create multiple FLSTOCKS to show stock developemnt in each OM
stks1 <- FLStocks(OM1_03=OM1, OM1_05=OM2)
stks2 <- FLStocks(OM2_03=OM3, OM2_05=OM4)
stks3 <- FLStocks(OM3_03=OM5, OM3_05=OM6)

png("./Input/Piera/OMs_StockDev_%3d.png", width=2000, height=1800, res=200)
plot(window(stks1,1,25)) + theme_bw(22) + 
  theme(legend.position = "top", legend.title = element_blank()) + xlab("Year") 
plot(window(stks2,1,25)) + theme_bw(22) + 
  theme(legend.position = "top", legend.title = element_blank()) + xlab("Year") 
plot(window(stks3,1,25)) + theme_bw(22) + 
  theme(legend.position = "top", legend.title = element_blank()) + xlab("Year") 
dev.off()

allOMs <- FLStocks(OM1=OM1, OM2=OM2,OM3=OM3, OM4=OM4,OM5=OM5, OM6=OM6)
png("./Input/Piera/allOMs_StockDev.png", width=2000, height=1800, res=200)
plot(window(allOMs,1,25)) + theme_bw(22) + 
  theme(legend.position = "top", legend.title = element_blank()) + xlab("Year") 
dev.off()


# ================================================
# CHECK CORRELATION BETWEEN SSB AND IDX.BIO 
# ================================================
IDX.OM1 <- as.data.frame(quantSums(window(index(OM_list[[1]]$observations$idx[[1]]),start=20, end=25) * stock.wt(window(OM_list[[1]]$stk, end=25))))
SSB.OM1 <- as.data.frame(ssb(window(OM1,start=20,end=25)))
corIdxSSB <- data.table(cbind(IDX.OM1, SSB.OM1$data))
names(corIdxSSB)[7:8] <- c("IDX", "SSB")
myCor <- corIdxSSB[,cor(IDX,SSB),by=iter]

names(myCor)[2] <- "cor"
# calculate summary statistics
myCor_sum <- summarySE(myCor, measurevar="cor")

# Plot with CI
# Use 95% confidence intervals
P.COR <- ggplot(myCor_sum, aes(x=1, y=cor)) + 
  geom_errorbar(aes(ymin=cor-ci, ymax=cor+ci), size=1.2,width=.2) +
  geom_point(size=3) + theme_light(25) + 
  theme(legend.position = "none")



# ================================================
# CHECK INTERNAL CONSISTENCY IN SURVEY
# ================================================
myIDX <- window(OM_list[[1]]$observations$idx[[1]], start=20,end=25)
plot(myIDX[,,,,,50], type="internal")

myIDX <- window(OM_list[[2]]$observations$idx[[1]], start=20,end=25)
plot(myIDX, type="internal")

myIDX <- window(OM_list[[3]]$observations$idx[[1]], start=20,end=25)
plot(myIDX[,,,,,1], type="internal")

myIDX <- window(OM_list[[4]]$observations$idx[[1]], start=20,end=25)
plot(myIDX, type="internal")

myIDX <- window(OM_list[[5]]$observations$idx[[1]], start=20,end=25)
plot(myIDX, type="internal")


# *************************************
# Plot other biological characteristics
# MATURITY
STK_mat <- data.table(getQuant(allOMs[c(1,3)], mat))
P.MAT <- ggplot(STK_mat[iter==10 & year==1,], aes(age, data, col=OM)) + geom_line(size=1.2) +
  theme_classic(25) + ylab("% mat") +
  scale_x_continuous(breaks=seq(0,7,1)) +
  scale_color_manual(labels = c("OM1-2 & OM5-6", "OM3-4"), values = c("salmon", "deepskyblue4")) +
  theme(legend.position = c(0.7,0.3), legend.title = element_blank())
# NATURAL MORTALITY
STK_m <- data.table(getQuant(allOMs[c(1,3,5)], m))
P.M <- ggplot(STK_m[iter==1 & year==1,], aes(age, data, col=OM)) + geom_line(size=1.2) +
  theme_classic(25) + ylab("M") + scale_y_continuous(breaks=seq(0.5,1.8,0.2)) +
  scale_x_continuous(breaks=seq(0,7,1)) +
  scale_color_manual(labels = c("OM1-2", "OM3-4", "OM5-6"), values=c("salmon", "deepskyblue4", "seagreen3")) +
  theme(legend.position = c(0.7,0.7), legend.title = element_blank())
# WEIGHT AT AGE
STK_wt <- data.table(getQuant(allOMs[c(1,3,5)], catch.wt))
P.WT <- ggplot(STK_wt[iter==1 & year==1,], aes(age, data, col=OM)) + geom_line(size=1.2) +
  theme_classic(25) + ylab("kg") + scale_y_continuous(breaks=seq(0.004,0.032,0.003)) +
  scale_x_continuous(breaks=seq(0,7,1)) +
  scale_color_manual(labels = c("OM1-2", "OM3-4", "OM5-6"), values=c("salmon", "deepskyblue4", "seagreen3")) +
  theme(legend.position = c(0.7,0.2), legend.title = element_blank())
# SELECTIVITY AT AGE
STK_f <- data.table(getQuant(allOMs[c(1,3,5)], harvest))
STK_sel <- STK_f[iter==1 & year==1,]
refAge <- STK_f[iter==1 & year==1 & age==2,data]
refAge <- c(rep(refAge[1], 6), rep(refAge[2],8),rep(refAge[3],6))
STK_sel[,sel:=data/refAge]#,rep(0.04036131,6)]
P.SEL <- ggplot(STK_sel, aes(age, sel, col=OM)) + geom_line(size=1.2) +
  theme_classic(25) + ylab("sel") + scale_y_continuous(breaks=seq(0.2,1,0.1)) +
  scale_x_continuous(breaks=seq(0,7,1)) +
  scale_color_manual(labels = c("OM1-2", "OM3-4", "OM5-6"), values=c("salmon", "deepskyblue4", "seagreen3")) +
  theme(legend.position = c(0.7,0.2), legend.title = element_blank())

# Print all biological assumptions
png("./Input/Piera/OMs_BIOL_%3d.png", width=2000, height=1500, res=300)
P.MAT
P.M
P.WT
P.SEL
dev.off()


