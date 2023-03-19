
# Prepare data for the HMP IDB

library(tidyverse)
library(data.table)
library(lmerTest)
library(refund)
library(fdapace)


setwd("D:/PC backup 02142023/Longitudinal/FDA/IBD")
load("ibddata.Rdata")


# get the data with more than 5 samples per subject
subID5 <- names(which(table(selmeta_rmna$subject_id) >= 5))

meta1 <- selmeta_rmna[selmeta_rmna$subject_id %in% subID5, ]
taxa1 <- taxadatafinal[selmeta_rmna$subject_id %in% subID5, ]
taxa1 <- taxa1[, names(colMeans(taxa1) %>% sort(decreasing = T))]


dim(selmeta_rmna)

length(subID5)

# --------------- merge the samples from the same sample window -----------------------#
alltime <- meta1$days_from_first_collection
nwindow <- 50

allwindow <- cut(alltime, breaks = nwindow) 
uwindow   <- unique(allwindow);length(uwindow)
newtime   <- round(seq(1, max(alltime), length = length(uwindow)))



meta12 <- c()
taxa12 <- c()
for(i in 1:length(subID5)){
  for(j in 1:length(uwindow)){
    
    
    tempdmeta <- meta1[meta1$subject_id == subID5[i] & allwindow == uwindow[j], ,drop = F]
    tempdtaxa <- taxa1[meta1$subject_id == subID5[i] & allwindow == uwindow[j], ,drop = F]
    
    if(nrow(tempdmeta) == 1){
      tempdmeta$window <- newtime[j]
      meta12 <- rbind(meta12, tempdmeta[, c(-5), drop = F])
      taxa12 <- rbind(taxa12, tempdtaxa)
    }else if(nrow(tempdmeta) > 1){
      tempdmeta1 <- tempdmeta[1,-5, drop = F]
      tempdmeta1$window <- newtime[j]
      tempdmeta1$fecalcal <- mean(tempdmeta$fecalcal, na.rm = T)
      meta12 <- rbind(meta12, tempdmeta1)
      taxa12 <- rbind(taxa12, colMeans(tempdtaxa))
      
    }
  }
}

rownames(taxa12) <- meta12$sample_id

taxa12 <- taxa1
meta12 <- meta1
# filter the taxa based on the prevalence
dim(taxa12)
pre <- colSums(taxa12 != 0)/nrow(taxa12)

taxa12sub <- taxa12[, pre > 0.05]
dim(taxa12sub)

taxa12sub <- taxa12sub/rowSums(taxa12sub)
rowSums(taxa12sub)

# --------------- merge the samples from the same sample window -----------------------#

getclr <- function(x, pseudocount = 1e-6, ifall= TRUE, method = c("pseduo"), denominator = NULL){
  
  if(method == "pseduo"){
    if(ifall){x <- x + pseudocount}else{x[x == 0] <- pseudocount}}
  if(method == "zcomp_czm"){
    x <- zCompositions::cmultRepl(x, output = "p-counts", method = "CZM")
  }
  if(is.null(denominator)){for(i in 1:nrow(x)){x[i,] <- log(x[i,])-(1/ncol(x))*sum(log(x[i,]))}}else{
    x   <- x/rowSums(x)
    ind <- rank(-1*colSums(x)) %in% 1:denominator
    for(i in 1:nrow(x)){x[i,] <- log(x[i,])-(1/denominator)*sum(log(x[i,ind]))}
  }
  x
}
taxa1clr <- getclr(taxa12)



# --------------- Fit LMM to select the taxa --------------------  #

# use either proportional data or CLR data
modeldata <- cbind(meta12, taxa12sub)

tname     <- colnames(taxa12sub)
myformula <- paste0("fecalcal ~ 1  + disease + age + ", tname," + (1 | subject_id)")



allp <- c()
for(i in 1:length(tname)){
  
  fit.m   <- lmerTest::lmer(as.formula(myformula[i]), data = modeldata)
  fit.res <- summary(fit.m)
  allp[i] <- fit.res$coefficients[4,5]
  
}
allpadj <- p.adjust(allp, method = "BH")

signame <- tname[which(allpadj < 0.06)];signame
pre[signame]
allpadj[which(allpadj < 0.06)]


mygenus <- signame[2]
mygenus

modeldata$window <- modeldata$days_from_first_collection


ffr_res <- list()
st_res <- list()
pfr_res <- list()

signamesel <- signame[c(1)]

for(i in 1:length(signamesel)){
  
  
  otum <- modeldata[,c("subject_id","window", signamesel[i])]
  otum <- otum[!duplicated(otum[,c(1,2)]),]
  otum <- otum %>% spread(key = "window", value = signamesel[i]) 
  subID <-  otum$subject_id
  rownames(otum) <- subID
  otum <- otum[, -1]
  
  sum(is.na(otum))/c(nrow(otum)*ncol(otum))
  
  nn <- nrow(otum)
  tt <- ncol(otum)
  ttime <- colnames(otum) %>% as.numeric()
  
  
  FPAC_input <- MakeFPCAInputs(IDs = rep(1:nn, each=tt), tVec=rep(ttime, nn), t(otum))
  FPCAres    <- FPCA(FPAC_input$Ly, FPAC_input$Lt, list(FVEthreshold = 0.95, nRegGrid = tt))
  otum_fpca  <- predict(FPCAres, newLy = FPAC_input$Ly, newLt = FPAC_input$Lt, 
                        K = FPCAres$selectK, xiMethod = "CE")$predCurves
  
  colnames(otum_fpca) <- colnames(otum)
  
  
  # get the outcome matrix
  yym <- modeldata[,c("subject_id","window","fecalcal")]
  yym <- yym[!duplicated(yym[,c(1,2)]),]
  yym <- yym %>% spread(key = "window", value = "fecalcal") #%>% dplyr::select(-subject_id)
  rownames(yym) <- yym$subject_id
  yym <- as.matrix(yym[,-1])
  
  
  yymcov <- modeldata[,c("subject_id","window","fecalcal", "disease", "age")]
  yymcov <- yymcov[!duplicated(yymcov[,c(1,2)]),]
  yymcov <- yymcov %>% spread(key = "window", value = "fecalcal") #%>% dplyr::select(-subject_id)
  rownames(yymcov) <- yymcov$subject_id
  mydisease <- yymcov$disease
  myage <- yymcov$age 
  myage <- scale(myage) %>% as.vector()
  
  # imput yy
  FPAC_input_yy <- MakeFPCAInputs(IDs = rep(1:nn, each=tt), tVec=rep(ttime, nn), t(yym))
  FPCAres_yy    <- FPCA(FPAC_input_yy$Ly, FPAC_input_yy$Lt, list(FVEthreshold = 0.95, nRegGrid = tt, dataType = "Sparse"))
  yym_imput     <- predict(FPCAres_yy, newLy = FPAC_input_yy$Ly, newLt = FPAC_input_yy$Lt, 
                           K = FPCAres_yy$selectK, xiMethod = "CE")$predCurves
  
  
  colnames(yym_imput) <- colnames(yym)
  
  otum_fpca <- scale(otum_fpca)
  
  fit.ff <- pffr(yym_imput ~ ff(otum_fpca, xind=round(FPCAres$workGrid), limits = "s<=t") + c(myage) + c(mydisease), yind=round(FPCAres_yy$workGrid))
  
  
  tgrid <- round(155)
  coef.ff      <- coef(fit.ff, n2 = tgrid)
  pffr_res     <- coef.ff[[2]][[2]]$coef
  
  pffr_res$lc     <- pffr_res$value - 1.96*pffr_res$se
  pffr_res$uc     <- pffr_res$value + 1.96*pffr_res$se
  pffr_res$pvalue <- (1 - pnorm(abs(pffr_res$value/pffr_res$se)))*2
  
  
  pvalue     <- matrix(pffr_res$pvalue, nrow = tgrid, dimnames = list(unique(pffr_res$otum_imput.smat), unique(pffr_res$otum_imput.smat)))
  coefficent <- matrix(pffr_res$value, nrow = tgrid, dimnames = list(unique(pffr_res$otum_imput.smat), unique(pffr_res$otum_imput.smat)))
  
  ffr_res[[i]] <- list(pvalue = t(pvalue), coef = t(coefficent))
  
  s <- 1:dim(yym_imput)[2]
  t <- s
  ##############################
  
  ST_pvalue  <- matrix(NA, nrow = length(s), ncol = length(t))
  ST_esimate <- matrix(NA, nrow = length(s), ncol = length(t))

  for(iii in 1:length(s)){
    for(jjj in 1:length(s)){

      mycoef          <- summary(glm(yym_imput[,jjj]~otum_fpca[,iii] + myage + mydisease))$coef
      ST_pvalue[iii,jjj]  <- mycoef[2,4]
      ST_esimate[iii,jjj] <- mycoef[2,1]

    }
  }

  st_res[[i]] <- list(pvalue = ST_pvalue, coef = ST_esimate)
  
  ##########################################
  
  pfr.res <- c()
  for(j in 1:length(s)){
    
    yy  <- yym_imput[,j]
    fit.lf     <- refund::pfr(yy ~ lf(otum_fpca, k=30, bs="ps", argvals = as.numeric(colnames(otum_fpca))) + myage + mydisease)
    pfr.pffr_res      <- coef(fit.lf, n = tt)
    
    pfr.pffr_res$lc     <- pfr.pffr_res$value - 1.96*pfr.pffr_res$se
    pfr.pffr_res$uc     <- pfr.pffr_res$value + 1.96*pfr.pffr_res$se
    pfr.pffr_res$pvalue <- (1 - pnorm(abs(pfr.pffr_res$value/pfr.pffr_res$se)))*2
    
    pfr.res <- rbind(pfr.res, pfr.pffr_res)
    
  }
  
  pfrcoeff <- matrix(pfr.res$value, nrow = 155)
  pfrpvalue <- matrix(pfr.res$pvalue, nrow = 155)
  
  
  pfr_res[[i]] <- list(pvalue = t(pfrpvalue), coef = t(pfrcoeff))
  
}


pffr_resplot <- ffr_res
st_resplot   <- st_res
pfr_resplot  <- pfr_res


for(i in 1:1){
  pffr_resplot[[i]]$coef[upper.tri(pffr_resplot[[i]]$coef, diag = F)] <- NA
  pffr_resplot[[i]]$pvalue[upper.tri(pffr_resplot[[i]]$pvalue, diag = F)] <- NA
  pffr_resplot[[i]]$pvalue[pffr_resplot[[i]]$pvalue > 0.05] <- NA
  pffr_resplot[[i]]$coef[is.na( pffr_resplot[[i]]$pvalue)] <- NA
  
  
  st_resplot[[i]]$coef[upper.tri(st_resplot[[i]]$coef, diag = F)] <- NA
  st_resplot[[i]]$pvalue[upper.tri(st_resplot[[i]]$pvalue, diag = F)] <- NA
  st_resplot[[i]]$pvalue[st_resplot[[i]]$pvalue > 0.05] <- NA
  st_resplot[[i]]$coef[is.na( st_resplot[[i]]$pvalue)] <- NA
  
  # 
  pfr_resplot[[i]]$coef[upper.tri(pfr_resplot[[i]]$coef, diag = F)] <- NA
  pfr_resplot[[i]]$pvalue[upper.tri(pfr_resplot[[i]]$pvalue, diag = F)] <- NA
  pfr_resplot[[i]]$pvalue[pfr_resplot[[i]]$pvalue > 0.05] <- NA
  pfr_resplot[[i]]$coef[is.na( pfr_resplot[[i]]$pvalue)] <- NA
  
}


library(GA)

x1 <- round(unique(pffr_res$otum_fpca.smat))
x2 <- x1



i=1
png(filename = "D:/FDA/realres_pvalue2.png", width = 7, height = 4, res = 1000, units = "in")
par(mfrow = c(1, 2))
persp3D(x1, x2, pffr_resplot[[i]]$coef, xlab = "", ylab = "", zlab="", 
        theta = -50, phi = 30, expand = 0.7, shade = 0.1, d = 2)

persp3D(x1, x2, pffr_resplot[[i]]$pvalue, xlab = "", ylab = "", zlab="", 
        theta = -50, phi = 30, expand = 0.7, shade = 0.1, d = 2)

persp3D(x1, x2, st_resplot[[i]]$pvalue, xlab = "", ylab = "", zlab="", 
        theta = -50, phi = 30, expand = 0.7, shade = 0.1, d = 2)


dev.off()


persp3D(x1, x2, coefficent, xlab = "Outcome", ylab = "OTU", zlab="", 
        theta = -70, phi = 30, expand = 0.7, shade = 0.1, d = 2)

png(filename = "g_Blautia_coef.png", width = 8, height = 7, res = 1000, units = "in")
persp3D(x1, x2, coefficent, xlab = "Outcome", ylab = "OTU", zlab="", 
        theta = -70, phi = 30, expand = 0.7, shade = 0.1, d = 2)
dev.off()




















p2 <- ggplot(modeldata, aes(x = window, y = g__Blautia, group = subject_id)) + 
  geom_line() + 
  geom_point(shape=18) +
  xlab("Days") + ylab(mygenus) + 
  theme_bw()
p2

ggsave(p2, filename = "g_Blautia.png", width = 5, height = 4)
# get the otu matrix for selected geneus



otum <- modeldata[,c("subject_id","window", mygenus)]
otum <- otum[!duplicated(otum[,c(1,2)]),]
otum <- otum %>% spread(key = "window", value = "g__Blautia") 
subID <-  otum$subject_id
rownames(otum) <- subID
otum <- otum[, -1]

sum(is.na(otum))/c(nrow(otum)*ncol(otum))

nn <- nrow(otum)
tt <- ncol(otum)
ttime <- colnames(otum) %>% as.numeric()

# check the correlation
cor(otum, use = "pairwise.complete.obs") %>% as.vector() %>% summary()

# check missing proportion
sum(is.na(otum))/(nn*tt)


#source("C:/Users/lymna/OneDrive/Desktop/Longitudinal/FDA/sparse_fun.R")

# use two method to impute missing
#otum_pfr <- imputx(as.matrix(otum), tlength = tt, by = 1, N.smooth = round(tt/2), N.smooth.single = round(tt/2))

FPAC_input <- MakeFPCAInputs(IDs = rep(1:nn, each=tt), tVec=rep(ttime, nn), t(otum))
FPCAres    <- FPCA(FPAC_input$Ly, FPAC_input$Lt, list(FVEthreshold = 0.95, nRegGrid = tt))
otum_fpca  <- predict(FPCAres, newLy = FPAC_input$Ly, newLt = FPAC_input$Lt, 
                      K = FPCAres$selectK, xiMethod = "CE")$predCurves



#colnames(otum_pfr) <- colnames(otum)
colnames(otum_fpca) <- colnames(otum)


plototum <- otum_fpca %>% as.data.frame() %>% 
  mutate(subject_id = unique(modeldata$subject_id)) %>% 
  gather(key = "time", value = "taxa", 1:tt) %>% 
  mutate(time = as.numeric(time))

plototum <- plototum %>% filter(!is.na(taxa))
p3 <- ggplot(plototum, aes(x = time, y = taxa, group = subject_id)) + 
  geom_line() + 
  geom_point(shape=18) +
  xlab("Days") + ylab(mygenus) + 
  theme_bw()
p3
ggsave(p3, filename = "g_Blautia_PFR.png", width = 5, height = 4)


# looks FPCA is better


subIDm <- cbind(subID, 1:length(subID))
rownames(otum_fpca) <- 1:length(subID)


# --------------- fit pffr --------------------- #

# get the outcome matrix
yym <- modeldata[,c("subject_id","window","fecalcal")]
yym <- yym[!duplicated(yym[,c(1,2)]),]
yym <- yym %>% spread(key = "window", value = "fecalcal") #%>% dplyr::select(-subject_id)
rownames(yym) <- yym$subject_id
yym <- as.matrix(yym[,-1])


yymcov <- modeldata[,c("subject_id","window","fecalcal", "disease", "age")]
yymcov <- yymcov[!duplicated(yymcov[,c(1,2)]),]
yymcov <- yymcov %>% spread(key = "window", value = "fecalcal") #%>% dplyr::select(-subject_id)
rownames(yymcov) <- yymcov$subject_id
mydisease <- yymcov$disease
myage <- yymcov$age 


# imput yy
FPAC_input_yy <- MakeFPCAInputs(IDs = rep(1:nn, each=tt), tVec=rep(ttime, nn), t(yym))
FPCAres_yy    <- FPCA(FPAC_input_yy$Ly, FPAC_input_yy$Lt, list(FVEthreshold = 0.95, nRegGrid = tt))
yym_imput     <- predict(FPCAres_yy, newLy = FPAC_input_yy$Ly, newLt = FPAC_input_yy$Lt, 
                         K = FPCAres_yy$selectK, xiMethod = "CE")$predCurves


colnames(yym_imput) <- colnames(yym)



subIDm_IBD <- subIDm[mydisease != "IBD", ]
yym_imput_IBD <- yym_imput[mydisease != "IBD", ]
otum_fpca_IBD <- otum_fpca[mydisease != "IBD", ]
myage_IBD   <- myage[mydisease != "IBD"]


subIDm_IBD <- subIDm[mydisease == "IBD",]
yym_imput_IBD <- yym_imput[mydisease == "IBD", ]
otum_fpca_IBD <- otum_fpca[mydisease == "IBD", ]
myage_IBD   <- myage[mydisease == "IBD"]


subIDm_IBD[,2] <- 1:nrow(subIDm_IBD)

plotyym <- yym_imput_IBD %>% as.data.frame() %>%
  mutate(subject_id = subIDm_IBD[,1]) %>% 
  gather(key = "time", value = "outcome", 1:ncol(yym)) %>% 
  mutate(time = as.numeric(time))

plotyym <- plotyym %>% filter(!is.na(outcome))
ydata   <- plotyym
colnames(ydata) <- c(".obs",   ".index", ".value")

for(i in 1:nrow(ydata)){
  ydata$.obs[i] <- subIDm_IBD[which(subIDm_IBD[,1] == ydata$.obs[i]), 2]
}

ydata$.obs <- as.numeric(ydata$.obs)

p4 <- ggplot(plotyym, aes(x = time, y = outcome, group = subject_id)) + 
  geom_line() + 
  geom_point(shape=18) +
  xlab("Days") + ylab("fecalcal") + 
  theme_bw()
p4

ggsave(p4, filename = "fecalcal.png", width = 5, height = 4)




all(ydata$.obs %in% rownames(otum_fpca))




yym_imput_IBD <- yym_imput[mydisease != "IBD", ]
otum_fpca_IBD <- otum_fpca[mydisease != "IBD", ]
myage_IBD   <- myage[mydisease != "IBD"]


yym_imput_IBD <- yym_imput[mydisease == "IBD", ]
otum_fpca_IBD <- otum_fpca[mydisease == "IBD", ]
myage_IBD   <- myage[mydisease == "IBD"]


# fit pffr
#fit.ff <- pffr(yym_imput ~ ff(otum_fpca, xind=round(FPCAres$workGrid), limits = "s<=t") + c(mydisease) + c(myage), yind=round(FPCAres_yy$workGrid))

fit.ff <- pffr(yym_imput_IBD ~ ff(otum_fpca_IBD, xind=round(FPCAres$workGrid), limits = "s<=t") + c(myage_IBD), yind=round(FPCAres_yy$workGrid))



#rownames(otum_fpca_IBD) <- 1:nrow(otum_fpca_IBD)
#fit.ff <- pffr(yym_imput_IBD ~ ff(otum_fpca_IBD, xind=round(FPCAres$workGrid) , limits = "s<=t") + c(myage_IBD), data = as.data.frame(otum_fpca_IBD), ydata = ydata, yind=ttime)
#plot(fit.ff, pers=TRUE, pages=1)


tgrid <- round(355/2)
coef.ff      <- coef(fit.ff, n2 = tgrid)
pffr_res     <- coef.ff[[2]][[2]]$coef

pffr_res$lc     <- pffr_res$value - 1.96*pffr_res$se
pffr_res$uc     <- pffr_res$value + 1.96*pffr_res$se
pffr_res$pvalue <- (1 - pnorm(abs(pffr_res$value/pffr_res$se)))*2


pvalue     <- matrix(pffr_res$pvalue, nrow = tgrid, dimnames = list(unique(pffr_res$otum_imput.smat), unique(pffr_res$otum_imput.smat)))
coefficent <- matrix(pffr_res$value, nrow = tgrid, dimnames = list(unique(pffr_res$otum_imput.smat), unique(pffr_res$otum_imput.smat)))

#plotly::plot_ly(z=~coefficent, type = "surface")
#plotly::plot_ly(z=~pvalue, type = "surface")

range(coefficent)

coefficent[upper.tri(coefficent, diag = F)] <- NA
pvalue[upper.tri(pvalue, diag = F)] <- NA

library(GA)

x1 <- round(unique(pffr_res$otum_fpca.smat))
x2 <- x1

png(filename = "g_Blautia_coef.png", width = 8, height = 7, res = 1000, units = "in")
persp3D(x1, x2, coefficent, xlab = "Outcome", ylab = "OTU", zlab="", 
        theta = -70, phi = 30, expand = 0.7, shade = 0.1, d = 2)
dev.off()

png(filename = "g_Blautia_pvalue.png", width = 8, height = 7, res = 1000, units = "in")
persp3D(x1, x2, pvalue, xlab = "Outcome", ylab = "OTU", zlab="", 
        theta = -70, phi = 30, expand = 0.7, shade = 0.1, d = 2)
dev.off()

sum(pvalue < 0.05)
# use mean of y to do one time point analysis
yymean <- apply(yym, 1, mean, na.rm = T)
STres <- apply(otum_imput, 2, function(x){
  
  mycoef          <- summary(lm(yymean~x))$coef
  c(mycoef[2,4],mycoef[2,1])

})

ST_p <- STres[1,]
ST_coef <- STres[2,]

plot(names(ST_p), ST_p, type = "b", pch = 19)
plot(names(ST_coef), ST_coef, type = "b", pch = 19)

plot(colnames(otum))

#one time point analysis use yym_imput

workleng <- length(FPCAres_yy$workGrid)
ST_pvalue  <- matrix(NA, nrow = workleng, ncol = workleng)
ST_esimate <- matrix(NA, nrow = workleng, ncol = workleng)
ST_se      <- matrix(NA, nrow = workleng, ncol = workleng)

for(i in 1:workleng){
  for(j in 1:workleng){
    
    mycoef          <- summary(lm(yym_imput[,j]~otum_imput[,i]))$coef
    ST_pvalue[i,j]  <- mycoef[2,4]
    ST_esimate[i,j] <- mycoef[2,1]
    ST_se[i,j]      <- mycoef[2,2]
    
  }
}

plotly::plot_ly(z=~ST_pvalue,type = "surface")
plotly::plot_ly(z=~ST_esimate,type = "surface")


range(ST_esimate)


