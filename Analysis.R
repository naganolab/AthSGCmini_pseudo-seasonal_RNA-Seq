# This is the R script used for RNA-Seq analysis in the following paper.
# "Artificial mimicry of seasonal transcriptome dynamics in Arabidopsis thaliana revealed short- and long-term response to the environmental condition."
# Kurita et al. 

oldpar <- par(no.readonly = TRUE)

# read input data ------------------------------------------------------------------
## Ahg
desAhg <- read.csv("input/GeneDescription_desAhg.csv", row.names = 1, stringsAsFactors = F)
atAhg5.5 <- read.csv("input/SampleAttributes_Ahg5.5.csv", row.names = 1, stringsAsFactors = F)
load("input/cpmAhg5.5")

## Ath
desAth <- read.csv("input/GeneDescription_desAth.csv", row.names = 1, stringsAsFactors = F)
desAth <- desAth[desAth$NormalizationGroup == "data",]
atAth5.5 <- read.csv("input/SampleAttributes_Ath5.5_all.csv", row.names = 1, stringsAsFactors = F)
load("input/cpmAth5.5")

## reciprocalBLASTresults
reci <- read.csv("input/reciprocalBLASTresults_AhgAth.csv", row.names = 1, stringsAsFactors = F)


# caluculate log2(mean cpm + 1) and draw hist, take ex gene (log2(mean cpm + 1) > 2)---------------------
## Ahg,  Fig. S2B
tmp <- apply(cpmAhg5.5, 1, mean)
tmp2 <- log2(tmp+1)
exgene_Ahg <- names(tmp2[tmp2 > 2]) 
length(exgene_Ahg)   # expressed gene in Ahg, 16158 genes

dir.create("plot_Fig")
tiff("plot_Fig/expressed_gene_hist.tif", pointsize = 13, width = 300, height = 300)
par(mar = c(4, 5, 3, 1))
hist(tmp2, breaks = seq(0, 16, 0.5), main = sprintf("Histogram of mean (Ahg)"), 
     xlab = "log2(rpm mean +1)", ylab = "", ylim = c(0,13000), las = 1)
mtext(side=2, text = "Genes", line = 3.5)
ex <- sum(tmp2 > 2)
mtext(sprintf("Expressed gene: %s", ex), side = 3, line = -0.8, at = 8.2)
abline(v = 2, col = "violetred1", lwd = 2)
dev.off()  

par(oldpar)


## Ath 
days <- c(1,3,7) # number of culture days
cols <- c("firebrick1", "#00FFFF50", "#1E90FF50", "#0000FF50")

for(i in 1:3){  
  d <- days[i]
  tmpcpm <- cpmAth5.5[,(atAth5.5$day == d)&(atAth5.5$set == "seasonal")]  
  tmp <- apply(tmpcpm, 1, mean)
  tmp2 <- log2(tmp+1)
  hist(tmp2, breaks = seq(0, 16, 0.5), main = sprintf("Histogram of mean (d%s culture)" , d),
       col = cols[(i+1)] , xlab = "log2(rpm mean +1)",ylim = c(0,15000))
  abline(v = 2, col = cols[1], lwd = 2)
  ex <- sum(tmp2 > 2)
  mtext(sprintf("Expressed gene: %s", ex), side = 3, line = -0.8, at = 8.2)
  
  assign(sprintf("exgene_d%d", days[i]), names(tmp2[tmp2 > 2]))
}

exgene_list <- list(Ahg = exgene_Ahg, Ath_d1 = exgene_d1, Ath_d3 = exgene_d3, Ath_d7 = exgene_d7)
rm(exgene_Ahg, exgene_d1, exgene_d3, exgene_d7)
summary(exgene_list)  # expressed gene in Ahg and Ath


# cos curve fitting -----------------------------------------------------------------------------
## set cos curve function (return: list(sinfit, predmx))
library(lubridate)

fn.cosfit1 <- function(xtime, tmpl2cpm){
  # xtime treat
  xtime <- as.POSIXct(xtime, format="%Y-%m-%d")
  xtime.yday <- yday(xtime)                       # convert day into number
  xtime2.yday <- c(xtime.yday, xtime.yday + 365, xtime.yday + 365*2)
  
  # tripled l2cpm
  tmpl2cpm3 <- cbind(tmpl2cpm,tmpl2cpm,tmpl2cpm)
  
  # make result matrix
  tmpresult <- matrix(NA, nrow=nrow(tmpl2cpm3), ncol=3)
  colnames(tmpresult) <- c("alpha", "phi", "modi.phi")
  
  # make prediction result matrix
  predmx <- matrix(NA, nrow = nrow(tmpl2cpm3), ncol = 365)
  
  # fit sin curve by nls  
  for(i in 1:nrow(tmpl2cpm3)){
    y <- tmpl2cpm3[i, ]
    C <- mean(y)
    fit.nls <- nls(y~C+(1/2)*alpha*cos(2*pi*(xtime2.yday-phi)/365),
                   start = list(alpha = 0.01, phi = 0), control=nls.control(warnOnly=TRUE))
    
    p <- coef(fit.nls)    # take coefficients
    
    if (p["alpha"] < 0) {                 # treatment for alpha < 0 
      p["alpha"] <- p["alpha"]*(-1)
      p["phi"] <- p["phi"]+365/2
    }
    
    modi.phi <- p["phi"] %% 365     # set phi 0~365
    p2 <- c(p, modi.phi)
    names(p2) <- c("alpha", "phi", "modi.phi")
    
    tmpresult[i,] <- p2
    
    rg <- 1:365
    predmx[i,] <- predict(fit.nls, newdata = list(xtime2.yday=rg))
    
    if(i%%1000==0){cat(sprintf("%s: %s\n", Sys.time(), i))}
  }
  
  rownames(tmpresult) <- rownames(tmpl2cpm)
  tmpresult <- as.data.frame(tmpresult)
  rownames(predmx) <- rownames(tmpl2cpm)
  
  list(sinFit = tmpresult, predmx = predmx)
} 

## cos curve fitting by using above function
### Ahg 
tmpl2cpm <- log2(cpmAhg5.5+1)
xtime <- sprintf("%s-%s-%s", 2012, atAhg5.5$month, atAhg5.5$day)

dir.create("Robject")
Ahg_cosfit1 <- fn.cosfit1(xtime = xtime, tmpl2cpm = tmpl2cpm) # It will take a while.
save(Ahg_cosfit1, file = "Robject/Ahg_cosfit1")

### Ath
# choose seasonal samples
s.at <- atAth5.5[atAth5.5$set == "seasonal",]
nrow(s.at)          # [1] 189

# cos curve fitting for each day sample.
for(j in c(1,3,7)){
  # take only ex gene  
  tmpat <- s.at[s.at$day==j,]
  tmpcpm <- cpmAth5.5[,rownames(tmpat)]
  tmpl2cpm <- log2(tmpcpm+1)
  
  # set time
  xtime <- sprintf("%s-%s-%s", rep(2017,nrow(tmpat)), tmpat$num.month, rep(15,nrow(tmpat)))
  
  # fitting by using above function
  cosfit1 <- fn.cosfit1(xtime = xtime, tmpl2cpm = tmpl2cpm)
  assign(sprintf("Athd%d_cosfit1", j), cosfit1)
}

save(Athd1_cosfit1, file = "Robject/Athd1_cosfit1")
save(Athd3_cosfit1, file = "Robject/Athd3_cosfit1")
save(Athd7_cosfit1, file = "Robject/Athd7_cosfit1")

# Take the exgene of Ahg, reciprocal orthologue with Ath (exbase) -------------------------
sum(reci$reciprocal_level == 2)  # 26789
tmp <- reci[exgene_list$Ahg,]
exbase <-  tmp[tmp$reciprocal_level == 2,]
nrow(exbase)         # 14587
nrow(exbase)/nrow(tmp)*100  # 90.27% 


# Histgram of alpha and phi of exbase gene, Fig. 3BC, S5AB, S6AB--------------------------------------------------------
Fits <- list(Ahg = Ahg_cosfit1$sinFit[exbase$Ahg_ID,], 
             Ath_d1 = Athd1_cosfit1$sinFit[exbase$Ath_ID,], 
             Ath_d3 = Athd3_cosfit1$sinFit[exbase$Ath_ID,],
             Ath_d7 = Athd7_cosfit1$sinFit[exbase$Ath_ID,])
tx <- c("Ahg", "Ath d1", "Ath d3", "Ath d7")
sp <- c("Ahg", "Ath", "Ath", "Ath")
d <- c("", "d1", "d3", "d7")
cols <- c("#FFCC99", "#00FFFF50", "#1E90FF50", "#0000FF50")

for(i in 1:4){
  tiff(sprintf("plot_Fig/Histgrams_%s.tif", tx[i]), pointsize = 23, width = 930, height = 480)
  par(lwd = 2, mfcol = c(1,2))
  tmp <- parse(text = paste0("bold('Histogram of amplitude '( ","bolditalic('", sp[i], "')~'", d[i], "'))")) 
  hist(Fits[[i]]$alpha, breaks = seq(0,12,0.25), main =tmp, ylim = c(0,4000),
       xlab = expression('Amplitude ' (italic('\u03B1'))), ylab = "Genes", cex.main = 1.1, lwd = 2, las = 1)
  hist(Fits[[i]]$alpha[Fits[[i]]$alpha > 1], breaks = seq(0,12,0.25), col = cols[i], add = TRUE)
  abline(v = 1, lwd = 2, lty = 3)
  n <- sum(Fits[[i]]$alpha > 1)
  tmp <-  parse(text = paste0("italic('\u03B1') > '1, '~", n, " ~'/' ~", nrow(Fits[[i]]), " ~gene")) 
  mtext(tmp, side = 3)
  
  tmp <- parse(text = paste0("bold('Histogram of phase ' (","bolditalic('", sp[i], "')~'", d[i], "'))")) 
  hist(Fits[[i]]$modi.phi, breaks = seq(0,365,14.6), ylim = c(0,3000), xlim = c(0,400), cex.main = 1.1,
       xlab = "Days from Jan 1", ylab = "Genes", main = tmp, lwd = 2, las =1)
  hist(Fits[[i]]$modi.phi[Fits[[i]]$alpha > 1], xaxt= "n", breaks = seq(0,365,14.6), col = cols[i], add = TRUE)
  tmp <-  parse(text = paste0("italic('\u03B1') > '1, '~", n, " / ", nrow(Fits[[i]]), " ~gene")) 
  mtext(tmp, side = 3)
  dev.off()
}
par(oldpar)


# take alpha and phi of exbase genes -------------------------------------------
## alpha 
alpha.all <- cbind(Fits$Ahg[exbase$Ahg_ID,]$alpha, Fits$Ath_d1[exbase$Ath_ID,]$alpha,
                   Fits$Ath_d3[exbase$Ath_ID,]$alpha, Fits$Ath_d7[exbase$Ath_ID,]$alpha)
colnames(alpha.all) <- c("Ahg", "Ath_d1", "Ath_d3", "Ath_d7")
alpha.all <- as.data.frame(alpha.all)
rownames(alpha.all) <- rownames(exbase)

## phi 
phi.all <- cbind(Fits$Ahg[exbase$Ahg_ID,]$modi.phi, Fits$Ath_d1[exbase$Ath_ID,]$modi.phi,
                 Fits$Ath_d3[exbase$Ath_ID,]$modi.phi, Fits$Ath_d7[exbase$Ath_ID,]$modi.phi)
colnames(phi.all) <- c("Ahg", "Ath_d1", "Ath_d3", "Ath_d7")
phi.all <- as.data.frame(phi.all)
rownames(phi.all) <- rownames(exbase)


# density plot of phi 1 vs 3 vs 7day, Fig. 6B -------------------------------------------
cols2 <- c("darkorange1", "cyan2", "dodgerblue", "blue3")
tmp <- hist(phi.all[alpha.all[,2]>1,2], breaks = seq(0,365,14.6), plot = F)
tmp2 <- c(16,17,15,18)

tiff("plot_Fig/phaseDistribution_plot.tif", pointsize = 18, width = 600, height = 400)
plot(tmp$breaks[-1]-7.3, tmp$counts, type="o",col = cols2[2], pch = tmp2[2],
     ylim = c(0,1600), lwd=2, main ="Distribution of phases", xlab = "Days from Jan 1", ylab = "Genes")

for (i in c(4,3,1)){
  tmp3 <- hist(phi.all[alpha.all[,i]>1,i], breaks = seq(0,365,14.6), plot = F)
  points(tmp3$breaks[-1]-7.3, tmp3$counts, col=cols2[i], pch = tmp2[i])
  lines(tmp3$breaks[-1]-7.3, tmp3$counts, col=cols2[i],lwd=2)
}
legend("topright", colnames(phi.all)[c(1,4,3,2)], text.col=cols2[c(1,4,3,2)], col=cols2[c(1,4,3,2)],
       pch = tmp2[c(1,4,3,2)], bty="n")
dev.off()


# define Ahg SOgene --------------------------------------------------------------------
alpha.SO <- alpha.all[alpha.all$Ahg > 1,]
nrow(alpha.SO)  #  4312 Ahg SO genes
SOgenes <- exbase[rownames(alpha.SO),]
phi.SO <- phi.all[alpha.all$Ahg > 1,]

# common with Ath_d7, Fig. 3D
sum(alpha.all$Ath_d7 > 1)  # 4016 Ath 7d SO genes
sum(alpha.SO$Ath_d7 > 1)   # 1757 common SO genes with Ahg
sum(alpha.all$Ath_d7 > 1)-sum(alpha.SO$Ath_d7 > 1)   # 2259, only Ath 
nrow(alpha.SO)-sum(alpha.SO$Ath_d7 > 1)   # 2555, only Ahg

# common with Ath_d3, Fig. S5E
sum(alpha.all$Ath_d3 > 1)  # 3975 Ath 3d SO genes
sum(alpha.SO$Ath_d3 > 1)   # 1766 common SO genes with Ahg
sum(alpha.all$Ath_d3 > 1)-sum(alpha.SO$Ath_d3 > 1)   # 2209, only Ath 
nrow(alpha.SO)-sum(alpha.SO$Ath_d3 > 1)   # 2546, only Ahg 

# common with Ath_d1, Fig. S6E
sum(alpha.all$Ath_d1 > 1)  # 4575, Ath 1d SO genes
sum(alpha.SO$Ath_d1 > 1)   # 1886, common SO genes with Ahg
sum(alpha.all$Ath_d1 > 1)-sum(alpha.SO$Ath_d1 > 1)   # 2689, only Ath 
nrow(alpha.SO)-sum(alpha.SO$Ath_d1 > 1)   # 2426, only Ahg 


# count phi in winter or summer ------------------------------------------------------------------
## Ahg
sum(phi.SO$Ahg < 59) + sum(phi.SO$Ahg > 335)   # 2099, Dec to Feb
(sum(phi.SO$Ahg < 59) + sum(phi.SO$Ahg > 335))/nrow(phi.SO)*100   # 48.67811 %

sum((phi.SO$Ahg > 152)&(phi.SO$Ahg < 243))     #  1445, Jun to Aug
sum((phi.SO$Ahg > 152)&(phi.SO$Ahg < 243))/nrow(phi.SO)*100   # 33.51113 %

## Ath_d7
tmp <- alpha.all$Ath_d7 > 1
sum(tmp)  # 4016
sum(phi.all$Ath_d7[tmp] < 59) + sum(phi.all$Ath_d7[tmp] > 335)   # 2125, Dec to Feb
(sum(phi.all$Ath_d7[tmp] < 59) + sum(phi.all$Ath_d7[tmp] > 335))/sum(tmp)*100  # 52.91335 %

sum((phi.all$Ath_d7[tmp] > 152)&(phi.all$Ath_d7[tmp] < 243))     #  1693, Jun to Aug
sum((phi.all$Ath_d7[tmp] > 152)&(phi.all$Ath_d7[tmp] < 243))/sum(tmp)*100   # 42.15637 %

## Ath_d3
tmp <- alpha.all$Ath_d3 > 1
sum(tmp)  # 3975
sum(phi.all$Ath_d3[tmp] < 59) + sum(phi.all$Ath_d3[tmp] > 335)   # 2403, Dec to Feb
(sum(phi.all$Ath_d3[tmp] < 59) + sum(phi.all$Ath_d3[tmp] > 335))/sum(tmp)*100  # 60.45283 %

sum((phi.all$Ath_d3[tmp] > 152)&(phi.all$Ath_d3[tmp] < 243))     #  1523, Jun to Aug
sum((phi.all$Ath_d3[tmp] > 152)&(phi.all$Ath_d3[tmp] < 243))/sum(tmp)*100   # 38.31447 %

## Ath_d1
tmp <- alpha.all$Ath_d1 > 1
sum(tmp)  # 4575
sum(phi.all$Ath_d1[tmp] < 59) + sum(phi.all$Ath_d1[tmp] > 335)   # 2476, Dec to Feb
(sum(phi.all$Ath_d1[tmp] < 59) + sum(phi.all$Ath_d1[tmp] > 335))/sum(tmp)*100  # 54.12022 %

sum((phi.all$Ath_d1[tmp] > 152)&(phi.all$Ath_d1[tmp] < 243))     #  2058, Jun to Aug
sum((phi.all$Ath_d1[tmp] > 152)&(phi.all$Ath_d1[tmp] < 243))/sum(tmp)*100   # 44.98361 %


# plot the number of SOgene 1 vs 3 vs 7 days, Fig. 6A -------------------------------------------
library(ggplot2)

tmp <- c(sum(alpha.all$Ath_d7 > 1), sum(alpha.all$Ath_d3 > 1), sum(alpha.all$Ath_d1 > 1),
         sum(alpha.SO$Ath_d7 > 1), sum(alpha.SO$Ath_d3 > 1), sum(alpha.SO$Ath_d1 > 1))
tmp2 <- rep(c("7d", "3d", "1d"),2)
tmp3 <- c(rep("SOgene",3),rep("common",3))
tmp <- as.data.frame(cbind(tmp, tmp2, tmp3))
colnames(tmp) <- c("Genes", "day", "group")
tmp$Genes <- as.numeric(as.character(tmp$Genes))
tmp$day <- factor(tmp$day, levels = c("7d", "3d", "1d"))
tmp$group <- factor(tmp$group, levels = c("SOgene", "common"))

tiff("plot_Fig/Histgrams_SOgene.tif", pointsize = 23, width = 350, height = 400)
g <- ggplot(tmp[1:3,], aes(x = day, y = Genes, fill = day)) +
  geom_bar(stat = "identity", width = 0.7, colour = "gray25", alpha = 0.4) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,6000)) +
  scale_fill_manual(values = c("royalblue1","skyblue1", "cadetblue1")) +
  theme_classic(base_line_size = 0.1, base_size = 18) +
  theme(panel.border = element_rect(fill = NA, size=1),
        axis.title.y = element_text(vjust = 1),
        axis.text=element_text(size = 16, color = "black")) +
  geom_bar(data = tmp[4:6,],aes(x = day, y = Genes, fill = day),
           stat = "identity", width = 0.7, colour = "gray25") 
plot(g)
dev.off()


# schatter plot of alpha alpha and phi, Fig. 3EF, S5CD, S6CD -----------------------------------------    
set <- t(combn(c(1,2,3,4),2))

## alpha
dir.create("plot_Fig/Scatterplots")
for(i in 1:nrow(set)){
  tiff(sprintf("plot_Fig/Scatterplots/Scatter_alpha_%s_%s.tif", 
               tx[set[i,1]],tx[set[i,2]]), pointsize = 19, width = 480, height = 480)
  par(pty = "s", lwd = 2)
  a <- set[i,1]
  b <- set[i,2]
  
  x1 <- alpha.all[,a][alpha.all[,a] > 1 | alpha.all[,b] > 1] 
  y1 <- alpha.all[,b][alpha.all[,a] > 1 | alpha.all[,b] > 1] 
  x2 <- alpha.all[,a][alpha.all[,a] > 1 & alpha.all[,b] > 1] # for colored points
  y2 <- alpha.all[,b][alpha.all[,a] > 1 & alpha.all[,b] > 1] # for colored points
  
  coltmp <- densCols(x2, y2, colramp = colorRampPalette(c("blue2","cyan","yellow", "white"))) 
  
  plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 10), xlab = tx[a], ylab = tx[b], cex.lab  = 1.1,
       main = "Scatterplot of amplitude", lwd = 2, las = 1)
  points(x1, y1, col = "#b0c4de", pch = 1, cex = 0.9, lwd= 1) 
  points(x2, y2, col = coltmp, pch = 1, cex = 0.9, lwd= 1)
  abline(0,1, col="maroon1", lwd= 1.5) 
  tmp <-  parse(text = paste0("italic('\u03B1') > '1, '~", length(x2), " / ", length(x1), " ~gene")) 
  mtext(tmp, line = 0.2)
  dev.off()
}

# phi
for(i in 1:nrow(set)){
  tiff(sprintf("plot_Fig/Scatterplots/Scatter_phi_%s_%s.tif", 
               tx[set[i,1]],tx[set[i,2]]), pointsize = 19, width = 480, height = 480)
  par(pty = "s", lwd = 2)
  a <- set[i,1]
  b <- set[i,2]
  
  x1 <- phi.all[,a][alpha.all[,a] > 1 | alpha.all[,b] > 1]
  y1 <- phi.all[,b][alpha.all[,a] > 1 | alpha.all[,b] > 1]
  x2 <- phi.all[,a][alpha.all[,a] > 1 & alpha.all[,b] > 1]
  y2 <- phi.all[,b][alpha.all[,a] > 1 & alpha.all[,b] > 1]
  
  coltmp <- densCols(x2, y2, colramp = colorRampPalette(c("blue3","blue","cyan","yellow", "white")))
  
  plot(0, 0, type = "n", xlim = c(1, 365), ylim = c(1, 365), xlab = tx[a], ylab = tx[b], cex.lab  = 1.1, 
       main = "Scatterplot of phase", lwd = 2, las = 1)
  points(x1, y1, col = "#b0c4de", pch = 1, cex = 0.9, lwd= 1)
  points(x2, y2, col = coltmp, pch = 1, cex = 0.9, lwd= 1)
  abline(0,1, col="maroon1", lwd= 1.5)
  tmp <-  parse(text = paste0("italic('\u03B1') > '1, '~", length(x2), " / ", length(x1), " ~gene")) 
  mtext(tmp, line = 0.2)
  dev.off()
}
par(oldpar)

# color scale for schatter plots
coltmp <- colorRampPalette(c("blue3","blue","cyan","yellow", "white"))
col <- coltmp(20)

plot(NULL, xlim=c(0,length(col)), ylim=c(0,10),
     axes=FALSE, xlab="", ylab="")
rect(0:(length(col)-1), 0, 1:length(col), 1, col=col,lty="blank")


# calculate the delta of alpha and phi between Ahg and Ath ----------------------------------
## delta alpha
delta.alpha <- matrix(NA, nrow=nrow(alpha.SO), ncol=3)
rownames(delta.alpha) <- rownames(alpha.SO)
colnames(delta.alpha) <- c("Ath_d1","Ath_d3","Ath_d7")

for(i in 2:4){
  tmp <- alpha.SO[,i]-alpha.SO$Ahg
  delta.alpha[,(i-1)] <- tmp
  
  hist(tmp, xlab="the delta alpha", xlim = c(-12,12), ylim = c(0,1500), breaks = seq(-12,12,0.5),
       main=sprintf("The delta alpha of %s vs Ahg", colnames(alpha.SO)[i]))
  abline(v=0, col = "deeppink", lwd = 2)
  
  tmp <- alpha.SO[alpha.SO[,i]>1,]
  tmp2 <- tmp[,i]-tmp$Ahg
  
  hist(tmp2, col = cols[i], breaks = seq(-8,8,0.5), add = TRUE)
  mtext(sprintf("alpha > 1,  %s / 4589 gene", nrow(tmp)), side = 3, line = 0.3)
}

delta.alpha <- as.data.frame(delta.alpha)

## delta phi 
delta.phi <- matrix(NA, nrow=nrow(phi.SO), ncol=3)
rownames(delta.phi) <- rownames(phi.SO)
colnames(delta.phi) <- c("Ath_d1","Ath_d3","Ath_d7")

for(i in 2:4){
  tmp <- phi.SO[,i]-phi.SO$Ahg
  for(j in 1:length(tmp)){
    if(tmp[j] < -365/2){tmp[j] <- 365+tmp[j]}else{tmp[j] <- tmp[j]}
    if(tmp[j] > 365/2){tmp[j] <- -(365-tmp[j])}else{tmp[j] <- tmp[j]}
  }
  delta.phi[,(i-1)] <- tmp
  
  hist(tmp, breaks = seq(-190,190,10), xlab="the delta phi", ylim = c(0,500),
       main=sprintf("The delta phi of %s vs Ahg", colnames(phi.SO)[i]))
  abline(v=0, col = "deeppink", lwd = 2)
  
  tmp2 <- tmp[alpha.SO[,i]>1]
  hist(tmp2, breaks = seq(-190,190,10), col = cols[i], add = TRUE)
  mtext(sprintf("alpha > 1,  %s / 4589 gene", length(tmp2)), side = 3, line = 0.3)
}

delta.phi <- as.data.frame(delta.phi)


# define mimicked genes in Ahg SOgene -----------------------------------------------------
mic.alpha <- abs(delta.alpha)<1
mic.phi <- abs(delta.phi) < 45
Ath.alpha.th <- alpha.SO[,2:4] > 1

mic.gene <- mic.alpha&mic.phi&Ath.alpha.th
mic.gene <- as.data.frame(mic.gene)
nrow(mic.gene) # 4312 mimicked genes were defined, T ot F

# mimicked gene list
tmp <- reci[rownames(alpha.SO),]
tmp$mic_d7 <- mic.gene$Ath_d7
tmp$mic_d3 <- mic.gene$Ath_d3
tmp$mic_d1 <- mic.gene$Ath_d1

tmp <- cbind(tmp, alpha.SO, phi.SO)
write.csv(tmp, file = "mimickedGeneList.csv")

## the number of mimicked genes, Fig. 3G, S5F, S6F ----------------------------------------------
sum(mic.gene$Ath_d7) # 946
sum(mic.gene$Ath_d3) # 1027
sum(mic.gene$Ath_d1) # 1018

sum(mic.gene$Ath_d1 & mic.gene$Ath_d3 & mic.gene$Ath_d7) # 412
sum(mic.gene$Ath_d1 | mic.gene$Ath_d3 | mic.gene$Ath_d7) # 1627

## delta alpha and alpha > 1 
sum(mic.alpha[,3]&Ath.alpha.th[,3]) # d7 1357
sum(mic.alpha[,2]&Ath.alpha.th[,2]) # d3 1386
sum(mic.alpha[,1]&Ath.alpha.th[,1]) # d1 1410

## delta phi and alpha > 1
sum(mic.phi[,3]&Ath.alpha.th[,3]) # d7 1215
sum(mic.phi[,2]&Ath.alpha.th[,2]) # d3 1282
sum(mic.phi[,1]&Ath.alpha.th[,1]) # d1 1363

### venn diagram of 7,3,1 day conditions, Fig. 6C -----------------------------------------------
library(VennDiagram)
vennlist <- list(d1=rownames(alpha.SO[mic.gene$Ath_d1,]),
                 d3=rownames(alpha.SO[mic.gene$Ath_d3,]),
                 d7=rownames(alpha.SO[mic.gene$Ath_d7,]))
venn.diagram(vennlist, filename="plot_Fig/micGene_venn_d137.tif", 
             fill=1:3, scaled=T, height = 1500, width = 1500)


# schatter plot of alpha of mimicked genes   d1 vs d7, Fig. 6D -----------------------------------
# scatter plot of alpha
tiff(file = "plot_Fig/Scatterplots/Schatterplot of micgene d1vsd7.tif", 
     pointsize = 19, width = 480, height = 480)
par(pty = "s", lwd = 2)
com <- rownames(alpha.SO[mic.gene[[1]]|mic.gene[[3]],])
x <- alpha.SO[com,2] 
y <- alpha.SO[com,4]

coltmp <- densCols(x, y, colramp = colorRampPalette(c("blue2","cyan","yellow", "white")))

plot(0, 0, type = "n", xlim = c(0, 10), ylim = c(0, 10), xlab = colnames(alpha.all)[2], ylab = colnames(alpha.all)[4],
     main = "Scatterplot of alpha", lwd = 2, las = 1)
points(x, y, col = coltmp, pch = 1, cex = 0.9, lwd= 1)
points(alpha.SO["Ahg470177",2], alpha.SO["Ahg470177",4], col = "maroon1", pch = 16) # LHY
points(alpha.SO["Ahg941635",2], alpha.SO["Ahg941635",4], col = "maroon1", pch = 16) # SUS1
points(alpha.SO["Ahg927111",2], alpha.SO["Ahg927111",4], col = "maroon1", pch = 16) # ADH1
points(alpha.SO["Ahg914363",2], alpha.SO["Ahg914363",4], col = "maroon1", pch = 16) # JAL5
points(alpha.SO["Ahg493293",2], alpha.SO["Ahg493293",4], col = "maroon1", pch = 16) # SWEET17 
points(alpha.SO["Ahg489834",2], alpha.SO["Ahg489834",4], col = "maroon1", pch = 16) # TT8
abline(0,1, col="lightslategrey", lwd= 1.5) 
dev.off()
par(oldpar)


# define non-mimicked gene ----------------------------------------------------------------------
sum(!(mic.gene$Ath_d1 | mic.gene$Ath_d3 | mic.gene$Ath_d7)) # 2685
tmp <- mic.gene[!(mic.gene$Ath_d1 | mic.gene$Ath_d3 | mic.gene$Ath_d7),]
sum(tmp == TRUE)

nonmic.gene <- rownames(tmp)
hist(delta.alpha[nonmic.gene,1], breaks = 50)


# GO analysis ----------------------------------------------------------------------------------
dir.create("GOanalysis")
source("input/GOanalysis_functions.R")
load("input/ulg.TAIR_190205") # objyect "ulg"

## make exgene subset of ulg
tmp <- substr(exbase$Ath_ID, 1,9)
exbase$Ath.locus <- tmp
ulg2 <- ulg[is.element(ulg[,"locus"], exbase$Ath.locus),]
dim(ulg2)                 # 830195      2
length(unique(ulg2[,1]))  # 14472

## set function
fn.GOanalysis <- function(tmp, fn, ulg2){
  # fisher's exact test 
  tage <- exbase[tmp,]$Ath.locus    
  result <- ng.mft(cgt=ulg2, gn.test=tage)
  cat(min(result[,1]))               # check minimum value
  result.sort <- result[order(result[,1]),]
  
  # output results as a csv file
  result2 <- ng.prepGOtestOutTable(result.sort, alpha = 0.05)
  write.csv(result2, file = fn)
}

## d7 mimicked gene, Supplementary Table 3
tmp <- rownames(alpha.SO[mic.gene$Ath_d7,])
fn <- "GOanalysis/GoOut_d7mic.csv"
fn.GOanalysis(tmp = tmp, fn = fn, ulg2 = ulg2)

## total mimicked gene, Supplementary Table 9
tmp <- rownames(alpha.SO[mic.gene$Ath_d1 | mic.gene$Ath_d3 | mic.gene$Ath_d7,])
fn <- "GOanalysis/GoOut_totalmic.csv"
fn.GOanalysis(tmp = tmp, fn = fn, ulg2 = ulg2)

## non-mimicked gene, Supplementary Table 10 
tmp <- rownames(alpha.SO[!(mic.gene$Ath_d1 | mic.gene$Ath_d3 | mic.gene$Ath_d7),])
fn <- "GOanalysis/GoOut_all_non-mim.csv"
fn.GOanalysis(tmp = tmp, fn = fn, ulg2 = ulg2)

## non-mimicked gene(phase), Supplementary Table 11
tmp <- mic.alpha&!(mic.phi)&Ath.alpha.th
tmp2 <- apply(tmp, 1, any)
sum(tmp2) # 747

tmp <- names(tmp2[tmp2])
fn <- "GOanalysis/GoOut_non_phese.csv"
fn.GOanalysis(tmp = tmp, fn = fn, ulg2 = ulg2)


# GO plot from REVIGO file, Fig. 6E -------------------------------------------------------
# The above GO analysis output files were summarized using REVIGO (http://revigo.irb.hr) (Supek et al., 2011) 
# allowed similarity = 0.5,  Others are default.
# The REVIGO output files(only biological process) were placed in a folder "GOanalysis/REVIGO_Biological_0.5"
# This script was written with reference to the previous protocol, Bonnot T, Gillard M, Nagel D. 2019. "A Simple Protocol for Informative Visualization of Enriched Gene Ontology Terms." Bio-Protocol e3429. 

## read REVIGO file
fn <- c("totalmic", "all_non-mim", "non_phese")

GOlist <- NULL
for(i in fn){
  n <- sprintf("GOanalysis/GoOut_%s.csv", i)
  tmp <- read.csv(n, row.names = 1,stringsAsFactors = F) 
  n <- sprintf("GOanalysis/REVIGO_Biological_0.5/%s_REVIGO.csv", i)
  tmp2 <- read.csv(n, stringsAsFactors = F)
  if(i == "totalmic"){
  }else{
  tmp2[tmp2$term_ID == "GO:0016143",]$eliminated <- 1
  tmp2[tmp2$term_ID == "GO:0019760",]$eliminated <- 0
  }
  
  tmp2 <- tmp2[tmp2$eliminated == 0,]
  tmp <- tmp[tmp2$term_ID,]
  tmp$kind <- rep(i, nrow(tmp))
  rownames(tmp) <- NULL
  tmp <- tmp[order(tmp$Adjusted.P.value),]
  
  GOlist <- rbind(GOlist, tmp)
}

colnames(GOlist) <- c("Adjusted.P.value", "ID", "Description", "Gene_number",
                      "A", "B", "U", "kind" )

GOlist$mlog10P <- -log10(GOlist$Adjusted.P.value)
GOlist$ylabel <- paste(GOlist$Description, ":", GOlist$ID)
GOlist$ylabel <- as.factor(GOlist$ylabel)
GOlist$kind <- factor(GOlist$kind, levels = c("totalmic", "all_non-mim", "non_phese"))
xlabel <- c("total SI gnene","total non-SI gene","non-SI gnene(phase)")
GOlist$ylabel <- factor(GOlist$ylabel, levels = unique(GOlist$ylabel))

## draw GO plot (only biologicl process)
pdf("GOanalysis/REVIGO_Biological_0.5/GOplot_Biological.pdf", width = 5.5, height = 6)
ggplot(GOlist, aes(x = ylabel, y = kind)) +
  geom_point(data=GOlist,aes(x=ylabel, y=kind, size = Gene_number, colour = mlog10P), alpha=.7)+
  scale_y_discrete(labels =xlabel)+
  scale_color_gradient(low = "cyan", high = "blue", limits=c(0, NA))+
  coord_flip()+
  theme_bw()+
  theme(axis.ticks.length=unit(-0.1, "cm"),
        axis.text.x = element_text(margin=margin(5,5,0,5,"pt"),angle = 45, hjust = 1),
        axis.text.y = element_text(margin=margin(5,5,5,5,"pt")),
        axis.text = element_text(color = "black"),
        axis.title.x=element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8))+
  xlab("GO biological processes")+
  labs(color=expression(atop("-log"[10],paste("(Adj. p)"))), size="Number\nof genes")
dev.off()


# set plot function for raw data points and fitted curve plot   -----------------------------------
## predicted lines plot
fun.predline.plot <- function(tg.ahg, tg){
  plot(1:365, Ahg_cosfit1$predmx[tg.ahg,], type = "l", main = "Fitted curve plot", ylim = c(0,15),
       xlab = "days from Jan 1", ylab = "log2(rpm+1)", lwd =3, col = cols2[1], las = 1)
  lines(1:365, Athd1_cosfit1$predmx[tg,], col=cols2[2], lwd =3)
  lines(1:365, Athd3_cosfit1$predmx[tg,], col=cols2[3], lwd =3)
  lines(1:365, Athd7_cosfit1$predmx[tg,], col=cols2[4], lwd =3)
  abline(v = phi.all[tg.ahg,]$Ahg, lwd =3, col = cols2[1], lty = 4)
  abline(v = phi.all[tg.ahg,]$Ath_d1, lwd =3, col=cols2[2], lty = 3)
  abline(v = phi.all[tg.ahg,]$Ath_d3, lwd =3, col=cols2[3], lty = 3)
  abline(v = phi.all[tg.ahg,]$Ath_d7, lwd =3, col=cols2[4], lty = 3)
  legend("topright", legend=c("Ahg", "Ath d7", "Ath d3", "Ath d1"), text.col=cols2[c(1,4,3,2)], 
         bty="n", cex = 0.8)
} 

## predicted lines plot (Ahg and 7d only)
fun.predline7.plot <- function(tg.ahg, tg){
  plot(1:365, Ahg_cosfit1$predmx[tg.ahg,], type = "l", main = "Fitted curve plot", ylim = c(0,15),
       xlab = "days from Jan 1", ylab = "log2(rpm+1)", lwd =3, col = cols2[1],las = 1)
  lines(1:365, Athd7_cosfit1$predmx[tg,], col=cols2[4], lwd =3)
  abline(v = phi.all[tg.ahg,]$Ahg, lwd =3, col = cols2[1], lty = 4)
  abline(v = phi.all[tg.ahg,]$Ath_d7, lwd =3, col=cols2[4], lty = 3)
  legend("topright", legend=c("Ahg", "Ath d7"), text.col = cols2[c(1,4)], bty="n", cex = 0.8)
} 

## Ahg raw points and predicted line plot
fun.Ahgplot <- function(tg.ahg){
  l2cpmAhg5.5 <- log2(cpmAhg5.5+1)      # set data
  xtime <- sprintf("%s-%s-%s", 2012, atAhg5.5$month, atAhg5.5$day)    # set time
  xtime <- as.POSIXct(xtime, format="%Y-%m-%d")
  xtime.yday <- yday(xtime)      # convert day into number
  
  # plot raw data and fitted curve
  plot(xtime.yday, l2cpmAhg5.5[tg.ahg, ], main = sprintf("%s", tg.ahg), ylim = c(0,15),
       xlab = "days from Jan 1", ylab = "log2(rpm+1)", lwd = 1, cex = 0.7, col = "gray15", las = 1)
  lines(1:365, Ahg_cosfit1$predmx[tg.ahg,], col=cols2[1], lwd = 3)
  abline(v = phi.all[tg.ahg,]$Ahg, lwd =2, col = cols2[1], lty = 4)
  tmp <- sprintf("%.2f", Ahg_cosfit1$sinFit[tg.ahg,]$alpha)
  tmp2 <- sprintf("%.0f", Ahg_cosfit1$sinFit[tg.ahg,]$modi.phi)
  tmp3 <- parse(text = paste0("list(italic('\u03B1') == ", tmp, ", italic('  \u03C6') == ", tmp2, ")")) 
  mtext(tmp3, side = 3, line = 0, cex = 0.75)
}

## Ath raw points and predicted line plot
cosfitAth <- list(d1 = Athd1_cosfit1, d3 = Athd3_cosfit1, d7 = Athd7_cosfit1)  # bind fit

fun.Athplot <- function(n,tg){
  d <- days[n]　　　　　　　　　  　# set data
  tmpat <- s.at[s.at$day==d,]
  tmpcpm <- cpmAth5.5[, rownames(tmpat)]　　
  tmpl2cpm <- log2(tmpcpm+1)
  
  xtime <- sprintf("%s-%s-%s", rep(2017,nrow(tmpat)), tmpat$num.month,       # set time
                   rep(15,nrow(tmpat)))
  xtime <- as.POSIXct(xtime, format="%Y-%m-%d")
  xtime.yday <- yday(xtime)      # convert day into number
  
  # plot raw data and fitted curve
  plot(xtime.yday, tmpl2cpm[tg, ], main = sprintf("%s (Ath d%s)", tg, d), ylim = c(0,15),
       xlab = "days from Jan 1", ylab = "log2(rpm+1)", lwd = 1, col = "gray15", cex = 0.7, las = 1)
  lines(1:365, cosfitAth[[n]]$predmx[tg,], col=cols2[n+1], lwd = 3)
  abline(v = cosfitAth[[n]]$sinFit[tg,]$modi.phi, lwd =2, col = cols2[n+1], lty = 3)
  tmp <- sprintf("%.2f", cosfitAth[[n]]$sinFit[tg,]$alpha)
  tmp2 <- sprintf("%.0f", cosfitAth[[n]]$sinFit[tg,]$modi.phi)
  tmp3 <- parse(text = paste0("list(italic('\u03B1') == ", tmp, ", italic('  \u03C6') == ", tmp2, ")")) 
  mtext(tmp3, side = 3, line = 0, cex = 0.75)
}

## Ath HT raw points and predicted line plot
HTat <- atAth5.5[atAth5.5$set == "HT",]      # set HT data
HTcpm <- cpmAth5.5[, rownames(HTat)]
HTl2cpm <- log2(HTcpm+1)

HT.yday <- yday(as.POSIXct("2017-8-15", format="%Y-%m-%d")) # set HT time
HT.xtime <- rep(HT.yday, 5)

fun.AthHTplot <- function(n,tg){ 
  d <- days[n]                            # set data
  tmpat <- s.at[s.at$day==d,]
  tmpcpm <- cpmAth5.5[, rownames(tmpat)]
  tmpl2cpm <- log2(tmpcpm+1)
  
  # set time
  xtime <- sprintf("%s-%s-%s", rep(2017,nrow(tmpat)), tmpat$num.month, rep(15,nrow(tmpat)))
  xtime <- as.POSIXct(xtime, format="%Y-%m-%d")
  xtime.yday <- yday(xtime)      # convert day into number
  
  # plot raw data and fitted curve
  plot(xtime.yday, tmpl2cpm[tg, ], main = sprintf("%s (Ath d%s)", tg, d), ylim = c(0,15),
       xlab = "days from Jan 1", ylab = "log2(rpm+1)", lwd = 1, cex = 0.7, col = "gray15", las = 1)
  lines(1:365, cosfitAth[[n]]$predmx[tg,], col=cols2[n+1], lwd = 3)
  abline(v = cosfitAth[[n]]$sinFit[tg,]$modi.phi, lwd =2, col = cols2[n+1], lty = 3)
  tmp <- sprintf("%.2f", cosfitAth[[n]]$sinFit[tg,]$alpha)
  tmp2 <- sprintf("%.0f", cosfitAth[[n]]$sinFit[tg,]$modi.phi)
  tmp3 <- parse(text = paste0("list(italic('\u03B1') == ", tmp, ", italic('  \u03C6') == ", tmp2, ")")) 
  mtext(tmp3, side = 3, line = 0, cex = 0.75)
  points(HT.xtime, HTl2cpm[tg, HTat$month == "HT32"], col = "deeppink1", cex = 1, pch = 6)
  points(HT.xtime, HTl2cpm[tg, HTat$month == "HT35"], col = "red", cex = 1, pch = 0)
}

## Temperature vs expression plot 
fun.heatExplot <- function(n,tg){ 
  d <- days[n]                     # set data
  tmpat <- s.at[s.at$day==d,]
  tmpcpm <- cpmAth5.5[, rownames(tmpat)]
  tmpl2cpm <- log2(tmpcpm+1)
  
  # draw temperature vs expression plot
  plot(as.numeric(as.character(tmpat$day_temp)), tmpl2cpm[tg, ], main = "", 
       ylim = c(0,15), xlim = c(0,35), xlab = expression("Temperature"~~(degree*"C")), ylab = "log2(rpm+1)",
       col = "gray15", cex = 0.8, las = 1)
  points(rep(32,5), HTl2cpm[tg, HTat$month == "HT32"], col = "deeppink1", cex = 1, pch = 6)
  points(rep(35,5), HTl2cpm[tg, HTat$month == "HT35"], col = "red", cex = 1, pch = 0)
}

# Draw raw data points and fitted curve plot-------------------------------------------------------
dir.create("plot_Fig/RawdataPlot")

# read target gene list
tgl <- read.csv("input/TargetGeneList.csv", stringsAsFactors = F)
rownames(tgl) <- tgl$Ahg_ID

## example plot, Fig. 3A
tg.ahg <- "Ahg927111"
tg <- reci[tg.ahg,]$Ath_ID

tiff("plot_Fig/RawdataPlot/example_plot.tif", pointsize = 26, width = 480*2, height = 380)
par(lwd=2, mfrow=c(1,3))
fun.predline7.plot(tg.ahg = tg.ahg, tg = tg)
fun.Ahgplot(tg.ahg = tg.ahg) 
fun.Athplot(n = 3, tg =tg)
dev.off()

## LHY plot, Fig. 3H
tg.ahg <- tgl[tgl$name == "LHY",]$Ahg_ID
tg <- reci[tg.ahg,]$Ath_ID

tiff("plot_Fig/RawdataPlot/d7_LHY_Ahg470177.tif", pointsize = 26, width = 480*2, height = 400)
par(lwd=1.5, oma=c(0,0,1,0), mfrow=c(1,3))
fun.predline7.plot(tg.ahg = tg.ahg, tg = tg)
fun.Ahgplot(tg.ahg = tg.ahg) 
fun.Athplot(n = 3, tg =tg)
mtext(sprintf("%s,  LHY,  %s", tg, desAth[tg,]$Note), 
      outer=T, cex=0.9, line=-0.8, adj=0, at = 0.01) 
dev.off()
par(oldpar)

## HT plot, Fig. 4C-F
tmp <- tgl[tgl$group == "HT",]
for (i in 1:nrow(tmp)) {
  tg.ahg <- tmp$Ahg_ID[i]
  tg <- reci[tg.ahg,]$Ath_ID
  
  fn <- sprintf("plot_Fig/RawdataPlot/HT_%s_%s.tif", tmp$name[i],tmp$Ahg_ID[i])
  tiff(fn, pointsize = 26, width = 480*2.5, height = 400)
  par(lwd=1.5, oma=c(0,0,1,0), mfrow=c(1,4))
  fun.predline7.plot(tg.ahg = tg.ahg, tg = tg)
  fun.Ahgplot(tg.ahg = tg.ahg) 
  fun.AthHTplot(n = 3, tg = tg)
  fun.heatExplot(n = 3, tg = tg)
  mtext(sprintf("%s,  %s,  %s", tg, tmp$name[i], desAth[tg,]$Note), 
        outer=T, cex=0.9, line=-0.8, adj=0, at = 0.01) 
  dev.off()
}
par(oldpar)

## Fig.7 plot
tmp <- tgl[tgl$group == "Fig7",]
for (i in 1:nrow(tmp)) {
  tg.ahg <- tmp$Ahg_ID[i]
  tg <- reci[tg.ahg,]$Ath_ID
  
  tmp2 <- sprintf("plot_Fig/RawdataPlot/Fig7_%s_%s.tif", tmp$name[i], tmp$Ahg_ID[i])
  tiff(tmp2, pointsize = 24, width = 480*3, height = 380)
  par(lwd=1.5, oma=c(0,0,1,0),mfrow=c(1,5))
  fun.predline.plot(tg.ahg = tg.ahg, tg = tg)
  fun.Ahgplot(tg.ahg = tg.ahg) 
  fun.Athplot(n = 3, tg =tg)
  fun.Athplot(n = 2, tg =tg)
  fun.Athplot(n = 1, tg =tg)
  mtext(sprintf("%s,  %s,  %s", tg, tmp$name[i], desAth[tg,]$Note), 
        outer=T, cex=0.9, line=-0.8, adj=0, at = 0.01) 
  dev.off()
}

## Fig.S7 plot
tmp <- tgl[tgl$group == "FigS7",]
tmp <- rbind(tmp, tgl[tgl$name == "TT8",])
for (i in 1:nrow(tmp)) {
  tg.ahg <- tmp$Ahg_ID[i]
  tg <- reci[tg.ahg,]$Ath_ID
  
  tmp2 <- sprintf("plot_Fig/RawdataPlot/FigS7_%s_%s.tif", tmp$name[i], tmp$Ahg_ID[i])
  tiff(tmp2, pointsize = 18, width = 290, height = 380)
  par(lwd=1.5, oma=c(0,0,1,0), mar=c(3, 3, 3.4, 1), mgp=c(2.0, 0.7, 0))
  fun.predline.plot(tg.ahg = tg.ahg, tg = tg)
  mtext(sprintf("%s,  %s", tg, tmp$name[i]), 
        outer=T, cex=1.3, line=-0.8) 
  dev.off()
}

## Fig.S8-10 plot
tmp <- tgl[tgl$group == "FigS7",]
for (i in 1:nrow(tmp)) {
  tg.ahg <- tmp$Ahg_ID[i]
  tg <- reci[tg.ahg,]$Ath_ID
  
  tmp2 <- sprintf("plot_Fig/RawdataPlot/FigS8-10_%s_%s.tif", tmp$name[i], tmp$Ahg_ID[i])
  tiff(tmp2, pointsize = 24, width = 480*3, height = 380)
  par(lwd=1.5, oma=c(0,0,1,0),mfrow=c(1,5))
  fun.predline.plot(tg.ahg = tg.ahg, tg = tg)
  fun.Ahgplot(tg.ahg = tg.ahg) 
  fun.Athplot(n = 3, tg =tg)
  fun.Athplot(n = 2, tg =tg)
  fun.Athplot(n = 1, tg =tg)
  mtext(sprintf("%s,  %s,  %s", tg, tmp$name[i], desAth[tg,]$Note), 
        outer=T, cex=0.9, line=-0.8, adj=0, at = 0.01) 
  dev.off()
}

###### END ####################################################################################
# Package versions
sessionInfo()
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.6
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] VennDiagram_1.6.20  futile.logger_1.4.3 lubridate_1.7.4    
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.3           rstudioapi_0.11      magrittr_1.5         tidyselect_1.0.0    
# [5] munsell_0.5.0        colorspace_1.4-1     R6_2.4.1             rlang_0.4.5         
# [9] stringr_1.4.0        dplyr_0.8.5          tools_3.6.1          gtable_0.3.0        
# [13] lambda.r_1.2.4       assertthat_0.2.1     tibble_2.1.3         lifecycle_0.2.0     
# [17] crayon_1.3.4         purrr_0.3.3          ggplot2_3.3.0        formatR_1.7         
# [21] futile.options_1.0.1 glue_1.3.1           stringi_1.4.6        compiler_3.6.1      
# [25] pillar_1.4.3         scales_1.1.0         pkgconfig_2.0.3   
