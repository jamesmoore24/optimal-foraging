rm(list=ls())
load("onlinedata.Rdata")
################################################ Figure 1
dccc <- entrydata
dccc <- data.frame(dccc$sid, dccc$entry)
names(dccc) <- c("sid", "entry")
simback1 <- rep(NA, nrow(dccc))
simback2 <- rep(NA, nrow(dccc))
simback3 <- rep(NA, nrow(dccc))
simback4 <- rep(NA, nrow(dccc))
simback5 <- rep(NA, nrow(dccc))

for(i in 2:nrow(dccc)){
  if (dccc$entry[i] %in% rownames(ancos) && dccc$entry[i-1] %in% rownames(ancos)){
	if (dccc$sid[i] == dccc$sid[i-1]){
		simback1[i] <- ancos[toString(dccc$entry[i]), toString(dccc$entry[i-1])]
		}
		}
	}

for(i in 3:nrow(dccc)){
  if (dccc$entry[i] %in% rownames(ancos) && dccc$entry[i-2] %in% rownames(ancos)){
  if (dccc$sid[i] == dccc$sid[i-2]){
		simback2[i] <- ancos[toString(dccc$entry[i]), toString(dccc$entry[i-2])]
		}
		}
	}
	
for(i in 4:nrow(dccc)){
  if (dccc$entry[i] %in% rownames(ancos) && dccc$entry[i-3] %in% rownames(ancos)){
  if (dccc$sid[i] == dccc$sid[i-3]){
		simback3[i] <- ancos[toString(dccc$entry[i]), toString(dccc$entry[i-3])]
		}
		}
	}
	
for(i in 5:nrow(dccc)){
  if (dccc$entry[i] %in% rownames(ancos) && dccc$entry[i-4] %in% rownames(ancos)){
  if (dccc$sid[i] == dccc$sid[i-4]){
		simback4[i] <- ancos[toString(dccc$entry[i]), toString(dccc$entry[i-4])]
		}
		}
	}
	
for(i in 6:nrow(dccc)){
  if (dccc$entry[i] %in% rownames(ancos) && dccc$entry[i-5] %in% rownames(ancos)){
  if (dccc$sid[i] == dccc$sid[i-5]){
		simback5[i] <- ancos[toString(dccc$entry[i]), toString(dccc$entry[i-5])]
		}
		}
	}
	
dddd <- data.frame(dccc, simback1, simback2, simback3, simback4, simback5)

sb1 <- with(dddd, tapply(simback1, sid, mean, na.rm=T))
sb2 <- with(dddd, tapply(simback2, sid, mean, na.rm=T))
sb3 <- with(dddd, tapply(simback3, sid, mean, na.rm=T))
sb4 <- with(dddd, tapply(simback4, sid, mean, na.rm=T))
sb5 <- with(dddd, tapply(simback5, sid, mean, na.rm=T))

x <- rep(1:5, each=141)

a_table <- data.frame(x, sim = c(sb1, sb2, sb3, sb4, sb5))

anova(lm(sim~factor(x), data=a_table))  # hypothesis 1

ms <- c(mean(sb1, na.rm=T),
mean(sb2, na.rm=T),
mean(sb3, na.rm=T),
 mean(sb4, na.rm=T),
mean(sb5, na.rm=T)
        )
sdd <- c(sd(sb1, na.rm=T),
         sd(sb2, na.rm=T),
         sd(sb3, na.rm=T),
         sd(sb4, na.rm=T),
         sd(sb5, na.rm=T)
         )
see <- sdd / sqrt(141)

par(mar=c(5,5,2,2))
x <- barplot(ms[5:1], ylab="BEAGLE similarity", xlab="Item's position preceding most recent item", cex.lab = 1.5, ylim=c(0,0.5),names=c("5", "4", "3", "2", "1"), col = "white")
arrows(x, ms[5:1]+see[5:1], x, ms[5:1]-see[5:1], angle=90, code=3, length=.1, lwd = 1.4)


################################################ Figure 2

db <- entrydata
rowsim <- rep(0, nrow(db))
lastid <- 0
for(i in 1:nrow(db)){
  if (lastid != toString(db$sid[i])){
		lastid <- db$sid[i]
		subcos <- ancos
		}
	subsub <- subcos[toString(db$entry[i]),]
	subsub <- subsub[subsub < 1]
	rowsim[i] <- mean(subsub)
	if( toString(db$entry[i] %in% colnames(subcos))){
		x <- 1:ncol(subcos)
		y <- x[toString(db$entry[i]) == colnames(subcos)]
		subcos <- subcos[,-y]
	}	
}

dtest <- data.frame(entrydata, rowsim)

#save(dtest, file = "dtest")




db <- dtest

hits <- matrix(999, nrow(db), ncol = 8)

for( i in 1:(nrow(db))){
	if(db$fpatchnum[i] > 0){
	if (db$fpatchitem[i] < 5) hits[i,4+db$fpatchitem[i]] <- db$rowsim[i]
	if (db$fitemsfromend[i] < 5 && db$fpatchitem[i] != 1) hits[i,5-db$fitemsfromend[i]] <- db$rowsim[i]
	}
	}


msss <- c(
mean(hits[hits[,1] < 999, 1], na.rm = T),
mean(hits[hits[,2] < 999, 2], na.rm = T),
mean(hits[hits[,3] < 999, 3], na.rm = T),
mean(hits[hits[,4] < 999, 4], na.rm = T),
mean(hits[hits[,5] < 999,5], na.rm = T),
mean(hits[hits[,6] < 999, 6], na.rm = T),
mean(hits[hits[,7] < 999, 7], na.rm = T),
mean(hits[hits[,8] < 999, 8], na.rm = T))

msse <- c(
sd(hits[hits[,1] < 999,1], na.rm = T) / sqrt(length(hits[hits[,1] < 999,1])-sum(as.numeric(is.na(hits[,1])))),
sd(hits[hits[,2] < 999,2], na.rm = T) / sqrt(length(hits[hits[,2] < 999,2])-sum(as.numeric(is.na(hits[,2])))),
sd(hits[hits[,3] < 999,3], na.rm = T) / sqrt(length(hits[hits[,3] < 999,3])-sum(as.numeric(is.na(hits[,3])))),
sd(hits[hits[,4] < 999,4], na.rm = T) / sqrt(length(hits[hits[,4] < 999,4])-sum(as.numeric(is.na(hits[,4])))),
sd(hits[hits[,5] < 999,5], na.rm = T) / sqrt(length(hits[hits[,5] < 999,5])-sum(as.numeric(is.na(hits[,5])))),
sd(hits[hits[,6] < 999,6], na.rm = T) / sqrt(length(hits[hits[,6] < 999,6])-sum(as.numeric(is.na(hits[,6])))),
sd(hits[hits[,7] < 999,7], na.rm = T) / sqrt(length(hits[hits[,7] < 999,7])-sum(as.numeric(is.na(hits[,7])))),
sd(hits[hits[,8] < 999,8], na.rm = T) / sqrt(length(hits[hits[,8] < 999,8])-sum(as.numeric(is.na(hits[,8])))))

 par(mar=c(5,5,2,2))
x <- barplot(msss[c(3,4,5,6,7)], names.arg = c("-2", "-1", "1", "2", "3"), xlab = "Order of entry relative to patch switch", ylab = "Global Similarity", cex.lab = 1.7, ylim = c(0.20, 0.212), cex.axis = 1.3, cex.names = 1.3, col="grey20")
arrows(x, msss[c(3,4,5,6,7)]+msse[c(3,4,5,6,7)], x, msss[c(3,4,5,6,7)]-msse[c(3,4,5,6,7)], code = 3, angle = 90, length = .1, lwd = 2)
abline(h= 1, lty = 2)




################################################ Figure 3

dccc <- entrydata

cosim <- rep(0, nrow(dccc))
for(i in 2:nrow(dccc)){
  if (dccc$entry[i] %in% rownames(ancos) && dccc$entry[i-1] %in% rownames(ancos)){
  if (dccc$sid[i] == dccc$sid[i-1]){
		cosim[i] <- ancos[toString(dccc$entry[i]), toString(dccc$entry[i-1])]
		}
		}
	}

db <- data.frame(dccc, cosim)

## ASSIGN MEAN SIMILARITIES

meansims <- tapply(db$cosim, db$sid, mean, na.rm = T)

## assign mean to each subject
meansim <- rep(NA, nrow(db))
for (i in 1:nrow(db)){
  meansim[i] <- meansims[toString(db$sid[i])]
}

dbb <- data.frame(db, meansim)

######

dbb <- subset(dbb, catitem > 1)
hits <- matrix(999, nrow(dbb), ncol = 8)

for( i in 1:(nrow(dbb))){
	if (dbb$fpatchitem[i] < 5) hits[i,4+dbb$fpatchitem[i]] <- dbb$cosim[i] / dbb$meansim[i]
	if (dbb$fitemsfromend[i] < 5) hits[i,5-dbb$fitemsfromend[i]] <- dbb$cosim[i] / dbb$meansim[i]
	}


msss <- c(
mean(hits[hits[,1] < 999, 1], na.rm = T),
mean(hits[hits[,2] < 999, 2], na.rm = T),
mean(hits[hits[,3] < 999, 3], na.rm = T),
mean(hits[hits[,4] < 999, 4], na.rm = T),
mean(hits[hits[,5] < 999,5], na.rm = T),
mean(hits[hits[,6] < 999, 6], na.rm = T),
mean(hits[hits[,7] < 999, 7], na.rm = T),
mean(hits[hits[,8] < 999, 8], na.rm = T))

msse <- c(
sd(hits[hits[,1] < 999,1], na.rm = T) / sqrt(length(hits[hits[,1] < 999,1])-sum(as.numeric(is.na(hits[,1])))),
sd(hits[hits[,2] < 999,2], na.rm = T) / sqrt(length(hits[hits[,2] < 999,2])-sum(as.numeric(is.na(hits[,2])))),
sd(hits[hits[,3] < 999,3], na.rm = T) / sqrt(length(hits[hits[,3] < 999,3])-sum(as.numeric(is.na(hits[,3])))),
sd(hits[hits[,4] < 999,4], na.rm = T) / sqrt(length(hits[hits[,4] < 999,4])-sum(as.numeric(is.na(hits[,4])))),
sd(hits[hits[,5] < 999,5], na.rm = T) / sqrt(length(hits[hits[,5] < 999,5])-sum(as.numeric(is.na(hits[,5])))),
sd(hits[hits[,6] < 999,6], na.rm = T) / sqrt(length(hits[hits[,6] < 999,6])-sum(as.numeric(is.na(hits[,6])))),
sd(hits[hits[,7] < 999,7], na.rm = T) / sqrt(length(hits[hits[,7] < 999,7])-sum(as.numeric(is.na(hits[,7])))),
sd(hits[hits[,8] < 999,8], na.rm = T) / sqrt(length(hits[hits[,8] < 999,8])-sum(as.numeric(is.na(hits[,8])))))

msss <- msss[3:7]
msse <- msse[3:7]
par(mar=c(5,5,2,2))
x <- barplot(msss, names.arg = c( "-2", "-1", "1", "2", "3"), xlab = "Order of entry relative to patch switch", ylab = "Ratio of pairwise similirity over subject's mean similarity", cex.lab = 1.7, ylim = c(0, 1.4), col = 0, cex.axis = 1.3, cex.names = 1.3)
arrows(x, msss+msse, x, msss-msse, code = 3, angle = 90, length = .1, lwd = 2)
abline(h= 1, lty = 2)


################################################ Figure 4

ssd <- entrydata
hits <- matrix(999, nrow(ssd), ncol = 8)
for( i in 1:(nrow(ssd))){
  if (ssd$fpatchitem[i] < 5) hits[i,4+ssd$fpatchitem[i]] <- ssd$irt[i] / ssd$meanirt[i]  
	if (ssd$fitemsfromend[i] < 5 ) hits[i,5-ssd$fitemsfromend[i]] <- ssd$irt[i] / ssd$meanirt[i]
	}


msss <- c(
mean(hits[hits[,1] < 999, 1], na.rm = T),
mean(hits[hits[,2] < 999, 2], na.rm = T),
mean(hits[hits[,3] < 999, 3], na.rm = T),
mean(hits[hits[,4] < 999, 4], na.rm = T),
mean(hits[hits[,5] < 999,5], na.rm = T),
mean(hits[hits[,6] < 999, 6], na.rm = T),
mean(hits[hits[,7] < 999, 7], na.rm = T),
mean(hits[hits[,8] < 999, 8], na.rm = T))

msse <- c(
sd(hits[hits[,1] < 999,1], na.rm = T) / sqrt(length(hits[hits[,1] < 999,1])-sum(as.numeric(is.na(hits[,1])))),
sd(hits[hits[,2] < 999,2], na.rm = T) / sqrt(length(hits[hits[,2] < 999,2])-sum(as.numeric(is.na(hits[,2])))),
sd(hits[hits[,3] < 999,3], na.rm = T) / sqrt(length(hits[hits[,3] < 999,3])-sum(as.numeric(is.na(hits[,3])))),
sd(hits[hits[,4] < 999,4], na.rm = T) / sqrt(length(hits[hits[,4] < 999,4])-sum(as.numeric(is.na(hits[,4])))),
sd(hits[hits[,5] < 999,5], na.rm = T) / sqrt(length(hits[hits[,5] < 999,5])-sum(as.numeric(is.na(hits[,5])))),
sd(hits[hits[,6] < 999,6], na.rm = T) / sqrt(length(hits[hits[,6] < 999,6])-sum(as.numeric(is.na(hits[,6])))),
sd(hits[hits[,7] < 999,7], na.rm = T) / sqrt(length(hits[hits[,7] < 999,7])-sum(as.numeric(is.na(hits[,7])))),
sd(hits[hits[,8] < 999,8], na.rm = T) / sqrt(length(hits[hits[,8] < 999,8])-sum(as.numeric(is.na(hits[,8])))))

msss <- msss[3:7]
msse <- msse[3:7]
par(mar=c(5,5,2,2))
x <- barplot(msss, names.arg = c( "-2", "-1", "1", "2", "3"), xlab = "Order of entry relative to patch switch", ylab = "Item IRT / Mean IRT", cex.lab = 1.7, ylim = c(0, 1.4), cex.axis = 1.3, cex.names = 1.3, col = "white")
arrows(x, msss+msse, x, msss-msse, code = 3, angle = 90, length = .1, lwd = 2)
abline(h= 1, lty = 2)

sidmeanirtswitch <- with(subset(ssd, fpatchitem == 1), tapply(irt,sid,mean ))
sidmeanirt2nd<- with(subset(ssd, fpatchitem == 2 ), tapply(irt,sid,mean ))
sidmeanirt <- with(subset(ssd, fpatchitem > 1), tapply(meanirt,sid,mean ))

t.test(sidmeanirtswitch, sidmeanirt, paired = T)
t.test(sidmeanirt2nd, sidmeanirt, paired = T)

################################################ Figure 5

# compute mean irt before switches and mean irt overall and compare with total productions

dccc <- entrydata
subs <- unique(dccc$sid)
meanswi.irt <- rep(0, length(subs))
meanoverall.irt <- rep(0, length(subs))
prods <- rep(0, length(subs))

for (i in 1:length(subs)){
  datsub <- subset(dccc, sid == subs[i])
  datsubflast <- subset(datsub, flastitem == 1)
  meanswi.irt[i] <- mean(datsubflast$irt)
  meanoverall.irt[i] <- sum(datsub$irt)/nrow(datsub)
  prods[i] <- nrow(datsub)
  }


plot(abs(meanswi.irt-meanoverall.irt), prods, ylab = "Number of words produced", xlab = "Absolute difference between mean last item IRT and mean overall IRT (sec)", cex.lab = 1.2)
summary(at <- lm(prods~abs(meanswi.irt-meanoverall.irt)))
abline(at)
