##################################################################
######
######    TEST FOR 'stigene' LIBRARY
######
######
##################################################################

library(network)
library(ggplot2)
library(plyr)
library(stigene,lib.loc = "./lib")

t0 <- Sys.time()
# path to model input files:
folder_inputs = "./test-inputs/"
folder_calib = "../calibration/"

# founder population
founder.size <- 300
founder.fem.prop <- 0.5
founder.csw.prop <- 0.02

# scenario 
# (includes intervention filenames 
#  that will be done during this simualtion)
scen.file <- "in_scenario_baseline.csv"

# simply run the model:
x <- stiagent_runsim(params = list(folder_inputs   = folder_inputs,
                                   folder_calib    = folder_calib,
								   founder_size    = founder.size,
								   founder_femprop = founder.fem.prop,
								   founder_cswprop = founder.csw.prop,
                                   scenario_file   = scen.file,
                                   MC_id = 1,  # <-- MC id (*not* number of MC iterations!!!)
								   displayProgress = 11)
                     )

par(mfrow=c(2,3))
with(x$df_sim, 
	 plot(time, nAlive, typ='l'))
with(x$df_sim, 
	 plot(time, nSp, typ='l'))
with(x$df_sim, 
	 plot(time, HIV, typ='l'))
with(x$df_sim, 
	 plot(time, Tp, typ='l'))
with(x$df_sim, 
	 plot(time, nNewBorn, typ='l'))


# secondary cases:
sec <- x$sec_cases
z <- sec[[1]]

sec.count <- numeric(length = length(z))
for(i in 1:length(z)) {
	if(length(z[[i]]>0)) {
		sec.count[i] <- length(z[[i]])
	}
}
sec.count <- sec.count[sec.count>0]
hist(sec.count,breaks = 50)
abline(v=mean(sec.count))

# Data frame of WIW
infector <- vector()
infectee <- vector()
sti      <- vector()

cnt <- 1
nsti <- length(sec)
for(s in 1:nsti){
	tmp <- sec[[s]]
	for(i in 1:length(tmp)){
		ns <- length(tmp[i][[1]])
		if(ns>0){
			for(j in 1:ns){
			infector[cnt] <- names(tmp[i])
			infectee[cnt] <- tmp[i][[1]][j]
			sti[cnt]      <- names(sec)[s]	
			cnt <- cnt+1
			}
		}
	}
}
WIW <- data.frame(sti, infector=as.character(infector), infectee = as.character(infectee))


wiw.hiv <- subset(WIW, sti=='HIV')
wiw.hiv <- wiw.hiv[,2:3]

net     <- network(wiw.hiv, directed = TRUE,loops = F)
par(mfrow=c(1,1))
plot(net, 
	 displaylabels = TRUE, 
	 label.cex = 0.5,
	 vertex.cex= 0.4, 
	 edge.lwd=0.5, edge.col=rgb(0,0,0,0.2),
	 arrowhead.cex=0.5)


# infection (acquisition) times
acq <- x$acquisition_times
acq$Tp

infectee <- vector()
acqTime <- vector()
sti <- vector()
cnt <- 1
nsti <- length(acq)
for(s in 1:nsti){
	tmp <- acq[[s]]
	for(i in 1:length(tmp)){
		ns <- length(tmp[i][[1]])
		if(ns>0){
			for(j in 1:ns){
				infectee[cnt] <- names(tmp[i])
				acqTime[cnt] <- tmp[i][[1]][j]
				sti[cnt]      <- names(acq)[s]	
				cnt <- cnt+1
			}
		}
	}
}
acqTimes <- data.frame(sti, infectee=as.character(infectee), acqTime = as.numeric(acqTime))


# Merge infector, infectee with acquisition times:
chain.transm <- join(WIW,acqTimes, by = c('infectee','sti'))

chain.transm$round.acqTime <- round(chain.transm$acqTime/5,2)*5

z <- ddply(chain.transm,c('round.acqTime','sti'),summarize, inc=length(infectee) )
g <- ggplot(z) + geom_step(aes(x=round.acqTime, y=inc, colour=sti)) #+geom_point()
g <- g + geom_line(data = as.data.frame(x$df_sim), aes(x=time, y=Tp))
plot(g)


# Genetic information:
g <- x$genomes
guid <- x$genomes_UID


# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
t1 <- Sys.time()
message(paste("time elapsed:",round(t1-t0,1),"sec"))
