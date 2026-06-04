# Simulation Study for Div. & Distr.
# Key purpose is to evaluate accuracy of AOO estimates using:
# 1) iSDM following Landau et al. 2020
#   -includes some amount of repeated detection/non-detection data
#   -also some presence-only data
# 2) standard presence-background SDM. here, prevalence is unknown/inestimable. 
#    We thus we consider a couple common "threshold" rules.
#    -Mean prediction (above this value, assumed occupied)
#    -Max sensitivity + specificity (above this threshold, assumed occupied)
# 
#    

# Notation:
# nsites = number of simulated sites (grid cells) per simulation replicate
# nsims = number of simulation replicates
# noccasions = number of replicates for detection/non-detection data
# nsamp1 = number of sites with detection/non-detection data
# z = true (simulated) occupancy state per site * sim replicate
# x = covariate value per site * sim replicate
# psi = probability of occupancy
# beta0 = intercept ("prevalence" term)
# beta1 = log-odds effect of x on psi
# p = probability of detection for y1
# y1 = simulated detection/non-detection data
# y2 = simulated presence-only data
# b = "thinning" parameter--i.e., probability that y2 for a given site =1 given z= 1 for that site

# We imagine here a small data situation with 200 sites (2x2 km cells) 
# encompassing a species range. We consider 300 simulation replicates.
# We imagine nsamp1= 3 rather than varying this, as the key consideration is p-star
# or 1-(1-p)^nsamp1.
# We allow other simulation inputs to vary. This is very much a toy exploration.

nsites = 200
nsims = 300
noccasions = 3
nsamp1 = round(runif(nsims, 15, 50)) ##basically, have between 15 and 50 sites with occupancy data.
p = runif(nsims, .1, .9) ## p-star ranges between 0.271 and 0.999
b = runif(nsims, .1, .5) ## presence-only data captures between 10 and 50% of species occurrences.
beta0 = rnorm(nsims, 0, .5) ## psi when x=0 ranges between roughly 0.2 and 0.8.
beta1 = rnorm(nsims, 0, 1)

###random variables
x=matrix(NA, nsites, nsims)
z=matrix(NA, nsites, nsims)

for (i in 1:nsites){
  for (s in 1:nsims){
    x[i, s]=rnorm(1)
    z[i, s]=rbinom(1, 1, plogis(beta0[s]+beta1[s]*x[i, s]))
  }
}

TrueAOO<-colSums(z) 
###how many km^2 "occupied". Here, between roughly 200 and 600. 
###between 40 and 160 sites are occupied. We will just focus on the number
###of sites occupied rather than multiplying this by 4.

###simulate presence-only (or background) data
y2<-matrix(NA, nsites, nsims)

for (i in 1:nsites){
  for (s in 1:nsims){
    y2[i, s]<-rbinom(1, 1, z[i, s]*b[s])
  }
}
###between 4 and ~60 cells have observed presences--hist(colSums(y2))

###simulate detection/non-detection data
y1<-array(NA, dim=c(max(nsamp1), nsims, noccasions))

for (s in 1:nsims){
  for (i in 1:nsamp1[s]){
    for (j in 1:noccasions){
    y1[i, s, j]<-rbinom(1, 1, z[i, s]*p[s])
    }
  }
}

y1occobserved<-rep(NA, nsims)
for (s in 1:nsims){
  y1occobserved[s]<-sum(apply(y1[, s, ], 1, max), na.rm=T)}
###here, species is observed at between 1 and 35 sites in detection/non-detection data
### pretty tiny data.

###now we will set up a nimble model to fit the simple iSDM.
library(nimble)

code<-nimbleCode({
  for (s in 1:nsims){
    beta0[s]~dnorm(0, 1)
    beta1[s]~dnorm(0, 1)
    p[s]~dunif(0, 1)
    b[s]~dunif(0, 1)
    for (i in 1:nsites){
      logit(psi[i, s])<-beta0[s]+beta1[s]*x[i, s]
      z[i, s]~dbern(psi[i, s])
      y2[i, s]~dbern(z[i, s]*b[s])
    }
    for (i in 1:nsamp1[s]){
      for (k in 1:3){
        y1[i, s, k]~dbern(p[s]*z[i, s])
      }
    }
  nocc[s]<-sum(z[1:nsites, s])
  }
})

Consts=list(x=x, nsims=nsims, nsites=nsites, nsamp1=nsamp1)

Dat=list(y1=y1, y2=y2)

###here, just going to initialize z=1 at any site where observed
zIn<-y2

##now have to deal with potential discrepancies in occupancy data
tmpz<-matrix(NA, max(nsamp1), nsims)
for (s in 1:nsims){
  tmpz[1:nsamp1[s], s]<-apply(y1[1:nsamp1[s], s, ], 1, max, na.rm=T)
  }

for (s in 1:nsims){
  for (i in 1:nsamp1[s]){
    if (tmpz[i, s]==1){
      zIn[i, s]<-1
    }
  }
}

##again, between roughly 10 & roughly 80 of the 200 sites show some detection.

###easiest to initalize as very overdispersed with low zi
Inits=list(z=zIn, b=runif(300, .1, .5), p=runif(300, .1, .9), beta0=rnorm(300, 0, .5),
           beta1=rnorm(300, 0, .5))

Model<-nimbleModel(code = code, name = "AOO", constants = Consts,
                  data = Dat, inits = Inits, calculate = FALSE)


AOOConf <- configureMCMC(Model, monitors=c("nocc"))
AOOMCMC <- buildMCMC(AOOConf)
CAOOMCMC <- compileNimble(AOOMCMC, Model)

samps<-runMCMC(CAOOMCMC$AOOMCMC, niter=20000, nburnin = 10000, thin=3, nchains=3)

samps2<-rbind(samps[[1]], samps[[2]], samps[[3]])

get_mode <- function(x) {
  uniq_x <- unique(x)
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}


df<-data.frame(AOO=TrueAOO, AOO_hat=apply(samps2, 2, get_mode), nOccsamps=nsamp1, p=p, b=b,
               beta0=beta0, beta1=beta1)
df$rse<-sqrt(((df$AOO - df$AOO_hat)^2))
df$percentrse<-df$rse/df$AOO*100


PAOCI<-t(apply(samps2, 2, HDInterval::hdi))
mean(PAOCI2[, 1]<= df$AOO & df$AOO <=PAOCI2[, 2])
PAOCI2<-t(apply(samps2, 2, quantile, probs=c(.025, .975)))


###Now sdm
library(maxnet)
library(ROCR)
### need to reduce/collapse y1 and supplement presences to y2 if needed.
### This is actually exactly the same as zIN above

ymaxent<-zIn
zmaxent1<-matrix(0, nsites, nsims)
zmaxent2<-matrix(0, nsites, nsims)

for (i in 1:nsims){
presences = zIn[, i]
covs=data.frame(X=x[, i], X2=rep(0, nsites))
formula=maxnet.formula(presences, covs, classes="l")

tmp <- maxnet(presences, covs, formula, regmult=1)
tmppredict<-predict(tmp, newdata = covs)
zmaxent1[which(tmppredict>=mean(tmppredict)), i]<-1
predROCR<-prediction(tmppredict, presences)
ss<-performance(predROCR, "sens", "spec")
thresh_idx<-which.max(ss@x.values[[1]]+ss@y.values[[1]])
thresh<-ss@alpha.values[[1]][thresh_idx]
zmaxent2[which(tmppredict>=thresh), i]<-1
}

plot(colSums(zmaxent2), df$AOO)
cor(colSums(zmaxent1), df$AOO)
cor(df$AOO_hat, df$AOO)

df$AOO_hat_maxent1<-colSums(zmaxent1)
df$AOO_hat_maxent2<-colSums(zmaxent2)

###Not really the focus, but if interested, checks to see
###that the sites predicted as occupied using the threshold tend to be 
###the sites actually occupied.
meancor1<-rep(NA, nsims)
meancor2<-rep(NA, nsims)

for (s in 1:nsims){
  meancor1[s]<-cor(zmaxent1[, s], z[, s]) ###mean correlation =0.246, -.038 - 0.575
  meancor2[s]<-cor(zmaxent2[, s], z[, s]) ###mean correlation =0.267, -0.074 - 0.568
}

df$AOO_hatlci<-PAOCI[, 1]
df$AOO_hatuci<-PAOCI[, 2]

df$AOO_hat_mean=apply(samps2, 2, mean)

###plotting

pan1a<-ggplot(df, aes(AOO, AOO_hat))+geom_errorbar(ymin=df$AOO_hatlci,ymax=df$AOO_hatuci, col='lightgray')+
  geom_point()+geom_abline(slope=1, intercept=0)+
  ylab(expression(widehat(AOO)~~iSDM(mode, HDI)))+theme_bw()+stat_smooth(method='lm')+
  xlim(0, 200)+ylim(0, 200)

pan1b<-ggplot(df, aes(AOO, AOO_hat_mean))+geom_errorbar(ymin=df$AOO_hatlci2,ymax=df$AOO_hatuci2, col='lightgray')+
  geom_point()+geom_abline(slope=1, intercept=0)+
  ylab(expression(widehat(AOO)~~iSDM(mean, Central)))+theme_bw()+stat_smooth(method='lm')+
  xlim(0, 200)+ylim(0, 200)


pan1c<-ggplot(df[df$AOO_hat_maxent1!=200, ], aes(AOO, AOO_hat_maxent1))+geom_point()+geom_abline(slope=1, intercept=0)+
  ylab(expression(widehat(AOO)~~Maxent--Mean~Threshold))+theme_bw()+
  stat_smooth(data=df[df$AOO_hat_maxent1!=200, ], aes(AOO,AOO_hat_maxent1), method='lm')+
  xlim(0, 200)+ylim(0, 200)

pan1d<-ggplot(df[df$AOO_hat_maxent2!=0, ], aes(AOO, AOO_hat_maxent2))+geom_point()+geom_abline(slope=1, intercept=0)+
  ylab(expression(widehat(AOO)~~Maxent--SSS~Threshold))+theme_bw()+
  stat_smooth(data=df[df$AOO_hat_maxent2!=0, ], aes(AOO, AOO_hat_maxent2), method='lm')+
  xlim(0, 200)+ylim(0, 200)


df$rse_maxent1<-sqrt(((df$AOO - df$AOO_hat_maxent1)^2))
df$percentrse_maxent1<-df$rse_maxent1/df$AOO*100

df$rse_maxent2<-sqrt(((df$AOO - df$AOO_hat_maxent2)^2))
df$percentrse_maxent2<-df$rse_maxent2/df$AOO*100

###more plotting
pan3a<-ggplot(df, aes(Consts$nsamp1, percentrse))+geom_point()+geom_smooth()+theme_bw()+
  xlab("# Sites with Occupancy Surveys")+ylab("Root Squared Error (%)")

pan3b<-ggplot(df, aes(p, percentrse))+geom_point()+geom_smooth()+theme_bw()+
  xlab("p (Detection Probability)")+ylab("Root Squared Error (%)")  

pan3c<-ggplot(df, aes(beta0, percentrse))+geom_point()+geom_smooth()+theme_bw()+
  xlab(expression(beta[0]))+ylab("Root Squared Error (%)")    


ggplot(df, aes(b, percentrse))+geom_point()+geom_smooth(method='gam', method.args = list(family = gaussian(link = "log")))
ggplot(df, aes(beta0, percentrse))+geom_point()+stat_smooth(method='gam', method.args = list(family = Gamma(link = "log")))
ggplot(df, aes(beta1, percentrse))+geom_point()+geom_smooth()


##save this is as desired
#save(df, samps2 ,Dat, Consts, Inits, b, p, beta0, beta1,
#     file='PatRap.RData')
