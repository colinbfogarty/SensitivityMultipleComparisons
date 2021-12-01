################################
#This script reproduces the results of the paper. There are 
#three sections. We first reproduce the inference performed
#and the sensitivity analysis for 1-naphthol and 2-naphthol.
#We then recreate the match described in the paper. Finally,
#we give code for conducting the simulation studies of section 5.
################################

###You will need to download Gurobi 9.0 or higher and install its R package. Instructions are available at www.gurobi.com; Gurobi is free for academics.
#You will also need the following R packages:
# slam, Matrix, MASS, optmatch, as well as their dependencies

####First, download the code and place everything in the working
#directory 


#load the data set and the required functions
dat = read.csv("naphthalene.csv")
source("multiCompareFunctions.R")
#treatment vector
treatment = dat[,1]
#indices resulting from full match. Individuals in the same matched set have the same index number
index = dat[,2]
#the  responses, log 1-naphthol and 2-naphthol concentrations
lnaph = log(dat[,3:4])
#the covariates
covariates = dat[,-(1:4)]

###########################
#1) Inference and sensitivity analysis
######################


#Get HL estimates for multiplicative effect (additive on log scale)
HL = rep(0,2)
for(aa in 1:2)
{
outcome = lnaph[,aa]
taugrid=seq(1, 3, .01);
pvalgrid=rep(0,length(taugrid));
for(i in 1:length(taugrid)){
adjusted.response=outcome-treatment*taugrid[i];
pvalgrid[i]=alignedranktest(adjusted.response,index,treatment);
}
HL[aa] = taugrid[which.max(pvalgrid)]	
}

###to get the multiplicative factors in the paper, take exp{} of answer
estimates = exp(HL)

#Now for confidence intervals:

cibound.func=function(treatment.effect,outcome,matchedset,treated,alternative="greater"){
adjusted.outcome=outcome-treatment.effect*treated;
alignedranktest(adjusted.outcome,matchedset,treated, alternative)-.025;
}


CI = matrix(0,2,2)
#invert a series of tests to get the confidence intervals on log scale, then exponentiate
for(aa in 1:2)
{
# Confidence bounds
lci.bound=uniroot(cibound.func,c(-2,5),lnaph[,aa],index, treatment, alternative="greater")
uci.bound=uniroot(cibound.func,c(-2,5),lnaph[,aa],index, treatment, alternative="less")
l = lci.bound$root
u = uci.bound$root
CI[aa,] = exp(c(l, u))
}

CI



################
#Now, the sensitivity analysis. We use the function multipleComparisons for all of these. It requires the gurobi optimization suite, which is freely available for academics
###############

#First, define the aligned ranks:
PO = matrix(0, nrow(lnaph), ncol(lnaph))
nostratum = length(unique(index))
for(i in 1:ncol(PO))
{
	alignedval = rep(0, nostratum)
	out = lnaph[,i]
	for(j in 1:nostratum)
	{
		ind = which(index == j)
		alignedval[ind] = out[ind] - mean(out[ind])
	}
	ranks = rank(alignedval)
	PO[,i] = ranks	
}
##The uncorrected maximal Gammas
G105 = uniroot(multipleComparisonsRoot, c(1, 9), index, PO[,1], treatment, alpha = 0.05, alternative = "G")$root
G205 = uniroot(multipleComparisonsRoot, c(1, 9), index, PO[,2], treatment, alpha = 0.05, alternative = "G")$root
#Now correcting
G1025 = uniroot(multipleComparisonsRoot, c(1, 9), index, PO[,1], treatment, alpha = 0.025, alternative = "G")$root
G2025 = uniroot(multipleComparisonsRoot, c(1, 9), index, PO[,2], treatment, alpha = 0.025, alternative = "G")$root

#Now finding the maximal Gamma for overall null through the minimax procedure

G12 = uniroot(multipleComparisonsRoot, c(1, 11), index, PO, treatment, alpha = 0.05, alternative = "G")$root


#Now, let's show the results in Table 3


res1 = multipleComparisons(index, PO[,1], treatment, alternative = "G", Gamma = 10)
res2 = multipleComparisons(index, PO[,2], treatment, alternative = "G", Gamma = 10)
res12 = multipleComparisons(index, PO, treatment, alternative = "G", Gamma = 10)

#u gives the worst case confounder, and rho gives the conditional probabilities under the worst case u at the Gamma being tested rho_{ij} = (exp(Gamma*u_{ij})/sum_{i=1}^{n_i}(exp(Gamma*u_{ij}))
rho1pair = res1$rho[index==348]
u1pair = res1$u[index==348]
rho2pair = res2$rho[index==348]
u2pair = res2$u[index==348]
POpair = PO[index==348,]
rho12pair = res12$rho[index==348]
u12pair = res12$u[index==348]

#worst case U's
u1pair
u2pair
u12pair

#worst expectations - individ
sum(rho1pair*POpair[,1])
sum(rho2pair*POpair[,2])
#worst expectations - minimax
sum(rho12pair*POpair[,1])
sum(rho12pair*POpair[,2])



################
#2) The following code returns the indices found in column 2 of the data set naphthalene.csv by conducting a full match. It also allows the user to create Figure 1 of the manuscript
################



###########
library(optmatch)
#Get the covaraites from the data set
covariates = dat[,-(1:4)]
X = covariates

names.vars = colnames(X)
dat.new = as.data.frame(X)
dat.new.miss = dat.new
names.miss = c()
missing.indicators = c()
####Make the indicators for missing values, and fill in the overall mean for the NA's. This does not affect the propensity score fits for the variables themselves, as discussed in the references in this part of the manuscript
MISS1 = matrix(0, nrow(dat.new), ncol(dat.new))
for(i in 1:ncol(dat.new))
{
  MISS1[,i] = is.na(dat.new[,i])*1
  if(any(is.na(dat.new[,i])))
  {
    miss = is.na(dat.new[,i])*1
    temp = dat.new[,i]
    
    temp[is.na(dat.new[,i])] = mean(dat.new[,i], na.rm = T)
    dat.new.miss[,i] = temp
    missing.indicators = cbind(missing.indicators, miss)
    names.miss = c(names.miss, paste(names.vars[i], "MISS", sep = ""))
      }
  
  
}
#Since we are treating the missingness indicators as covariates, we also include a block of ``missingness indicators" for the missingness indicators (all zeroes of course). This gets used down the line for calculating standardized differences
MISS1a = cbind(MISS1, matrix(0, nrow(MISS1), ncol(missing.indicators)))
treat1 = 1*treatment
colnames(missing.indicators) = names.miss
dat.temp= cbind(dat.new.miss, missing.indicators)#[!row.remove,]
dat.final = lm(treat1~., data = dat.temp, x=T)$x[,-1]
dat.final = as.data.frame(dat.final)
model = glm(treat1~., data = dat.final, family = binomial, x = T)
#remove the race category "other" to get the final covariate matrix for matching (for singularity reasons). Also remove the intercept
 X = model$x[,-c(1,23)]
 
 write.csv(cbind(treat1, X), file = "Xmatrix.csv" )
 
 logit.propscore1 = predict(model)
 logit.propscore = logit.propscore1
 propscore = exp(logit.propscore1)/(1+exp(logit.propscore1))


subject.index=seq(1,length(treat1),1)
rownames(X) = subject.index
#Make the rank based mahalanobis distances
distmat = smahal(treat1, X)
#Add the propensity score caliper
distmat2 = addcaliper(distmat, treat1, logit.propscore, calipersd = .08)
rownames(distmat2)=subject.index[treat1==1]
colnames(distmat2)=subject.index[treat1==0]
treated=(treat1==1)
#conduct the full match
matchvec=fullmatch(distmat2,min.controls = 1/9, max.controls = 9)

############This code strips the resulting match vector, defines indices for the matched sets, and associates the proper indices with the corresponding row in the original dataset
matchedset.index=substr(matchvec,start=3,stop=10)
matchedset.index.numeric=as.numeric(matchedset.index)
subjects.match.order=as.numeric(names(matchvec)) 

# Create a numeric variable for which stratum each unit belongs to
# 0 denotes that the unit was not matched
stratum.short=substr(matchvec,start=3,stop=10);
stratum.numeric=as.numeric(stratum.short);
# Reassign numbers to each stratum that go from 1,..., no. of stratum
sort.unique.stratum=sort(unique(stratum.numeric));
stratum.myindex.matchvecorder=rep(0,length(stratum.numeric));
for(i in 1:length(sort.unique.stratum)){
  stratum.myindex.matchvecorder[stratum.numeric==sort.unique.stratum[i]]=i;
}
##makes stratum.myindex in the proper order

stratum.myindex=rep(0,length(stratum.myindex.matchvecorder))
stratum.myindex[subjects.match.order]=stratum.myindex.matchvecorder
####check that this aligns with index used for inference in the first bit of this code
all(stratum.myindex==dat[,2])

#add "other" back into the X matrix (it got removed to avoid singularity issues)
Xfinal = cbind(X[,1:21], covariates[,22], X[,22:31])
MISSfinal = MISS1a

std.diff.before=rep(0,ncol(Xfinal));
std.diff.after=rep(0,ncol(Xfinal));
names(std.diff.before)=colnames(Xfinal);
names(std.diff.after)=colnames(Xfinal);
#define the standardized differences
for(i in 1:ncol(Xfinal)){ temp.stand.diff=standardized.diff.func(Xfinal[,i],treatment,stratum.myindex, MISSfinal[,i]);
  std.diff.before[i]=temp.stand.diff$std.diff.before.matching;
  std.diff.after[i]=temp.stand.diff$std.diff.after.matching;
  print(i)
}
#Check the max absolute standardized difference. It is quite low (below 0.1), so we declare the match successful in reducing overt biases
max(abs(std.diff.after))
covname = c("Age", "Gender", "Poverty:Income Ratio (PIR)", "Education", "Weight", "Height", "Recreational Exercise", "Regular Walking", "Moderate Workplace Exertion", "Drinks per Day", "Any Drinks this Year?", "Urinary Creatine", "Mineral Dust Exposure", "Organic Dust Exposure", "Exhaust Fume Exposure", "Other Fume Exposure", "Charred Meat Consumption", "Mexican American", "Other Hispanic", "White", "Black", "Other Race", "PIR MISS",  "Weight MISS",    "Height MISS",      "Drinks per Day MISS", "Any Drinks this Year? MISS",     "Mineral Dust MISS",     "Organic Dust MISS",  "Exhaust Fumes MISS",  "Other Fumes MISS", "Charred Meats MISS")


#####Make the balance plot (Figure 1 of the manuscript)
plotBalancesign(std.diff.before, std.diff.after, covname, titleOfPlot = "Standardized Differences Before and After Matching" )



##############
# 3) Simulation Study
################
 
#1: Five outcomes, look at overall null
library(MASS)
nsim = 10000
rbonf = rep(0, nsim)
RBONF = rep(0, nsim)
RMULT = rep(0, nsim)
rmulti = rep(0, nsim)
timecode = rep(0, nsim)
K=5
mu4 = c(0.3, 0, 0, 0, 0)
mu3 = c(0.3, 0.3, 0, 0, 0)
mu1 = c(.25, .25, .25, .25, .25)
mu2 = c(.25, .25, .25, .25, 0)
Sig1 = diag(5)
Sig2 = matrix(.5, 5, 5)
diag(Sig2) = 1
temp = rep(0, 4)

mlist = list(mu1, mu1, mu2, mu2, mu3, mu3, mu4, mu4, mu1, mu1, mu2, mu2, mu3, mu3, mu4, mu4, mu1, mu1, mu2, mu2, mu3, mu3, mu4, mu4)
Slist = list(Sig1, Sig2,Sig1, Sig2,Sig1, Sig2,Sig1, Sig2,Sig1, Sig2,Sig1, Sig2,Sig1, Sig2,Sig1, Sig2,Sig1, Sig2, Sig1, Sig2, Sig1, Sig2, Sig1, Sig2)
Gvec = c(rep(1.25, 8), rep(1.5,8), rep(1.75, 8))
RBlist = list(RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF, RBONF,RBONF,RBONF,RBONF,RBONF,RBONF)
RMlist = RBlist
for(aa in 1:24)
{
  RBONF = rep(0, nsim)
  RMULT = rep(0, nsim)
  Gam = Gvec[aa]

  I = 250
  Sig = Slist[[aa]]
   
  mu = mlist[[aa]]
  
  for(mm in 1:nsim)
  {
  	#simulate pairs
    mat = mvrnorm(I, mu, Sig)
    PO = rbind(mat, -mat)   
    treatment = rep(c(T,F), each = I)
    index = rep(1:I, 2)
    OUT = PO
    QQ = OUT
    DIFF = QQ
    ns = rep(2, I)
    nostratum = I
   #Define Huber's M statistic
      for(j in 1:nostratum)
      {
        ind = which(index==j)
        
        DIFF[ind[1],] = OUT[ind[1],] - OUT[ind[2],]
        DIFF[ind[2],] = -DIFF[ind[1],]
      }
      sk = apply(abs(DIFF), 2, median)
      SK = matrix(sk, 2, K, byrow = T)
      for(j in 1:nostratum)
      {
        ind = which(index==j)
        for(jj in 1:length(ind))
        {
          sg = sign(DIFF[ind,])
          ad = abs(DIFF[ind,]/SK)
          ad[ad > 2.5] = 2.5
          QQ[ind,] = sg*ad  
        }
        
      }
      
    PO = QQ
    MG = rep(0, K)
    #Figure out if we can reject the overall null by separately conducting sensitivity analysis
    for(k in 1:K)
    {
    
      MG[k] = multipleComparisons(index, PO[,k], treatment, alpha = .05/K, Gamma = Gam)$Reject[1]
    }
    if(any(MG ==1))
    {
      RBONF[mm] = 1
    }
    #While timing code, figure out if we can reject the overall null by considering all outcomes simultaneously
    ptm = proc.time()
    rm = multipleComparisons(index, PO, treatment, Gamma = Gam)$Reject[1]
    if(rm == 1)
    {
      RMULT[mm] = 1
     
    }
    timecode[mm] = (proc.time() - ptm)[3]
    print(c(aa,mm))
  }
  RBlist[[aa]] = mean(RBONF)
  RMlist[[aa]] = mean(RMULT)
}


###RBlist and RMlist Give the proportions of rejections for both procedures in the order given in Table 1

############Simulation 2: Individual Tests, K=3 outcomes

nsim = 10000
rbonf = rep(0, nsim)
RBONF = matrix(0, nsim, 4)
RMULT = matrix(0, nsim, 4)
rmulti = rep(0, nsim)
timecode = rep(0, nsim)
K = 3

mu4 = c(0.15, .2, .35)
mu3 = c(0.2, .25, .35)
mu2 = c(.25, .3, .35)
mu1 = c(.2, .225, .25)
Sig1 = diag(3)
Sig2 = matrix(.5, 3, 3)
diag(Sig2) = 1

mlist = list(mu1, mu1, mu2, mu2, mu3, mu3, mu4, mu4, mu1, mu1, mu2, mu2, mu3, mu3, mu4, mu4, mu1, mu1, mu2, mu2, mu3, mu3, mu4, mu4)
Slist = list(Sig1, Sig2,Sig1, Sig2,Sig1, Sig2,Sig1, Sig2,Sig1, Sig2,Sig1, Sig2,Sig1, Sig2,Sig1, Sig2,Sig1, Sig2, Sig1, Sig2, Sig1, Sig2, Sig1, Sig2)
Gvec = c(rep(1.25, 8), rep(1.375,8), rep(1.5, 8))
RBlist = list(RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF,RBONF, RBONF,RBONF,RBONF,RBONF,RBONF,RBONF)
RMlist = RBlist
for(aa in 1:24)
{
  RBONF = matrix(0, nsim, 4)
  RMULT = matrix(0, nsim, 4)
  Gam = Gvec[aa]
  I = 250
  Sig = Slist[[aa]]
  mu = mlist[[aa]]
  for(mm in 1:nsim)
  {
  mat = mvrnorm(I, mu, Sig)
  PO = rbind(mat, -mat)
  
  treatment = rep(c(T,F), each = I)
  index = rep(1:I, 2)
  OUT = PO

    QQ = OUT
    DIFF = QQ
    ns = rep(2, I)
    nostratum = I
   #Define Huber's M statistic
      for(j in 1:nostratum)
      {
        ind = which(index==j)
        
        DIFF[ind[1],] = OUT[ind[1],] - OUT[ind[2],]
        DIFF[ind[2],] = -DIFF[ind[1],]
      }
      sk = apply(abs(DIFF), 2, median)
      SK = matrix(sk, 2, K, byrow = T)
      for(j in 1:nostratum)
      {
        ind = which(index==j)
        for(jj in 1:length(ind))
        {
          sg = sign(DIFF[ind,])
          ad = abs(DIFF[ind,]/SK)
          ad[ad > 2.5] = 2.5
          QQ[ind,] = sg*ad  
        }
        
      }
      
    PO = QQ

  MG = rep(0, K)
  ####Here we use a different approach for a sensitivity analysis, described in a submitted paper which is included as supplementary material. It allows us to get p-values for the sensitivity analysis without using asymptotic separability. Equivalently, one could have coded this based on our function "multipleComparison" as a closed testing procedure.
  for(k in 1:K)
  {
    MG[k] = sensitivity(index, PO[,k], treatment, alpha = .05/K, alternative = "two.sided", Gamma.vec = Gam, calculate.pval = T)$pval
  }
  MGsort = sort(MG)
  ord = order(MG)
  if(MGsort[1] < .05/3)
  {
    RBONF[mm,ord[1]] = 1
    RBONF[mm, 4] = 1
    if(MGsort[2] < .05/2)
    {
      RBONF[mm,ord[2]] = 1
      if(MGsort[3] < .05)
      {
        RBONF[mm,ord[3]] = 1
        
      }
      
    }
  }
  ptm = proc.time()
  rm = multipleComparisons(index, PO, treatment, Gamma = Gam)$Reject[1]
  if(rm == 1)
  {
    RMULT[mm,4] = 1
    rm12 = multipleComparisons(index, PO[,-3], treatment, Gamma = Gam)$Reject[1]
    rm23 = multipleComparisons(index, PO[,-1], treatment, Gamma = Gam)$Reject[1]
    rm13 = multipleComparisons(index, PO[,-2], treatment, Gamma = Gam)$Reject[1]
    if(rm12 == 1 & rm13 == 1)
    {
      RMULT[mm,1] = (MG[1] < .05)
    }
    if(rm12 == 1 & rm23 == 1)
    {
      RMULT[mm,2] = (MG[2] < .05)
    }
    if(rm13 == 1 & rm23 == 1)
    {
      RMULT[mm,3] = (MG[3] < .05)
    }
    
  }
  timecode[mm] = (proc.time() - ptm)[3]
  print(c(aa,mm))
}
RBlist[[aa]] = colMeans(RBONF)
RMlist[[aa]] = colMeans(RMULT)
}

###RBlist and RMlist contain the estimated power for each situation, in the order of table 2




