# ======== ANCOM LOG LINEAR ========
library(parallel)
library(doSNOW)

sim.data<-function(n=100,n.taxa=250,
                   FC=c(1,3,2), phi=0.5,
                   SF0=0.8,SF1=1.2,
                   phi.var=1) {
  # NB data is simulated
  # FC is vector with fold changes for the n.taxa
  # phi is the overdispersion factor. A very small value of phi
  # approximates the Poisson distribution
  # phi.var makes phi randomly vary between subjects,
  # between phi/phi.var and phi*phi.var
  # SF0 and SF1 are the sampling fractions in group 1 and 2
  db<-data.frame(group=rep(c(0,1),c(n,n)*n.taxa),
                 subject=rep(1:(2*n),rep(n.taxa,2*n)),
                 taxon=rep(1:n.taxa,2*n))
  db$group<-as.factor(db$group)
  db$subject<-as.factor(db$subject)
  db$taxon<-as.factor(db$taxon) 
  mu0<-seq(1,100,length.out=n.taxa)
  FC<-c(FC,rep(1,n.taxa-length(FC)))
  mu1<-mu0*FC
  SF0k<-rep(SF0*exp(rnorm(n,sd=0.2)),rep(n.taxa,n))
  SF1k<-rep(SF1*exp(rnorm(n,sd=0.2)),rep(n.taxa,n))
  phi2<-runif(n.taxa,min=phi/phi.var,max=phi*phi.var)
  db$O[db$group==0]<-rnbinom(n*n.taxa,
                             mu=rep(mu0,n)*SF0k,
                             size=1/rep(phi2,n))
  db$O[db$group==1]<-rnbinom(n*n.taxa,
                             mu=rep(mu1,n)*SF1k,
                             size=1/rep(phi2,n))
  return(list(
    db=db,
    mu0=mu0,
    mu1=mu1))
}

#data1=sim.data(n=100,n.taxa=250,FC=c(1,3,2), phi=0.5,SF0=0.8,SF1=1.2,phi.var=1)
#db=data1$db

#### Parameter estimation ####

est.par2<-function(db) {
  D0<-numeric(length(unique(db$subject[db$group==0])))
  D1<-numeric(length(unique(db$subject[db$group==1])))
  mu1<-mu0<-rep(0,length(unique(db$taxon)))
  cnt<-1
  for(i in unique(db$subject[db$group==1])) { 
    D1[cnt]<-log(sum(db$O[(db$subject==i)&(db$group==1)]*exp(mu1)))-
      log(sum(exp(2*mu1)))
    cnt<-cnt+1
  }
  cnt<-1
  for(i in unique(db$subject[db$group==0])) { 
    D0[cnt]<-log(sum(db$O[(db$subject==i)&(db$group==0)]*exp(mu0)))-
      log(sum(exp(2*mu0)))
    cnt<-cnt+1
  }
  #D0<-D0[(D0!=-Inf)&(!is.na(D0))]
  #D1<-D1[(D1!=-Inf)&(!is.na(D1))]
  cnt<-1
  for(j in unique(db$taxon)) {
    mu1[cnt]<-log(sum(db$O[(db$taxon==j)&(db$group==1)]*exp(D1)))-
      log(sum(exp(2*D1)))
    mu0[cnt]<-log(sum(db$O[(db$taxon==j)&(db$group==0)]*exp(D0)))-
      log(sum(exp(2*D0)))
    cnt<-cnt+1
  }  
  cnt<-1
  for(i in unique(db$subject[db$group==1])) { 
    D1[cnt]<-log(sum(db$O[(db$subject==i)&(db$group==1)]*exp(mu1)))-
      log(sum(exp(2*mu1)))
    cnt<-cnt+1
  }
  cnt<-1
  for(i in unique(db$subject[db$group==0])) { 
    D0[cnt]<-log(sum(db$O[(db$subject==i)&(db$group==0)]*exp(mu0)))-
      log(sum(exp(2*mu0)))
    cnt<-cnt+1
  }
  
  return(list(
    mu0=mu0,mu1=mu1,D0=D0,D1=D1,
    all=c(mu0,mu1,D0,D1)
  ))
}


#### Getting reference taxa ####
# from text
ref = function(db){
  parest <- est.par2(db)
  mu1 <- parest$mu1
  mu0 <- parest$mu0
  mu_diff <- mu1-mu0
  # step 1
  #V_j = c()
  #mu_diff <- c()
  V_j = vector(mode="numeric", length=length(mu_diff))
  for (j in 1:length(mu_diff)) {
    delta_j = mu_diff[j]
    gamma_hat = delta_j - mu_diff[-j]
    #V_j = c(V_j, var(gamma_hat)) #gets slower with larger data
    V_j[j] = var(gamma_hat)
  }
  
  # step 2
  q = quantile(V_j, probs = 0.10)
  
  # step 3
  J = c(1:length(V_j))[which(V_j < q)]
  
  # step 4
  mu1_J = mu1[J]
  mu0_J = mu0[J]
  
  gamma_hat1 = mean(mu1_J)
  gamma_hat0 = mean(mu0_J)
  gamma_hatdiff = gamma_hat1 - gamma_hat0
  
  return(list(
    V_j=V_j,
    mu0=parest$mu0,
    mu1=parest$mu1,
    SF0=parest$D0,
    SF1=parest$D1,
    mu_diff=mu_diff,
    gamma_hatdiff=gamma_hatdiff
  ))
}

ANCOMLL = function(db, B=1000){
  
  refs = ref(db)
  V_j=refs$V_j
  
  #### Variance ####
  n0 = nrow(subset(db, group==0))/length(V_j)
  n1 = nrow(subset(db, group==1))/length(V_j)
  
  boot_results=matrix(nrow=B, ncol=length(V_j))
  for (b in 1:B) {
    set.seed(b)
    subject0 = sample(x=(1:n0), size=n0, replace = TRUE)
    subject1 = sample(x=((n0+1):(n0+n1)), size=n1, replace = TRUE)
    db_boot = data.frame(matrix(ncol=ncol(db), nrow = 0))
    colnames(db_boot) = colnames(db)
    for (i in 1:length(subject0)){
      sub = subset(db, (group==0)&(subject==subject0[i]))
      sub$subject=i
      db_boot = rbind(db_boot, sub)
    }
    for (i in 1:length(subject1)){
      sub = subset(db, (group==1)&(subject==subject1[i]))
      sub$subject=n0+i
      db_boot = rbind(db_boot, sub)         
    }
    
    refs_boot = ref(db_boot)
    boot_results[b,] = refs_boot$mu_diff - refs_boot$gamma_hatdiff
    #if(b%%100==0)print(b)
    print(b)
  }
  
  var_j = apply(boot_results, 2, var)
  
  teststat_j = (refs$mu_diff-refs$gamma_hatdiff)/sqrt(var_j)
  
  pvalues_j = 2*(1 - pnorm(abs(teststat_j)))
  
  padj_j = p.adjust(pvalues_j, method = "BH")
  
  #which(padj_j<0.05)
  
  return(list(
    parameter_estimates=list(mu0=refs$mu0,
                             mu1=refs$mu1,
                             SF0=refs$D0,
                             SF1=refs$D1),
    mu_diff=refs$mu_diff,
    gamma_hatdiff=refs$gamma_hatdiff,
    variance_estimates=var_j,
    test_statistic=teststat_j,
    unadjusted_p=pvalues_j,
    adjusted_p=padj_j
  ))
}


doSIM <- function(i){
  set.seed(i)
  mib_data <- sim.data(n=100,n.taxa=250,
                       FC=c(1,3,2), phi=0.5,
                       SF0=0.8,SF1=1.2,
                       phi.var=1)
  res = ANCOMLL(mib_data$db, B=1000)
  
  return(list(mu0=res$parameter_estimates$mu0,
              mu1=res$parameter_estimates$mu1,
            adjusted_p=res$adjusted_p))
}

#detectCores() # how many cores does the pc have, 8

cl <- makeCluster(detectCores()-1)
registerDoSNOW(cl)
iterations <- 100
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
output2 <- foreach(i=1:iterations,.packages=c("MASS"), 
                  .options.snow=opts, .combine = rbind, .verbose=T) %dopar% {
                    result = doSIM(i)
                    mu0 = result$mu0 
                    mu1=result$mu1
                    adjusted_p=result$adjusted_p
              return(c(mu0,mu1,adjusted_p))
                  }
save(output2, file= '~/output2.RData')

close(pb)
stopCluster(cl)
save(output, file= '~/output.RData')

pvalues<-output2[,501:750]

FDR = function(x){
y=which(x<0.05)
  return(length(which(!(y%in%c(2,3))))/length(y))
}

Sensitivity = function(x){
  y=which(x<0.05)
  return(length(which(y%in%c(2,3)))/length(c(2,3)))
}

result_FDR<- apply(pvalues, 1, FDR)
mean(result_FDR)

result_sens<- apply(pvalues, 1, Sensitivity)
mean(result_sens)
