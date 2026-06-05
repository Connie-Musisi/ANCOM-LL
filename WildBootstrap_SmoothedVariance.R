# Load data
load("G:/My Drive/PhD UHasselt Connie.zip (Unzipped Files)/ANCOM-BC/Results/Results Leyla/simsB_2.2_p.RData")

# Function that identifies the taxa that are not differentially abundant taxa and also removes them from the data 
ref <- function(x){
  if(taxa_are_rows(x))
  {
    librarysize <- rowSums(otu_table(x))
  }else{
    librarysize <- rowSums(otu_table(x))
  }
  librarysize <- as.data.frame(librarysize)
  newtax_table <- cbind(tax_table(x)[,"DE.ind"],librarysize)
  # if(newtax_table[,"DE.ind"]=="FALSE")
  #   {refTaxa1 <- newtax_table[sample(nrow(newtax_table), 5), ]}
  # ifelse(newtax_table[,"DE.ind"]=="FALSE", refTaxa1 <- newtax_table[sample(nrow(newtax_table), 5), ])
  for(i in 1:nrow(newtax_table)) {                     # for-loop & if-statement
    if(newtax_table$DE.ind[i] =="FALSE"){refTaxa <- newtax_table[sample(nrow(newtax_table), 5), ]}
  }
  return(refTaxa)
} 

sim.data = function(physeq){
  ref <- ref(physeq)
  refName <- row.names(ref)
  allTaxa <- taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% refName)]
  refTaxa <- allTaxa[(allTaxa %in% refName)]
  db_red <- prune_taxa(myTaxa, physeq)
  db_refTaxa <- prune_taxa(refTaxa, physeq)
  ref_phy <- as.matrix(apply(db_refTaxa@otu_table,2,median))
  colnames(ref_phy) <- "V1"
  db_otu <- otu_table(cbind(ref_phy,t(db_red@otu_table)), taxa_are_rows = FALSE)
  db_wide <- phyloseq(db_otu, tax_table(db_red@tax_table), sample_data(db_red@sam_data))
  V1 <- matrix(data = c("FALSE", 0), nrow = 1, ncol = 2)
  colnames(V1) <- c("DE.ind","source.ID")
  rownames(V1) <- "V1"
  db_tax <- tax_table(rbind(V1,db_red@tax_table))
  db_wide <- phyloseq(db_otu, db_tax, sample_data(db_red@sam_data))
  
  # Transforming the data to long format 
  if(taxa_are_rows(db_wide))
  {
    db_otu <- data.frame(t(db_wide@otu_table))
  }else{
    db_otu <- data.frame(db_wide@otu_table)
  }
  db_sam <- data.frame(db_wide@sam_data)
  db_long <- gather(db_otu, key="taxon", value="O", factor_key=TRUE)
  db_long$group <- rep(db_sam[,], times=ncol(db_otu))
  db_long$subject <- rep(c(1:nrow(db_otu)),times=ncol(db_otu))
  db <- arrange(db_long, subject)
  n.taxa <- ncol(db_otu)
  db$taxname <- names(db_otu)
  db$taxon <- rep(1:ncol(db_otu))
  db$group<-as.factor(db$group)
  db$subject<-as.factor(db$subject)
  db$taxon<-as.factor(db$taxon)
  return(list(db=db,
              n.taxa=n.taxa))
}

## b. Parameter estimation
est.mu<-function(db,V) {
  db$S<-db$O/V
  mu0<-numeric(length(unique(db$taxon)))
  mu1<-M0<-mu0<-M1<-mu1<-S0<-S1<-mu0
  cnt<-1
  for(i in unique(db$taxon)) {
    S0[cnt]<-sum((db$S)[(db$taxon==i)&(db$group==0)])
    M0[cnt]<-sum((1/V)[(db$taxon==i)&(db$group==0)])
    mu0[cnt]<-log(S0[cnt]/M0[cnt])
    S1[cnt]<-sum((db$S)[(db$taxon==i)&(db$group==1)])
    M1[cnt]<-sum((1/V)[(db$taxon==i)&(db$group==1)])
    mu1[cnt]<-log(S1[cnt]/M1[cnt])
    cnt<-cnt+1
  }
  return(list(
    mu0=mu0[-1]-mu0[1],
    mu1=mu1[-1]-mu1[1],
    delta=c(mu0[1],mu1[1])
  ))
}

## c. Variance estimation
### 1-Wild bootstrap
ks<-function(x,EE,db) {
  u<-runif(length(EE),min=0.1,max=2)
  w<-x[1]*u+x[2]*u^2+x[3]*u^3+x[4]*u^4
  
  r<-mean((sort(w*EE)-sort(EE))^2)
  
  return(r)
}

wb<-function(db,V,est,EE,B=100) {
  # Wild bootstrap for Estimating Equations
  est.delta<-matrix(nrow=B,ncol=2)
  est.mu0<-est.mu1<-matrix(nrow=B,ncol=length(est$mu0))
  ks.opt<-optim(c(0.1,0.1,0.1,0.1),fn=ks,EE=EE,db=db)
  x<-ks.opt$par
  for(i in 1:B) {
    u<-runif(nrow(db),min=0.1,max=2)
    w<-x[1]*u+x[2]*u^2+x[3]*u^3+x[4]*u^4
    V2<-V/w
    est.wb<-est.mu(db,V2)
    est.delta[i,]<-est.wb$delta
    est.mu0[i,]<-est.wb$mu0
    est.mu1[i,]<-est.wb$mu1
  }
  return(list(
    var.delta=diag(var(est.delta)),
    var.mu0=diag(var(est.mu0)),
    var.mu1=diag(var(est.mu1)),
    wild.par=x
  ))
}

### 2-Smoothed variance
est.all<-function(db,V) {
  # with smoothed variance
  db$S<-db$O/V
  mu0<-mu1<-var0<-var1<-cov0<-cov1<-S0b<-S1b<-V0b<-V1b<-C0b<-C1b<-numeric(length(unique(db$taxon)))
  n0<-sum((db$group==0)&(db$taxon==1))
  n1<-sum((db$group==1)&(db$taxon==1))
  s0.1<-(db$S)[(db$taxon==1)&(db$group==0)]
  s1.1<-(db$S)[(db$taxon==1)&(db$group==1)]
  S0.1<-sum(s0.1)
  S1.1<-sum(s1.1)
  cnt<-1
  for(i in unique(db$taxon)) { 
    s0<-(db$S)[(db$taxon==i)&(db$group==0)]
    S0b[cnt]<-S0<-sum(s0)
    V0b[cnt]<-var(s0)
    C0b[cnt]<-cor(s0,s0.1)
    M0<-sum((1/V)[(db$taxon==i)&(db$group==0)])
    mu0[cnt]<-log(S0/M0)
    
    s1<-(db$S)[(db$taxon==i)&(db$group==1)]
    S1b[cnt]<-S1<-sum(s1)
    V1b[cnt]<-var(s1)
    C1b[cnt]<-cor(s1,s1.1)
    M1<-sum((1/V)[(db$taxon==i)&(db$group==1)])
    mu1[cnt]<-log(S1/M1)
    
    cnt<-cnt+1
  }
  
  m0<-glm(V0b~S0b+I(S0b^2),family = poisson())
  v0<-predict(m0, type="response")
  m1<-glm(V1b~S1b+I(S1b^2),family = poisson())
  v1<-predict(m1, type="response")
  
  cnt<-1
  for(i in unique(db$taxon)) {
    sim<-rmvnorm(1000, mean=c(S0.1,S0b[cnt]),
                 sigma=n0*matrix(c(
                   1*var(s0.1)+0*v0[1],
                   C0b[cnt]*sqrt((1*var(s0.1)+0*v0[1])*(0.5*V0b[cnt]+0.5*v0[cnt])),
                   C0b[cnt]*sqrt((1*var(s0.1)+0*v0[1])*(0.5*V0b[cnt]+0.5*v0[cnt])),
                   0.5*V0b[cnt]+0.5*v0[cnt]),
                   ncol=2,byrow = TRUE))
    
    sim<-sim[pmin(sim[,1],sim[,2])>1,]
    varcovar0<-var(log(sim))
    
    sim<-rmvnorm(1000,mean=c(S1.1,S1b[cnt]),
                 sigma=n1*matrix(c(
                   1*var(s1.1)+0*v1[1],
                   C1b[cnt]*sqrt((1*var(s1.1)+0*v1[1])*(0.5*v1[cnt]+0.5*V1b[cnt])),
                   C1b[cnt]*sqrt((1*var(s1.1)+0*v1[1])*(0.5*v1[cnt]+0.5*V1b[cnt])),
                   0.5*v1[cnt]+0.5*V1b[cnt]),
                   ncol=2,byrow = TRUE))
    
    sim<-sim[pmin(sim[,1],sim[,2])>1,]
    varcovar1<-var(log(sim))
    if(cnt>1) {
      var0[cnt]<-varcovar0[1,1]+varcovar0[2,2]-2*varcovar0[1,2]
      var1[cnt]<-varcovar1[1,1]+varcovar1[2,2]-2*varcovar1[1,2]
    }
    if(cnt==1) {
      var0[cnt]<-varcovar0[1,1]
      var1[cnt]<-varcovar1[1,1]
    }
    cnt<-cnt+1
  }
  
  return(list(
    mu0=mu0[-1]-mu0[1],
    mu1=mu1[-1]-mu1[1],
    delta=c(mu0[1],mu1[1]),
    var0=var0,
    var1=var1
  ))
}

predictions<-function(db,mu.hat) {
  # this functions calculates the predictions of db$O based on
  # the parameter estimates in mu.hat
  m.hat<-sapply(as.numeric(db$group),function(x) {
    exp(mu.hat$delta[x])
  })
  m.hat[(db$taxon!=1)&(db$group==0)]<-m.hat[(db$taxon!=1)&(db$group==0)]*
    sapply(as.numeric(db$taxon[(db$taxon!=1)&(db$group==0)]),function(x){
      exp(mu.hat$mu0[x-1])
    })
  m.hat[(db$taxon!=1)&(db$group==1)]<-m.hat[(db$taxon!=1)&(db$group==1)]*
    sapply(as.numeric(db$taxon[(db$taxon!=1)&(db$group==1)]),function(x){
      exp(mu.hat$mu1[x-1])
    })
  return(m.hat)
}

run.scenario <- function(physeq, V.method="none",
                     var.method="wild", B=100){
  tmp<-sim.data(physeq)
  n.taxa<-db<-tmp$n.taxa 
  
  par.estimates<-length(2*n.taxa)
  var.estimates<-length(2*n.taxa)
  p.values<-length(n.taxa-1)
  wild.pars<-length(4)
  
  db<-tmp$db
  if(V.method=="none") {
    V<-rep(1,nrow(db))
  }
  if(var.method=="wild") {
    mu.hat<-est.mu(db,V)
    m.hat<-predictions(db,mu.hat)
    var.wb<-wb(db,V=V,est=mu.hat,EE=db$O-m.hat,B=B)
    var.hat<-c(var.wb$var.delta,var.wb$var.mu0,var.wb$var.mu1)
    wild.pars<-var.wb$wild.par
  }
  if(var.method=="smoothed") {
    mu.hat<-est.all(db,V)
    var.hat<-c(mu.hat$var0[1],mu.hat$var1[1],mu.hat$var0[-1],mu.hat$var1[-1])
  }
  
  par.estimates<-c(mu.hat$delta,mu.hat$mu0,mu.hat$mu1)
  var.estimates<-var.hat
  p.values<-(1-pnorm(abs(mu.hat$mu1-mu.hat$mu0)/
                       sqrt(var.hat[3:(n.taxa+1)]+var.hat[(n.taxa+2):(2*n.taxa)])))*2
  p.adj<-p.adjust(p.values,method="BH")
  
  return(list(est.par=par.estimates,
              est.var=var.estimates,
              p.values=p.values,
              p.adj=p.adj))
}

ancom_LL <- function(physeq){
results_wb <- run.scenario(physeq, V.method = "none", var.method = "wild",B=100)
results_sv <- run.scenario(physeq, V.method = "none", var.method = "smoothed",B=100)
return(list(results_bootstrap = results_wb,
            results_smoothed = results_wb))
}







