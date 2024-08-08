library("nleqslv")
library("MASS")

#-----------------------
# Data feeding
#----------------------
data0 = read.csv("ada_gusek_plcbo_noBMI.csv")
(head(data0))
# STUDYID = {1, 2, 3}; L = [AGE : int, SEX : (1,2), PASI_BLN : float, PASI_OTC : float]; EXTRT = {PLacebo : 3, Guselkumab : 2, Adalimumab : 1}
# GUSEL vs ADA
data0 = data0[data0$EXTRT!=3]
data0$EXTRT[data0$EXTRT==1] = 0
data0$EXTRT[data0$EXTRT==2] = 1
data0$SEX = data0$SEX -1
# Calculate PASI 90
data0$y = (data0$PASI_OTC-data0$PASI_BLN)/data0$PASI_BLN*100
data0$y_f = -1
data0$y_f[data0$y>=(-90)] = 0
data0$y_f[data0$y<(-90)] = 1

#-----------------------
# Data description
#----------------------
for (s in unique(data0$STUDYID)) {
  for (ttm in unique(data0$EXTRT)) {
    df_proc = data0[data0STUDYID==s & data0$EXTRT==ttm,]
    print(paste("Study",s,"Treatment", ttm,rep=" "))
    print(paste("Age - Mean (SD)", mean(df_proc$AGE), sd(df_proc$AGE), rep=" "))
    print(paste("PASI week0 - Mean (SD)", mean(df_proc$PASI_BLN), sd(df_proc$PASI_BLN), rep=" "))
    print(paste("PASI week16 - Mean (SD)", mean(df_proc$PASI_OTC), sd(df_proc$PASI_OTC), rep=" "))
    print("Summary")
    print(summary(df_proc))
  }
}

#---------------
# Data final
#----------------


expit = function(x) exp(x)/(1+exp(x))


data = data.frame(int = 1, 
                  x = data0$EXTRT, 
                  l1 = data0$AGE, 
                  l2 = data0$PASI_BLN,
                  l3 = data0$SEX,
                  y = data0$y_f, 
                  s = data0$STUDYID)

#data = na.omit(data[which(data$s %in% (1:4)),])
#prop.table(table(data$y,data$s),margin = 2)
#boxplot(data$l3~data$s)
#table(data$s)
#------------------
# Data preparation 
#------------------
n.cov = 3
n = dim(data)[1]
K = 3
Z = 1

# list of mean of L in the aggregated trials
mu.ag = lapply(1:Z, function(ag) colMeans(data[which(data$s == ag),paste0(rep("l",n.cov),1:n.cov)]))
# list of variance of L in the aggregated trials
v.ag = lapply(1:Z, function(ag) sapply(data[which(data$s == ag),paste0(rep("l",n.cov),1:n.cov)], var))
# mean of y in the treatment group of aggregated trials
y1.ag = sapply(1:Z, function(ag) mean(data$y[which(data$x == 1 & data$s == ag)]))
# mean of y in the treatment group of aggregated trials
y0.ag = sapply(1:Z, function(ag) mean(data$y[which(data$x == 0 & data$s == ag)]))
# sample size of aggregated trials
n.ag = sapply(1:Z, function(ag) sum(data$s == ag))
# number of treated patients in aggregated trials
n.xag = sapply(1:Z, function(ag) sum(data$s == ag & data$x == 1))

# P(S = 1, ..., K): proportion of patients in each trial
ps = c(n.ag/n, sapply((Z+1):K, function(i) sum(data$s == i)/n))
# P(X = 1|S): randomization ratio in all trials
px = c(n.xag/n.ag, sapply((Z+1):K, function(i) mean(data$x[which(data$s == i)])))
# N(S=j, X=x): number of treated and control patients in each aggregated trials
nj = lapply(1:Z, function(i) c(n.ag[i]-n.xag[i], n.xag[i]))
# E(L|x,j); E(Y|x,j): mean of L and of Y in each treatment group of aggregated trials
mu = lapply(1:Z, function(i) sapply(0:1, function(ttm) 
  colMeans(data[which(data$s == i & data$x == ttm),c(paste0(rep("l",n.cov),1:n.cov),"y")]))) # mu[[j]][L,ttm] = E(L|X=ttm, S = j)

# E(L^2|x,j); E(Y^2|x,j): mean of L^2 and Y^2 in each treatment group of aggregated trials 
mu.q = lapply(1:Z, function(i) 
  sapply(0:1, function(ttm){
    dataq = data[,3:(n.cov+2)]^2
    colnames(dataq) = paste0(rep("l",n.cov),1:n.cov,rep("q",n.cov))
    data = cbind(data,dataq)
    colMeans(data[which(data$s == i & data$x == ttm),c(paste0(rep("l",n.cov),1:n.cov,rep("q",n.cov)),"y")])})) # mu.q[[j]][L,ttm] = E(L^2|X=ttm, S = j)

# trials with individual-leve data: int = 1; x; L = (l1, l2, ...); y; s = Z+1, ..., K
ori = data
data = data[which(data$s>Z),]
head(data)

# some parameters
n.int = c(3) # index of covariates that have interaction with the treatment. Rule: these covariates are always put at the end.
l.name = c("x", paste0(rep("l",n.cov),1:n.cov), paste0(rep("x:l",length(n.int)),n.int))
lp = paste(l.name, collapse = ' + ')

phi.all <- sapply(1:K,function(i)
  coef(glm(as.formula(paste("y~", lp)), 
           data = ori[which(ori$s == i),], family = binomial('logit'))))


# Step 1
phi.c <- sapply((Z+1):K,function(i)
  coef(glm(as.formula(paste("y~", lp)), 
           data = data[which(data$s == i),], family = binomial('logit'))))

# Step 2
beta.func = function(b, AG, S, data){
  l.cov = paste0(rep("l",n.cov),1:n.cov)
  L = as.matrix(cbind(int = 1, data[,l.cov]))
  st = (data$s == S)*exp(L%*% b)
  indi = sapply(1:ncol(L), function(i) st*L[,i])
  part1 = colMeans(indi)*(dim(data)[1])/n
  part2 = c(1,mu.ag[[AG]])*(n.ag[AG]/n)
  return(part1-part2)
}

bstart = rep(0, n.cov + 1)
b.est = lapply(1:Z, function(ag) sapply((Z+1):K, 
                                        function(st) nleqslv(bstart, beta.func, 
                                                             AG = ag, S = st, data = data)$x))
# Step 3
l.cov = paste0(rep("l",n.cov),1:n.cov)
cova = as.matrix(cbind(int = 1, data[,l.cov]))
w <- lapply(1:Z, function(ag) sapply(1:(K-Z), 
                                     function(i) as.vector(exp(cova %*% b.est[[ag]][,i]))))

k1 <- lapply(1:Z, function(ag) sapply((Z+1):K, function(i) 
  (px[ag]/px[i])))

k0 <- lapply(1:Z, function(ag) sapply((Z+1):K, function(i) 
  (1-px[ag])/(1-px[i])))

l.int = paste0(rep("l",length(n.int)),n.int)
data.int = sapply(1:length(l.int), function(i) data$x * data[,l.int[i]])
colnames(data.int) = paste0(rep("xl", length(n.int)),n.int)

data = cbind(data,data.int)
cova1 = as.matrix(data[,c(paste0(rep("l",n.cov),1:n.cov), colnames(data.int))])

est_eq <- function(p, ag, ipd, data = data) {
  
  # System for S1.1
  oc = as.vector(expit(p[1] + p[2]*data$x + cova1 %*% phi.c[3:length(phi.c[,1]), ipd]))
  
  # Using data of trial S = ipd
  eq <- sapply(0:1, 
               function(ttm) (sum((data$x==ttm) * (data$s == (ipd+Z)) * oc * w[[ag]][,ipd])/n) * (ttm*k1[[ag]][ipd] + (1-ttm)*k0[[ag]][ipd]) - 
                 (y1.ag[ag]*px[ag]*ttm + y0.ag[ag]*(1-px[ag])*(1-ttm))*ps[ag])
  
  return(eq)
}

phi.ag.e = lapply(1:Z, function(ag) 
  sapply(1:(K-Z), # list ag-th: study S=ag with aggregated data
         function(i) nleqslv(rep(0,2), 
                             fn = est_eq, 
                             ag = ag, 
                             ipd = i, 
                             data = data, 
                             method = 'Broyden')$x)) 
# column ith: estimates by using study S=i+Z

#------------------------
# Pseudo-data generation
#------------------------

new = data

# ----------- Estimate E(LaLb|x,j) --------------

# A function to estimate E(u[,1]*u[,2]|x,j): 
cova.func = function(u){
  mu.cov = lapply(1:Z, function(j){
  sapply(1:(K-Z), function(k){
    sapply(0:1, function(ttm){
      (sum((new$x == ttm) * (new$s == (k+Z)) * new[,u[1]] * new[,u[2]] * exp(cova %*% b.est[[j]][,k]))/n)/
        ((px[k+Z]*(ttm == 1) + (1-px[k+Z])*(ttm == 0)) * ps[j])
    })
  })})
  mu.cov = lapply(1:Z, function(j) rowMeans(mu.cov[[j]])) 
  mu.cov = do.call(rbind,mu.cov) #mu.l12[j,ttm] = E(L1L2|X = ttm, S = j)
  return(mu.cov)
}
ind2 = rep(1:n.cov,n.cov)
ind1 = ind2[order(ind2, decreasing = F)]
ind = data.frame(ind1 = ind1, ind2 = ind2)
ind$cov1 =  paste0("l",ind[,1])
ind$cov2 =  paste0("l",ind[,2])

# mu.lab$sj.x: E(LaLb|x,j) where La is ind[i,3] and Lb is ind[i,4]
mu.lab = lapply(1:dim(ind)[1], function(i) {
  u = c(ind[i,3],ind[i,4])
  cova.func(u = u)})
for (i in 1:length(mu.lab)) 
  mu.lab[[i]] = as.vector(t(mu.lab[[i]]))
mu.lab = as.data.frame(do.call(rbind, mu.lab))
ag.ind = rep(1:Z,2)
ag.ind = ag.ind[order(ag.ind, decreasing = F)]
ag.ind = paste0(rep("s",2*Z),ag.ind)
x.ind = rep(paste0(rep("x",2),0:1),Z)
colnames(mu.lab) = paste(ag.ind,x.ind, sep=".")
mu.lab = cbind(ind, mu.lab)

# ---------- Estimate E(YL|x,j) ----------------
cova3 = cova1[,1:n.cov]

mu.yl = lapply(1:n.cov, function(l){
  lapply(1:Z, function(j){
    sapply(1:(K-Z), function(k){
      sapply(0:1, function(ttm){
        (ps[j]*(px[k+Z]*(ttm == 1) + (1-px[k+Z])*(ttm == 0)))^(-1)*
          (sum(
            (new$x == ttm)*(new$s == (k+Z))*
              cova3[,l]*expit(phi.ag.e[[j]][1,k] + phi.ag.e[[j]][2,k]*(new$x) + cova1%*%phi.c[3:length(phi.c[,1]),k])*
              exp(cova %*% b.est[[j]][,k]))/n)})
    })})})
for(l in 1:ncol(cova3)) mu.yl[[l]] = lapply(1:Z, function(j) rowMeans(mu.yl[[l]][[j]])) 
for(l in 1:ncol(cova3)) mu.yl[[l]] = do.call(rbind,mu.yl[[l]]) # mu.yl[[l]][j,ttm] = E(YL_l|S=j,X=ttm)

# ---------- Re-estimate E(Y|x,j) ----------------
mu.yq1 = lapply(1:Z, function(j){
  sapply(1:(K-Z), function(k){
    sapply(0:1, function(ttm){
      new.int = sapply(1:length(l.int), function(i) ttm * new[,l.int[i]])
      colnames(new.int) = paste0(rep("ttm.l", length(n.int)),n.int)
      
      new = cbind(new,new.int)
      cova1.new = as.matrix(new[,c(paste0(rep("l",n.cov),1:n.cov), colnames(new.int))])
      
      (ps[j]^(-1))*sum((new$s == (k+Z)) * 
                                     expit(phi.ag.e[[j]][1,k] + phi.ag.e[[j]][2,k]*ttm + cova1.new%*%phi.c[3:length(phi.c[,1]),k]) * 
                                     exp(cova %*% b.est[[j]][,k]))/n
    })
  })})
mu.yq1 = lapply(1:Z, function(j) rowMeans(mu.yq1[[j]]))

# ----------- Regenerate L1, L2 and Y given ttm and j ----------
psd = lapply(1:Z, function(j){
  lapply(0:1, function(ttm){
    m.yl = mu[[j]][,(ttm+1)]
    col.var = paste0("s",j,".x",ttm)
    mat.yl.cov = matrix(mu.lab[,col.var], nrow = n.cov, ncol = n.cov) 
    mat.yl = rbind(mat.yl.cov, 
                  sapply(1:n.cov, function(i) mu.yl[[i]][j,ttm+1]))
    mat.yl = cbind(mat.yl, 
                  c(sapply(1:n.cov, function(i) mu.yl[[i]][j,ttm+1]),mu.yq1[[j]][(ttm+1)]))
    mat.yl = mat.yl - m.yl %*% t(m.yl)
    
    psd.j = as.data.frame(mvrnorm(n = nj[[j]][(ttm+1)],
                                  mu = m.yl,
                                  Sigma = mat.yl,
                                  empirical = TRUE)) # Pseudo-data
    psd.j$x = ttm
    psd.j$s = j
    psd.j$int = 1
    psd.int = sapply(1:length(l.int), function(i) psd.j$x * psd.j[,l.int[i]])
    colnames(psd.int) = paste0(rep("xl", length(n.int)),n.int)
    psd.j = cbind(psd.j,psd.int)
    psd.j = psd.j[,c('int','x',l.cov,'y','s',colnames(psd.int))]
    return(psd.j)})
})
psd = lapply(1:Z, function(j) do.call(rbind, psd[[j]]))
psd = do.call(rbind, psd)
new = rbind(new, psd)

#----------------
# M-estimation 
#----------------

#------------ Estimating function ------------------

# Estimating function
func = function(data, ps, px, phi.c, b.est, phi.ag.e){
  
  # ps and px
  f0 = cbind(sapply(1:K, function(i) (data$s==i) - ps[i]), # P(S = j) # (1 to 5)
             sapply(1:K, function(i) (data$s==i)*data$x/ps[i] - px[i])) # P(X=1|S) # (6 to 10)
  
  # phi.c
  cova = as.matrix(data[,c("int", "x", l.cov, paste0(rep("x",length(l.int)), l.int))])
  resi.y = lapply(1:(K-Z), 
                  function(i) (data$s==(i+Z))*as.vector(data$y - expit(cova%*%phi.c[,i])))
  f1 = lapply(1:(K-Z),
              function(ipd) t(sapply(1:ncol(cova), function(i) resi.y[[ipd]] * cova[,i])))
  f1.out = do.call(rbind, f1) # 11 to 20
  
  # beta.jk
  cova2 = cova[,c("int",l.cov)]
  resi.s = lapply(1:Z, 
                  function(ag)
                    sapply(1:(K-Z), function(ipd) 
                      (data$s==(ipd+Z))*exp(cova2%*%b.est[[ag]][,ipd]))) # resi.s[[j]][,k+Z]: study j with AG, study k+Z with IPD
  
  f2 = lapply(1:Z,
              function(ag)
                lapply(1:(K-Z), function(ipd)
                  t(sapply(1:(n.cov + 1), function(i) resi.s[[ag]][,ipd] * cova2[,i] - (data$s == ag)*cova2[,i]))))  #f2[[j]][[k+Z]][,c(est.1, est.l1,est.l2,...)]
  f2.out = lapply(1:Z, function(i) do.call(rbind,f2[[i]]))
  f2.out = do.call(rbind, f2.out) # b.1(Z+1), ..., b.1K, ..., b.Z(Z+1), ..., b.ZK
  
  # phi.ag
  cova3 = cova[,c(l.cov,paste0(rep("x",length(l.int)), l.int))]
  f3 = lapply(1:Z,
              function(ag)
                lapply(1:(K-Z),
                       function(ipd) t(sapply(0:1, function(i) 
                         (data$s==(ipd+Z))*(data$x==i)*
                           expit(phi.ag.e[[ag]][1,ipd] + phi.ag.e[[ag]][2,ipd]*(data$x) + cova3%*%phi.c[3:length(phi.c[,1]),ipd])*
                           as.vector(exp(cova2%*%b.est[[ag]][,ipd]))*
                           (px[ag]*(i==1)/px[ipd+Z] + (i==0)*(1-px[ag])/(1-px[ipd+Z]))- 
                           data$y*(data$s==ag)*(data$x==i))))) # f3[[j]][[k+Z]][,c(phi_0j,phi_1j)]
  
  f3.out = lapply(1:Z, function(i) do.call(rbind,f3[[i]]))
  f3.out = do.call(rbind, f3.out) # phi.1(Z+1), ..., phi.1K, ..., phi.Z(Z+1), ..., phi.ZK
  
  # summarize
  return(rbind(t(f0), f1.out, f2.out, f3.out))
}

# ------ Matrix B ------
phi1 = func(data = new, 
            ps=ps, 
            px=px, 
            phi.c=phi.c, 
            b.est=b.est, 
            phi.ag.e=phi.ag.e)
mat.b1 = var(t(phi1))

# ------ Matrix A ------

# ps & px part
h = 1e-7
identity = diag(length(c(ps,px)))
dev = lapply(1:length(c(ps,px)), function(u){
  p. = complex(real = c(ps,px), imaginary = h*identity[u,])
  out = rowSums(Im(func(data = new, 
                        ps=p.[1:K], 
                        px=p.[(K+1):(2*K)], 
                        phi.c=phi.c, 
                        b.est=b.est, 
                        phi.ag.e=phi.ag.e))/h)
  return(out)})
mat.a.psx = -do.call(cbind,dev)/n

# phi.c
identity = diag(length(phi.c[,1]))
dev = lapply(1:(K-Z), function(ipd){
  sapply(1:length(phi.c[,ipd]), function(i){
    phi.cx= phi.c
    phi.cx[,ipd] = complex(real = phi.c[,ipd], imaginary = h*identity[i,])
    out = rowSums(Im(func(data = new, 
                          ps=ps, 
                          px=px, 
                          phi.c=phi.cx, 
                          b.est=b.est, 
                          phi.ag.e=phi.ag.e))/h)
    return(out)})
})
mat.a.phic = -do.call(cbind,dev)/n

# b.jk
identity = diag(length(b.est[[1]][,1]))
dev = lapply(1:Z, function(ag){
  lapply(1:(K-Z), function(ipd){
    sapply(1:length(b.est[[ag]][,ipd]), function(i){
      b.estx= b.est
      b.estx[[ag]][,ipd] = complex(real = b.est[[ag]][,ipd], imaginary = h*identity[i,])
      out = rowSums(Im(func(data = new, 
                            ps=ps, 
                            px=px, 
                            phi.c=phi.c, 
                            b.est=b.estx, 
                            phi.ag.e=phi.ag.e))/h)
      return(out)})
  })
})
dev = lapply(1:Z, function(i) do.call(cbind,dev[[i]]))
mat.a.best = -do.call(cbind,dev)/n

# phi.ag.e
identity = diag(length(phi.ag.e[[1]][,1]))
dev = lapply(1:Z, function(ag){
  lapply(1:(K-Z), function(ipd){
    sapply(1:length(phi.ag.e[[ag]][,ipd]), function(i){
      phi.ag.ex= phi.ag.e
      phi.ag.ex[[ag]][,ipd] = complex(real = phi.ag.e[[ag]][,ipd], imaginary = h*identity[i,])
      out = rowSums(Im(func(data = new, 
                            ps=ps, 
                            px=px, 
                            phi.c=phi.c, 
                            b.est=b.est, 
                            phi.ag.e=phi.ag.ex))/h)
      return(out)})
  })
})
dev = lapply(1:Z, function(i) do.call(cbind,dev[[i]]))
mat.a.phiagt = -do.call(cbind,dev)/n

mat.a = cbind(mat.a.psx, mat.a.phic, mat.a.best, mat.a.phiagt)

# ----------- Sandwich estimator -------------

var.p1= solve(mat.a) %*% mat.b1 %*% t(solve(mat.a))/n

#---------------
# Delta method
#---------------

# Point estimates of c(phi0s, phi1s)
phi.ag = lapply(1:Z, function(i) as.vector(phi.ag.e[[i]]))
phi = c(as.vector(phi.c[1:2,]),  unlist(phi.ag)) # phi.1(Z+1), ..., phi.1K, ..., phi.Z(Z+1), ..., phi.ZK

# Covariance matrix of c(phi0s, phi1s)
ind.phic = 2*K + as.vector(sapply(0:(K-Z-1), function(i) c(1,2) + i*(ncol(data) - 2)))
ind.phiag = (1 + 2*K + ncol(mat.a.phic) + ncol(mat.a.best)) : nrow(mat.b1)
ind = c(ind.phic, ind.phiag)
var.phi1 = var.p1[ind,ind] # partial IPD

# Linear operator matrix
i1 = cbind(rep(c(1/(K-Z),0),K-Z), rep(c(0,1/(K-Z)),K-Z))
i2 = matrix(0, ncol = 2, nrow = 2*(K-Z))
Tr = lapply(1:Z, function(ag) cbind(do.call(cbind, replicate(ag-1, i2,simplify = F)),
                                    i1, 
                                    do.call(cbind, replicate(Z-ag, i2,simplify = F))))
Tr = do.call(rbind, Tr)
Tr.ag = rbind(matrix(0, ncol = ncol(Tr), nrow = length(ind) - nrow(Tr)),
              Tr)
Tr.ipd = rbind(diag((K-Z)*2),
               matrix(0, nrow = length(ind) - (K-Z)*2, ncol = (K-Z)*2))
Tr.full = t(cbind(Tr.ipd, Tr.ag))

# Final estimates of c(phi0s, phi1s) after aggregation
phi.fn = as.vector(Tr.full %*% phi)

# Covariance matrix of aggregated estimates: partial IPD
v1 = Tr.full %*% var.phi1 %*% t(Tr.full)
v.phiag1 = diag(v1)[(2*(K-Z)+1):(2*K)]

#----------------------------
# Meta-analysis: partial IPD
#----------------------------

phi1.ma = phi.fn[2*(1:K)]
v.phi1ma = v1[2*(1:K),2*(1:K)]
Tr.ma = t(rep(1/K,K))

# Summary estimate
(te = mean(phi1.ma))
# Var(te)
(v.te = Tr.ma %*% v.phi1ma %*% t(Tr.ma))
# 95% CI
(ci = te + qnorm(c(.025,.975))*sqrt(as.vector(v.te)))

phi.all
phi1.ma
phi.ag.e

rowMeans(phi.all)
rowMeans(phi.c)

lapply(1:K,function(i)
  summary(glm(as.formula(paste("y~", lp)), 
           data = ori[which(ori$s == i),], family = binomial('logit'))))

