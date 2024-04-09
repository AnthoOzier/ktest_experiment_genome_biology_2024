# Data generation

truth <- c(rep("DE",250),rep("DM",250),rep("DP",250),rep("DB",250),rep("EE",4500),rep("EP",4500))
nameRow <- c(paste0("DE_",1:250), paste0("DM_",1:250), paste0("DP_",1:250), paste0("DB_",1:250), paste0("EE_",1:4500),
             paste0("EP_",1:4500))

sample_mat_NBP <- function(n_G=10000, n=100, s = sample(1:500, size = 1),pzmin=.7,pzmax=.9){

  stopifnot(n_G == 10000)
  set.seed(s)  
  # set.seed(2101 + slar_taskid + n)
  seeds <- runif(n = n_G, min = 1, max = 2^30)

  X <- c(rep(1,n/2),rep(2,n/2))
  Y <- matrix(NA, nrow = n_G, ncol = n)

  propNonZeros <- runif(n = n_G, min = pzmin, max = pzmax)
  p <- floor(rep(n/4,n_G))
  p1 <- floor(rep(n/6,n_G))
  p2 <- (n/2)-p1
  moy0 <- runif(n_G,0,0.5)
  moy1 <- runif(n_G,10,20)
  moy2 <- 3*moy1
  moy3 <- (moy1+moy2)/2

  for (i in 1:n_G) {
    set.seed(seeds[i])
    zeros <- rbinom(n = n, size = 1, p=propNonZeros[i])

    if (i<=250){ # DE
      val <- c(rnbinom(n/2,prob=.5,size=moy1[i]),
               rnbinom(n/2,prob=.5,size=moy3[i]))
    }else if (i<=500 & i>250){ # DM
      bimod <- rep(NA, n/2)
      unimod <- rep(NA, n/2)
      bimod[1:p[i]] <- rnbinom(floor(length(1:p[i])), prob=.5, size=moy1[i])
      bimod[(p[i]+1):(n/2)] <- rnbinom(floor(length((p[i]+1):(n/2))),prob=.5,size=moy2[i])
      unimod <- rnbinom(n/2,prob=.5,size=moy2[i])

      val <- c(unimod,bimod)
    }else if (i<=750 & i>500){ # DP
      bimod1 <- rep(NA, n/2)
      bimod2 <- rep(NA, n/2)
      bimod1[1:p1[i]] <- rnbinom(length(1:p1[i]),prob=.5,size=moy1[i])
      bimod1[(p1[i]+1):(n/2)] <- rnbinom(length((p1[i]+1):(n/2)),prob=.5,size=moy2[i])
      bimod2[1:p2[i]] <- rnbinom(length(1:p2[i]),prob=.5,size=moy1[i])
      bimod2[(p2[i]+1):(n/2)] <- rnbinom(length((p2[i]+1):(n/2)),prob=.5,size=moy2[i])

      val <- c(bimod1,bimod2)
    }else if (i<=1000 & i>750){ # DB
      bimod <- rep(NA, n/2)
      unimod <- rep(NA, n/2)
      bimod[1:p[i]] <- rnbinom(length(1:p[i]),prob=.5,size=moy1[i])
      bimod[(p[i]+1):(n/2)] <-rnbinom(length((p[i]+1):(n/2)),prob=.5,size=moy2[i])
      unimod <- rnbinom(n/2,prob=.5,size=moy3[i])

      val <- c(unimod,bimod)
    }else if (i<=5500 & i>1000){# EE
      val <- rnbinom(n,size=moy1[i],mu=10)
    }else if (i>5500){# EP
      val <- sample(c(rnbinom(2*p[i],size=moy1[i],mu=30),
                    rnbinom(n-2*p[i], size=moy2[i],mu=10)))
    }

    Y[i,] <- val * zeros
  }

  rownames(Y) <- nameRow
  return(list(Y=Y, X=X,propNonZeros=propNonZeros))
}



seeds=0:500
size <- c(20,40,60,80,100,160,200)
path='/home/data/experimental_data/'
pzs = c()
for (s in seeds){
  for (n in size){
    pzmin=.7,pzmax=.9
    temp <- sample_mat_NBP(n_G = 10000, n = n, s = s)
    Y <- temp$Y
    X <- temp$X
    propNonZeros <- temp$propNonZeros
    write.csv(Y,file=paste(path,'Y',n,'_s',s,'.csv',sep=''))
    write.csv(X,file=paste(path,'X',n,'_s',s,'.csv',sep=''))
    write.csv(propNonZeros,file=paste(path,'propNonZeros',n,'_s',s,'.csv',sep=''))
    
  }
}


# Zero inflation
seeds=0:100
size <- c(200)
pzs = list('a'=c(.7,1),'b'=c(.8,1),'c'=c(.9,1),'d'=c(.95,1))
# pzs = list('d'=c(.95,1))
print(pzs)
for (pz in names(pzs)){
    pzmin = pzs[[pz]][1]
    pzmax = pzs[[pz]][2]
}

for (pz in names(pzs)){
    pzmin = pzs[[pz]][1]
    pzmax = pzs[[pz]][2]
    pzstr=paste0('pz_',pzmin*100,'_',pzmax*100)
    path=paste0('/home/data/experimental_data_',pzstr,'/')
    dir.create(path)
    for (s in seeds){
        for (n in size){
            print(paste('pz=',pz,'pzmin=',pzmin,'pzmax=',pzmax))
            temp <- sample_mat_NBP(n_G = 10000, n = n, s = s ,pzmin = pzmin , pzmax = pzmax)
            Y <- temp$Y
            X <- temp$X
            propNonZeros <- temp$propNonZeros
            write.csv(Y,file=paste(path,'Y',n,'_s',s,'.csv',sep=''))
            write.csv(X,file=paste(path,'X',n,'_s',s,'.csv',sep=''))
            write.csv(propNonZeros,file=paste(path,'propNonZeros',n,'_s',s,'.csv',sep='')) 
        }   
    }
}


