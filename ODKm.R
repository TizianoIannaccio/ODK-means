odkm<-function(data,K,margin=0,restart=10,mult=1){
  
  nr<-nrow(data) # size of dataset
  
  epsilon<-.0000001 # convergence constant
  
  # output containers
  c_final<-NA 
  f_final<-Inf
  lab_final<-NA

  # distance matrix
  di<-as.matrix(dist(data))
  
  # as many loops as inputed
  for (s in 1:restart){
    
    # randomly sample starting centroids
    idx<-sample(1:nr,K)
    c<-data[idx,]
    
    # initialize convergence as not achieved
    convergence<-FALSE
    f0<-Inf
    
    W<-numeric(nr) # weights container
    
    # this variable will manage exceptions in which an empty cluster
    # is formed (beta, working)
    restart = F
    
    while (!convergence){
      
      lab<-new_assign(data,c) # evaluate clustering matrix
      
      # if there are missing integers in the range 1,...,K it means
      # that an empty cluster has been formed
      if (length(unique(lab))!=K){
        restart = T
        break
      }
      
      # fill W
      for (k in 1:K){
        u_k<-which(lab==k)
        W[u_k]<-outdet(u_k,di,margin,mult)$w

        
        if (length(u_k)==1){
          c[k,]<-data[u_k,]
        } else {
          c[k,]<-colSums(W[u_k]*data[u_k,])
        }
        
      }
      
      # compute within deviance
      f<-0
      
      for (i in 1:nr){
        f=f+sum((data[i,]-c[lab[i],])^2)
      }
      
      # check for convergence
      convergence<-ifelse(f0-f<epsilon,TRUE,FALSE)
      
      # store new objective 
      f0<-f
      
    }
    
    # if an empty cluster is detected at iteration i, said iteration is
    # aborted and restarted
    if (restart){
      i=i-1
      warning("Found an inizialization with an empty cluster, 
              repeating iteration.")
      next
    }
    
    # update results if a lower within deviance has been achieved
    if (f0<f_final){
      f_final<-f0
      c_final<-c
      lab_final<-lab
    }
  }
  
  # build the dummy matrix with cluster assignments
  U<-matrix(0,nr,K)
  for (i in 1:nr){
    U[i,lab_final[i]]=1
  }
  
  # detect final iteration outliers
  outliers<-numeric(0)
  for ( k in 1:K){
    lab<-which(U[,k]==1)
    outliers<-c(outliers,outdet(lab,di,margin,mult)$out)
  }
  
  # this piece of code is needed for compatibility issues
  if (identical(numeric(0),outliers)){
    outliers = NULL
  }
  
  out<-list("U"=U,"c"=c_final,"f"=f_final,"out"=sort(outliers))
  
  return(out)
}

new_assign<-function(data,c){
  
  # evaluate needed dimensions
  k <- nrow(c)
  n <- nrow(data)
  j <- ncol(data)
  
  dif<-matrix(NA,k,n)
  
  # evaluate distances between data points and centroids
  for (k in 1:k){
    dif[k,]<-rowSums(abs(matrix(c[k,],n,j,TRUE)-data))
  }
  
  # find minimal distances and assign data points
  clusters<-numeric(n)
  for (i in 1:n){
    clusters[i]<-which.min(dif[,i])
  }
  
  return(clusters)
  
}

penalty<-function(sums,outliers,mult){
  
  # vector that stores the penalties
  pen<-rep(1,length(sums))
  
  # mean and sd of the distances of non-outliers
  non_outliers<-sums[-outliers]
  m<-mean(non_outliers)
  se<-sd(non_outliers)
  
  # penalty evaluation
  if (length(outliers)>0){
    pen[outliers]<- mult*(se/(sums[outliers]-m))^2
  } 
  
  w<-pen/sum(pen)
  
  return(w)
}

outdet<-function(labels,distance,margin=0,mult){
  
  # we are only interested in a submatrix of the distance matrix
  dist<-distance[labels,labels]
  
  l<-length(labels)

  # managing exceptions
  if (l==1){
    out<-list("k"=1,
              "sums"=numeric(0),
              "out"=numeric(0),
              "w"=1)
    return(out)
  }

  if (l==2){
    out<-list("k"=1,
              "sums"=numeric(0),
              "out"=numeric(0),
              "w"=c(0.5,0.5))
    return(out)
  }

  for (i in 1:l){
    dist[i,]<-sort(dist[i,])
  }

  # index for knn distances
  k<-ifelse(l>2,floor(log(l)),1)

  # core of the outdet function
  sums_k<-sort(rowSums(dist[,(k+1):(2*k+1)]))

  mr<-median(sums_k)

  if (margin==0){
    margin=2*mr # by default, outliers register twice the median isolation
  } else {
    margin = margin*mr
  }

  outl<-which(sums_k>margin)

  # apply penalty function as a weight
  w<-penalty(sums_k,outl,mult)

  out<-list("k"=k,
            "sums"=sums_k,
            "out"=as.numeric(names(sums_k[outl])),
            "w"=w)

  return(out)
}

metrics <- function(n,true,est){

  true_m<-est_m<-matrix(0,n,2)
  
  # build nx2 matrices which tells if a unit (row) is an outlier or not
  if (!is.null(true)){
    true_m[true,2]<-1
    true_m[-true,1]<-1
  } else {
    true_m[,1]<-1
  }
  

  if (!is.null(est)){
    est_m[est,2]<-1
    est_m[-est,1]<-1
  } else {
    est_m[,1]<-1
  }

  # build the confusion matrix
  conf_mat <- t(est_m)%*%true_m # confusion matrix
  
  # build performance measures on the confusion matrix
  accuracy <- (conf_mat[1,1]+conf_mat[2,2])/n
  precision <- conf_mat[1,1]/(conf_mat[1,1]+conf_mat[1,2])
  sensitivity <- conf_mat[1,1]/(conf_mat[1,1]+conf_mat[2,1])
  F1 <- (2*precision*sensitivity)/(precision+sensitivity)
  
  out <- list("measures" = c(accuracy,precision,sensitivity,F1),
              "conf" = conf_mat)
  
  return(out)
}

performance <- function(comb, iter){
  
  mats <- array(NA,dim = c(2,2,nrow(comb)))
  results <- results2 <- matrix(NA,nrow(comb),4)
  colnames(comb)<-c("perc_out","var_fluc","side","r")
  
  for ( i in 1:nrow(comb)){
    mom_results<-mom_results2<-numeric(4)
    message(paste("Exploring combination", i,": ", 
                " perc_out = ", as.numeric(combinations[i,1]),
                " var_fluc = ", as.numeric(combinations[i,2]),
                " side = ", as.numeric(combinations[i,3]),
                " radius = ", as.numeric(combinations[i,4])))
    
    for (j in 1:iter){
      message(paste("   Iteration",j))
      
      # data generation
      s <- tethraset(nk,combinations[i,1],combinations[i,2],combinations[i,3],combinations[i,4])
      data <- s$data
      true_outliers <- s$outliers
      
      kmj <- kmeans_jump(data,4,0,10)
      trim<-trimkmeans(data,4,comb[i,1])
      # accuracy analysis
      out_jump <- kmj$out
      out_trim <- which(trim$classification==5)
      
      met <- metrics(4*nk,true_outliers,out_jump)
      mom_results<-mom_results+met$measures
      
      met2 <- metrics(4*nk,true_outliers,out_trim)
      mom_results2<-mom_results2+met2$measures
    }
    results[i,]<-mom_results/iter
    results2[i,]<-mom_results2/iter
    mats[,,i]<-met$conf
  }
  colnames(results) <- c("Accuracy","Precision","Sensitivity","F1")
  
  out <- list("results" = cbind(comb,results,results2),
              "conf" = mats)
  
  return(out)
}

smallcube <- function(dev,coord){
  
  # simply build a small cube around the clusters
  mat <- matrix(NA,3,2)
  for (i in 1:3){
    mat[i,1]<-coord[i]-dev
    mat[i,2]<-coord[i]+dev
  }
  return(mat)
}

bigcube <- function(ag,tbg,s,r,cube){

  ng<-matrix(runif(3*tbg,-3*r,s+3*r),tbg,3)

  index<-rep(0,tbg)

  # check if generated data points are in the small cubes
  # and in case delete them
  for (i in 1:tbg){
    for (c in 1:4){
      if (ng[1]>cube[1,1,c] & ng[1]<cube[1,2,c] &
          ng[2]>cube[2,1,c] & ng[2]<cube[2,2,c] &
          ng[3]>cube[3,1,c] & ng[3]<cube[3,2,c]){
        index[i]<-1
        break
      }
    }
  }

  if (sum(index)>0){
    ng=ng[-which(index==1),]
  }

  # update already generated data points size
  ag<-rbind(ag,ng)
  
  # update to be generated data points size
  tbg<-sum(index)

  if (tbg>0){
    return(bigcube(ag,tbg,s,r,cube))
  } else {
    return(ag)
  }
}

tethraset<-function(nk,perc_out,var_fluc,side,r){
  
  # impose limit on the outliers percentage
  if (perc_out>0.05){
    perc_out <- 0.05
  } else if (perc_out<0){
    perc_out <- 0
  }

  k<-3 # three dimensions
  
  n<-4*nk # size of the dataset
  
  # build data uniformly distributed in the k-sphere
  Z<-matrix(rnorm(nk*k),nk,k)
  X<-t(apply(Z,1,function(x){x/(sqrt(sum(x*x)))}))

  X<-X*r # adjust radius of each cluster

  # evaluate coordinates of the clusters
  centers<-matrix(c(0,0,0,
                    side,0,0,
                    side/2,side,0,
                    side/2,side/2,sqrt(side^2-2*(side/2)^2)),4,3,T)
  
  # first cluster
  dev<-matrix(rnorm(n*k,0,var_fluc),n,k)
  X_dev<-X+matrix(dev[1:nk,],nk,k,TRUE)

  # second cluster
  X2_dev<-X+matrix(dev[(nk+1):(2*nk),],nk,k,TRUE)
  X2_dev<-X2_dev+matrix(centers[2,],nk,k,TRUE)

  # third cluster
  X3_dev<-X+matrix(dev[(2*nk+1):(3*nk),],nk,k,TRUE)
  X3_dev<-X3_dev+matrix(centers[3,],nk,k,TRUE)

  # fourth cluster
  X4_dev<-X+matrix(dev[(3*nk+1):(4*nk),],nk,k,TRUE)
  X4_dev<-X4_dev+matrix(centers[4,],nk,k,TRUE)

  # merge the data together
  data<-rbind(X_dev,X2_dev,X3_dev,X4_dev)


  # at this point, we want to remove some "regular" data points
  # and replace them with some outliers/inliers
  if (perc_out>0){
    
    # remove data
    toberemoved<-floor(perc_out*n)
    idx<-sample(1:n,toberemoved)
    data<-data[-idx,]
    
    # evaluate true centroids
    centroids <- matrix(NA,4,3)
    lab<-new_assign(data,centers)
    for (v in 1:4){
      centroids[v,]<-colMeans(data[which(lab==v),])
    }
    
    # inliers are 30%
    m = floor(toberemoved*0.3)
    inl_id <- sample(1:4,m,replace = T)
    dev <- matrix(rnorm(3*m,0,var_fluc),m,3)
    inl <- centers[inl_id,]+ dev

    # the rest is outliers
    cubes<-array(NA,c(3,2,4))

    for (c in 1:4){
      cubes[,,c]<-smallcube(var_fluc+3*r,centers[c,])
    }

    outl<-bigcube(ag = numeric(0),
                  tbg = toberemoved-m,
                  s = side,
                  r = r,
                  cube = cubes)

    data <- rbind(data,inl,outl)
    o <- (4*nk-toberemoved+1):(4*nk)
  } else {
    idx <- NULL
    o <- NULL
    centroids <- matrix(NA,4,3)
    lab<-new_assign(data,centers)
    for (v in 1:4){
      centroids[v,]<-colMeans(data[which(lab==v),])
    }
  }


  out<-list("data"=data, "centroids"=centroids, "outliers"=o, "removed"=idx)

  return(out)

}
