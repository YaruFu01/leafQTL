djg_first <- 1;
target    <- "ck";

dat <- read.table("biaoxing_changkuan.txt", head=FALSE);
changdu <- as.vector(dat[,2]);
kuandu  <- as.vector(dat[,3]);
dat <- read.table("biaoxing_mianji.txt", head=FALSE);
mianji  <- as.vector(dat[,2]);
dat <- as.vector(dat[,1]);

n <- length(dat);
raw_info <- matrix(0, n, 5);
for(i in 1 : n){
     res <- dat[i];
     res <- gsub(".txt", "", res);
     res <- strsplit(res, '_');
     res <- res[[1]];   #0145_CK_P1_T1_rep1.txt
     raw_info[i,1:5] <- res[1:5]; 
}

raw_info  <- cbind(raw_info, changdu);
colnames(raw_info) <- c("xihao", "chuli", "buwei", "shijian", "rplct", "changdu");
rownames(raw_info) <- seq(1 : dim(raw_info)[1]);
raw_info <- data.frame(raw_info);
raw_info$changdu <- as.numeric(changdu);


get_subdata <- function(raw_info, factor1, factor2){
  id1 <- which(raw_info[,2] == factor1);
  id2 <- which(raw_info[,4] == factor2);
  id  <- intersect(id1, id2);
  res <- raw_info[id,];
  return( res );
}

raw_ck_t1 <- get_subdata(raw_info, "CK", "T1"); 
raw_ck_t2 <- get_subdata(raw_info, "CK", "T2"); 
raw_ck_t3 <- get_subdata(raw_info, "CK", "T3"); 
raw_ck_t4 <- get_subdata(raw_info, "CK", "T4");
raw_ck_t5 <- get_subdata(raw_info, "CK", "T5");
raw_ck_t6 <- get_subdata(raw_info, "CK", "T6");
raw_salt_t1 <- get_subdata(raw_info, "Salt", "T1");
raw_salt_t2 <- get_subdata(raw_info, "Salt", "T2");
raw_salt_t3 <- get_subdata(raw_info, "Salt", "T3");
raw_salt_t4 <- get_subdata(raw_info, "Salt", "T4");
raw_salt_t5 <- get_subdata(raw_info, "Salt", "T5");
raw_salt_t6 <- get_subdata(raw_info, "Salt", "T6");


get_longitudinal_data <- function(sub_data){
    chuli   <- as.vector(sub_data[1,2]);
    shijian <- as.vector(sub_data[1,4]);

    res <- NULL;
    xihao_uniq    <- unique(sub_data[,1]);
    for(xh in xihao_uniq){
      id1 <- which(sub_data[,1] == xh);
      target_data <- sub_data[id1, ];
      
      chongfu_uniq  <- unique(target_data[,5]);
      buwei_uniq    <- unique(target_data[,3]);
      
      row_num <- length(chongfu_uniq);
      res_tmp <- matrix(NA, row_num, 8);
      rownames(res_tmp) <- rep(xh, row_num);

      jishu <- 0;
      for(cf in chongfu_uniq){
      	jishu <- jishu + 1;
      	for(bw in buwei_uniq){
      	   id2 <- which(target_data[,5] == cf);
      	   id3 <- which(target_data[,3] == bw);
           id123 <- intersect(id2, id3);
           if(length(id123) == 1){
           	  buwei_id <- as.numeric( gsub("P","",bw) );
           	  res_tmp[jishu,buwei_id] <- target_data[id123,6];
           }
      	}
      }
      res <- rbind(res, res_tmp);  
    }
    colnames(res) <- c("P1","P2","P3","P4","P5","P6","P7","P8");
    return(res);
}

data_ck_t1 <- get_longitudinal_data(raw_ck_t1);
data_ck_t2 <- get_longitudinal_data(raw_ck_t2);
data_ck_t3 <- get_longitudinal_data(raw_ck_t3);
data_ck_t4 <- get_longitudinal_data(raw_ck_t4);
data_ck_t5 <- get_longitudinal_data(raw_ck_t5);
data_ck_t6 <- get_longitudinal_data(raw_ck_t6);
data_salt_t1 <- get_longitudinal_data(raw_salt_t1);
data_salt_t2 <- get_longitudinal_data(raw_salt_t2);
data_salt_t3 <- get_longitudinal_data(raw_salt_t3);
data_salt_t4 <- get_longitudinal_data(raw_salt_t4);
data_salt_t5 <- get_longitudinal_data(raw_salt_t5);
data_salt_t6 <- get_longitudinal_data(raw_salt_t6);

EI_ck_t1 <- as.vector(apply(data_ck_t1, 2, function(X){mean(X,na.rm=TRUE);}))
EI_ck_t2 <- as.vector(apply(data_ck_t2, 2, function(X){mean(X,na.rm=TRUE);}))
EI_ck_t3 <- as.vector(apply(data_ck_t3, 2, function(X){mean(X,na.rm=TRUE);}))
EI_ck_t4 <- as.vector(apply(data_ck_t4, 2, function(X){mean(X,na.rm=TRUE);}))
EI_ck_t5 <- as.vector(apply(data_ck_t5, 2, function(X){mean(X,na.rm=TRUE);}))
EI_ck_t6 <- as.vector(apply(data_ck_t6, 2, function(X){mean(X,na.rm=TRUE);}))
EI_salt_t1 <- as.vector(apply(data_salt_t1, 2, function(X){mean(X,na.rm=TRUE);}))
EI_salt_t2 <- as.vector(apply(data_salt_t2, 2, function(X){mean(X,na.rm=TRUE);}))
EI_salt_t3 <- as.vector(apply(data_salt_t3, 2, function(X){mean(X,na.rm=TRUE);}))
EI_salt_t4 <- as.vector(apply(data_salt_t4, 2, function(X){mean(X,na.rm=TRUE);}))
EI_salt_t5 <- as.vector(apply(data_salt_t5, 2, function(X){mean(X,na.rm=TRUE);}))
EI_salt_t6 <- as.vector(apply(data_salt_t6, 2, function(X){mean(X,na.rm=TRUE);}))


processing_data <- function(data, EI){
   rk    <- rank(EI);
   id1 <- which(rk == 1);
   id2 <- which(rk == 2);
   id3 <- which(rk == 3);
   id4 <- which(rk == 4);
   id5 <- which(rk == 5);
   id6 <- which(rk == 6);
   id7 <- which(rk == 7);
   id8 <- which(rk == 8);
   xuhao <- c(id1, id2, id3, id4,id5,id6,id7,id8);

   new_EI   <- EI[xuhao];
   
   new_y  <- data[,xuhao];
   n <- dim(new_y)[1];
   new_x  <- matrix( rep(new_EI,n), n, 8, byrow=TRUE);

   res <- list();
   res[[1]] <- new_x;
   res[[2]] <- new_y;
   return( res );
}

ck_t1 <- processing_data(data_ck_t1, EI_ck_t1);
ck_t2 <- processing_data(data_ck_t2, EI_ck_t2);
ck_t3 <- processing_data(data_ck_t3, EI_ck_t3);
ck_t4 <- processing_data(data_ck_t4, EI_ck_t4);
ck_t5 <- processing_data(data_ck_t5, EI_ck_t5);
ck_t6 <- processing_data(data_ck_t6, EI_ck_t6);
salt_t1 <- processing_data(data_salt_t1, EI_salt_t1);
salt_t2 <- processing_data(data_salt_t2, EI_salt_t2);
salt_t3 <- processing_data(data_salt_t3, EI_salt_t3);
salt_t4 <- processing_data(data_salt_t4, EI_salt_t4);
salt_t5 <- processing_data(data_salt_t5, EI_salt_t5);
salt_t6 <- processing_data(data_salt_t6, EI_salt_t6);

ck_x   <- NULL;  ck_y   <- NULL;
salt_x <- NULL;  salt_y <- NULL;
ck_x  <- rbind(ck_x, ck_t1[[1]]);  ck_y  <- rbind(ck_y, ck_t1[[2]]);
ck_x  <- rbind(ck_x, ck_t2[[1]]);  ck_y  <- rbind(ck_y, ck_t2[[2]]);
ck_x  <- rbind(ck_x, ck_t3[[1]]);  ck_y  <- rbind(ck_y, ck_t3[[2]]);
ck_x  <- rbind(ck_x, ck_t4[[1]]);  ck_y  <- rbind(ck_y, ck_t4[[2]]);
ck_x  <- rbind(ck_x, ck_t5[[1]]);  ck_y  <- rbind(ck_y, ck_t5[[2]]);
ck_x  <- rbind(ck_x, ck_t6[[1]]);  ck_y  <- rbind(ck_y, ck_t6[[2]]);
salt_x  <- rbind(salt_x, salt_t1[[1]]);  salt_y  <- rbind(salt_y, salt_t1[[2]]);
salt_x  <- rbind(salt_x, salt_t2[[1]]);  salt_y  <- rbind(salt_y, salt_t2[[2]]);
salt_x  <- rbind(salt_x, salt_t3[[1]]);  salt_y  <- rbind(salt_y, salt_t3[[2]]);
salt_x  <- rbind(salt_x, salt_t4[[1]]);  salt_y  <- rbind(salt_y, salt_t4[[2]]);
salt_x  <- rbind(salt_x, salt_t5[[1]]);  salt_y  <- rbind(salt_y, salt_t5[[2]]);
salt_x  <- rbind(salt_x, salt_t6[[1]]);  salt_y  <- rbind(salt_y, salt_t6[[2]]);

zhengli <- function(ck_x, ck_y){
	new_ckx <- matrix(0, dim(ck_x)[1], dim(ck_x)[2]);
	new_cky <- matrix(0, dim(ck_x)[1], dim(ck_x)[2]);
	rownames(new_ckx) <- rownames(ck_y);
	rownames(new_cky) <- rownames(ck_y);

	for(i in 1:dim(ck_y)[1]){
       tmp_id <- which( !is.na(ck_y[i,]) );
       if(length(tmp_id) > 0){
       	  new_ckx[i,] <- c(ck_x[i,tmp_id], rep(0,dim(ck_x)[2]-length(tmp_id))); 
       	  new_cky[i,] <- c(ck_y[i,tmp_id], rep(0,dim(ck_y)[2]-length(tmp_id))); 
       }
	}

	res <- list();
	res[[1]] <- new_ckx;
	res[[2]] <- new_cky;
	return( res );
}

res_ck   <- zhengli(ck_x, ck_y);
res_salt <- zhengli(salt_x, salt_y);
ck_x   <- res_ck[[1]];      #####****
ck_y   <- res_ck[[2]];      #####****
salt_x <- res_salt[[1]];    #####****
salt_y <- res_salt[[2]];    #####****
res_ck  <- NULL;
res_salt <- NULL;



##################################################################################
estimate_h0 <- function(par_all, geti_dianshu, pheno, datee){
    if( any(par_all <= 0) | par_all[3] >= 1 ){return(NaN);}
    parin <- par_all[1:2];
    para  <- par_all[3:4]; 
    logL <- .C("optimize4", parin, para, geti_dianshu, pheno, datee)[[1]][1];
    loglike <- 0 - logL;
    return( loglike );
}

estimate_h0_ab <- function(parin, para, geti_dianshu, pheno, datee){
    if( any(parin <= 0) ){return(NaN);}
    logL <- .C("optimize4", parin, para, geti_dianshu, pheno, datee)[[1]][1];
    loglike <- 0 - logL;
    return( loglike );
}

estimate_h0_psiv2 <- function(para, parin, geti_dianshu, pheno, datee){
    if( any(para <= 0) | para[1] >= 1 ){return(NaN);}
    logL <- .C("optimize4", parin, para, geti_dianshu, pheno, datee)[[1]][1];
    loglike <- 0 - logL;
    return( loglike );
}

estimate_H0 <- function(PARS, geti_dianshu, pheno, datee){
   min.val <- Inf;
   min.par <- c();

   parin <- PARS[1:2];
   para  <- PARS[3:4];
   for(ii in 1:1){
     out <- try( optim(para, estimate_h0_psiv2, parin = parin, 
                geti_dianshu = geti_dianshu, pheno = pheno, datee = datee, method="Nelder-Mead"), TRUE);
     if(class(out) == "try-error"){
         para <- c(0.6365, 296) * runif(2, 0, 2);
         next;
     }else{
         if(out$val < min.val){
            min.val <- out$val;
            min.par <- c(parin, out$par);
            para <- out$par * runif(2, 0.8, 1.2); 
         }else{ para <- para * runif(2, 0.8, 1.2); }
         cat("loop_h0_psiv2", ii, min.val, min.par, "  ", out$val, parin, out$par, "\n" );
     }
   }


   par_all <- min.par;
   count <- 1;
   while(count <= 3){
   if(count == 1){ loops <- 10; }
   if(count == 2){ loops <- 2; }
   if(count == 3){ loops <- 1; }
   for(ii in 1:loops){
     out <- try( optim(par_all, estimate_h0,
                geti_dianshu = geti_dianshu, pheno = pheno, datee = datee, method="Nelder-Mead"), TRUE);
     if(class(out) == "try-error"){
         par_all <- c(0.713, 1.10, 0.6365, 296) * runif(4, 0, 1);
         if(par_all[3] > 1){ par_all[3] <- runif(1, 0, 1); }
         next;
     }else{
         if(out$val < min.val){
            min.val <- out$val;
            min.par <- out$par;
            par_all <- out$par * runif(4, 0, 1); 
         }else{ par_all <- par_all * runif(4, 0, 1); }
         cat("loop_h0", ii, min.val, min.par, "  ", out$val, out$par, "\n" );
     }
   }
   
   parin <- min.par[1:2];
   para  <- min.par[3:4];
   for(ii in 1:loops){
     out <- try( optim(parin, estimate_h0_ab, para = para,
                geti_dianshu = geti_dianshu, pheno = pheno, datee = datee, method="Nelder-Mead"), TRUE);
     if(class(out) == "try-error"){
         parin <- c(0.713, 1.10) * runif(2, 0.9, 1.1);
         next;
     }else{
         if(out$val < min.val){
            min.val <- out$val;
            min.par <- c(out$par, para);
            parin <- out$par * runif(2, 0.9, 1.1); 
         }else{ parin <- parin * runif(2, 0.9, 1.1); }
         cat("loop_h0_ab", ii, min.val, min.par, "  ", out$val, out$par, para, "\n" );
     }
   }
   
   parin <- min.par[1:2];
   para  <- min.par[3:4];
   for(ii in 1:1){
     out <- try( optim(para, estimate_h0_psiv2, parin = parin, 
                geti_dianshu = geti_dianshu, pheno = pheno, datee = datee, method="Nelder-Mead"), TRUE);
     if(class(out) == "try-error"){
         para <- c(0.6365, 296) * runif(2, 0, 1);
         next;
     }else{
         if(out$val < min.val){
            min.val <- out$val;
            min.par <- c(parin, out$par);
            para <- out$par * runif(2, 0, 1); 
         }else{ para <- para * runif(2, 0, 1); }
         cat("loop_h0_psiv2", ii, min.val, min.par, "  ", out$val, parin, out$par, "\n" );
     }
   }
   
   count <- count + 1;
   par_all <- min.par;
   }
   
   res <- c(min.val, min.par);
   return( res );
}

#########
estimate_h1_g2 <- function(par_all, pheno0, datee0, geti_dianshu0, pheno1, datee1, geti_dianshu1){
    if(any(par_all <= 0)){return(NaN);}
    if(any(par_all[5] > 1)){return(NaN);}
    
    parin0 <- par_all[1:2];
    parin1 <- par_all[3:4];
    para   <- par_all[5:6];
        
    logL0 <- .C("optimize4", parin0, para, geti_dianshu0, pheno0, datee0)[[1]][1];
    logL1 <- .C("optimize4", parin1, para, geti_dianshu1, pheno1, datee1)[[1]][1];
    loglike <- 0 - logL0 - logL1;
    return( loglike );
}

estimate_h1_g2_ab <- function(par_all, para, pheno0, datee0, geti_dianshu0, pheno1, datee1, geti_dianshu1){
    if(any(par_all <= 0)){return(NaN);}
    
    parin0 <- par_all[1:2];
    parin1 <- par_all[3:4];
    
    logL0 <- .C("optimize4", parin0, para, geti_dianshu0, pheno0, datee0)[[1]][1];
    logL1 <- .C("optimize4", parin1, para, geti_dianshu1, pheno1, datee1)[[1]][1];
    loglike <- 0 - logL0 - logL1;
    return( loglike );
}

estimate_h1_g2_psiv2 <- function(para, par_all, pheno0, datee0, geti_dianshu0, pheno1, datee1, geti_dianshu1){
    if(any(para <= 0)){return(NaN);}
    if(para[1]  >= 1){ return(NaN);}
    
    parin0 <- par_all[1:2];
    parin1 <- par_all[3:4];
    
    logL0 <- .C("optimize4", parin0, para, geti_dianshu0, pheno0, datee0)[[1]][1];
    logL1 <- .C("optimize4", parin1, para, geti_dianshu1, pheno1, datee1)[[1]][1];
    loglike <- 0 - logL0 - logL1;
    return( loglike );  
}

estimate_H1_G2 <- function(PARS, pheno0, datee0, geti_dianshu0, pheno1, datee1, geti_dianshu1){
   min.val <- Inf;
   min.par <- PARS;
   
   par_all <- PARS;
   count <- 1;
   while(count <= 1){
   for(ii in 1:30){
     out <- try( optim(par_all, estimate_h1_g2,
                pheno0 = pheno0, datee0 = datee0, geti_dianshu0 = geti_dianshu0, 
                pheno1 = pheno1, datee1 = datee1, geti_dianshu1 = geti_dianshu1, method="Nelder-Mead"), TRUE);
     if(class(out) == "try-error"){
         par_all <- c(0.713, 1.10, 0.713, 1.10, 0.6365, 296) * runif(6, 0.5, 1.5);
         if(par_all[5] > 1){ par_all[5] <- runif(1, 0, 1); }
         next;
     }else{
         if(out$val <= min.val){
            min.val <- out$val;
            min.par <- out$par;
            par_all <- out$par * runif(6, 0.5, 1.5); 
         }else{ par_all <- par_all * runif(6, 0.5, 1.5); }
         cat("loop_h1_g2", ii, min.val, min.par, "  ", out$val, out$par, "\n" );
     }
   }
   
   par_all <- min.par[1:4];
   para    <- min.par[5:6];
   for(ii in 1:20){
     out <- try( optim(par_all, estimate_h1_g2_ab, para=para,
                pheno0 = pheno0, datee0 = datee0, geti_dianshu0 = geti_dianshu0, 
                pheno1 = pheno1, datee1 = datee1, geti_dianshu1 = geti_dianshu1, method="Nelder-Mead"), TRUE);
     if(class(out) == "try-error"){
         par_all <- c(0.713, 1.10, 0.713, 1.10) * runif(4, 0.5, 1.5);
         next;
     }else{
         if(out$val <= min.val){
            min.val <- out$val;
            min.par <- out$par;
            par_all <- out$par * runif(4, 0.5, 1.5); 
         }else{ par_all <- par_all * runif(4, 0.5, 1.5); }
         cat("loop_h1_g2_ab", ii, min.val, min.par, para, "  ", out$val, out$par, para, "\n" );
     }
   }
   min.par <- c(min.par, para);
   
   par_all <- min.par[1:4];
   para  <- min.par[5:6];
   for(ii in 1:1){
     out <- try( optim(para, estimate_h1_g2_psiv2, par_all = par_all,
                pheno0 = pheno0, datee0 = datee0, geti_dianshu0 = geti_dianshu0, 
                pheno1 = pheno1, datee1 = datee1, geti_dianshu1 = geti_dianshu1, method="Nelder-Mead"), TRUE);
     if(class(out) == "try-error"){
         para <- c(0.6365, 296) * runif(2, 0.5, 1.5);
         next;
     }else{
         if(out$val <= min.val){
            min.val <- out$val;
            min.par <- out$par;
            para <- out$par * runif(2, 0.5, 1.5); 
         }else{ para <- para * runif(2, 0.5, 1.5); }
         cat("loop_h1_g2_psiv2", ii, min.val, par_all, min.par,"  ", out$val, par_all, out$par, "\n" );
     }
   }
   min.par <- c(par_all, min.par);
   
   count <- count + 1;
   par_all <- min.par;
   }

   res <- c(min.val, min.par);
   return( res );
}

#########
calc_LR2 <- function(geno, new_D, new_H){
    min.val <- Inf;
    min.par <- c();
    
    id  <- which(geno >= 0);
    id0 <- which(geno == 0);
    id1 <- which(geno == 1);
    
    if(length(id) > 0){
      datee   <- as.vector( t(new_D[id,]) );
      pheno   <- as.vector( t(new_H[id,]) );
      id_na  <- c(which(is.na(pheno)), which(is.na(datee)));
      if( length(id_na) > 0 ){ pheno[id_na] <- 0.0; datee[id_na] <- 0.0; }
      geti_dianshu <- as.numeric(c(dim(new_H)[1], dim(new_H)[2]));
    }
    
    if(length(id0) > 0){
      datee0   <- as.vector( t(new_D[id0,]) );
      pheno0   <- as.vector( t(new_H[id0,]) );
      id_na0  <- c(which(is.na(pheno0)), which(is.na(datee0)));
      if( length(id_na0) > 0 ){ pheno0[id_na0] <- 0.0; datee0[id_na0] <- 0.0; }
      geti_dianshu0 <- as.numeric(c(dim(new_H[id0,])[1], dim(new_H[id0,])[2]));
    }
    
    if(length(id1) > 0){
      datee1   <- as.vector( t(new_D[id1,]) );
      pheno1   <- as.vector( t(new_H[id1,]) );
      id_na1  <- c(which(is.na(pheno1)), which(is.na(datee1)));
      if( length(id_na1) > 0 ){ pheno1[id_na1] <- 0.0; datee1[id_na1] <- 0.0; }
      geti_dianshu1 <- as.numeric(c(dim(new_H[id1,])[1], dim(new_H[id1,])[2]));
    }
    
    PARS <- c(0.713, 1.10, 0.6365, 296);
    res_H0 <- estimate_H0(PARS, geti_dianshu, pheno, datee);
    
    PARS1 <- c(res_H0[2:3], res_H0[2:3], res_H0[4:5]);
    res_H1 <- estimate_H1_G2(PARS1, pheno0, datee0, geti_dianshu0, pheno1, datee1, geti_dianshu1);
    res <- c(2*(res_H0[1]-res_H1[1]), res_H0[1], res_H1[1], res_H1[2:7]);
    
    cat("snp:", res, "\n");
    return( res );   
}



##########
get.snp <- function( filename1, filename2, common_xihao ){
  snp  <- read.csv(filename1, head=FALSE);
  indv <- read.csv(filename2, head=FALSE);
  indv <- as.vector(indv[,1]);

  ids <- c();
  for(id in common_xihao){
     tmp <- which(indv == id);
     if(length(tmp) == 1){ ids <- c(ids, tmp); }
  }
  snp <- snp[,ids];

  return( snp );
}

processing_snp_pheno <- function(xuhao, snp, dt_x, dt_y, common_xihao){
   snp_tmp <- as.vector(snp[xuhao, ]);

   geti <- as.vector( rownames(dt_y) );
   id_snp <- c();
   id_ind <- c();
   for(i in 1:length(geti)){
       tmp_id <- which(common_xihao == geti[i]);
       if(length(tmp_id)==1){
         id_snp <- c(id_snp, tmp_id);
         id_ind <- c(id_ind, i);
       }
   }

   new_snp <- snp_tmp[id_snp];
   new_dtx <- dt_x[id_ind, ];
   new_dty <- dt_y[id_ind, ];

   res <- list();
   res[[1]] <- new_snp;
   res[[2]] <- new_dtx;
   res[[3]] <- new_dty;
   return( res );
}

##########--part 1---#################
if(target == "ck"){
   dt_x  <- ck_x;
   dt_y  <- ck_y;
}else{
   dt_x  <- salt_x;
   dt_y  <- salt_y;
}

cexu_nm <- read.csv("huyang_indv.csv", head=FALSE);
cexu_nm <- as.vector(cexu_nm[,1]);
common_xihao <- intersect(unique(rownames(ck_y)), cexu_nm);

snp   <- get.snp( "huyang_012.csv", "huyang_indv.csv", common_xihao );
n_snp <- dim(snp)[1];
jieguo <- matrix(0, n_snp, (1+1+1+1+6) ); ###-- xuhao, LR2, min.val, min.val2, min.par2)
dyn.load("mvn_optim_HD_derived.so");


djg <- djg_first;
each <- ceiling(n_snp/7);
qishi <- (djg - 1) * each + 1;
zhongzhi <- djg*each;
if(zhongzhi > dim(snp)[1]){ zhongzhi <- dim(snp)[1]; } 

savef_nm <- paste("qtl", target, djg, sep = "_");
savef_nm <- paste(savef_nm, ".Rdata", sep = "");

pdf_nm <- paste("qtl", target, djg, sep = "_");
pdf_nm <- paste(savef_nm, ".pdf", sep = "");
pdf( pdf_nm );

for( xuhao in seq(qishi, zhongzhi, 1) ){
    cat("xuhao =", xuhao, "\n");  
    res_tmp <- processing_snp_pheno(xuhao, snp, dt_x, dt_y, common_xihao);
    
    marker  <- as.numeric( res_tmp[[1]] );
    dtx_tmp <- res_tmp[[2]];
    dty_tmp <- res_tmp[[3]];

    id <- which( !is.na(marker) );
    if(length(id) > 0){
       marker <- marker[id];
       dtx_tmp <- dt_x[id,];
       dty_tmp <- dt_y[id,];
    }
    
    if(length(unique(marker)) == 1){
      continue;
    }else{
      res0 <- calc_LR2(marker, dtx_tmp, dty_tmp);
      res <- c(xuhao, res0);
      jieguo[xuhao, ] <- res;
      save(jieguo, file = savef_nm);
      cat("res:", res,"\n\n");

      plot(c(-1,1),c(-1,1), type="n", xlab="EI (mm)", ylab="Leaf length (mm)",
      	   xlim=c(30,max(dt_x,na.rm=TRUE)), ylim=c(min(dt_y,na.rm=TRUE),1.2*max(dt_y,na.rm=TRUE)));
      legend( "topleft", legend=xuhao, bty="n", cex=1.2 );
      for(i in 1:dim(dtx_tmp)[1]){ 
        nn <- length(which(dtx_tmp[i,] > 0));
        lines(dtx_tmp[i,1:nn], dty_tmp[i,1:nn], type="l", col="lightgreen", pch=16); 
      }
      heng <- seq(30, max(dt_x,na.rm=TRUE), 0.5);
      alpha1 <- res0[4];  beta1  <- res0[5];
      alpha2 <- res0[6];  beta2  <- res0[7];
      zong1 <- alpha1 * heng^beta1;
      zong2 <- alpha2 * heng^beta2;
      lines(heng, zong1, type="l", lwd=1, lty=1, col="red");
      lines(heng, zong2, type="l", lwd=1, lty=1, col="blue");
      #Sys.sleep(5);
    }
}
dev.off();
