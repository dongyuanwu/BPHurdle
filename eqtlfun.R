library(dplyr)
library(tidyselect)
library(fastDummies)
library(Matrix)
library(glmnet)
library(boot)
library(car)

source("obtain_gene_info.R")
genepos <- readRDS("genepos.rds")

eqtlfun <- function(exprmat, snpmat, metadata, donor_id, covariate=NULL,
                    snppos, genepos, genelist=NULL, cisrange=1000000,
                    family="poisson", nboot=100, sig.level=0.05, alpha=1,
                    parallel=FALSE, BPPARAM=bpparam()) {
    
    snppos <- as.data.frame(snppos)
    
    ##### check snps #####
    
    n1 <- ncol(snpmat)
    n2 <- nrow(snppos)
    colnames(snpmat) <- gsub("[[:punct:]]+", "_", colnames(snpmat))
    colnames(snpmat) <- gsub("^(\\d)", "chr\\1", colnames(snpmat))
    rownames(snppos) <- gsub("[[:punct:]]+", "_", rownames(snppos))
    rownames(snppos) <- gsub("^(\\d)", "chr\\1", rownames(snppos))
    name1 <- colnames(snpmat)
    name2 <- rownames(snppos)
    if (n1 != n2) {
        stop(sprintf("There are %s SNPs in snpmat, while %s SNPs in snppos. Please check the data.",
                     n1, n2))
    } else if (!all(name1 %in% name2)){
        stop("The SNP names in snpmat cannot be matched to the SNP names in snppos. Please check the data.")
    } else {
        snppos <- snppos[order(match(name2, name1)), ]
    }
    snppos[, 1] <- gsub("chr", "", tolower(snppos[, 1]))
    
    ##### check covariates #####
    
    newmeta <- metadata[, -c(1:ncol(metadata))]
    partformu <- ""
    if (!is.null(covariate)) {
        if (!(any(covariate %in% colnames(metadata)))) {
            stop("None of the covariates exist in metadata. Please check the data.")
        } else if (!all(covariate %in% colnames(metadata))) {
            covariate <- covariate[covariate %in% colnames(metadata)]
            warning(sprintf("Only %s exist in metadata.", paste(covariate, collapse=", ")))
        }
        newmeta <- cbind(newmeta, metadata[, covariate])
        
        ## make dummy variables
        newmeta <- dummy_cols(newmeta, remove_first_dummy = TRUE, 
                              remove_selected_columns = TRUE)
        rownames(newmeta) <- rownames(metadata)
    }
    
    ##### check donor_id #####
    
    if (length(donor_id) > 1) {
        donor_id <- donor_id[1]
        warning("More than one element exists in donor_id, only the first one will be used.")
    }
    if (!(donor_id %in% colnames(metadata))) {
        stop(sprintf("%s does not exist in metadata. Please check the data."))
    }
    newmeta[, donor_id] <- metadata[, donor_id]
    
    ##### check cells #####
    
    n1 <- ncol(exprmat)
    n2 <- nrow(snpmat)
    n3 <- nrow(newmeta)
    name1 <- colnames(exprmat)
    name2 <- rownames(snpmat)
    name3 <- rownames(newmeta)
    if (n1 != n3) {
        stop(sprintf("There are %s cells in exprmat, while %s cells in metadata. Please check the data.",
                     n1, n3))
    } else if (!all(name1 %in% name3)){
        stop("The cell names in exprmat cannot be matched to the cell names in metadata. Please check the data.")
    } else {
        if(ncol(newmeta) == 1) {
            newmeta[, 1] <- newmeta[order(match(name3, name1)), 1]
        } else {
            newmeta <- newmeta[order(match(name3, name1)), ]
        }
    }
    
    if (!any(name2 %in% newmeta[, donor_id])) {
        stop("The donor ids in snpmat cannot be matched to the donor ids in metadata. Please check the data.")
    } else if (!all(name2 %in% newmeta[, donor_id]) | !all(newmeta[, donor_id] %in% name2)) {
        snpmat <- snpmat[name2 %in% newmeta[, donor_id], ]
        exprmat <- exprmat[, newmeta[, donor_id] %in% name2]
        newmeta <- newmeta[newmeta[, donor_id] %in% name2, ]
    }
    
    ##### check genelist #####
    
    rownames(exprmat) <- gsub("^(\\d)", "gene\\1", rownames(exprmat))
    if (!is.null(genelist)) {
        genelist <- gsub("^(\\d)", "gene\\1", genelist)
        if (!(any(genelist %in% rownames(exprmat)))) {
            stop("None of the genes exist in exprmat. Please check the data.")
        } else if (!all(genelist %in% rownames(exprmat))) {
            genelist <- genelist[genelist %in% rownames(exprmat)]
            warning(sprintf("Only %s genes exist in exprmat.", length(genelist)))
        }
    } else {
        genelist <- rownames(exprmat)
    }
    
    geneinfo <- obtain_gene_info(genelist, genepos)
    geneinfo <- match_gene_snp(geneinfo, snppos, cisrange)
    if (nrow(geneinfo) < length(genelist)) {
        genelist <- genelist[genelist %in% rownames(geneinfo)]
    }

    ##### prepare and fit the model #####
    
    tempfun <- function(k) {
        cispos <- which(snppos[, 1] == geneinfo$chromosome[k] &
                            snppos[, 2] > (geneinfo$start[k] - cisrange) &
                            snppos[, 2] < (geneinfo$start[k] + cisrange))
        snpmat_temp <- as.data.frame(snpmat[, cispos])
        snpmat_temp[, donor_id] <- rownames(snpmat_temp)
        
        exprmat_temp <- exprmat[rownames(exprmat) == genelist[k], ]
        
        dat_temp <- cbind(data.frame(expr=c(t(exprmat_temp))), newmeta)
        dat_temp <- left_join(dat_temp, snpmat_temp, by=donor_id) %>% 
            dplyr::select(-all_of(donor_id))
        
        Y <- as.vector(dat_temp[, 1])
        X <- Matrix(as.matrix(dat_temp[, -1]), sparse = TRUE)
        
        snp_start_in_X <- ncol(newmeta)
        penalfac <- rep(1, ncol(X))
        penalfac[1:(snp_start_in_X-1)] <- 0
        
        start <- Sys.time()
        res_temp <- bphurdle(X, Y, family, nboot, sig.level, alpha, penalfac=penalfac)
        usedtime <- Sys.time() - start
        
        if (!is.null(res_temp$logi)) {
          logires <- res_temp$logi[(ncol(newmeta):ncol(X)) + 1, ]
          rownames(logires) <- colnames(snpmat_temp)[-ncol(snpmat_temp)]
          logiposi <- which(logires[, 5] > 0 | logires[, 6] < 0)
          if (length(logiposi) > 1) {
            logires <- as.data.frame(logires[logiposi, ])
            colnames(logires) <- c("Estimate", "SE", "CILower", "CIUpper", "CILowerP", "CIUpperP")
            logires$Part <- "Logistic"
            logires$Gene <- genelist[k]
            logires$SNP <- rownames(logires)
            logires <- logires[, c(8:9, 1:7)]
          } else if (length(logiposi) == 1) {
            logires <- as.data.frame(t(logires[logiposi, ]))
            colnames(logires) <- c("Estimate", "SE", "CILower", "CIUpper", "CILowerP", "CIUpperP")
            logires$Part <- "Logistic"
            logires$Gene <- genelist[k]
            logires$SNP <- names(logiposi)
            logires <- logires[, c(8:9, 1:7)]
          } else {
            logires <- NULL
          }
        } else {
          logires <- NULL
        }
        
        if (!is.null(res_temp$trun)) {
          trunres <- res_temp$trun[(ncol(newmeta):ncol(X)) + 1, ]
          rownames(trunres) <- colnames(snpmat_temp)[-ncol(snpmat_temp)]
          trunposi <- which(trunres[, 5] > 0 | trunres[, 6] < 0)
          if (length(trunposi) > 1) {
            trunres <- as.data.frame(trunres[trunposi, ])
            colnames(trunres) <- c("Estimate", "SE", "CILower", "CIUpper", "CILowerP", "CIUpperP")
            trunres$Part <- "Zero-truncated"
            trunres$Gene <- genelist[k]
            trunres$SNP <- rownames(trunres)
            trunres <- trunres[, c(8:9, 1:7)]
          } else if (length(trunposi) == 1) {
            trunres <- as.data.frame(t(trunres[trunposi, ]))
            colnames(trunres) <- c("Estimate", "SE", "CILower", "CIUpper", "CILowerP", "CIUpperP")
            trunres$Part <- "Zero-truncated"
            trunres$Gene <- genelist[k]
            trunres$SNP <- names(trunposi)
            trunres <- trunres[, c(8:9, 1:7)]
          } else {
            trunres <- NULL
          }
        } else {
          trunres <- NULL
        }
        eqtlres <- rbind(logires, trunres)
        if (is.null(eqtlres)) {
          eqtlres <- NA
        }
        
        return (eqtlres)
    }
    
    if (parallel) {
        
        res <- bplapply(1:length(genelist), tempfun, BPPARAM = BPPARAM)
        
    } else {
        
        res <- vector(mode="list", length=length(genelist))
        for (k in 1:length(genelist)) {
            message(sprintf("%s: %s", k, genelist[k]))
            res[[k]] <- tempfun(k)
        }
        
    }
    
    names(res) <- genelist
    
    return (res)
    
}

bphurdle <- function(X, Y, family = "poisson", nboot = 100, sig.level = 0.05,
                     alpha = 1, penalfac = NULL) {
  
    if (is.null(penalfac)) {
      penalfac <- rep(1, ncol(X))
    }
    
    Ylogi <- ifelse(Y > 0, 1, 0)
    Xlogi <- X
    Ytrun <- Y[Y > 0]
    Xtrun <- X[Y > 0, ]
    
    quant <- qnorm(1 - sig.level / 2)
    
    zeroprop <- mean(Ylogi == 0)
    zerocnt <- sum(Ylogi == 0)
    nonzerocnt <- sum(Ylogi == 1)
    
    if ((zeroprop < 0.001) | (zerocnt < 5)) {
        logi <- NULL
        cvfittrun <- tryCatch({
          cv.glmnet(x = Xtrun, y = Ytrun, family = family, alpha = alpha, nfolds = 5,
                    penalty.factor = penalfac)
        }, error=function(e) {
          cv.glmnet(x = Xtrun, y = Ytrun, family = family, alpha = alpha, nfolds = 5, 
                    lambda=get_lambdas(Xtrun[, penalfac != 0], Ytrun),
                    penalty.factor = penalfac)
        })
        trunfun <- function(data, ids) {
          dat <- data[ids, ]
          fit <- glmnet(x = as.matrix(dat[, -1]), y = dat[, 1], family = family,
                        alpha = alpha,
                        penalty.factor = penalfac)
          as.vector(coef(fit, s=cvfittrun$lambda.min))
        }
        trunres <- tryCatch({
          boot(data=cbind(Ytrun, Xtrun), statistic=trunfun, R=nboot)
        }, error=function(e) NULL)
        if (is.null(trunres)) {
          trun <- NULL
        } else {
          trunres.est <- summary(trunres)$original
          trunres.SE <- summary(trunres)$bootSE
          trunres.lci <- trunres.est - quant * trunres.SE
          trunres.uci <- trunres.est + quant * trunres.SE
          trunres.lcip <- apply(trunres$t, 2, function(x) quantile(x, sig.level/2))
          trunres.ucip <- apply(trunres$t, 2, function(x) quantile(x, 1-sig.level/2))
          trun <- cbind(trunres.est, trunres.SE, trunres.lci, trunres.uci,
                        trunres.lcip, trunres.ucip)
        }
    } else if ((zeroprop > 0.999) | nonzerocnt < 5) {
      logi <- NULL  
      trun <- NULL
    } else {
      cvfitlogi <- tryCatch({
        cv.glmnet(x = Xlogi, y = Ylogi, family = "binomial", alpha = alpha, nfolds = 5,
                  penalty.factor = penalfac)
      }, error=function(e) {
        cv.glmnet(x = Xlogi, y = Ylogi, family = "binomial", alpha = alpha, nfolds = 5, 
                  lambda=get_lambdas(Xlogi[, penalfac != 0], Ylogi),
                  penalty.factor = penalfac)
      })
      logifun <- function(data, ids) {
        dat <- data[ids, ]
        fit <- glmnet(x = as.matrix(dat[, -1]), y = dat[, 1], family = "binomial",
                      alpha = alpha,
                      penalty.factor = penalfac)
        as.vector(coef(fit, s=cvfitlogi$lambda.min))
      }
      logires <- tryCatch({
        boot(data=cbind(Ylogi, Xlogi), statistic=logifun, R=nboot, strata=Ylogi)
      }, error=function(e) {
        NULL
      })
      if (is.null(logires)) {
        logi <- NULL
      } else {
        logires.est <- summary(logires)$original
        logires.SE <- summary(logires)$bootSE
        logires.lci <- logires.est - quant * logires.SE
        logires.uci <- logires.est + quant * logires.SE
        logires.lcip <- apply(logires$t, 2, function(x) quantile(x, sig.level/2))
        logires.ucip <- apply(logires$t, 2, function(x) quantile(x, 1-sig.level/2))
        logi <- cbind(logires.est, logires.SE, logires.lci, logires.uci,
                      logires.lcip, logires.ucip)
      }

      if (length(unique(Ytrun)) == 1) {
        trun <- NULL
      } else {
        cvfittrun <- tryCatch({
          cv.glmnet(x = Xtrun, y = Ytrun, family = family, alpha = alpha, nfolds = 5,
                    penalty.factor = penalfac)
        }, error=function(e) {
          cv.glmnet(x = Xtrun, y = Ytrun, family = family, alpha = alpha, nfolds = 5, 
                    lambda=get_lambdas(Xtrun[, penalfac != 0], Ytrun),
                    penalty.factor = penalfac)
        })
        if (is.null(cvfittrun)) {
          
        }
        
        trunfun <- function(data, ids) {
          dat <- data[ids, ]
          fit <- glmnet(x = as.matrix(dat[, -1]), y = dat[, 1], family = family,
                        alpha = alpha,
                        penalty.factor = penalfac)
          as.vector(coef(fit, s=cvfittrun$lambda.min))
        }
        trunres <- tryCatch({
          boot(data=cbind(Ytrun, Xtrun), statistic=trunfun, R=nboot)
        }, error=function(e) NULL)
        if (is.null(trunres)) {
          trun <- NULL
        } else {
          trunres.est <- summary(trunres)$original
          trunres.SE <- summary(trunres)$bootSE
          trunres.lci <- trunres.est - quant * trunres.SE
          trunres.uci <- trunres.est + quant * trunres.SE
          trunres.lcip <- apply(trunres$t, 2, function(x) quantile(x, sig.level/2))
          trunres.ucip <- apply(trunres$t, 2, function(x) quantile(x, 1-sig.level/2))
          trun <- cbind(trunres.est, trunres.SE, trunres.lci, trunres.uci,
                        trunres.lcip, trunres.ucip)
        }
      }
    }
    
    return(list(logi = logi, trun = trun))
    
}

get_lambdas <- function(X, Y, K=50, epsilon=0.0001) {
  X_sd <- apply(X, 2, function(z) sqrt(sum((z-mean(z))^2)/length(z)))
  nonzeroposi <- which(X_sd != 0)
  newX <- X
  newX[, nonzeroposi] <- scale(X[, nonzeroposi], scale = X_sd[nonzeroposi])
  lambda_max <- max(abs(colSums(newX*Y)))/nrow(newX)
  lambdapath <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon), 
                              length.out = K)), digits = 10)
  lambdapath
}

