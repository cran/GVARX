#AAA=vars::VAR(data.frame(y=rnorm(10),x=rnorm(10)),p=1,type="none")

.GVAR0 <- function (y, p=1, exogen=NULL,type=c("const","trend","both","none"),season = NULL, lag.max = NULL, ic = c("AIC", "HQ", "SC", "FPE"))
{
type=type
exogen=exogen
p=p
  y <- as.matrix(y)
    if (any(is.na(y)))
        stop("\nNAs in y.\n")
    if (ncol(y) < 2)
        stop("The matrix 'y' should contain at least two variables. For univariate analysis consider ar() and arima() in package stats.\n")
    if (is.null(colnames(y))) {
        colnames(y) <- paste("y", 1:ncol(y), sep = "")
        warning(paste("No column names supplied in y, using:",
            paste(colnames(y), collapse = ", "), ", instead.\n"))
    }
    colnames(y) <- make.names(colnames(y))
    y.orig <- y
    type <- match.arg(type)
    obs <- dim(y)[1]
    K <- dim(y)[2]
    if(!is.null(lag.max)){
      lag.max <- abs(as.integer(lag.max))
      ic <- paste(match.arg(ic), "(n)", sep = "")
      p <- .VARselect(y, lag.max = lag.max, type = type, season = season, exogen = exogen)$selection[ic]
    }
    sample <- obs - p
    ylags <- embed(y, dimension = p + 1)[, -(1:K)]
    temp1 <- NULL
    for (i in 1:p) {
        temp <- paste(colnames(y), ".L", i, sep = "")
        temp1 <- c(temp1, temp)
    }
    colnames(ylags) <- temp1
    yend <- y[-c(1:p), ]
    if (type == "const") {
        rhs <- cbind(ylags, rep(1, sample))
        colnames(rhs) <- c(colnames(ylags), "const")
    }
    else if (type == "trend") {
        rhs <- cbind(ylags, seq(p + 1, length = sample))
        colnames(rhs) <- c(colnames(ylags), "trend")
    }
    else if (type == "both") {
        rhs <- cbind(ylags, rep(1, sample), seq(p + 1, length = sample))
        colnames(rhs) <- c(colnames(ylags), "const", "trend")
    }
    else if (type == "none") {
        rhs <- ylags
        colnames(rhs) <- colnames(ylags)
    }
    if (!(is.null(season))) {
        season <- abs(as.integer(season))
        dum <- (diag(season) - 1/season)[, -season]
        dums <- dum
        while (nrow(dums) < obs) {
            dums <- rbind(dums, dum)
        }
        dums <- dums[1:obs, ]
        colnames(dums) <- paste("sd", 1:ncol(dums), sep = "")
        rhs <- cbind(rhs, dums[-c(1:p), ])
    }
    if (!(is.null(exogen))) {
        exogen <- as.matrix(exogen)
        if (!identical(nrow(exogen), nrow(y))) {
            stop("\nDifferent row size of y and exogen.\n")
        }
        if (is.null(colnames(exogen))) {
            colnames(exogen) <- paste("exo", 1:ncol(exogen),
                sep = "")
            warning(paste("No column names supplied in exogen, using:",
                paste(colnames(exogen), collapse = ", "), ", instead.\n"))
        }
        colnames(exogen) <- make.names(colnames(exogen))
        tmp <- colnames(rhs)
        rhs <- cbind(rhs, exogen[-c(1:p), ])
        colnames(rhs) <- c(tmp, colnames(exogen))
    }
    datamat <- as.data.frame(rhs)
    colnames(datamat) <- colnames(rhs)
    equation <- list()
    NW<-list()
    for (i in 1:K) {
        y <- yend[, i]
        equation[[colnames(yend)[i]]] <- lm(y ~ -1 + ., data = datamat)
        NW[[colnames(yend)[i]]]<- sandwich::NeweyWest(lm(y ~ -1 + ., data = datamat))
        if(any(c("const", "both") %in% type)){
         attr(equation[[colnames(yend)[i]]]$terms, "intercept") <- 1
        }
    }
  call <- match.call()
  if("season" %in% names(call)) call$season <- eval(season)
result <- list(varresult = equation, NWHAC=NW, datamat = data.frame(cbind(yend, rhs)), y = y.orig, type = type, p = p, K = K, obs = sample, totobs = sample + p, restrictions = NULL, call = call)
    class(result) <- "varest"
    return(result)
}



.VARselect <-
function (y, lag.max = 10, type = c("const", "trend", "both",
    "none"), season = NULL, exogen = NULL)
{
    y <- as.matrix(y)
    if (any(is.na(y)))
        stop("\nNAs in y.\n")
    colnames(y) <- make.names(colnames(y))
    K <- ncol(y)
    lag.max <- abs(as.integer(lag.max))
    type <- match.arg(type)
    lag <- abs(as.integer(lag.max + 1))
    ylagged <- embed(y, lag)[, -c(1:K)]
    yendog <- y[-c(1:lag.max), ]
    sample <- nrow(ylagged)
    rhs <- switch(type, const = rep(1, sample), trend = seq(lag.max + 1,
        length = sample), both = cbind(rep(1, sample), seq(lag.max + 1, length = sample)), none = NULL)
    if (!(is.null(season))) {
        season <- abs(as.integer(season))
        dum <- (diag(season) - 1/season)[, -season]
        dums <- dum
        while (nrow(dums) < sample) {
            dums <- rbind(dums, dum)
        }
        dums <- dums[1:sample, ]
        rhs <- cbind(rhs, dums)
    }
    if (!(is.null(exogen))) {
        exogen <- as.matrix(exogen)
        if (!identical(nrow(exogen), nrow(y))) {
            stop("\nDifferent row size of y and exogen.\n")
        }
        if (is.null(colnames(exogen))) {
            colnames(exogen) <- paste("exo", 1:ncol(exogen),
                sep = "")
            warning(paste("No column names supplied in exogen, using:",
                paste(colnames(exogen), collapse = ", "), ", instead.\n"))
        }
        colnames(exogen) <- make.names(colnames(exogen))
        rhs <- cbind(rhs, exogen[-c(1:lag.max), ])
    }
    idx <- seq(K, K * lag.max, K)
    if(!is.null(rhs)){
      detint <- ncol(as.matrix(rhs))
    } else {
      detint <- 0
    }
    criteria <- matrix(NA, nrow = 4, ncol = lag.max)
    rownames(criteria) <- c("AIC(n)", "HQ(n)", "SC(n)", "FPE(n)")
    colnames(criteria) <- paste(seq(1:lag.max))
    for (i in 1:lag.max) {
        ys.lagged <- cbind(ylagged[, c(1:idx[i])], rhs)
        sampletot <- nrow(y)
        nstar <- ncol(ys.lagged)
        resids <- resid(lm(yendog ~ -1 + ys.lagged))
        sigma.det <- det(crossprod(resids)/sample)
        criteria[1, i] <- log(sigma.det) + (2/sample) * (i * K^2 + K * detint)
        criteria[2, i] <- log(sigma.det) + (2 * log(log(sample))/sample) * (i * K^2 + K * detint)
        criteria[3, i] <- log(sigma.det) + (log(sample)/sample) * (i * K^2 + K * detint)
        criteria[4, i] <- ((sample + nstar)/(sample - nstar))^K * sigma.det
    }
    order <- apply(criteria, 1, which.min)
    return(list(selection = order, criteria = criteria))
}

