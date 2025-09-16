# We define the fitting function for gamlss2
copula.fit <- function(fname1, frm1, lnk1, fname2, frm2, lnk2, fname3, frm3, lnk3, name){
	# We get the relevant functions for the first marginal
	dfun1 <- eval(parse(text = paste0('d', fname1)))
	pfun1 <- eval(parse(text = paste0('p', fname1)))
	qfun1 <- eval(parse(text = paste0('q', fname1)))

	# We get the relevant functions for the second marginal
	dfun2 <- eval(parse(text = paste0('d', fname2)))
	pfun2 <- eval(parse(text = paste0('p', fname2)))
	qfun2 <- eval(parse(text = paste0('q', fname2)))

	# We get the density of the copula
	dcop1 <- eval(parse(text = paste0('d.', name)))
	rcop1 <- eval(parse(text = paste0('r.', name)))
	p2c1cop1 <- eval(parse(text = paste0('p2c1.', name)))

	lnks <- c(lnk1, lnk2, lnk3)
	frms <- c(frm1, frm2, frm3)
	names(lnks) <- frms

	# Construct function body by substituting necessary variables
	fun_body <- substitute({
		fam <- list(
			"family" = paste0(fname1, ',', fname2, ',', fname3),
			"names" = frms,
			"links" = lnks,
			"pdf" = function(y, par, log = FALSE) {
				# Some sanity checks
				if(length(y) %% 2 != 0){
					stop('Problem with the number of observations: we should have an even number')
				}

				set.seed(1) # This helps for randomized residuals and does not affect non-randomized residuals
				n <- length(y)/2

				# Retrieves the arguments of the first marginal
				if (L1 > 0) {
					args1 <- lapply(F1, function(i) par[[i]][1:n])
					names(args1) <- F1
				} else {
					stop('Problem with the marginal 1 provided')
				}
				y1 <- list(y1 = y[1:n]) # Observed value of the first marginal

				# Retrieves the arguments of the second marginal
				if (L2 > 0) {
					args2 <- lapply(F2, function(i) par[[i]][1:n])
					names(args2) <- F2
				} else {
					stop('Problem with the marginal 2 provided')
				}
				y2 <- list(y2 = y[n + (1:n)]) # Observed value of the second marginal

				# Retrieve the arguments for the copula
				if (L3 > 0) {
					args_cop <- lapply(F3, function(i) par[[i]][1:n])
					names(args_cop) <- F3
				} else {
					stop('Problem with the copula provided')
				}

				# Compute marginal densities
				args.all <- c(y1, y2, args1, args2, args_cop, list(log = FALSE))
				dk <- do.call(DCOP1, args.all)

				if(log) dk <- log(dk)
				dk <- rep(dk/2, 2)

				return(dk)
			}, "residuals" = function(object, type = 'quantile', id = 5){
				# id: 
				## 0 = all
				## 1 = marginal 1
				## 2 = marginal 2
				## 3 = marginal 2 given marginal 1
				## 4 = chi-square distribution with df=2 induced by p = res(id=1)^2 + res(id=3)^2
				## 5 = joint residuals, according to Kalliovirta (2008) "Quantile Residuals for Multivariate Models", Equation 5

				# Some sanity checks
				y <- object$y
				if(length(y) %% 2 != 0){
					stop('Problem with the number of observations: we should have an even number')
				}

				n <- length(y)/2
				pr <- predict(object)[1:n,,drop = FALSE]

				# Retrieves the arguments of the first marginal
				if (L1 > 0) {
					args1 <- pr[,F1,drop = FALSE]
					names(args1) <- strp.prf(F1)

				} else {
					stop('Problem with the marginal 1 provided')
				}
				y1 <- y[1:n] # Observed value of the first marginal
				r.id1 <- qnorm(do.call(PFUN1, c(list(q = y1), args1))) # res(id = 1)

				# Retrieves the arguments of the second marginal
				if (L2 > 0) {
					args2 <- pr[,F2,drop = FALSE]
					names(args2) <- strp.prf(F2)
				} else {
					stop('Problem with the marginal 2 provided')
				}
				y2 <- y[n + (1:n)] # Observed value of the second marginal
				r.id2 <- qnorm(do.call(PFUN2, c(list(q = y2), args2))) # res(id = 1)

				# Retrieve the arguments for the copula
				if (L3 > 0) {
					args_cop <- pr[,F3,drop = FALSE]
					names(args_cop) <- F3
				} else {
					stop('Problem with the copula provided')
				}

				# Compute marginal densities
				names(args1) <- F1
				names(args2) <- F2
				args.all <- c(list(y1 = y1, y2 = y2), args1, args2, args_cop)
				r.id3 <- qnorm(do.call(P2C1COP1, args.all)) # res(id = 3)

				r.id4 <- qnorm(pchisq(r.id1^2 + r.id3^2, df = 2)) # res(id = 4)

				X.id5 <- pnorm(r.id1) * pnorm(r.id3)
				r.id5 <- qnorm(X.id5 * (1 - log(X.id5))) # res(id == 5)

				if(id == 0){
					res <- cbind(r.id1, r.id2, r.id3, r.id4, r.id5) # res(id == 0)
					colnames(res) <- c('qres_m1', 'qres_m2', 'qres_m2c1', 'qres_chisq', 'qres_joint')
				} else if(id == 1){
					res <- r.id1
				} else if(id == 2){
					res <- r.id2
				} else if(id == 3){
					res <- r.id3
				} else if(id == 4){
					res <- r.id4
				} else if(id == 5){
					res <- r.id5
				}

				return(res)
			}, "quantile" = function(p, par, id = 0){
				# id: 
				## 0 = all
				## 1 = marginal 1
				## 2 = marginal 2

				# Retrieves the arguments of the first marginal
				if (L1 > 0) {
					args1 <- lapply(F1, function(i) par[[i]][1:n])
					names(args1) <- strp.prf(F1)
				} else {
					stop('Problem with the marginal 1 provided')
				}
				q.id1 <- do.call(QFUN1, c(list(p = p), args1)) # qtl(id = 1)

				# Retrieves the arguments of the second marginal
				if (L2 > 0) {
					args2 <- lapply(F2, function(i) par[[i]][1:n])
					names(args2) <- strp.prf(F2)
				} else {
					stop('Problem with the marginal 2 provided')
				}
				q.id2 <- do.call(QFUN2, c(list(p = p), args2)) # qtl(id = 2)
				
				if(id == 0){
					res <- c(q.id1, q.id2)
				} else if(id == 1){
					res <- q.id1
				} else if(id == 2){
					res <- q.id2
				}
				
				return(res)
			}, "random" = function(n, par, id = 0){
				# id:
				## 0 = joint
				## 1 = marginal 1
				## 2 = marginal 2

				# Retrieves the arguments of the first marginal
				if (L1 > 0) {
					args1 <- par[F1]
					names(args1) <- strp.prf(F1)
				} else {
					stop('Problem with the marginal 1 provided')
				}

				# Retrieves the arguments of the second marginal
				if (L2 > 0) {
					args2 <- par[F2]
					names(args2) <- strp.prf(F2)
				} else {
					stop('Problem with the marginal 2 provided')
				}

				# Retrieve the arguments for the copula
				if (L3 > 0) {
					args_cop <- par[F3]
					names(args_cop) <- F3
				} else {
					stop('Problem with the copula provided')
				}

				if(id == 0){
					# Compute marginal densities
					names(args1) <- F1
					names(args2) <- F2
					args.all <- c(list(n = n), args1, args2, args_cop)

					res <- do.call(RCOP1, args.all) # rand(id = 0)
				} else if(id == 1){
					p <- runif(n)
					res <- do.call(QFUN1, c(list(p = p), args1)) # rand(id = 1)
				} else if(id == 2){
					p <- runif(n)
					res <- do.call(QFUN2, c(list(p = p), args2)) # rand(id = 2)
				}

				return(res)
			}
		)
		class(fam) <- "gamlss2.family"
		return(fam)

		return(p)
	}, list(
		F1 = frm1, L1 = length(frm1),
		F2 = frm2, L2 = length(frm2),
		F3 = frm3, L3 = length(frm3),
		PFUN1 = pfun1, PFUN2 = pfun2,
		QFUN1 = qfun1, QFUN2 = qfun2,
		DCOP1 = dcop1, RCOP1 = rcop1,
		P2C1COP1 = p2c1cop1
	))

	# Create function
	fit.fun <- function(...){}
	body(fit.fun) <- fun_body

	return(fit.fun)
}
