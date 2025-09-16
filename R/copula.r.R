# We define the random number generator
copula.r <- function(fname1, frm1, lnk1, fname2, frm2, lnk2, fname3, frm3, lnk3, name){
	frms <- c('n', frm1, frm2, frm3)
	vls <- rep(list(NULL), length(frms))
	names(vls) <- frms

	# We get the relevant functions for the first marginal
	dfun1 <- eval(parse(text = paste0('d', fname1)))
	pfun1 <- eval(parse(text = paste0('p', fname1)))
	qfun1 <- eval(parse(text = paste0('q', fname1)))

	# We get the relevant functions for the second marginal
	dfun2 <- eval(parse(text = paste0('d', fname2)))
	pfun2 <- eval(parse(text = paste0('p', fname2)))
	qfun2 <- eval(parse(text = paste0('q', fname2)))

	# We get the constructor of the copula
	q2c1.cop1 <- eval(parse(text = paste0('q2c1.', name)))

	# Construct function body by substituting necessary variables
	fun_body <- substitute({
		# Retrieves the arguments of the first marginal
		if (L1 > 0) {
			args1 <- mget(F1)
			names(args1) <- F1
		} else {
			stop('Problem with the marginal 1 provided')
		}

		# Retrieves the arguments of the second marginal
		if (L2 > 0) {
			args2 <- mget(F2)
			names(args2) <- F2
		} else {
			stop('Problem with the marginal 2 provided')
		}

		# Retrieve the arguments for the copula
		if (L3 > 0) {
			args_cop <- mget(F3)
			names(args_cop) <- F3
		} else {
			stop('Problem with the copula provided')
		}

		# Ok - tested
		u1 <- runif(n)
		w2 <- runif(n)

		## Uses the inverse cdf of the first marginal density
		args1.q <- args1
		names(args1.q) <- strp.prf(F1)
		y1 <- do.call(QFUN1, c(list(u1), args1.q))

		## Computes the inverse conditional cdf of the copula
		y2 <- do.call(QCOP1, c(list(p2 = w2, y1 = y1), args1, args2, args_cop))

		rn <- cbind(y1, y2)
		return(rn)
	}, list(
		F1 = frm1, L1 = length(frm1),
		F2 = frm2, L2 = length(frm2),
		F3 = frm3, L3 = length(frm3),
		QFUN1 = qfun1, QFUN2 = qfun2,
		QCOP1 = q2c1.cop1
	))

	# Create function
	p.fun <- function(){}
	formals(p.fun) <- vls
	body(p.fun) <- fun_body

	return(p.fun)
}

