# We define the joint density
copula.d <- function(fname1, frm1, lnk1, fname2, frm2, lnk2, fname3, frm3, lnk3){
	frms <- c('y1', 'y2', frm1, frm2, frm3, 'log')
	vls <- rep(list(NULL), length(frms))
	names(vls) <- frms
	vls['log'] <- FALSE

	# We get the relevant functions for the first marginal
	dfun1 <- eval(parse(text = paste0('d', fname1)))
	pfun1 <- eval(parse(text = paste0('p', fname1)))
	qfun1 <- eval(parse(text = paste0('q', fname1)))

	# We get the relevant functions for the second marginal
	dfun2 <- eval(parse(text = paste0('d', fname2)))
	pfun2 <- eval(parse(text = paste0('p', fname2)))
	qfun2 <- eval(parse(text = paste0('q', fname2)))

	# We get the constructor of the copula
	cop1 <- fname3

	# Construct function body by substituting necessary variables
	fun_body <- substitute({
		# Retrieves the arguments of the first marginal
		if (L1 > 0) {
			args1 <- mget(F1)
			names(args1) <- strp.prf(F1)
		} else {
			stop('Problem with the marginal 1 provided')
		}
		u1 <- do.call(PFUN1, c(list(y1), args1)) # Calls the first cdf

		# Retrieves the arguments of the second marginal
		if (L2 > 0) {
			args2 <- mget(F2)
			names(args2) <- strp.prf(F2)
		} else {
			stop('Problem with the marginal 2 provided')
		}
		u2 <- do.call(PFUN2, c(list(y2), args2)) # Calls the second cdf

		# Retrieve the arguments for the copula
		if (L3 > 0) {
			args_cop <- mget(F3)
			names(args_cop) <- strp.prf(F3)
		} else {
			stop('Problem with the copula provided')
		}

		# All formulas come from ISBN 978-1-4665-8323-8
		n <- length(u1)
		theta <- args_cop$theta
		if(length(theta) > 1 & length(theta) < n){
			stop('The copula parameter theta has the wrong length.')
		}
		if(length(theta) == 1){
			theta <- rep(theta, n)
		}

		if(COP1 == 'NO'){
			# OK - tested
			# Page 163
			if(any(theta < -1 | theta > 1)) stop('Theta parameter out of range')

			x <- qnorm(u1)
			y <- qnorm(u2)
			cden <- 1/sqrt(1 - theta^2) * exp(- (theta^2 * (x^2 + y^2) - 2 * theta * x * y)/(2 * (1 - theta^2)))
		} else if(COP1 == 'PL'){
			# OK - tested
			# Page 164
			if(any(theta <= 0)) stop('Theta parameter out of range')

			num <- theta * (1 + (theta - 1) * (u1 + u2 - 2 * u1 * u2))
			den <- ((1 + (theta - 1) * (u1 + u2))^2 - 4 * theta * (theta - 1) * u1 * u2)^1.5
			cden <- num/den
		} else if(COP1 == 'FR'){
			# OK - tested
			# Page 165
			if(any(theta <= 0)) stop('Theta parameter out of range')

			num <- theta * (1 - exp(- theta)) * exp(- theta * (u1 + u2))
			den <- (1 - exp(- theta) - (1 - exp(- theta * u1))*(1 - exp(- theta * u2)))^2
			cden <- num/den
		} else if(COP1 == 'MTCJ'){
			# OK - tested
			# Page 168
			if(any(theta <= 0)) stop('Theta parameter out of range')

			t1 <- (1 + theta)
			t2 <- (u1 * u2)^(- theta - 1)
			t3 <- (u1^(-theta) + u2^(- theta) - 1)^(-2 - 1/theta)
			cden <- t1 * t2 * t3
		} else if(COP1 == 'JOE'){
			# OK - tested
			# Page 170
			if(any(theta <= 1)) stop('Theta parameter out of range')

			u. <- 1 - u1
			v. <- 1 - u2
			
			t1 <- (u.^theta + v.^theta - u.^theta * v.^theta)^(-2 + 1/theta)
			t2 <- u.^(theta - 1)
			t3 <- v.^(theta-1)
			t4 <- (theta - 1 + u.^theta + v.^theta - u.^theta * v.^theta)
			cden <- t1 * t2 * t3 * t4
		} else if(COP1 == 'GU'){
			# OK - tested
			# Page 172
			if(any(theta <= 1)) stop('Theta parameter out of range')

			x <- - log(u1)
			y <- - log(u2)
			t1 <- exp(- (x^theta + y^theta)^(1/theta) )
			t2 <- ((x^theta + y^theta)^(1/theta) + theta - 1)
			t3 <- (x^theta + y^theta)^(1/theta - 2)
			t4 <- (x*y)^(theta - 1)
			t5 <- (u1*u2)^(-1)
			cden <- t1 * t2 * t3 * t4 * t5
		} else if(COP1 == 'GAL'){
			# OK - tested
			# Page 174
			if(any(theta <= 0)) stop('Theta parameter out of range')

			x <- - log(u1)
			y <- - log(u2)
			t1 <- exp((x^(-theta) + y^(-theta))^(-1/theta))
			t2 <- (1 - (x^(-theta) + y^(-theta))^(-1 - 1/theta) * (x^(-theta-1) + y^(-theta-1)) + (x^(-theta) + y^(-theta))^(-2 - 1/theta) * (x*y)^(-theta - 1) * (1 + theta + (x^(-theta) + y^(-theta))^(-1/theta)))
			cden <- t1 * t2
		} else if(COP1 == 'HR'){
			# OK - tested
			# Page 213
			if(any(theta <= 0)) stop('Theta parameter out of range')

			x <- - log(u1)
			y <- - log(u2)

			st1 <- pnorm(theta^(-1) + 1/2 * theta * log(x/y))
			st2 <- pnorm(theta^(-1) + 1/2 * theta * log(y/x))

			t1 <- exp(- x * st1 - y * st2)/(u1 * u2)
			t2 <- st2 * st1 + 1/2 * theta * y^(-1) * dnorm(qnorm(st1))

			cden <- t1 * t2
		} else if(COP1 == 'FGM'){
			# OK - tested
			# Page 176
			if(any(theta < -1 | theta > 1)) stop('Theta parameter out of range')

			cden <- 1 + theta * (1 - 2 * u1) * (1 - 2 * u2)
		} else if(COP1 == 'T'){
			# OK - tested
			# Page 181
			kappa <- args_cop$kappa
			if(length(kappa) > 1 & length(kappa) < n){
				stop('The copula parameter kappa has the wrong length.')
			}
			if(length(kappa) == 1){
				kappa <- rep(kappa, n)
			}

			if(any(theta < -1 | theta > 1)) stop('Theta parameter out of range')
			if(any(kappa <= 0)) stop('Theta parameter out of range')

			x <- qt(u1, df = kappa)
			y <- qt(u2, df = kappa)

			t1 <- 1/sqrt(1 - theta^2)
			t2 <- gamma((kappa + 2)/2) * gamma(kappa/2)/gamma((kappa + 1)/2)^2
			t3.num <- (1 + (x^2 + y^2 - 2 * theta * x * y)/(kappa * (1 - theta^2)))^(- (kappa + 2)/2)
			t3.den <- (1 + x^2/kappa)^(-(kappa+1)/2) * (1 + y^2/kappa)^(-(kappa+1)/2)

			cden <- t1 * t2 * t3.num/t3.den
		} else if(COP1 == 'BB1'){
			# OK - tested
			# Page 191
			kappa <- args_cop$kappa
			if(length(kappa) > 1 & length(kappa) < n){
				stop('The copula parameter kappa has the wrong length.')
			}
			if(length(kappa) == 1){
				kappa <- rep(kappa, n)
			}

			if(any(theta <= 0)) stop('Theta parameter out of range')
			if(any(kappa <= 1)) stop('Theta parameter out of range')

			x <- (u1^(-theta) - 1)^kappa
			y <- (u2^(-theta) - 1)^kappa

			t1 <- (1 + (x+y)^(1/kappa))^(-1/theta - 2)
			t2 <- (x + y)^(1/kappa - 2)
			t3 <- (theta * (kappa - 1) + (theta*kappa + 1)*(x+y)^(1/kappa))
			t4 <- (x*y)^(1 - 1/kappa)
			t5 <- (u1*u2)^(-theta-1)
			cden <- t1 * t2 * t3 * t4 * t5
		}
		
		# Compute marginal densities
		f1 <- do.call(DFUN1, c(list(y1), args1))
		f2 <- do.call(DFUN2, c(list(y2), args2))

		# Computes the joint density (i.e., copula density * marginal densities)
		res <- cden * f1 * f2
		if (log) res <- log(res)
		return(res)
	}, list(
		F1 = frm1, L1 = length(frm1),
		F2 = frm2, L2 = length(frm2),
		F3 = frm3, L3 = length(frm3),
		PFUN1 = pfun1, PFUN2 = pfun2,
		DFUN1 = dfun1, DFUN2 = dfun2,
		COP1 = cop1
	))

	# Create function
	d.fun <- function(){}
	formals(d.fun) <- vls
	body(d.fun) <- fun_body

	return(d.fun)
}
