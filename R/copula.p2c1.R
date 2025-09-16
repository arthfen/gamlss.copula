# We define the conditional CDF of the copula (i.e., Pr(Y <= y, X = x))
copula.p2c1 <- function(fname1, frm1, lnk1, fname2, frm2, lnk2, fname3, frm3, lnk3){
	frms <- c('y1', 'y2', frm1, frm2, frm3)
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

			p <- pnorm((qnorm(u2) - theta * qnorm(u1))/sqrt(1 - theta^2))
		} else if(COP1 == 'PL'){
			# Tested - not implemented in 'copula'
			# Page 164
			if(any(theta <= 0)) stop('Theta parameter out of range')

			num <- (theta - 1) * u1 + 1 - ((theta - 1) + 2) * u2
			den <- ((1 + (theta - 1) * (u1 + u2))^2 - 4 * theta * (theta - 1) * u1 * u2)^(1/2)
			p <- 1/2 - 1/2 * num/den
		} else if(COP1 == 'FR'){
			# OK - tested
			# Page 166
			if(any(theta <= 0)) stop('Theta parameter out of range')

			p <- exp(-theta * u1) * ((1 - exp(-theta)) * (1 - exp(-theta * u2))^(-1) - (1 - exp(-theta * u1)))^(-1)
		} else if(COP1 == 'MTCJ'){
			# Untested - not implemented in 'copula'
			# Page 168
			if(any(theta <= 0)) stop('Theta parameter out of range')

			p <- (1 + u1^theta * (u2^(-theta) -1))^(- 1 - 1/theta)
		} else if(COP1 == 'JOE'){
			# OK - tested
			# Page 170
			if(any(theta <= 1)) stop('Theta parameter out of range')

			p <- (1 + (1 - u2)^theta * (1 - u1)^(-theta) - (1 - u2)^theta)^(-1 + 1/theta) * (1 - (1 - u2)^theta)
		} else if(COP1 == 'GU'){
			# OK - tested
			# Page 172
			if(any(theta <= 1)) stop('Theta parameter out of range')

			x <- -log(u1)
			y <- -log(u2)
			p <- u1^(-1) * exp(-(x^theta + y^theta)^(1/theta)) * (1 + (y/x)^theta)^(1/theta - 1)
		} else if(COP1 == 'GAL'){
			# Untested - not implemented in 'copula'
			# Page 174
			if(any(theta <= 0)) stop('Theta parameter out of range')

			x <- -log(u1)
			y <- -log(u2)
			p <- u2 * exp((x^(-theta) + y^(-theta))^(-1/theta)) * (1 - (1 + (x/y)^theta)^(-1 - 1/theta))
		} else if(COP1 == 'HR'){
			# OK - tested
			# Page 176
			if(any(theta <= 0)) stop('Theta parameter out of range')

			x <- -log(u1)
			y <- -log(u2)

			st1 <- pnorm(theta^(-1) + 1/2 * theta * log(x/y))
			st2 <- pnorm(theta^(-1) + 1/2 * theta * log(y/x))
			
			t1 <- exp(- x * st1 - y * st2)
			p <- t1 * u1^(-1) * st1
		} else if(COP1 == 'FGM'){
			# OK - tested
			# Page 213
			if(any(theta < -1 | theta > 1)) stop('Theta parameter out of range')

			p <- u2 + theta * u2 * (1 - u2) * (1 - 2 * u1)
		} else if(COP1 == 'T'){
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

			p <- pt((y - theta * x)/sqrt((1 - theta^2) * (kappa + x^2)/(kappa + 1)), df = kappa + 1)
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

			p <- (1 + (x + y)^(1/kappa))^(-1/theta-1) * (x + y)^(1/kappa-1) * x^(1-1/kappa) * u1^(-theta-1)
		}

		return(p)
	}, list(
		F1 = frm1, L1 = length(frm1),
		F2 = frm2, L2 = length(frm2),
		F3 = frm3, L3 = length(frm3),
		PFUN1 = pfun1, PFUN2 = pfun2,
		COP1 = cop1
	))

	# Create function
	p.fun <- function(){}
	formals(p.fun) <- vls
	body(p.fun) <- fun_body

	return(p.fun)
}
