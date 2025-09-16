# We define the inverse of the conditional CDF of the copula (i.e., y: Pr(Y <= y, X = x) = p1)
copula.q2c1 <- function(fname1, frm1, lnk1, fname2, frm2, lnk2, fname3, frm3, lnk3){
	frms <- c('p2', 'y1', frm1, frm2, frm3)
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

		# Retrieve the arguments for the copula
		if (L3 > 0) {
			args_cop <- mget(F3)
			names(args_cop) <- strp.prf(F3)
		} else {
			stop('Problem with the copula provided')
		}

		# All formulas come from ISBN 978-1-4665-8323-8, except when mentioned otherwise
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

			u2 <- pnorm(qnorm(p2) * sqrt(1 - theta^2) + theta * qnorm(u1))
		} else if(COP1 == 'PL'){
			# Tested - not implemented in 'copula'
			# Taken from: https://civilengineerkey.com/6-plackett-copula/
			if(any(theta <= 0)) stop('Theta parameter out of range')

			S <- p2 * (1 - p2)
			b <- theta + S * (theta - 1)^2
			c <- 2 * S * (u1 * theta^2 + 1 - u1) + theta * (1 - 2 * S)
			d <- theta^(1/2) * (theta + 4 * S * u1 * (1 - u1) * (1 - theta)^2)^(1/2)
			u2 <- (c - (1 - 2 * p2) * d)/(2 * b)
		} else if(COP1 == 'FR'){
			# OK - tested
			# Page 166
			if(any(theta <= 0)) stop('Theta parameter out of range')

			u2 <- - theta^(-1) * log(1 - (1 - exp(-theta))/((p2^(-1) - 1) * exp(-theta * u1) + 1))
		} else if(COP1 == 'MTCJ'){
			# Tested - Page 168
			if(any(theta <= 0)) stop('Theta parameter out of range')

			u2 <- ((p2^(-theta/(1+theta)) - 1) * u1^(-theta) + 1)^(-1/theta)
		} else if(COP1 == 'JOE'){
			# OK - tested
			# Page 170
			if(any(theta <= 1)) stop('Theta parameter out of range')

			f <- function(v, u1, theta, p2, i){
				c2_1 <- (1 + (1 - v)^theta[i] * (1 - u1[i])^(-theta[i]) - (1 - v)^theta[i])^(-1 + 1/theta[i]) * (1 - (1 - v)^theta[i])
				return(c2_1 - p2[i])
			}

			u2 <- rep(NA, n)
			for (i in 1:n) {
				sv <- uniroot(f,
					interval = c(0 + 1e-6, 1 - 1e-6),
					tol = 1e-5,
					u1 = u1,
					theta = theta,
					p2 = p2,
					i = i)
				u2[i] <- sv$root
			}
		} else if(COP1 == 'GU'){
			# OK - tested
			# Page 172
			if(any(theta <= 1)) stop('Theta parameter out of range')

			x <- - log(u1)
			f <- function(v, u1, x, theta, p2, i){
				y <- - log(v)
				c2_1 <- u1[i]^(-1) * exp(-(x[i]^theta[i] + y^theta[i])^(1/theta[i])) * (1 + (y/x[i])^theta[i])^(1/theta[i] - 1)
				return(c2_1 - p2[i])
			}

			u2 <- rep(NA, n)
			for (i in 1:n) {
				sv <- uniroot(f,
					interval = c(0 + 1e-6, 1 - 1e-6),
					tol = 1e-5,
					u1 = u1,
					x = x,
					theta = theta,
					p2 = p2,
					i = i)
				u2[i] <- sv$root
			}
		} else if(COP1 == 'GAL'){
			# Untested - not implemented in 'copula'
			# Page 174
			if(any(theta <= 0)) stop('Theta parameter out of range')

			x <- - log(u1)
			f <- function(v, u1, x, theta, p2, i){
				y <- - log(v)
				c2_1 <- v * exp((x[i]^(-theta[i]) + y^(-theta[i]))^(-1/theta[i])) * (1 - (1 + (x[i]/y)^theta[i])^(-1 - 1/theta[i]))
				return(c2_1 - p2[i])
			}

			u2 <- rep(NA, n)
			for (i in 1:n) {
				sv <- uniroot(f,
					interval = c(0 + 1e-6, 1 - 1e-6),
					tol = 1e-5,
					u1 = u1,
					x = x,
					theta = theta,
					p2 = p2,
					i = i)
				u2[i] <- sv$root
			}
		} else if(COP1 == 'HR'){
			# OK - tested
			# Page 176
			if(any(theta <= 0)) stop('Theta parameter out of range')

			x <- -log(u1)

			f <- function(v, u1, x, theta, p2, i){
				y <- -log(v)
				st1 <- pnorm(theta[i]^(-1) + 1/2 * theta[i] * log(x[i]/y))
				st2 <- pnorm(theta[i]^(-1) + 1/2 * theta[i] * log(y/x[i]))
				
				t1 <- exp(- x[i] * st1 - y * st2)
				c2_1 <- t1 * u1[i]^(-1) * st1

				return(c2_1 - p2[i])
			}

			u2 <- rep(NA, n)
			for (i in 1:n) {
				sv <- uniroot(f,
					interval = c(0 + 1e-6, 1 - 1e-6),
					tol = 1e-5,
					u1 = u1,
					x = x,
					theta = theta,
					p2 = p2,
					i = i)
				u2[i] <- sv$root
			}
		} else if(COP1 == 'FGM'){
			# OK - tested
			# Page 213 - I believe there might be a typo, so I rewrote it
			if(any(theta < -1 | theta > 1)) stop('Theta parameter out of range')

			u2 <- (- theta + sqrt((theta * (2 * u1 - 1) - 1)^2 - 4 * p2 * (theta - 2 * theta * u1)) + 2 * theta * u1 -1)/(2 * theta * (2 * u1 - 1))
			nk <- which(u1 == 0.5)
			if(length(nk) > 0) u2[nk] <- p2[nk]
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

			f <- function(v, u1, x, theta, kappa, p2, i){
				y <- qt(v, df = kappa[i])
				c2_1 <- pt((y - theta[i] * x[i])/sqrt((1 - theta[i]^2) * (kappa[i] + x[i]^2)/(kappa[i] + 1)), df = kappa[i] + 1)

				return(c2_1 - p2[i])
			}

			u2 <- rep(NA, n)
			for (i in 1:n) {
				sv <- uniroot(f,
					interval = c(0 + 1e-6, 1 - 1e-6),
					tol = 1e-5,
					u1 = u1,
					x = x,
					theta = theta,
					kappa = kappa,
					p2 = p2,
					i = i)
				u2[i] <- sv$root
			}
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

			f <- function(v, u1, x, theta, kappa, p2, i){
				y <- (v^(-theta[i]) - 1)^kappa[i]

				c2_1 <- (1 + (x[i] + y)^(1/kappa[i]))^(-1/theta[i]-1) * (x[i] + y)^(1/kappa[i]-1) * x[i]^(1-1/kappa[i]) * u1[i]^(-theta[i]-1)

				return(c2_1 - p2[i])
			}

			u2 <- rep(NA, n)
			for (i in 1:n) {
				sv <- uniroot(f,
					interval = c(0 + 1e-6, 1 - 1e-6),
					tol = 1e-5,
					u1 = u1,
					x = x,
					theta = theta,
					kappa = kappa,
					p2 = p2,
					i = i)
				u2[i] <- sv$root
			}
		}
		
		y <- do.call(QFUN2, c(list(p = u2), args2))

		return(y)
	}, list(
		F1 = frm1, L1 = length(frm1),
		F2 = frm2, L2 = length(frm2),
		F3 = frm3, L3 = length(frm3),
		PFUN1 = pfun1, QFUN2 = qfun2,
		COP1 = cop1
	))

	# Create function
	p.fun <- function(){}
	formals(p.fun) <- vls
	body(p.fun) <- fun_body

	return(p.fun)
}

