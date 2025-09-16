# We define a function to plot the density of the copula with given marginals
copula.plot <- function(fname1, frm1, lnk1, fname2, frm2, lnk2, fname3, frm3, lnk3, name){
	frms <- c(frm1, frm2, frm3, 'log', 'plot')
	vls <- rep(list(NULL), length(frms))
	names(vls) <- frms
	vls['log'] <- TRUE
	vls['plot'] <- TRUE

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
	d.cop1 <- eval(parse(text = paste0('d.', name)))

	# Construct function body by substituting necessary variables
	fun_body <- substitute({
		# Retrieves the arguments of the first marginal
		if (L1 > 0) {
			args1 <- mget(F1)
			names(args1) <- strp.prf(F1)
		} else {
			stop('Problem with the marginal 1 provided')
		}

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
			names(args_cop) <- F3
		} else {
			stop('Problem with the copula provided')
		}

		# Ok - tested
		u1 <- u2 <- seq(0, 1, length.out = 200)
		u1 <- u1[-c(1, 200)]
		u2 <- u2[-c(1, 200)]
		y1 <- do.call(QFUN1, c(list(u1), args1)) # Calls the first inverse cdf
		y2 <- do.call(QFUN2, c(list(u2), args2)) # Calls the second inverse cdf

		# Calls the copula density function
		names(args1) <- F1
		names(args2) <- F2
		z <- outer(y1, y2, FUN = function(x, y){
			do.call(DCOP1, c(list(y1 = x, y2 = y, log = log), args1, args2, args_cop))
		})

		## plots the log-density
		tt <- ifelse(log, 'log density', 'density')
		if(plot){
			filled.contour(x = y1, y = y2, z = z, nlevels = 20,
				plot.title = {
					## ... adds the title
					title(xlab = 'First marginal', ylab = 'Second marginal', main = tt)
				})
		}

		if(!plot){
			return(list(x = y1, y = y2, z = z))
		}
	}, list(
		F1 = frm1, L1 = length(frm1),
		F2 = frm2, L2 = length(frm2),
		F3 = frm3, L3 = length(frm3),
		QFUN1 = qfun1, QFUN2 = qfun2,
		DCOP1 = d.cop1
	))

	# Create function
	plot.fun <- function(){}
	formals(plot.fun) <- vls
	body(plot.fun) <- fun_body

	return(plot.fun)
}
