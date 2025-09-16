gen.copula <- function(m1 = "NO",
	m2 = "NO", 
	copula = "NO",
	name = NULL,
	print = TRUE){
		# Checks if the copula is currently supported, currently:
		## NO: Gaussian (normal) copula
		## PL: Plackett copula
		## FR: Frank copula
		## MTCJ: Mardia-Takahasi-Clayton-Cook-Johnson copula
		## JOE: Joe copula
		## GU: Gumbel copula
		## GAL: Galambos copula
		## HR: Husler-Reiss copula
		## FGM: Farlie-Gumbel-Morgenstern copula
		## T: Student copula
		## BB1: BB1 copula
		spt <- c('NO', 'PL', 'FR', 'MTCJ', 'JOE', 'GU', 'GAL', 'HR', 'FGM', 'T', 'BB1')
		if(!(copula %in% spt)){
			stop('The copula provided is currently not supported')
		}

		# If a name is not provided, we set a default name
		if(is.null(name)){
			name <- paste0(m1, '_', m2, '_', copula, 'cop')
		}

		# We get the first marginal
		fam1  <- as.gamlss.family(m1)
		np1 <- fam1$nopar
		fname1 <- fam1$family[[1]]
		frm1 <- names(fam1$parameters)[unlist(fam1$parameters)]
		lnk1 <- sapply(paste0(frm1, '.link'), function(i) fam1[[i]])
		frm1 <- paste0('m1.', frm1)

		# We get the second marginal
		fam2 <- as.gamlss.family(m2)
		np2 <- fam2$nopar
		fname2 <- fam2$family[[1]]
		frm2 <- names(fam2$parameters)[unlist(fam2$parameters)]
		lnk2 <- sapply(paste0(frm2, '.link'), function(i) fam2[[i]])
		frm2 <- paste0('m2.', frm2)

		# We define the link functions and number of parameters of the copula
		fname3 <- copula
		if(fname3 == 'NO'){
			lnk3 <- '[-1,1]'
		} else if(fname3 == 'PL'){
			lnk3 <- 'log'
		} else if(fname3 == 'FR'){
			lnk3 <- 'log'
		} else if(fname3 == 'MTCJ'){
			lnk3 <- 'log'
		} else if(fname3 == 'JOE'){
			lnk3 <- 'logshiftto1'
		} else if(fname3 == 'GU'){
			lnk3 <- 'logshiftto1'
		} else if(fname3 == 'GAL'){
			lnk3 <- 'log'
		} else if(fname3 == 'HR'){
			lnk3 <- 'log'
		} else if(fname3 == 'FGM'){
			lnk3 <- '[-1,1]'
		} else if(fname3 == 'T'){
			lnk3 <- c('[-1,1]', 'log')
		} else if(fname3 == 'BB1'){
			lnk3 <- c('log', 'logshiftto1')
		}
		np3 <- length(lnk3)
		frm3 <- c('theta', 'kappa')[1:np3]

		# Generates d
		dummy <- copula.d(fname1, frm1, lnk1, fname2, frm2, lnk2, fname3, frm3, lnk3)
		dfun <- paste0('d.', name)
		eval(call('<-', as.name(dfun), dummy), envir = parent.frame(n = 1))

		# Generates p2c1
		dummy <- copula.p2c1(fname1, frm1, lnk1, fname2, frm2, lnk2, fname3, frm3, lnk3)
		p2c1fun <- paste0('p2c1.', name)
		eval(call('<-', as.name(p2c1fun), dummy), envir = parent.frame(n = 1))

		# Generates q2c1
		dummy <- copula.q2c1(fname1, frm1, lnk1, fname2, frm2, lnk2, fname3, frm3, lnk3)
		q2c1fun <- paste0('q2c1.', name)
		eval(call('<-', as.name(q2c1fun), dummy), envir = parent.frame(n = 1))

		# Generates p
		dummy <- copula.p(fname1, frm1, lnk1, fname2, frm2, lnk2, fname3, frm3, lnk3)
		pfun <- paste0('p.', name)
		eval(call('<-', as.name(pfun), dummy), envir = parent.frame(n = 1))

		# Generates r
		dummy <- copula.r(fname1, frm1, lnk1, fname2, frm2, lnk2, fname3, frm3, lnk3, name)
		rfun <- paste0('r.', name)
		eval(call('<-', as.name(rfun), dummy), envir = parent.frame(n = 1))

		# Generates plot
		dummy <- copula.plot(fname1, frm1, lnk1, fname2, frm2, lnk2, fname3, frm3, lnk3, name)
		plotfun <- paste0('plot.', name)
		eval(call('<-', as.name(plotfun), dummy), envir = parent.frame(n = 1))

		# Generates the fitting distribution
		dummy <- copula.fit(fname1, frm1, lnk1, fname2, frm2, lnk2, fname3, frm3, lnk3, name)
		fun <- name
		eval(call('<-', as.name(fun), dummy), envir = parent.frame(n = 1))

		alldislist <- c(dfun, p2c1fun, q2c1fun, pfun, rfun, plotfun)

		if (print){
			cat("A ", copula, " copula of ",  fname1, " and ", fname2, "distributions has been generated \n", "and saved under the names: ", "\n",paste(alldislist, sep=","),"\n")
		}  
}
