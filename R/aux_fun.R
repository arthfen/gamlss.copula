# Auxiliary function to remove the prefix from a string vector
strp.prf <- function(x){
	gsub("m2\\.", "", gsub("m1\\.", "", x))
}
