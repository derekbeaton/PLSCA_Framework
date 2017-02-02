rplsca <- function(X,Y,lambda.x=NULL,lambda.y=NULL,projmat.X=NULL,projmat.Y=NULL,k=5){
#rplsca <- function(X,Y,lambda.x=NULL,lambda.y=NULL,as.is=F,projmat.X=NULL,projmat.Y=NULL,k=5){


	if( is.null(lambda.x) & is.null(lambda.y) ){
		lambda.x <- lambda.y <- 0
	}else if( !is.null(lambda.x) & is.null(lambda.y) ){
		lambda.y <- lambda.x
	}else if( is.null(lambda.x) & !is.null(lambda.y) ){
		lambda.x <- lambda.y
	}
	
	X.preproc <- takanePreProc(X,lambda.x,projmat.X)
	Y.preproc <- takanePreProc(Y,lambda.y,projmat.Y)
	X.dat <- sqrt.mat(X.preproc$LW) %*% X.preproc$dat
	Y.dat <- sqrt.mat(Y.preproc$LW) %*% Y.preproc$dat
	
	##test classes of matrices here...
	
	R.dat <- t( X.dat %*% svd.pinv_sq(X.preproc$RW) ) %*% ( Y.dat %*% svd.pinv_sq(Y.preproc$RW) )
	rownames(R.dat) <- colnames(X.dat)
	colnames(R.dat) <- colnames(Y.dat)	
		
	gres.out <- gsvd(R.dat,X.preproc$RW,Y.preproc$RW,k=k)
	gres.out$LX <- X.dat %*% gres.out$p
	gres.out$LY <- Y.dat %*% gres.out$q
	gres.out$fi <- gres.out$p * matrix(gres.out$Dv,nrow(gres.out$p),ncol(gres.out$p),byrow=T)
	gres.out$fj <- gres.out$q * matrix(gres.out$Dv,nrow(gres.out$q),ncol(gres.out$q),byrow=T)
	
	gres.out$X.preproc <- X.preproc
	gres.out$Y.preproc <- Y.preproc	
	
	##Takes up precious time.
	# gres.out$fi <- gres.out$p %*% diag(gres.out$Dv)
	# gres.out$fj <- gres.out$q %*% diag(gres.out$Dv)

	##Takes up precious space.	
	# opt.arr <- array(NA,dim=c(length(gres.out$Dv),length(gres.out$Dv),4))
	# opt.arr[,,1] <- t(gres.out$fi) %*% X.preproc$RW %*% gres.out$fi
	# opt.arr[,,2] <- t(gres.out$fj) %*% Y.preproc$RW %*% gres.out$fj
	# opt.arr[,,3] <- t(gres.out$LX) %*% gres.out$LY
	# opt.arr[,,4] <- cor(gres.out$LX,gres.out$LY)
	# gres.out$optimizations <- opt.arr

	##Takes up precious space.
	#gres.out$preprocs <- list(X=X.preproc, Y=Y.preproc)

	return(gres.out)
	
}
