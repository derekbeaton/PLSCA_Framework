rplsca <- function(X,Y,lambda.x=NULL,lambda.y=NULL,projmat.X=NULL,projmat.Y=NULL,k=5){

	if( is.null(lambda.x) & is.null(lambda.y) ){
		lambda.x <- lambda.y <- 0
	}else if( !is.null(lambda.x) & is.null(lambda.y) ){
		lambda.y <- lambda.x
	}else if( is.null(lambda.x) & !is.null(lambda.y) ){
		lambda.x <- lambda.y
	}
	
	X.preproc <- takanePreProc(X,lambda.x,projmat.X)
	Y.preproc <- takanePreProc(Y,lambda.y,projmat.Y)

	R.dat <- t( X.preproc$dat * matrix(1/X.preproc$RW,nrow(X.preproc$dat),ncol(X.preproc$dat),byrow=T) ) %*% ( Y.preproc$dat * matrix(1/Y.preproc$RW,nrow(Y.preproc$dat),ncol(Y.preproc$dat),byrow=T) )
	
	rownames(R.dat) <- colnames(X.preproc$dat)
	colnames(R.dat) <- colnames(Y.preproc$dat)
				
	gres.out <- gsvd(R.dat,X.preproc$RW,Y.preproc$RW,k=k)
	
	gres.out$X.preproc <- X.preproc
	gres.out$Y.preproc <- Y.preproc	

	
	gres.out$LX <- as.matrix(X.preproc$dat %*% gres.out$p)
	gres.out$LY <- as.matrix(Y.preproc$dat %*% gres.out$q)
		rownames(gres.out$LX) <- rownames(gres.out$LY) <- rownames(X)
	gres.out$fi <- as.matrix(gres.out$p * matrix(gres.out$Dv,nrow(gres.out$p),ncol(gres.out$p),byrow=T))
		rownames(gres.out$fi) <- colnames(X)
	gres.out$fj <- as.matrix(gres.out$q * matrix(gres.out$Dv,nrow(gres.out$q),ncol(gres.out$q),byrow=T))
		rownames(gres.out$fj) <- colnames(Y)
	
	return(gres.out)
	
}
