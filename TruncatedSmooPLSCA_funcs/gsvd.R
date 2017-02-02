gsvd <- function(G.DAT,LW,RW,k=5){
	
	
	#dat.for.svd <- sqrt.mat(LW) %*% G.DAT %*% sqrt.mat(RW)
	dat.for.svd <- matrix(sqrt(LW),length(LW), ncol(G.DAT),byrow=F) * G.DAT * matrix(sqrt(RW), nrow(G.DAT), length(RW),byrow=T)
	
	if(k < 1){
		k <- min( nrow(dat.for.svd)-1, ncol(dat.for.svd)-1 )
	}
	res <- svd(dat.for.svd,nu=k,nv=k)
	d <- res$d
	precisionLimit <- 2 * .Machine$double.eps
	indToKeep <- which(d^2 > precisionLimit)
	d <- d[indToKeep]
	tau <- d^2/sum(d^2)

	comp.ret <- min(length(indToKeep),k)
	res$u <- as.matrix(res$u[,1:comp.ret])
	res$v <- as.matrix(res$v[,1:comp.ret])
#	P <- sqrt.mat(svd.pinv_sq(LW)) %*% res$u
	P <- matrix(1/sqrt(LW),nrow(res$u),ncol(res$u)) * res$u
	rownames(P) <- rownames(G.DAT)
#	Q <- sqrt.mat(svd.pinv_sq(RW)) %*% res$v	
	Q <- matrix(1/sqrt(RW),nrow(res$v),ncol(res$v)) * res$v
	rownames(Q) <- colnames(G.DAT)	
	
    return(list(p = P, q = Q, u=res$u, v=res$v, Dv = d[1:comp.ret], Dv.orig=d, tau = tau))
}