# Algebraically equivalent cumulative-incidence gradients, accumulated over
# interval log-hazards before multiplying by design matrices.
predict_joint_fast<-function(ff,ref,ss,pp='2015-2018'){
 p<-length(ff[[1]]$b);n<-nrow(ref);z<-copy(ref);z[,sex:=factor(ss,levels=levels(ref$sex))];z[,period:=factor(pp,levels=levels(ref$period))]
 B<-do.call(cbind,lapply(ff,`[[`,'b'));Xlist<-hlist<-inc<-Slist<-alist<-aplist<-vector('list',5);S<-rep(1,n)
 for(j in 1:5){
  z[,interval:=factor(j,levels=1:5)];X<-model.matrix(form,z);h<-exp(X%*%B);if(pp!='2015-2018'||j==1)h[,5]<-0
  H<-rowSums(h);q<-exp(-H);a<--expm1(-H)/H;ap<-(q*H+expm1(-H))/H^2
  Xlist[[j]]<-X;hlist[[j]]<-h;Slist[[j]]<-S;alist[[j]]<-a;aplist[[j]]<-ap;inc[[j]]<-h*(S*a);S<-S*q
 }
 F<-Reduce('+',inc);future<-matrix(0,n,7);G<-matrix(0,7,7*p);w<-ref$weight
 for(j in 5:1){
  h<-hlist[[j]];X<-Xlist[[j]];Sp<-Slist[[j]];a<-alist[[j]];ap<-aplist[[j]]
  for(l in 1:7){
   deriv<-(h*(Sp*ap)-future)*h[,l];deriv[,l]<-deriv[,l]+Sp*a*h[,l]
   G[,((l-1)*p+1):(l*p)]<-G[,((l-1)*p+1):(l*p)]+t(crossprod(X,deriv*w))
  }
  future<-future+inc[[j]]
 }
 rownames(G)<-causes;ga<-colSums(G);est<-c(all=sum((1-S)*w),setNames(colSums(F*w),causes));G<-rbind(all=ga,G)
 est<-c(est,other_total=sum(est[causes[2:6]]));G<-rbind(G,other_total=colSums(G[causes[2:6],,drop=FALSE]));stopifnot(abs(est[1]-sum(est[causes]))<1e-10,all(is.finite(G)));list(est=est,G=G)
}
