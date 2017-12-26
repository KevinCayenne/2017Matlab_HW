function fv=cpar(r0,a,b,x,xlimit,v)
%********************************************************************
%** This is a supporting function for trmcod.m                     **
%** Created by Min-Hua Huang, Dec 19, 2011                         **
%** r0: parameter estimates                                        **
%** a: lower boundary limit                                        **
%** b: upper boundary limit                                        **
%** x: covariate matrix with the constant                          **
%** xlimit: maximum and miminum vector of covariates               **
%** v: vectors of indicator variables                              **
%** output: constraints vector                                     **
%********************************************************************
m=size(x,2)-1; g=ones(2*m+6,1);
%maximum and minimum vectors of estimated y when xj is held at its min.
ynonj=ones(m+1,2)*-999;
ynonj(1,1)=(v(2:m+1,1).*r0(2:m+1))'*xlimit(2:m+1,1);
ynonj(1,2)=(v(2:m+1,2).*r0(2:m+1))'*xlimit(2:m+1,1);
for i=2:1:m+1
    ynonj(i,1)=r0(1)+(v(2:m+1,1).*r0(2:m+1))'*xlimit(2:m+1,1)-v(i,1)*r0(i)*xlimit(i,1);
    ynonj(i,2)=r0(1)+(v(2:m+1,2).*r0(2:m+1))'*xlimit(2:m+1,1)-v(i,2)*r0(i)*xlimit(i,1);
end
%constraints vector
g(1)=r0(1)+ynonj(1,1)-b;
g(2)=a-r0(1)-ynonj(1,2);
g(3)=r0(1)-b;
g(4)=-r0(1)+a;
for i=1:1:m
    g(2*i+3)=r0(i+1)-(b-ynonj(i+1,1))/xlimit(i+1,1);
    g(2*i+4)=-r0(i+1)+(a-ynonj(i+1,2))/xlimit(i+1,1);
end
g(2*m+5)=r0(m+2)-b+a;
g(2*m+6)=-r0(m+2)+0.001;
fv=g;   