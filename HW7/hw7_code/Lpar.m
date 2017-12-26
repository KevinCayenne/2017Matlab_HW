function fv=Lpar(r0,a,b,y,x,tau)
%********************************************************************
%** This is a supporting function for trmcod.m                     **
%** Created by Min-Hua Huang, Dec 19, 2011                         **
%** r0: parameter estimates                                        **
%** a: lower boundary limit                                        **
%** b: upper boundary limit                                        **
%** y: dependent variable                                          **
%** x: covariate matrix with the constant                          **
%** tau: step-size factor                                          **
%** output: Hessian matrix                                         **
%********************************************************************
m=size(x,2)-1; L=ones(m+2,m+2)*-999;
% d^2(logL)/d(bj)d(bk)
for i=1:1:m+1
    for j=1:1:i
        L(i,j)=d2Ldbb(r0,a,b,y,x,i,j,tau);
    end
end
for i=1:1:m+1
    for j=m+1:-1:i+1
        L(i,j)=L(j,i);
    end
end
% d^2(logL)/d(bj)d(sig)
for i=1:1:m+1
    L(i,m+2)=d2Ldbs(r0,a,b,y,x,i);
    L(m+2,i)=L(i,m+2);
end
% d^2(logL)/d(sig)d(sig)
L(m+2,m+2)=d2Ldss(r0,a,b,y,x);
fv=L; 