function fv=robuster(y,x,r0,code)   
%********************************************************************
%** This is a supporting function for model2.m                     **
%** Created by Min-Hua Huang, Dec 19, 2011                         **
%** y: dependent variable                                          **
%** x: covariate matrix with the constant                          **
%** r0: parameter estimates                                        **
%** code: group id                                                 **
%** output: robust standard error                                  **
%********************************************************************
n=size(y,1); m=size(r0,1)-1; er=y-x*r0(1:m); B=0;
for index=min(code):1:max(code) 
    count=0; dfi=0;
    for i=1:1:n 
        if code(i)==index;
           count=count+1;
           temp=x(i,:)*er(i)/r0(m+1)^2;
           dfi=dfi+temp;
        end   
    end
     B=B+dfi'*dfi;
end
sumL=0;
for i=1:1:n
    sumL=sumL-x(i,:)'*x(i,:)/r0(m+1)^2;
end
invnA=sumL\eye(m);
V=invnA*B*invnA;
fv=sqrt(diag(V));
