function fv=LL(y,x,r0,a,b)
%********************************************************************
%** This is a supporting function for model2.m                     **
%** Created by Min-Hua Huang, Dec 19, 2011                         **
%** y: dependent variable                                          **
%** x: covariate matrix with the constant                          **
%** r0: parameter estimates                                        **
%** a: lower boundary limit                                        **
%** b: upper boundary limit                                        **
%** output: Loglikelihood value                                    **
%********************************************************************
m=size(x,2)-1;
n=size(y,1); sum=0;
  for i=1:1:n
        if abs(b-x(i,:)*r0(1:m+1))<5*r0(m+2) || abs(a-x(i,:)*r0(1:m+1))<5*r0(m+2)
           D=normcdf((b-x(i,:)*r0(1:m+1))/r0(m+2))-normcdf((a-x(i,:)*r0(1:m+1))/r0(m+2));           
        else
           D=1; 
        end   
  sum=sum-log(D*sqrt(2*pi)*r0(m+2));
  end   
fv=sum-1/(2*r0(m+2)^2)*((y-x*r0(1:m+1))'*(y-x*r0(1:m+1)));


