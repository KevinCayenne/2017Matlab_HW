function fv=d2Ldbb(r0,a,b,y,x,j,k,tau)
%********************************************************************
%** This is a supporting function for Lpar.m                       **
%** Created by Min-Hua Huang, Dec 19, 2011                         **
%** r0: parameter estimates                                        **
%** a: lower boundary limit                                        **
%** b: upper boundary limit                                        **
%** y: dependent variable                                          **
%** x: covariate matrix with the constant                          **
%** jk: row and column indicators for  d^2(logL)/d(bj)d(bk)        **
%** tau: step-size factor                                          **
%** output: d^2(logL)/d(bj)d(bk) in the Hessian matrix             **
%********************************************************************
n=size(y,1); m=size(x,2)-1; sum=0;
   for i=1:1:n
         if abs(b-x(i,:)*r0(1:m+1))<5*r0(m+2) || abs(a-x(i,:)*r0(1:m+1))<5*r0(m+2)
           D=normcdf((b-x(i,:)*r0(1:m+1))/r0(m+2))-normcdf((a-x(i,:)*r0(1:m+1))/r0(m+2));           
        else
           D=1; 
         end        
      Di=D*sqrt(2*pi)*r0(m+2);  
      dDdbj=x(i,j)*(-exp(-(b-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2))+exp(-(a-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2)));
      dDdbk=x(i,k)*(-exp(-(b-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2))+exp(-(a-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2)));
      d2Ddbjbk_1=-x(i,j)*x(i,k)*(b-x(i,:)*r0(1:m+1))/r0(m+2)^2*exp(-(b-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2));
      d2Ddbjbk_2=x(i,j)*x(i,k)*(a-x(i,:)*r0(1:m+1))/r0(m+2)^2*exp(-(a-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2));
      sum=sum-1/Di^2*dDdbj*dDdbk+1/Di*(d2Ddbjbk_1+d2Ddbjbk_2);    
   end    
sum=sum+tau/r0(m+2)^2*(x(:,j)'*x(:,k));
fv=sum;