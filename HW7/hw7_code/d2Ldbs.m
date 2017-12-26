function fv=d2Ldbs(r0,a,b,y,x,j)
%********************************************************************
%** This is a supporting function for Lpar.m                       **
%** Created by Min-Hua Huang, Dec 19, 2011                         **
%** r0: parameter estimates                                        **
%** a: lower boundary limit                                        **
%** b: upper boundary limit                                        **
%** y: dependent variable                                          **
%** x: covariate matrix with the constant                          **
%** j: row or column indicator for  d^2(logL)/d(bj)d(sig)          **
%** output: d^2(logL)/d(bj)d(sig) in the Hessian matrix            **
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
      part1=-(b-x(i,:)*r0(1:m+1))/r0(m+2)*exp(-(b-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2));
      part2=(a-x(i,:)*r0(1:m+1))/r0(m+2)*exp(-(a-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2));
      dDds=(part1+part2+Di/r0(m+2));
      d2Ddbjs_1=-x(i,j)*(b-x(i,:)*r0(1:m+1))^2/r0(m+2)^3*exp(-(b-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2));
      d2Ddbjs_2=x(i,j)*(a-x(i,:)*r0(1:m+1))^2/r0(m+2)^3*exp(-(a-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2));
      sum=sum-1/Di^2*dDdbj*dDds+1/Di*(d2Ddbjs_1+d2Ddbjs_2);    
   end    
sum=sum+2/r0(m+2)^3*((y-x*r0(1:m+1))'*x(:,j));
fv=sum;


