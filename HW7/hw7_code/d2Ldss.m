function fv=d2Ldss(r0,a,b,y,x)
%********************************************************************
%** This is a supporting function for Lpar.m                       **
%** Created by Min-Hua Huang, Dec 19, 2011                         **
%** r0: parameter estimates                                        **
%** a: lower boundary limit                                        **
%** b: upper boundary limit                                        **
%** y: dependent variable                                          **
%** x: covariate matrix with the constant                          **
%** output: d^2(logL)/d(sig)d(sig) in the Hessian matrix           **
%********************************************************************
n=size(y,1); m=size(x,2)-1; sum=0;
   for i=1:1:n
         if abs(b-x(i,:)*r0(1:m+1))<5*r0(m+2) || abs(a-x(i,:)*r0(1:m+1))<5*r0(m+2)
           D=normcdf((b-x(i,:)*r0(1:m+1))/r0(m+2))-normcdf((a-x(i,:)*r0(1:m+1))/r0(m+2));           
        else
           D=1; 
         end   
      Di=D*sqrt(2*pi)*r0(m+2);          
      part1=-(b-x(i,:)*r0(1:m+1))/r0(m+2)*exp(-(b-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2));
      part2=(a-x(i,:)*r0(1:m+1))/r0(m+2)*exp(-(a-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2));
      dDds=(part1+part2+Di/r0(m+2));
      d2Ddss_1=-((b-x(i,:)*r0(1:m+1))^3/r0(m+2)^4)*exp(-(b-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2));
      d2Ddss_2=((a-x(i,:)*r0(1:m+1))^3/r0(m+2)^4)*exp(-(a-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2));
      sum=sum-1/Di^2*dDds^2+1/Di*(d2Ddss_1+d2Ddss_2);    
   end    
sum=sum+3/r0(m+2)^4*((y-x*r0(1:m+1))'*(y-x*r0(1:m+1)));
fv=sum;