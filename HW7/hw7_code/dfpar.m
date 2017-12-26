function fv=dfpar(r0,a,b,y,x)   
%********************************************************************
%** This is a supporting function for trmcod.m                     **
%** Created by Min-Hua Huang, Dec 19, 2011                         **
%** r0: parameter estimates                                        **
%** a: lower boundary limit                                        **
%** b: upper boundary limit                                        **
%** y: dependent variable                                          **
%** x: covariate matrix with the constant                          **
%** output: gradient vector                                        **
%********************************************************************
n=size(y,1); m=size(x,2)-1; df=ones(m+2,1)*-999; 
% d(logL)/dbj
for j=0:1:m
    sum=0; 
    for i=1:1:n
         if abs(b-x(i,:)*r0(1:m+1))<5*r0(m+2) || abs(a-x(i,:)*r0(1:m+1))<5*r0(m+2)
           D=normcdf((b-x(i,:)*r0(1:m+1))/r0(m+2))-normcdf((a-x(i,:)*r0(1:m+1))/r0(m+2));           
        else
           D=1; 
         end   
        Di=D*sqrt(2*pi)*r0(m+2);  
        % d(Di)/dbj
        dDdbj=x(i,j+1)*(-exp(-(b-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2))+exp(-(a-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2)));
        sum=sum+1/Di*dDdbj;
    end    
    df(j+1)=sum-1/(r0(m+2)^2)*((y-x*r0(1:m+1))'*x(:,j+1));       
end
% d(logL)/dsig
sum=0; 
for i=1:1:n
         if abs(b-x(i,:)*r0(1:m+1))<5*r0(m+2) || abs(a-x(i,:)*r0(1:m+1))<5*r0(m+2)
           D=normcdf((b-x(i,:)*r0(1:m+1))/r0(m+2))-normcdf((a-x(i,:)*r0(1:m+1))/r0(m+2));           
        else
           D=1; 
         end   
        Di=D*sqrt(2*pi)*r0(m+2);       
     % d(Di)/dsig
     part1=-(b-x(i,:)*r0(1:m+1))/r0(m+2)*exp(-(b-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2));
     part2=(a-x(i,:)*r0(1:m+1))/r0(m+2)*exp(-(a-x(i,:)*r0(1:m+1))^2/(2*r0(m+2)^2));
     dDdsig=part1+part2+Di/r0(m+2);
     sum=sum+1/Di*dDdsig;
 end    
df(m+2)=sum-1/(r0(m+2)^3)*((y-x*r0(1:m+1))'*(y-x*r0(1:m+1)));   
fv=df;