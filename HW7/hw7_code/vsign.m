function fv=vsign(r0)
%********************************************************************
%** This is a supporting function for trmcod.m                     **       
%** Created by Min-Hua Huang, Dec 19, 2011                         **
%** r0: parameter estimates                                        **
%** output: indicator variables (v^{+},v^{-}) for r0               **
%********************************************************************
m=size(r0,1); v=zeros(m,2); 
for i=1:1:m
    if r0(i)>0
          v(i,1)=1; 
       else
          v(i,2)=1; 
    end
    if abs(r0(i))<10^-6
          v(i,1)=0;
          v(i,2)=0;
    end
    if i==1
       v(1,:)=[1 1];
    end 
end    
fv=v;


