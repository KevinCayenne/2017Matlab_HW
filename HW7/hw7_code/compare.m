function fv=compare(r0)
%********************************************************************
%** This is a supporting function for aipar.m                      **
%** Created by Min-Hua Huang, Dec 19, 2011                         **
%** r0: parameter estimates                                        **
%** output: maximum and minimum of dummy variables                 **
%********************************************************************
k=size(r0,1); maxb=-10^10; minb=10^10;
for i=1:1:k
    if r0(i)>maxb
       maxb=r0(i); 
    end
    if r0(i)<minb
       minb=r0(i);
    end   
end
note=zeros(k,2);
for i=1:1:k
    if abs(r0(i)-maxb)<10^-10
       note(i,1)=1;
    end
    if abs(r0(i)-minb)<10^-10
       note(i,2)=1;
    end   
end
fv=note;