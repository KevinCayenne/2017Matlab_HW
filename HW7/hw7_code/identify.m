function fv=identify(r0)
%********************************************************************
%** This is a function to identify boundary violations             **
%** Created by Min-Hua Huang, Dec 19, 2011                         **
%** r0: parameter estimates                                        **
%** output: boundary violations                                    **
%********************************************************************
format long; load x.mat
n=size(x,1); xs=[ones(n,1) x]; dmy=4; m=size(xs,2)-1; 
xmax=zeros(m+1,1); xmin=zeros(m+1,1);
   for i=1:1:m+1
       xmax(i)=max(xs(:,i));
       xmin(i)=min(xs(:,i));
   end    
xlimit=[xmax xmin];  
v=vsign(r0);
m1=m+1-dmy; predict=ones(m+3,2);
ynonj=ones(m+1,2)*-999;
dmymax=max(r0(m1+1:m+1));
dmymin=min(r0(m1+1:m+1));
ynonj(1,1)=(v(2:m1,1).*r0(2:m1))'*xlimit(2:m1,1)+dmymax;
ynonj(1,2)=(v(2:m1,2).*r0(2:m1))'*xlimit(2:m1,1)+dmymin;
for i=2:1:m1
    ynonj(i,1)=r0(1)+(v(2:m1,1).*r0(2:m1))'*xlimit(2:m1,1)-v(i,1)*r0(i)*xlimit(i,1)+dmymax;
    ynonj(i,2)=r0(1)+(v(2:m1,2).*r0(2:m1))'*xlimit(2:m1,1)-v(i,2)*r0(i)*xlimit(i,1)+dmymin;
end
ynonj(m1+1:m+1,1)=r0(1)+(v(2:m1,1).*r0(2:m1))'*xlimit(2:m1,1);
ynonj(m1+1:m+1,2)=r0(1)+(v(2:m1,2).*r0(2:m1))'*xlimit(2:m1,1);
predict(1,1)=r0(1)+ynonj(1,1);
predict(1,2)=r0(1)+ynonj(1,2);
predict(2,1)=r0(1);
predict(2,2)=r0(1);
for i=1:1:m
    predict(i+2,1)=xlimit(i+1,1)*r0(i+1)+ynonj(i+1,1);
    predict(i+2,2)=xlimit(i+1,1)*r0(i+1)+ynonj(i+1,2);
end
predict(m+3,1)=r0(m+2);
predict(m+3,2)=r0(m+2);
fv=predict;
