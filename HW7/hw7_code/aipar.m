function fv=aipar(r0,x,xlimit,v)  
%********************************************************************
%** This is a supporting function for trmcod.m                     **
%** Created by Min-Hua Huang, Dec 19, 2011                         **
%** r0: parameter estimates                                        **
%** x: covariate matrix with the constant                          **
%** xlimit: maximum and miminum vector of covariates               **
%** v: vectors of indicator variables                              **
%** dmy: number of dummy variable                                  **
%** output: first derivative of the constraint vector              **
%********************************************************************
m=size(x,2)-1; ai=zeros(2*m+6,m+2); m1=m+1; v(1,:)=[1 1];
% derivative of ynonj
dynonj=ones(m1,2)*-999;
for i=1:1:m1
    dynonj(i,1)=v(i,1)*xlimit(i,1);
    dynonj(i,2)=v(i,2)*xlimit(i,1);
end
% *************************************************************
% yhat
ai(1,1)=1; ai(2,1)=-1; 
for i=1:1:m1
    ai(1,i)=dynonj(i,1); 
    ai(2,i)=-dynonj(i,2); 
end    
ridx=compare(r0(m1+1:m+1));
ai(1,m1+1:m+1)=ridx(:,1)';
ai(2,m1+1:m+1)=-ridx(:,2)';
% b0
ai(3,1)=1; ai(4,1)=-1; 
% b1 to bm1
for i=2:1:m+1
    for j=1:1:m1
        if i~=j
               ai(2*i+1,j)=dynonj(j,1)/xlimit(i,1);
               ai(2*i+2,j)=-dynonj(j,2)/xlimit(i,1);
        else
           ai(2*i+1,j)=1;
           ai(2*i+2,j)=-1;           
        end
    end
end
for i=2:1:m1
          ai(2*i+1,m1+1:m+1)=ridx(:,1)'/xlimit(i,1); 
          ai(2*i+2,m1+1:m+1)=-ridx(:,2)'/xlimit(i,1); 
end    
for i=m1+1:1:m+1
      for j=m1+1:1:m+1
          if i==j
              ai(2*i+1,j)=1;
              ai(2*i+2,j)=-1;  
          end
      end
end      
% sigma
ai(2*m+5,m+2)=1;
ai(2*m+6,m+2)=-1;
fv=ai;