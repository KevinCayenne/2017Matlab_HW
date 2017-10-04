function fv=ex1(X)

n=size(X,1);

Xs=ones(n,1)*-999;

idx=0;

 for i=1:1:n

     k=X(i);

     for j=i+1:1:n
         if k<X(j)
            k=X(j);
            idx=j;
         end
     end

     if k~=X(i)
        Xs(i)=k;
        X(idx)=[];
        X1=[k;X];
        X=X1;

     else
        Xs(i)=X(i); 
     end     
 end 

fv=Xs;
end