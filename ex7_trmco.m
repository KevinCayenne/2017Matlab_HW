
% 以下是資料的讀取與參數設定的部份，請大家務必自行設定
% 下面是先把資料讀取進來的語法，這裡是先預設大家把資料矩陣存成'final.mat'
load('final.mat');
% 這裡預設大家的final.mat第一個column是y % 
y=final(:,1);
% 這裡預設大家的final.mat第三個column是x1，然後中心化 % 
xs1 = final(:,3); 
mxs1 = log(final(:,3))-mean(log(final(:,3)));
mxs2 = final(:,4)-mean(final(:,4));
% 下面則是建構矩陣X %
xs=[mxs1 mxs2];
% 這裡設定起始值：r0(1)為beta0的起點，我們這裡預設是50，當然大家可以更改起始值試試看結果會不會不同！！ %
% 接著，r0(2)為beta1的起始值，r0(3)為beta2的起始值，r0(4)則為beta的變異數(sigma^2)的預測起始值 % 
r0=[50 0 0 15]';
% 下面的a和b，則是y的上下界，runmax是最多要疊代幾次，tau是step-size factor，tolv則是KKT條件能容忍的位數 %
a=1; b=7; runmax=1001; tau=10; tolv=10^-6;
% 下面的xmin和xmax，則是我們要定下的x0, x1, x2的最大最小值界限 %
xmin=[1 0 0]'; xmax=[1 100000 1]';
% 那麼下面的xmin 和 xmax反映的是什麼時候的x0, x1, x2最大最小值呢？大家可以自行發現 %
% xmin=[1 -3.82 0]'; xmax=[1 10.18 1]';

% 下面則是老師TRMCO主計算程式的部份
%******************************************************************
%** This is the core program of truncated regression model with  **
%** constrained optimization (TRMCO)                             **
%** Created by Min-Hua Huang, Dec 16, 2017                       **
%** y: dependent variable (n by 1)                               **
%** xs: covariate matrix without the constant (n by m)           **
%** r0: initial value of parameters (m+2 by 1)                   **
%** a: lower boundary limit (scalar)                             **
%** b: upper boundary limit (scalar)                             **
%** tau: step-size factor (>1, scalar)                           **
%** runmax: maximum number of iteration (scalar)                 **
%** tolv: tolerance value for KKT conditions                     **  
%** output: parameter estimates, negative inverse Hessian,       **
%** and related statistics                                       **  
%******************************************************************
format long; warning off all;
% basic setups
n=size(y,1); m=size(xs,2); x=[ones(n,1) xs]; record=ones(1,4*m+17); 
% initial values
lamda=ones(2*m+6,1)*10; kk=1; run=1; stop=1;
%****************** Record the initial results ***************************
sumD1=0;
  for i=1:1:n
        if abs(b-x(i,:)*r0(1:m+1))<5*r0(m+2) || abs(a-x(i,:)*r0(1:m+1))<5*r0(m+2)
           D=normcdf((b-x(i,:)*r0(1:m+1))/r0(m+2))-normcdf((a-x(i,:)*r0(1:m+1))/r0(m+2));           
        else
           D=1; 
        end   
  sumD1=sumD1-log(D*sqrt(2*pi)*r0(m+2));
  end   
record(1,1)=0;
record(1,2:m+3)=r0';
record(1,m+4:2*m+5)=zeros(1,m+2);
record(1,2*m+6:4*m+11)=lamda';
record(1,4*m+12)=sumD1-1/(2*r0(m+2)^2)*((y-x*r0(1:m+1))'*(y-x*r0(1:m+1)));
record(1,4*m+13:4*m+17)=ones(1,5)*-999;
% maximum and miminum vector of covariates 
% xmax=zeros(m+1,1); xmin=zeros(m+1,1);
%   for i=1:1:m+1
%       xmax(i)=max(x(:,i));
%       xmin(i)=min(x(:,i));
%   end    
xlimit=[xmax xmin];   
while stop>0 && run<runmax
% vectors of indicator variables
v=vsign(r0);
% constraint vector
ci=cpar(r0,a,b,x,xlimit,v);
% gradient vector
df=dfpar(r0,a,b,y,x);
% first derivative of the constraint vector
ai=aipar(r0,x,xlimit,v); 
%****************** Step 1:check KKT conditions ********************
% KKT
KKT1=df+ai'*lamda; check1=sum(abs(KKT1)>tolv);
KKT2=ci; check2=sum(KKT2>tolv);
KKT3=lamda; check3=sum(KKT3<-tolv);
KKT4=lamda'*ci; check4=sum(abs(KKT4)>tolv);
check=check1+check2+check3+check4;
   if check==0
      stop=0; 
   end
% If KKT conditions are not satisfied, then iteration continues.
%****************** Step 2:compute dk and lamdaQP ********************
% Hessian matrix
L=Lpar(r0,a,b,y,x,tau);   
% solving the QP problem  
d=quadprog(L,df,ai,-ci);  
lamdaQP=(ai')\(-df-L*d);
%****************** Step 3:set new x and lamda ***********************
r1=r0+d;
lamda=lamdaQP;
%****************** Record the updated results ***********************
sumD2=0;
  for i=1:1:n
        if abs(b-x(i,:)*r1(1:m+1))<5*r1(m+2) || abs(a-x(i,:)*r1(1:m+1))<5*r1(m+2)
           D=normcdf((b-x(i,:)*r1(1:m+1))/r1(m+2))-normcdf((a-x(i,:)*r1(1:m+1))/r1(m+2));           
        else
           D=1; 
        end   
  sumD2=sumD2-log(D*sqrt(2*pi)*r1(m+2));
  end  
record(kk+1,1)=kk;
record(kk+1,2:m+3)=r1';
record(kk+1,m+4:2*m+5)=d';
record(kk+1,2*m+6:4*m+11)=lamdaQP';
record(kk+1,4*m+12)=sumD2-1/(2*r1(m+2)^2)*((y-x*r1(1:m+1))'*(y-x*r1(1:m+1)));
record(kk+1,4*m+13)=record(kk+1,4*m+12)-record(kk,4*m+12);
record(kk+1,4*m+14:4*m+17)=[check1 check2 check3 check4];
%****************** Step 4 and 5: return to Step 1 *******************
r0=r1;
kk=kk+1;
run=run+1;
end
% check admissibility 
if run==runmax  
   maxLL=-10^20; maxt=-999;
   for t=1:1:runmax-1
        if record(t,4*m+12)>maxLL && record(t+1,4*m+15)==0
           maxLL=record(t,4*m+12);
           maxt=t;
        end
   end
   if maxt~=-999
      fs=record(maxt,2:m+3);
      note=record(maxt+1,4*m+15);
   end   
else
  if record(kk-1,m+3)<10^10
     fs=record(kk-1,2:m+3);   
     note=record(kk,4*m+15);
     maxt=kk-1;
  else
     maxLL=-10^20; maxt=-999; 
       for t=1:1:run-2
           if record(t,4*m+12)>maxLL && record(t+1,4*m+15)==0
              maxLL=record(t,4*m+12);
              maxt=t;
           end
       end
       if maxt~=-999 
          fs=record(maxt,2:m+3);
          note=record(maxt+1,4*m+15);
       end   
  end       
end
if maxt==-999
    fv=ones(m+6,m+2)*-999999999999999; 
else
    var=Lpar(fs',a,b,y,x,1)\eye(m+2);
    last=[record(2,4*m+15) record(1,4*m+12) zeros(1,m);note record(maxt,4*m+12) zeros(1,m); maxt run zeros(1,m)];
    save record;
	fv=record;
% 最後的結果，會存在record這個檔案裡面，最下面的一個row就是最後一個疊代的結果，如果收斂，則row數當然會小於1001 %
% 至於過程，請各位務必去看"record.mat"這個檔案！！ %
% 然後把結果print出來： %
fprintf('beta0:%6.3f, beta1:%6.3f, beta2:%6.3f\n',record(length(record),2),record(length(record),3),record(length(record),4));
end  

