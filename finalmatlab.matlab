record=ones(1,32)*-999; 

LT=[5 3 5 3 5]';

cont_a=-0.5; cont_b=-1; cont_v=-4;
a0=0; 
b0=0; 
a=11; 
b=-3;
r0=[1 0 1 0 1 -1 1 -2]'; % [x1 y1 x2 y2 x3 y3 x4 y4]
lamda=ones(9,1); 
kk=1; run=1; stop=1; 
runmax=1000; tolv=10^-4;

record(1,1)=0;
record(1,2:9)=r0';	
record(1,10:17)=zeros(1,8);
record(1,18:26)=lamda';
record(1,27)=0.5*(LT(1)*(b0+r0(2))+LT(2)*(r0(2)+r0(4))+LT(3)*(r0(4)+r0(6))+LT(4)*(r0(6)+r0(8))+LT(5)*(r0(8)+b));
record(1,28:32)=ones(1,5)*-999;


while stop>0 && run<runmax

	g1=(a0-r0(1))^2+(b0-r0(2))^2-LT(1)^2;
	g2=(r0(3)-r0(1))^2+(r0(4)-r0(2))^2-LT(2)^2;
	g3=(r0(5)-r0(3))^2+(r0(6)-r0(4))^2-LT(3)^2;
	g4=(r0(7)-r0(5))^2+(r0(8)-r0(6))^2-LT(4)^2;
	g5=(a-r0(7))^2+(b-r0(8))^2-LT(5)^2;
	g6=cont_a*r0(1)+cont_b*r0(2)+cont_v;
	g7=cont_a*r0(3)+cont_b*r0(4)+cont_v;
	g8=cont_a*r0(5)+cont_b*r0(6)+cont_v;
	g9=cont_a*r0(7)+cont_b*r0(8)+cont_v;
	ci=[g1 g2 g3 g4 g5 g6 g7 g8 g9]';

	df=[0 (LT(1)+LT(2))/2 0 (LT(2)+LT(3))/2 0 (LT(3)+LT(4))/2 0 (LT(4)+LT(5))/2]';

	ai1=[2*(r0(1)-a0) 2*(r0(2)-b0) 0 0 0 0 0 0]; 
	ai2=[-2*(r0(3)-r0(1)) -2*(r0(4)-r0(2))  2*(r0(3)-r0(1)) 2*(r0(4)-r0(2)) 0 0 0 0]; 
	ai3=[0 0 -2*(r0(5)-r0(3)) -2*(r0(6)-r0(4))  2*(r0(5)-r0(3)) 2*(r0(6)-r0(4)) 0 0]; 
	ai4=[0 0 0 0 -2*(r0(7)-r0(5)) -2*(r0(8)-r0(6))  2*(r0(7)-r0(5)) 2*(r0(8)-r0(6))]; 
	ai5=[0 0 0 0 0 0 -2*(a-r0(7)) 2*(b+r0(8))]; 
	ai6=[cont_a cont_b 0 0 0 0 0 0]; 
	ai7=[0 0 cont_a cont_b 0 0 0 0];
	ai8=[0 0 0 0 cont_a cont_b 0 0];
	ai9=[0 0 0 0 0 0 cont_a cont_b];
	ai=[ai1;ai2;ai3;ai4;ai5;ai6;ai7;ai8;ai9];

	L1=[2*lamda(1)+2*lamda(2) 0 -2*lamda(2)  0 0 0 0 0];
	L2=[0 2*lamda(1)+2*lamda(2) 0 -2*lamda(2) 0 0 0 0];
	L3=[-2*lamda(2) 0 2*lamda(2)+2*lamda(3) 0 -2*lamda(3) 0 0 0];
	L4=[0 -2*lamda(2) 0 2*lamda(2)+2*lamda(3) 0 -2*lamda(3) 0 0];
	L5=[0 0 -2*lamda(3) 0 2*lamda(3)+2*lamda(4) 0 -2*lamda(4) 0];
	L6=[0 0 0 -2*lamda(3) 0 2*lamda(3)+2*lamda(4) 0 -2*lamda(4)];
	L7=[0 0 0 0 -2*lamda(4) 0 2*lamda(4)+2*lamda(5) 0];
	L8=[0 0 0 0 0 -2*lamda(4) 0 2*lamda(4)+2*lamda(5)];
	L=[L1;L2;L3;L4;L5;L6;L7;L8];

	% KKT
	KKT1=df+ai'*lamda; check1=sum(abs(KKT1)>tolv);
	KKT2=ci; check2=sum(KKT2>tolv);
	KKT3=lamda; check3=sum(KKT3<-tolv);
	KKT4=lamda'*ci; check4=sum(abs(KKT4)>tolv);
	check=check1+check2+check3+check4;
	   if check==0
	      stop=0; 
	   end

	d=quadprog(L,df,ai,-ci);  
	lamdaQP=(ai')\(-df-L*d);

	r1=r0+d;
	lamda=lamdaQP;

	record(kk+1,1)=kk;
	record(kk+1,2:9)=r1';
	record(kk+1,10:17)=d';
	record(kk+1,18:26)=lamdaQP';
	record(kk+1,27)=0.5*(LT(1)*(b0+r1(2))+LT(2)*(r1(2)+r1(4))+LT(3)*(r1(4)+r1(6))+LT(4)*(r1(6)+r1(8))+LT(5)*(r1(8)+b));
	record(kk+1,28)=record(kk+1,21)-record(kk,21);
	record(kk+1,29:32)=[check1 check2 check3 check4];
	%****************** Step 4 and 5: return to Step 1 *******************
	r0=r1;
	kk=kk+1;
	run=run+1;
end