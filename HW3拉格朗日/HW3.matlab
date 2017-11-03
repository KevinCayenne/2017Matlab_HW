
%% 1.

function F=hw1fun(x)  
    F = [2*(x(1)-3)+2*x(1)*x(3);
         2*(x(2)-1)+2*x(2)*x(3);
         x(1)^2+x(2)^2-1];
end

x0=[1 0 0]';
options = optimset('Algorithm','trust-region-dogleg','Display','iter');

[x,fval] = fsolve(@hw1fun,x0,options)

%% 2.

function F=hw2fun(x)
	F =[(2/3)*200*x(1)^(-1/3)*x(2)(1/3) + 20*x(3);
		(1/3)*200*x(1)^(2/3)*x(2)^(-2/3) + 170*x(3);
		20*x(1) + 170*x(2) - 20000];
end

x0=[1 1 0]';
options = optimset('Algorithm','trust-region-dogleg','Display','iter');

[x,fval] = fsolve(@hw2fun,x0,options)