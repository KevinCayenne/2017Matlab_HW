%% fig 1

x = linspace(-3, 3);
y = 5*(x.^5)-3*x+6;
plot(x, y);

%% fig 2

bm=[-0.049 0.061 0.095 0.009 -0.011 0.057 0.083]';
xm=[1 0 1 0 0 0 1]; 
C=xm*bm;

x=0:0.5:6; 
y=0:5:40;
[X,Y] = meshgrid(x,y);

Z=X.*0.032+Y.*0.003-X.*(Y.*0.002)+C;

plot3(X,Y,Z);
xlabel('Regime Type');
ylabel('Limits on Content');
zlabel('Correlation of ICTs Use and Political Activism');

grid on
