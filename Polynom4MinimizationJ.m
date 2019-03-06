% function Polynom4MinimizationJ;
clc
clear all

% начальное и конечное положения системы, время 
x = [1,2,3,4]; % x = [x0, x'0, xt, x't]
y = [6,7,8,9]; % y = [y0, y'0, yt, y't] 

time = 3;
time_int = 0:0.001:time; 
maxVelocity = 50;
maxAcceleration = 50;
initApprox = [0,0,0];

% коэфф для многочленов Px и Py 4-ых степеней a2=a2(a4)...b3=b3(b4)
syms a4 b4
LeftSideX = [time^2, time^3; 2*time, 3*time^2];
RightSideX = [x(3)-x(1)-x(2)*time-a4*time^4; x(4)-x(2)-4*a4*time^3];
LeftSideY = [time^2, time^3; 2*time, 3*time^2];
RightSideY = [y(3)-y(1)-y(2)*time-b4*time^4; y(4)-y(2)-4*b4*time^3];

coeffsPolynomX = LeftSideX \ RightSideX;
coeffsPolynomY = LeftSideY \ RightSideY;

% многочлены Px и Py
syms t
polynomX =  x(1)+x(2)*t+t.^2*coeffsPolynomX(1)+t.^3*coeffsPolynomX(2)+t.^4*a4;
polynomY = y(1)+y(2)*t+t.^2*coeffsPolynomY(1)+t.^3*coeffsPolynomY(2)+t.^4*b4;

% нахождение Px', Px", Py', Py"
dpolynomX = diff(polynomX, t);
ddpolynomX = diff(dpolynomX, t);
dpolynomY = diff(polynomY, t);
ddpolynomY = diff(dpolynomY, t);

% ограничения на скорость и ускорение
velConstr = sqrt(dpolynomX^2+dpolynomY^2);
accConstr = sqrt(ddpolynomX^2+ddpolynomY^2);

% функционал потерь
intExpression = polynomX^2+polynomY^2+dpolynomX^2+dpolynomY^2+ddpolynomX^2+ddpolynomY^2;
J = 0.5*int(intExpression,0,time)-t*1e-9;
Jhandle = matlabFunction(J);


global c1handle;
global c2handle;
c1handle=matlabFunction(velConstr-maxVelocity);
c2handle=matlabFunction(accConstr-maxAcceleration);


nonlcon=@constraints;
 A=[0,0,1;0,0,-1];
 b=[time;0];

minCoeffs = fmincon(Jhandle,initApprox,A,b,[],[],[],[],nonlcon);

%построение многочлена, который является траекторий системы из (y0 y'0) в (0 0)
dotsPolynomX  = polynomX(time_int);
dotsPolynomY  = polynomY(time_int);

figure
hold on
plot(dotsPolynomX,dotsPolynomY,'Red', 'LineWidth', 1.5);
xlabel('x') 
ylabel('y') 
hold off

% end



function [c,ceq] = constraints(x,y,t)
global c1handle;
global c2handle;
c =[c1handle(x,y,t);c2handle(x,y,t)];
ceq = [];
end


 