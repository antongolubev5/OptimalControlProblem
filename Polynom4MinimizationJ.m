% function Polynom4MinimizationJ;
clc
clear all

% начальное и конечное положения системы, время 
x = [1,1,3,0]; % x = [x0, x'0, xt, x't]
y = [6,1,8,0]; % y = [y0, y'0, yt, y't] 

time = 10;
time_int = 0:0.1:time; 
maxVelocity = 1.6-0.05;
maxAcceleration = 2;
initApprox = [0,0];
optCoeffs = [-1,-1];

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
J = 0.5*int(intExpression,0,time);

% минимизация функционала
global Jhandle;
Jhandle = matlabFunction(J);

global c1handle;
global c2handle;
c1handle=matlabFunction(velConstr-maxVelocity);
c2handle=matlabFunction(accConstr-maxAcceleration);
global timemoment;

for tparam = 0:0.1:time
timemoment=tparam;
nonlcon=@constraints;
minCoeffs = fmincon(@fun,initApprox,[],[],[],[],[],[],nonlcon);
disturb1=c1handle(minCoeffs(1),minCoeffs(2),time_int) > 0.01;
disturb2=c2handle(minCoeffs(1),minCoeffs(2),time_int) > 0.01;

if not(any(disturb1) || any(disturb2))
   optCoeffs =vertcat(optCoeffs,minCoeffs);
end
end
[m,n]=size(optCoeffs);
currval=Jhandle(optCoeffs(2,1),optCoeffs(2,2));
imin=2;
if m>1
   for i = 2:1:m 
       if Jhandle(optCoeffs(i,1),optCoeffs(i,2))<currval
          currval= Jhandle(optCoeffs(i,1),optCoeffs(i,2));
          imin=i;
       end
   end

% траектории, скорость, ускорение через полученный коэфф
polynomXopt=matlabFunction(polynomX);
polynomYopt=matlabFunction(polynomY);

velocity=matlabFunction(velConstr);
acceleration=matlabFunction(accConstr);

dotsPolynomX=polynomXopt(optCoeffs(imin,1),time_int);
dotsPolynomY=polynomYopt(optCoeffs(imin,2),time_int);

dotsVel=velocity(optCoeffs(imin,1),optCoeffs(imin,2),time_int);
dotsAcc=acceleration(optCoeffs(imin,1),optCoeffs(imin,2),time_int);

% визуализация
figure
hold on
plot(dotsPolynomX,dotsPolynomY,'Red', 'LineWidth', 1.5);
xlabel('x') 
ylabel('y') 
hold off

figure
hold on
plot(time_int,dotsVel,'Blue', 'LineWidth', 1.5);
xlabel('t') 
ylabel('v') 
hold off

figure
hold on
plot(time_int,dotsAcc,'Blue', 'LineWidth', 1.5);
xlabel('t') 
ylabel('a') 
hold off
end

function [c,ceq] = constraints(x)
global c1handle;
global c2handle;
global timemoment;
c =[c1handle(x(1),x(2),timemoment);c2handle(x(1),x(2),timemoment)];
ceq = [];
end

function f = fun(x)
global Jhandle;
f =Jhandle(x(1),x(2));
end

 