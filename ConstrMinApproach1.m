% 1 подход - разбиваем задачу минимизации на несколькo
% последовательных задач
clc; clear all;
   
% начальное и конечное положения системы, время 
x = [0,0,1,0]; % x = [x0, x'0, xt, x't]
y = [0,0,1,0]; % y = [y0, y'0, yt, y't] 

% механические ограничения
centerA = [0.3 0.4]; 
rA = 0.08;
centerB = [0.6 0.5]; 
rB = 0.08;

time = 3;
maxVelocity = 0.6;
maxAcceleration = 2;
initApprox = [1,1];

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

% механические ограничения
mechConstrA = (polynomX-centerA(1))^2+(polynomY-centerA(2))^2;
mechConstrB = (polynomX-centerB(1))^2+(polynomY-centerB(2))^2;

% функционал потерь
intExpression = polynomX^2+polynomY^2+dpolynomX^2+dpolynomY^2+ddpolynomX^2+ddpolynomY^2;
J = 0.5*int(intExpression,0,time);

global Jhandle;
Jhandle = matlabFunction(J);

global c1handle;
global c2handle;
global c3handle;
global c4handle;
global timemoment;

c1handle=matlabFunction(-velConstr-maxVelocity);
c2handle=matlabFunction(-accConstr-maxAcceleration);
c3handle=matlabFunction(-mechConstrA+rA^2);
c4handle=matlabFunction(-mechConstrB+rB^2);

h = 0.1; time_int = 0:h:time; 
optCoeffs = zeros(round(time/h)+1,2); 

i = 1;
for tparam = 0:h:time
timemoment=tparam;
nonlcon=@constraintsMech;
rng default;
gs = GlobalSearch;
problem = createOptimProblem('fmincon', 'x0', initApprox, 'objective', @fun, 'nonlcon', nonlcon);
minCoeffs = run(gs,problem);
optCoeffs(i,1) = minCoeffs(1);
optCoeffs(i,2) = minCoeffs(2);
i=i+1;
end

 % траектории, скорость, ускорение через полученный коэфф
polynomXopt=matlabFunction(polynomX);
polynomYopt=matlabFunction(polynomY);
velocityOpt=matlabFunction(velConstr);
accelerationOpt=matlabFunction(accConstr);

dotsPolynomX = zeros(round(time/h)+1,1); 
dotsPolynomY = zeros(round(time/h)+1,1);
dotsVel = zeros(round(time/h)+1,1); 
dotsAcc = zeros(round(time/h)+1,1); 

cnt = 1;
for tparam = 0:h:time
 dotsPolynomX(cnt)=polynomXopt(optCoeffs(cnt,1),tparam);
 dotsPolynomY(cnt)=polynomYopt(optCoeffs(cnt,2),tparam);
 dotsVel(cnt)=velocityOpt(optCoeffs(cnt,1),optCoeffs(cnt,2),tparam);
 dotsAcc(cnt)=accelerationOpt(optCoeffs(cnt,1),optCoeffs(cnt,2),tparam);
 cnt=cnt+1;
end

% визуализация
figure
hold on
grid on;
plot(dotsPolynomX,dotsPolynomY,'Red', 'LineWidth', 1.5);
[x_circle,y_circle] = circle(centerA(1),centerA(2),rA); 
plot(x_circle,y_circle,'Black', 'linewidth',1.5); fill(x_circle,y_circle,rA);
[x_circle,y_circle] = circle(centerB(1),centerB(2),rB); 
plot(x_circle,y_circle,'Black', 'linewidth',1.5); fill(x_circle,y_circle,rB);
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


function [c,ceq] = constraints(x)
global c1handle;
global c2handle;
global c3handle;
global c4handle;
global timemoment;
c =[c1handle(x(1),x(2),timemoment);c2handle(x(1),x(2),timemoment);
    c3handle(x(1),x(2),timemoment);c4handle(x(1),x(2),timemoment)];
ceq = [];
end

function [c,ceq] = constraintsMech(x)
% global c1handle;
% global c2handle;
global c3handle;
global c4handle;
global timemoment;
c = [c3handle(x(1),x(2),timemoment);c4handle(x(1),x(2),timemoment)];
ceq = [];
end

function f = fun(x)
global Jhandle;
f =Jhandle(x(1),x(2));
end

function [x_circle,y_circle] = circle(x,y,r)
% построение механических ограничений
angle=0:0.01:2*pi; 
xp=r*cos(angle);
yp=r*sin(angle);
x_circle = x + xp; 
y_circle = y + yp; 
end

 
