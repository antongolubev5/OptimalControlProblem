%% исходные данные
clc; clear; tic;

xData = [0,0,1,0]; % x = [x0, x'0, xt, x't]
yData = [0,0,1,0]; % y = [y0, y'0, yt, y't] 

mechA = [0.3, 0.4, 0.08]; 
mechB = [0.6, 0.5, 0.08]; 
mechC = [0.7, 0.75, 0.08];

time = 3;
maxVelocity = 0.8;
maxAcceleration = 2;

h = 0.01; time_int = 0:h:time; 
optCoeffs = zeros(round(time/h)+1,2); 

maxVelocityList =  maxVelocity*ones(1, length(time_int));
maxAccelerationList =  maxAcceleration*ones(1, length(time_int));

%% построение многочленов и ограничений
syms a4 b4
LeftSideX = [time^2, time^3; 2*time, 3*time^2];
RightSideX = [xData(3)-xData(1)-xData(2)*time-a4*time^4; xData(4)-xData(2)-4*a4*time^3];
LeftSideY = [time^2, time^3; 2*time, 3*time^2];
RightSideY = [yData(3)-yData(1)-yData(2)*time-b4*time^4; yData(4)-yData(2)-4*b4*time^3];

coeffsPolynomX = LeftSideX \ RightSideX;
coeffsPolynomY = LeftSideY \ RightSideY;

% многочлены Px и Py
syms t
polynomX =  xData(1)+xData(2)*t+t.^2*coeffsPolynomX(1)+t.^3*coeffsPolynomX(2)+t.^4*a4;
polynomY = yData(1)+yData(2)*t+t.^2*coeffsPolynomY(1)+t.^3*coeffsPolynomY(2)+t.^4*b4;

% нахождение Px', Px", Py', Py"
dpolynomX = diff(polynomX, t);
ddpolynomX = diff(dpolynomX, t);
dpolynomY = diff(polynomY, t);
ddpolynomY = diff(dpolynomY, t);

% ограничения на скорость и ускорение
velConstr = sqrt(dpolynomX^2+dpolynomY^2);
accConstr = sqrt(ddpolynomX^2+ddpolynomY^2);

% механические ограничения
mechConstrA = (polynomX-mechA(1))^2+(polynomY-mechA(2))^2;
mechConstrB = (polynomX-mechB(1))^2+(polynomY-mechB(2))^2;
mechConstrC = (polynomX-mechC(1))^2+(polynomY-mechC(2))^2;

% целевой функционал 
intExpression = polynomX^2+polynomY^2+dpolynomX^2+dpolynomY^2+ddpolynomX^2+ddpolynomY^2;
J = 0.5*int(intExpression,0,time);

global timemoment;
global Jhandle; Jhandle = matlabFunction(J);
global c1handle; c1handle=matlabFunction(-velConstr-maxVelocity);
global c2handle; c2handle=matlabFunction(-accConstr-maxAcceleration);
global c3handle; c3handle=matlabFunction(-mechConstrA+mechA(3)^2+0.003);
global c4handle; c4handle=matlabFunction(-mechConstrB+mechB(3)^2+0.003);
global c5handle; c5handle=matlabFunction(-mechConstrC+mechC(3)^2+0.003);

i = 1;
for tparam = 0:h:time
timemoment=tparam;
nonlcon=@constraints;
rng default;
gs = GlobalSearch;
problem = createOptimProblem('fmincon', 'x0', [1,1], 'objective', @fun, 'nonlcon', nonlcon);
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

velocityXopt=matlabFunction(dpolynomX);
velocityYopt=matlabFunction(dpolynomY);
accelerationXopt=matlabFunction(ddpolynomX);
accelerationYopt=matlabFunction(ddpolynomY);

dotsPolynomX = zeros(round(time/h)+1,1); 
dotsPolynomY = zeros(round(time/h)+1,1);
dotsVel = [zeros(round(time/h)+1,1), zeros(round(time/h)+1,1), zeros(round(time/h)+1,1)]; % vx, vy, |v|
dotsAcc = [zeros(round(time/h)+1,1), zeros(round(time/h)+1,1), zeros(round(time/h)+1,1)]; % ax, ay, |a|

cnt = 1;
for tparam = 0:h:time
 dotsPolynomX(cnt)=polynomXopt(optCoeffs(cnt,1),tparam);
 dotsPolynomY(cnt)=polynomYopt(optCoeffs(cnt,2),tparam);
 
 dotsVel(cnt,1)=velocityXopt(optCoeffs(cnt,1),tparam);
 dotsVel(cnt,2)=velocityYopt(optCoeffs(cnt,2),tparam);
 dotsVel(cnt,3)=velocityOpt(optCoeffs(cnt,1),optCoeffs(cnt,2),tparam);
 
 dotsAcc(cnt,1)=accelerationXopt(optCoeffs(cnt,1),tparam);
 dotsAcc(cnt,2)=accelerationYopt(optCoeffs(cnt,2),tparam);
 dotsAcc(cnt,3)=accelerationOpt(optCoeffs(cnt,1),optCoeffs(cnt,2),tparam);
           
 cnt=cnt+1;
end

%% визуализация результатов  
figure
hold on
grid on;
plot(dotsPolynomX,dotsPolynomY,'Red', 'LineWidth', 1.5);
[x_circle,y_circle] = circle(mechA(1),mechA(2),mechA(3)); 
plot(x_circle,y_circle,'Black', 'linewidth',1.5); fill(x_circle,y_circle,mechA(3));
[x_circle,y_circle] = circle(mechB(1),mechB(2),mechB(3)); 
plot(x_circle,y_circle,'Black', 'linewidth',1.5); fill(x_circle,y_circle,mechB(3));
[x_circle,y_circle] = circle(mechC(1),mechC(2),mechC(3)); 
plot(x_circle,y_circle,'Black', 'linewidth',1.5); fill(x_circle,y_circle,mechC(3));
xlabel('x') 
ylabel('y') 
hold off

figure
hold on
plot(time_int,dotsPolynomX,'Blue', 'LineWidth', 1.5);
xlabel('t') 
ylabel('x') 
hold off

figure
hold on
plot(time_int,dotsPolynomY,'Blue', 'LineWidth', 1.5);
xlabel('t') 
ylabel('y') 
hold off

figure
hold on
grid on;
plot(time_int,dotsVel(:,3),'Blue', 'LineWidth', 1.5);
plot(time_int,maxVelocityList,'Red', 'LineWidth', 0.5);
ylim([0 maxVelocity+0.1])
xlabel('t') 
ylabel('v') 
hold off

figure
hold on
plot(time_int,dotsVel(:,1),'Blue', 'LineWidth', 1.5);
xlabel('t') 
ylabel('v_x') 
hold off

figure
hold on
plot(time_int,dotsVel(:,2),'Blue', 'LineWidth', 1.5);
xlabel('t') 
ylabel('v_y') 
hold off

figure
hold on
grid on;
plot(time_int,dotsAcc(:,3),'Blue', 'LineWidth', 1.5);
plot(time_int,maxAccelerationList,'Red', 'LineWidth', 0.5);
ylim([0 maxAcceleration+0.1])
xlabel('t') 
ylabel('a') 
hold off

figure
hold on
plot(time_int,dotsAcc(:,1),'Blue', 'LineWidth', 1.5);
xlabel('t') 
ylabel('a_x') 
hold off

figure
hold on
plot(time_int,dotsAcc(:,2),'Blue', 'LineWidth', 1.5);
xlabel('t') 
ylabel('a_y') 
hold off
toc;
%% 

function [c,ceq] = constraints(x)
global c1handle;
global c2handle;
global c3handle;
global c4handle;
global c5handle;
global timemoment;
c =[c1handle(x(1),x(2),timemoment);c2handle(x(1),x(2),timemoment);c3handle(x(1),x(2),timemoment);
                                   c4handle(x(1),x(2),timemoment);c5handle(x(1),x(2),timemoment)];
ceq = [];
end

function [c,ceq] = constraintsMech(x)
global c3handle;
global c4handle;
global c5handle;
global timemoment;
c = [c3handle(x(1),x(2),timemoment);c4handle(x(1),x(2),timemoment);c5handle(x(1),x(2),timemoment)];
ceq = [];
end

function [c,ceq] = constraintsPhys(x)
global c1handle;
global c2handle;
global timemoment;
c = [c1handle(x(1),x(2),timemoment);c2handle(x(1),x(2),timemoment)];
ceq = [];
end

function f = fun(x)
global Jhandle;
f =Jhandle(x(1),x(2));
end

function [x_circle,y_circle] = circle(x,y,r)
angle=0:0.01:2*pi; 
xp=r*cos(angle);
yp=r*sin(angle);
x_circle = x + xp; 
y_circle = y + yp; 
end

 
