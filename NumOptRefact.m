%% исходные данные
clc; clear; tic;
xData = [0,0,1,0]; % x = [x0, x'0, xt, x't]
yData = [0,0,1,0]; % y = [y0, y'0, yt, y't] 
T = 3;
mechA = [0.3, 0.4, 0.08]; 
mechB = [0.6, 0.5, 0.08]; 
mechC = [0.7, 0.75, 0.08];
tnet = 0:0.01:T;
tnet=tnet';

gs = GlobalSearch;
global tdegr; tdegr=[ones(size(tnet)),tnet,(tnet.^2),(tnet.^3),(tnet.^4)];
global tdegrVel; tdegrVel=[ones(size(tnet)),tnet,(tnet.^2),(tnet.^3),(tnet.^4), (tnet.^5), (tnet.^6)];
global tdegrAcc; tdegrAcc=[ones(size(tnet)),tnet,(tnet.^2),(tnet.^3),(tnet.^4)];

global maxVelocity2; maxVelocity2 = 0.8;
global maxAcceleration2; maxAcceleration2 = 2;
global xgu; xgu=xData;
global ygu; ygu=yData;
global Tgu; Tgu=T;
global AMech; AMech=mechA; 
global BMech; BMech=mechB;
global CMech; CMech=mechC;
global WeightConstr; WeightConstr=[1;1;1;1;1]; % ограничения на v, a, A, B, C

%% вычисление критических параметров и подготовка ограничений
% ищем последовательно наиболее жесткие ограничения
%{
% %ищем наиболее жесткое ограничение по скорости
% weights=[1, 0, 0, 0, 0];
% WeightConstr=[0;0;1;1;1];
% [vcrit, acrit, r2critA, r2critB, r2critC] = findOptParams2(xData, yData, T, mechA, mechB, mechC,weights);
% maxVelocity2 = 1.1*vcrit;
% 
% %ищем наиболее жесткое ограничение по ускорению
% weights=[0, 1, 0, 0, 0];
% WeightConstr=[1;0;1;1;1];
% [vcrit, acrit, r2critA, r2critB, r2critC] = findOptParams2(xData, yData, T, mechA, mechB, mechC,weights);
% maxAcceleration2=1.1*acrit;
% %ищем наиболее жесткое ограничение по радиусу A
% weights=[0, 0, 1, 0, 0];
% WeightConstr=[1;1;0;1;1];
% [vcrit, acrit, r2critA, r2critB, r2critC] = findOptParams2(xData, yData, T, mechA, mechB, mechC,weights);
% mechA(3)=sqrt(0.8*r2critA);
% AMech=mechA;BMech=mechB;CMech=mechC;
% %ищем наиболее жесткое ограничение по радиусу  B
% weights=[0, 0, 0, 1, 0];
% WeightConstr=[1;1;1;0;1];
% [vcrit, acrit, r2critA, r2critB, r2critC] = findOptParams2(xData, yData, T, mechA, mechB, mechC,weights);
% mechB(3)=sqrt(0.8*r2critB);
% AMech=mechA;BMech=mechB;CMech=mechC;
% %ищем наиболее жесткое ограничение по радиусу  C
% weights=[0, 0, 0, 0, 1];
% WeightConstr=[1;1;1;1;0];
% [vcrit, acrit, r2critA, r2critB, r2critC] = findOptParams2(xData, yData, T, mechA, mechB, mechC,weights);
% mechC(3)=sqrt(0.8*r2critC);
% AMech=mechA;BMech=mechB;CMech=mechC;
%}

% ищем единоразово 
weights=[0.1, 1, 6, 3, 6];
[vcrit, acrit, r2critA, r2critB, r2critC] = findOptParams2(xData, yData, T, mechA, mechB, mechC,weights);

% maxVelocity2 = vcrit*1.5;
% maxAcceleration2=acrit*1.5;
% mechA(3)=sqrt(r2critA);
% mechB(3)=sqrt(r2critB);
% mechC(3)=sqrt(r2critC);

maxVelocity2List =  maxVelocity2*ones(1, length(tnet));
maxAcceleration2List =  maxAcceleration2*ones(1, length(tnet));

%% минимизация целевого функционала
% lambda=1000;
% penalty= max([0;mechA(3)^2-min((trajectoryX(xData, T, a4)-mechA(1)).^2+(trajectoryY(yData, T,  b4)-mechA(2)).^2)]);
J = @(a,b)(trapz(tnet,(trajectoryX(xData, T, a)).^2+(trajectoryY(yData, T,  b)).^2+...
    velocity2(xData, yData, T, a, b)+acceleration2(xData, yData, T, a, b))...
    );%
%     +lambda*max([0;mechA(3)^2-min((trajectoryX(xData, T, a)-mechA(1)).^2+(trajectoryY(yData, T,  b)-mechA(2)).^2)])...
%     +lambda*max([0;mechB(3)^2-min((trajectoryX(xData, T, a)-mechB(1)).^2+(trajectoryY(yData, T,  b)-mechB(2)).^2)])...
% );

global Jhandle; Jhandle = J;
 
% problem = createOptimProblem('fmincon', 'x0', [-0.01,-0.01], 'objective', @costFun, 'lb', [-10;-10],'ub', [10;10],'nonlcon', @constrglobal);
problem = createOptimProblem('fmincon', 'x0', [-0.0,-0.0], 'objective', @costFun, 'lb', [-10;-10],'ub', [10;10],'nonlcon', @constrglobalrefact);
minCoeffs = run(gs,problem);
 
%% визуализация результатов
figure
hold on
grid on;
plot (trajectoryX(xData, T, minCoeffs(1)),trajectoryY(yData, T,  minCoeffs(2)),'Red', 'LineWidth', 1.5);
[x_circle,y_circle] = circle(mechA(1),mechA(2),mechA(3)); 
plot(x_circle,y_circle,'Black', 'linewidth',1.5); fill(x_circle,y_circle,mechA(3));
[x_circle,y_circle] = circle(mechB(1),mechB(2),mechB(3));
plot(x_circle,y_circle,'Black', 'linewidth',1.5); fill(x_circle,y_circle,mechB(3));
[x_circle,y_circle] = circle(mechC(1),mechC(2),mechC(3));
plot(x_circle,y_circle,'Black', 'linewidth',1.5); fill(x_circle,y_circle,mechC(3));
% xlim([0 1]) 
% ylim([0 1])
xlabel('x') 
ylabel('y') 
hold off

figure
grid on;
hold on;
plot(tnet,sqrt(velocity2(xData, yData, T, minCoeffs(1), minCoeffs(2))),'Blue', 'LineWidth', 1.5);
plot(tnet,sqrt(maxVelocity2List),'Red', 'LineWidth', 0.5);
ylim([0 sqrt(maxVelocity2)+0.2])
xlabel('t') 
ylabel('v') 
hold off;

figure
grid on;
hold on;
plot(tnet,sqrt(acceleration2(xData, yData, T, minCoeffs(1), minCoeffs(2))),'Blue', 'LineWidth', 1.5);
plot(tnet,sqrt(maxAcceleration2List),'Red', 'LineWidth', 0.5);   
ylim([0 sqrt(maxAcceleration2)+0.2])
xlabel('t') 
ylabel('a') 
hold off; 
toc;
%% 
function f = trajectoryX(dataX, T, a4)
a2 = -(3*dataX(1) - 3*dataX(3) + 3*T*dataX(2) + T*dataX(4) - T^4*a4 - T^2*dataX(2))/T^2;
a3 =  (2*dataX(1) - 2*dataX(3) + 2*T*dataX(2) + T*dataX(4) - 2*T^4*a4 - T^2*dataX(2))/T^3;
koeffs=[dataX(1);dataX(2);a2;a3;a4];
global tdegr;
f=tdegr*koeffs;
end

function f = trajectoryY(dataY, T,  b4)
b2 = -(3*dataY(1) - 3*dataY(3) + 3*T*dataY(2) + T*dataY(4) - T^4*b4 - T^2*dataY(2))/T^2;
b3 =  (2*dataY(1) - 2*dataY(3) + 2*T*dataY(2) + T*dataY(4) - 2*T^4*b4 - T^2*dataY(2))/T^3;
koeffs=[dataY(1);dataY(2);b2;b3;b4];
global tdegr;
f=tdegr*koeffs;
end

function f = velocity2(dataX, dataY, T, a4, b4)
a2 = -(3*dataX(1) - 3*dataX(3) + 3*T*dataX(2) + T*dataX(4) - T^4*a4 - T^2*dataX(2))/T^2;
a3 =  (2*dataX(1) - 2*dataX(3) + 2*T*dataX(2) + T*dataX(4) - 2*T^4*a4 - T^2*dataX(2))/T^3;
b2 = -(3*dataY(1) - 3*dataY(3) + 3*T*dataY(2) + T*dataY(4) - T^4*b4 - T^2*dataY(2))/T^2;
b3 =  (2*dataY(1) - 2*dataY(3) + 2*T*dataY(2) + T*dataY(4) - 2*T^4*b4 - T^2*dataY(2))/T^3;

koeffs=[dataX(2)^2 + dataY(2)^2; 4*a2*dataX(2) + 4*b2*dataY(2); 4*a2^2 + 4*b2^2 + 6*a3*dataX(2) + 6*b3*dataY(2);...
        12*a2*a3 + 12*b2*b3 + 8*a4*dataX(2) + 8*b4*dataY(2); 9*a3^2 + 9*b3^2 + 16*a2*a4 + 16*b2*b4;...
        24*a3*a4 + 24*b3*b4; 16*a4^2 + 16*b4^2];
global tdegrVel;
f=tdegrVel*koeffs;
end

function f = acceleration2(dataX, dataY, T, a4, b4)
a2 = -(3*dataX(1) - 3*dataX(3) + 3*T*dataX(2) + T*dataX(4) - T^4*a4 - T^2*dataX(2))/T^2;
a3 =  (2*dataX(1) - 2*dataX(3) + 2*T*dataX(2) + T*dataX(4) - 2*T^4*a4 - T^2*dataX(2))/T^3;
b2 = -(3*dataY(1) - 3*dataY(3) + 3*T*dataY(2) + T*dataY(4) - T^4*b4 - T^2*dataY(2))/T^2;
b3 =  (2*dataY(1) - 2*dataY(3) + 2*T*dataY(2) + T*dataY(4) - 2*T^4*b4 - T^2*dataY(2))/T^3;

koeffs=[4*a2^2 + 4*b2^2; 24*a2*a3 + 24*b2*b3; 36*a3^2 + 36*b3^2 + 48*a2*a4 + 48*b2*b4;...
        144*a3*a4 + 144*b3*b4; 144*a4^2 + 144*b4^2];
global tdegrAcc;
f=tdegrAcc*koeffs;
end

function [c,ceq] = constrglobalrefact(x)
global maxVelocity2;
global maxAcceleration2;
global xgu;
global ygu;
global Tgu;
global AMech; 
global BMech;
global CMech;
global WeightConstr;
lambda=10000;
c1=0;
if(abs(WeightConstr(1))>1e-20)
  c1=WeightConstr(1)*max(acceleration2(xgu, ygu, Tgu,x(1),x(2)))-maxAcceleration2;  
end
c2=0;
if(abs(WeightConstr(2))>1e-20)
  c2=WeightConstr(2)*max(velocity2(xgu, ygu, Tgu, x(1),x(2)))-maxVelocity2;  
end
c3=0;
if(abs(WeightConstr(3))>1e-20)
  c3=-WeightConstr(3)*lambda*(min((trajectoryX(xgu, Tgu, x(1))-AMech(1)).^2+(trajectoryY(ygu, Tgu,  x(2))-AMech(2)).^2)-AMech(3)^2-0.01);  
end
c4=0;
if(abs(WeightConstr(4))>1e-20)
  c4=-WeightConstr(4)*lambda*(min((trajectoryX(xgu, Tgu, x(1))-BMech(1)).^2+(trajectoryY(ygu, Tgu,  x(2))-BMech(2)).^2)-BMech(3)^2-0.01);  
end
c5=0;
if(abs(WeightConstr(5))>1e-20)
  c5=-WeightConstr(5)*lambda*(min((trajectoryX(xgu, Tgu, x(1))-CMech(1)).^2+(trajectoryY(ygu, Tgu,  x(2))-CMech(2)).^2)-CMech(3)^2-0.01);  
end
c =[c1;c2;c3;c4;c5];
ceq = [];
end

function [x_circle,y_circle] = circle(x,y,r)
angle=0:0.01:2*pi; 
xp=r*cos(angle);
yp=r*sin(angle);
x_circle = x + xp; 
y_circle = y + yp; 
end

function f = costFun(x)
global Jhandle;
f =Jhandle(x(1),x(2));
end

function [vcrit, acrit, r2critA, r2critB, r2critC] = findOptParams2(xData, yData, T, mechA, mechB, mechC,weights)
% boundCond, T, oA, oB, oC ---> v, a, rA, rB, rC
J = @(a,b)( weights(1)*max(velocity2(xData, yData, T, a, b))+weights(2)*max(acceleration2(xData, yData, T, a, b))...
             -weights(3)*min((trajectoryX(xData, T, a)-mechA(1)).^2+(trajectoryY(yData, T,  b)-mechA(2)).^2)...
             -weights(4)*min((trajectoryX(xData, T, a)-mechB(1)).^2+(trajectoryY(yData, T,  b)-mechB(2)).^2)...
             -weights(5)*min((trajectoryX(xData, T, a)-mechC(1)).^2+(trajectoryY(yData, T,  b)-mechC(2)).^2)...
 );

global Jhandle; Jhandle = J;
gs = GlobalSearch;
problem = createOptimProblem('fmincon', 'x0', [-0.01,-0.01], 'objective', @costFun, 'lb', [-10;-10],'ub', [10;10]);
minCoeffs = run(gs,problem);
vcrit=max(velocity2(xData, yData, T, minCoeffs(1), minCoeffs(2)));
acrit=max(acceleration2(xData, yData, T, minCoeffs(1), minCoeffs(2)));
r2critA=min((trajectoryX(xData, T, minCoeffs(1))-mechA(1)).^2+(trajectoryY(yData, T,  minCoeffs(2))-mechA(2)).^2);
r2critB=min((trajectoryX(xData, T, minCoeffs(1))-mechB(1)).^2+(trajectoryY(yData, T,  minCoeffs(2))-mechB(2)).^2);
r2critC=min((trajectoryX(xData, T, minCoeffs(1))-mechC(1)).^2+(trajectoryY(yData, T,  minCoeffs(2))-mechC(2)).^2);
end
