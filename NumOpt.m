%% исходные данные
clc; clear; tic;
xData = [0,0,1,0]; % x = [x0, x'0, xt, x't]
yData = [0,0,1,0]; % y = [y0, y'0, yt, y't] 
T = 3;
mechA = [0.3, 0.4, 0.08]; 
mechB = [0.6, 0.5, 0.08]; 
mechC = [0.7, 0.75, 0.08];
maxVelocity = 0.6;
maxAcceleration=2;
tnet = 0:0.01:T;
tnet=tnet';

global tdegr; tdegr=[ones(size(tnet)),tnet,(tnet.^2),(tnet.^3),(tnet.^4)];
global tdegrVel; tdegrVel=[ones(size(tnet)),tnet,(tnet.^2),(tnet.^3),(tnet.^4), (tnet.^5), (tnet.^6)];
global tdegrAcc; tdegrAcc=[ones(size(tnet)),tnet,(tnet.^2),(tnet.^3),(tnet.^4)];

% нахождение критических параметров 
weights=[1, 1, 1, 1, 1];
% [vcrit2, acrit2, r2critA, r2critB, r2critC] = findOptParams2(xData, yData, T, mechA, mechB, mechC, weights);
[vcrit2, acrit2, r2critA, r2critB, r2critC] = findOptParams(xData, yData, T, mechA, mechB, mechC);

%решаем задачу для своих данных
maxVelocity = sqrt(vcrit2)*1.1;
% maxAcceleration=sqrt(acrit2)*1.1;
% mechA(3)=sqrt(r2critA);
% mechB(3)=sqrt(r2critB);
% mechC(3)=sqrt(r2critC);

maxVelocityList =  maxVelocity*ones(1, length(tnet));
maxAccelerationList =  maxAcceleration*ones(1, length(tnet));

maxv=@(a,b)(max(velocity2(xData, yData, T, a, b)));
maxa=@(a,b)(max(acceleration2(xData, yData, T, a, b)));
minr2A=@(a,b)(min((trajectoryX(xData, T, a)-mechA(1)).^2+(trajectoryY(yData, T,  b)-mechA(2)).^2));
minr2B=@(a,b)(min((trajectoryX(xData, T, a)-mechB(1)).^2+(trajectoryY(yData, T,  b)-mechB(2)).^2));
minr2C=@(a,b)(min((trajectoryX(xData, T, a)-mechC(1)).^2+(trajectoryY(yData, T,  b)-mechC(2)).^2));

constrv=@(a,b)(maxv(a,b)-maxVelocity);
constra=@(a,b)(maxa(a,b)-maxAcceleration);
constrmexA=@(a,b)(minr2A(a,b)-mechA(3)^2-0.01);%-0.01
constrmexB=@(a,b)(minr2B(a,b)-mechB(3)^2-0.01);%-0.01
constrmexC=@(a,b)(minr2C(a,b)-mechC(3)^2-0.01);%-0.01

global maxvhandle; maxvhandle=maxv;
global maxahandle; maxahandle=maxa;
global minr2Ahandle; minr2Ahandle=minr2A;
global minr2Bhandle; minr2Bhandle=minr2B;
global minr2Chandle; minr2Chandle=minr2C;

global constrvhandle;  constrvhandle=constrv;
global constrahandle; constrahandle=constra;
global constrmexAhandle; constrmexAhandle=constrmexA;
global constrmexBhandle; constrmexBhandle=constrmexB;
global constrmexChandle; constrmexChandle=constrmexC;

%% функционал потерь
% lambda=1000;
% penalty= max([0;mechA(3)^2-min((trajectoryX(xData, T, a4)-mechA(1)).^2+(trajectoryY(yData, T,  b4)-mechA(2)).^2)]);
J = @(a,b)(trapz(tnet,(trajectoryX(xData, T, a)).^2+(trajectoryY(yData, T,  b)).^2+...
    velocity2(xData, yData, T, a, b)+acceleration2(xData, yData, T, a, b))...
    );%
%     +lambda*max([0;mechA(3)^2-min((trajectoryX(xData, T, a)-mechA(1)).^2+(trajectoryY(yData, T,  b)-mechA(2)).^2)])...
%     +lambda*max([0;mechB(3)^2-min((trajectoryX(xData, T, a)-mechB(1)).^2+(trajectoryY(yData, T,  b)-mechB(2)).^2)])...
% );

global Jhandle; Jhandle = J;

gs = GlobalSearch;
problem = createOptimProblem('fmincon', 'x0', [-0.01,-0.01], 'objective', @costFun, 'nonlcon', @constrglobal);
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
plot (tnet,velocity2(xData, yData, T, minCoeffs(1), minCoeffs(2)),'Red', 'LineWidth', 1.5);
plot(tnet,maxVelocityList,'Red', 'LineWidth', 0.5);
hold off;

figure
grid on;
hold on;
plot (tnet,acceleration2(xData, yData, T, minCoeffs(1), minCoeffs(2)),'Red', 'LineWidth', 1.5);
plot(tnet,maxAccelerationList,'Red', 'LineWidth', 0.5);
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

function f = velocityX(dataX, T, a4)
a2 = -(3*dataX(1) - 3*dataX(3) + 3*T*dataX(2) + T*dataX(4) - T^4*a4 - T^2*dataX(2))/T^2;
a3 =  (2*dataX(1) - 2*dataX(3) + 2*T*dataX(2) + T*dataX(4) - 2*T^4*a4 - T^2*dataX(2))/T^3;
koeffs=[dataX(2);2*a2;3*a3;4*a4;0];
global tdegr;
f=tdegr*koeffs;
end

function f = velocityY(dataY, T,  b4)
b2 = -(3*dataY(1) - 3*dataY(3) + 3*T*dataY(2) + T*dataY(4) - T^4*b4 - T^2*dataY(2))/T^2;
b3 =  (2*dataY(1) - 2*dataY(3) + 2*T*dataY(2) + T*dataY(4) - 2*T^4*b4 - T^2*dataY(2))/T^3;
koeffs=[dataY(2);2*b2;3*b3;4*b4;0];
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

function f = funVmax(x)
global maxvhandle;
f =maxvhandle(x(1),x(2));
end

function f = funAmax(x)
global maxahandle;
f =maxahandle(x(1),x(2));
end

function f = funr2Amax(x)
global minr2Ahandle;
f =-minr2Ahandle(x(1),x(2));
end

function f = funr2Bmax(x)
global minr2Bhandle;
f =-minr2Bhandle(x(1),x(2));
end

function f = funr2Cmax(x)
global minr2Bhandle;
f =-minr2Bhandle(x(1),x(2));
end

function [c,ceq] = constrglobV(x)
global constrvhandle;
c = double(constrvhandle(x(1),x(2)));
ceq = [];
end

function [c,ceq] = constrglobA(x)
global constrahandle;
global constrvhandle;
c =[double(constrahandle(x(1),x(2))); double(constrvhandle(x(1),x(2)))];
ceq = [];
end

function [c,ceq] = constrglobB(x)
global constrahandle;
global constrvhandle;
global constrmexAhandle;
lambda=10000;
c =[double(constrahandle(x(1),x(2))); double(constrvhandle(x(1),x(2)));-lambda*constrmexAhandle(x(1),x(2))];
ceq = [];
end

function [c,ceq] = constrglobC(x)
global constrahandle;
global constrvhandle;
global constrmexAhandle;
global constrmexBhandle;
lambda=10000;
c =[double(constrahandle(x(1),x(2))); double(constrvhandle(x(1),x(2)));...
    -lambda*constrmexAhandle(x(1),x(2)); -lambda*constrmexBhandle(x(1),x(2))];
ceq = [];
end

function [c,ceq] = constrglobal(x)
global constrahandle;
global constrvhandle;
global constrmexAhandle;
global constrmexBhandle;
global constrmexChandle;
lambda=10000;
c =[double(constrahandle(x(1),x(2))); double(constrvhandle(x(1),x(2)));...
    -lambda*constrmexAhandle(x(1),x(2)); -lambda*constrmexBhandle(x(1),x(2)); -lambda*constrmexChandle(x(1),x(2))];
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

function [vcrit2, acrit2, r2critA, r2critB, r2critC] = findOptParams(xData, yData, T, mechA, mechB, mechC)
maxv=@(a,b)(max(velocity2(xData, yData, T, a, b)));
maxa=@(a,b)(max(acceleration2(xData, yData, T, a, b)));
minr2A=@(a,b)(min((trajectoryX(xData, T, a)-mechA(1)).^2+(trajectoryY(yData, T,  b)-mechA(2)).^2));
minr2B=@(a,b)(min((trajectoryX(xData, T, a)-mechB(1)).^2+(trajectoryY(yData, T,  b)-mechB(2)).^2));
minr2C=@(a,b)(min((trajectoryX(xData, T, a)-mechC(1)).^2+(trajectoryY(yData, T,  b)-mechC(2)).^2));

global maxvhandle; maxvhandle=maxv;
global maxahandle; maxahandle=maxa;
global minr2Ahandle; minr2Ahandle=minr2A;
global minr2Bhandle; minr2Bhandle=minr2B;
global minr2Chandle; minr2Chandle=minr2C;

% определяем минимально возможную скорость при заданных гу
gs = GlobalSearch;
problem = createOptimProblem('fmincon', 'x0', [-0.01,-0.01], 'objective', @funVmax, 'lb', [-10;-10],'ub', [10;10]);
minCoeffs = run(gs,problem);
vcrit2=maxv(minCoeffs(1),minCoeffs(2));

constrv=@(a,b)(maxv(a,b)-1.2*vcrit2);
global constrvhandle; 
constrvhandle=constrv;

% определяем минимально возможное  ускорение при заданных v и гу
problem = createOptimProblem('fmincon', 'x0', [-0.01,-0.01], 'objective', @funAmax, 'lb', [-10;-10],'ub', [10;10],'nonlcon', @constrglobV);
minCoeffs = run(gs,problem);
acrit2=maxa(minCoeffs(1),minCoeffs(2));
constra=@(a,b)(maxa(a,b)-1.2*acrit2);

global constrahandle;
constrahandle=constra;

% максимальный размер препятствия А с учетом его расположения и предыдущих 3 ограничений + гу
problem = createOptimProblem('fmincon', 'x0', [-0.01,-0.01], 'objective', @funr2Amax, 'lb', [-10;-10],'ub', [10;10],'nonlcon', @constrglobA);
minCoeffs = run(gs,problem);
r2critA=minr2A(minCoeffs(1),minCoeffs(2));
% mechA(3)=0.8*sqrt(r2critA);
constrmexA=@(a,b)(minr2A(a,b)-mechA(3)^2-0.01);
% constrmexA=@(a,b)(minr2A(a,b)-mechA(3)^2);
global constrmexAhandle;
constrmexAhandle=constrmexA;

% максимальный размер препятствия В с учетом его расположения и предыдущих 3 ограничений + гу
problem = createOptimProblem('fmincon', 'x0', [-0.01,-0.01], 'objective', @funr2Bmax, 'lb', [-10;-10],'ub', [10;10],'nonlcon', @constrglobB);
minCoeffs = run(gs,problem);
r2critB=minr2B(minCoeffs(1),minCoeffs(2));
% mechB(3)=0.8*sqrt(r2critB);

constrmexB=@(a,b)(minr2B(a,b)-mechB(3)^2-0.01);

global constrmexBhandle;
constrmexBhandle=constrmexB;

%  максимальный размер препятствия С с учетом его расположения и предыдущих 4 ограничений + гу
problem = createOptimProblem('fmincon', 'x0', [-0.01,-0.01], 'objective', @funr2Cmax, 'lb', [-10;-10],'ub', [10;10],'nonlcon', @constrglobC);
minCoeffs = run(gs,problem);
r2critC=minr2C(minCoeffs(1),minCoeffs(2));
% mechC(3)=0.8*sqrt(r2critC);
constrmexC=@(a,b)(minr2C(a,b)-mechC(3)^2-0.01);

global constrmexChandle;
constrmexChandle=constrmexC;

end

function [vcrit2, acrit2, r2critA, r2critB, r2critC] = findOptParams2(xData, yData, T, mechA, mechB, mechC, weights)
maxv=@(a,b)(max(velocity2(xData, yData, T, a, b)));
maxa=@(a,b)(max(acceleration2(xData, yData, T, a, b)));
minr2A=@(a,b)(min((trajectoryX(xData, T, a)-mechA(1)).^2+(trajectoryY(yData, T,  b)-mechA(2)).^2));
minr2B=@(a,b)(min((trajectoryX(xData, T, a)-mechB(1)).^2+(trajectoryY(yData, T,  b)-mechB(2)).^2));
minr2C=@(a,b)(min((trajectoryX(xData, T, a)-mechC(1)).^2+(trajectoryY(yData, T,  b)-mechC(2)).^2));
global maxvhandle; maxvhandle=maxv;
global maxahandle; maxahandle=maxa;
global minr2Ahandle; minr2Ahandle=minr2A;
global minr2Bhandle; minr2Bhandle=minr2B;
global minr2Chandle; minr2Chandle=minr2C;

J = @(a,b)( weights(1)*max(velocity2(xData, yData, T, a, b))+weights(2)*max(acceleration2(xData, yData, T, a, b))...
             -weights(3)*min((trajectoryX(xData, T, a)-mechA(1)).^2+(trajectoryY(yData, T,  b)-mechA(2)).^2)...
             -weights(4)*min((trajectoryX(xData, T, a)-mechB(1)).^2+(trajectoryY(yData, T,  b)-mechB(2)).^2)...
             -weights(5)*min((trajectoryX(xData, T, a)-mechC(1)).^2+(trajectoryY(yData, T,  b)-mechC(2)).^2)...
 );

global Jhandle; Jhandle = J;
gs = GlobalSearch;
problem = createOptimProblem('fmincon', 'x0', [-0.01,-0.01], 'objective', @costFun, 'lb', [-10;-10],'ub', [10;10]);
minCoeffs = run(gs,problem);

vcrit2=max(velocity2(xData, yData, T, minCoeffs(1), minCoeffs(2)));
acrit2=max(acceleration2(xData, yData, T, minCoeffs(1), minCoeffs(2)));
r2critA=min((trajectoryX(xData, T, minCoeffs(1))-mechA(1)).^2+(trajectoryY(yData, T,  minCoeffs(2))-mechA(2)).^2);
r2critB=min((trajectoryX(xData, T, minCoeffs(1))-mechB(1)).^2+(trajectoryY(yData, T,  minCoeffs(2))-mechB(2)).^2);
r2critC=min((trajectoryX(xData, T, minCoeffs(1))-mechC(1)).^2+(trajectoryY(yData, T,  minCoeffs(2))-mechC(2)).^2);

end
