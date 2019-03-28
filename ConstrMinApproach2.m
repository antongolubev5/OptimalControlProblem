% 2 подход - решаем одну задачу оптимизации на всем интервале
clc; clear all;

% начальное и конечное положения системы, время 
x = [0,0,1,0]; % x = [x0, x'0, xt, x't]
y = [0,0,1,0]; % y = [y0, y'0, yt, y't] 

% механические ограничения
centerA = [0.3 0.4]; 
rA = 0.08;
centerB = [0.6 0.5]; 
rB = 0.08;

time = 10;
time_int = 0:0.01:time;  
maxVelocity = 0.6;
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
% 
% velConstrM = matlabFunction(velConstr);
% accConstrM =  matlabFunction(accConstr);

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
% global c3handle;
% global c4handle;
global timemoment;

c1handle=matlabFunction(-velConstr-maxVelocity);
c2handle=matlabFunction(-accConstr-maxAcceleration);
% c3handle=matlabFunction(-mechConstrA+rA^2);
% c4handle=matlabFunction(-mechConstrB+rB^2);

% минимизируем 4 функции ограничений
% ищем максимум по скорости как t*=t*(a4,b4) 
dvelConstr=diff(velConstr, t);
eqn = dvelConstr == 0;
solt = solve(eqn,t);
% получим 3 корня уравнения
r1=solt(1);
r2=solt(2);
r3=solt(3);

% скорости на левом и правом концах
leftVelocity=sqrt(x(2)^2+y(2)^2);
rightVelocity=sqrt(x(4)^2+y(4)^2);

% tstar=@(a,b) max([leftvelocity,rightvelocity,subs(solt(1),[a4,b4],[a, b]),subs(solt(2),[a4,b4],[a, b]),subs(solt(3),[a4,b4],[a, b])]);
% syms vstar a b;
vstar =@(a,b)( vpa(max([leftVelocity,rightVelocity, subs(velConstr,[a4,b4,t],[a,b,subs(solt(1),[a4,b4],[a, b])]),...
                                         subs(velConstr,[a4,b4,t],[a,b,subs(solt(2),[a4,b4],[a, b])]),...
                                         subs(velConstr,[a4,b4,t],[a,b,subs(solt(3),[a4,b4],[a, b])])]) ));
vstarvect =@(xx)( -maxVelocity+vpa(max([leftVelocity,rightVelocity, subs(velConstr,[a4,b4,t],[xx(1),xx(2),subs(solt(1),[a4,b4],[xx(1), xx(2)])]),...
                                         subs(velConstr,[a4,b4,t],[xx(1),xx(2),subs(solt(2),[a4,b4],[xx(1), xx(2)])]),...
                                         subs(velConstr,[a4,b4,t],[xx(1),xx(2),subs(solt(3),[a4,b4],[xx(1), xx(2)])])]) ));
                                     
                                     vstarvect([1,1])
     nonlcon=vstarvect;                                
minCoeffs = fmincon(@fun,[-1 -1],[],[],[],[],[],[],nonlcon);                                     
                                     
                                     
                                     
                                     
                                     
                                     
                                     
                                     
    Z=zeros(21,21);
    for X=-10:10
        for Y=-10:10
            if (not(X==0 || Y==0))
               Z(X+11,Y+11)=vstar(X,Y);
            end   
        end
    end

surf(X,Y,Z) ;                               

subs(velConstr,t,tstar);
% multi start - находит несколько локальных min
% global search - глобальный min
% https://www.mathworks.com/help/gads/how-globalsearch-and-multistart-work.html
% https://www.mathworks.com/help/gads/globalsearch.html
% https://www.mathworks.com/help/gads/createoptimproblem.html 


for tparam = 0:0.01:time
timemoment=tparam;
nonlcon=@constraints;
minCoeffs = fmincon(@fun,initApprox,[],[],[],[],[],[],nonlcon);

disturb1=c1handle(minCoeffs(1),minCoeffs(2),time_int) >0;% 0.01;
disturb2=c2handle(minCoeffs(1),minCoeffs(2),time_int) >0;% 0.01; 
disturb3=c3handle(minCoeffs(1),minCoeffs(2),time_int) > 0.0;
disturb4=c4handle(minCoeffs(1),minCoeffs(2),time_int) > 0.0;

if not(any(disturb1) || any(disturb2) || any(disturb3) || any(disturb4))
   optCoeffs =vertcat(optCoeffs,minCoeffs);
end
% if not(any(disturb3) || any(disturb4))
%    optCoeffs =vertcat(optCoeffs,minCoeffs);
% end
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

mechConstrAM = matlabFunction(mechConstrA);
mechConstrBM = matlabFunction(mechConstrB);

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
end

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

function f = fun(x)
global Jhandle;
f =Jhandle(x(1),x(2));
end

function f = funVelocity(x)
global c1handle;
f =c1handle(x(1),x(2),x(3));
end

function f = funAcceleration(x)
global c2handle;
f =c2handle(x(1),x(2),x(3));
end

function [x_circle,y_circle] = circle(x,y,r)
% построение механических ограничений
angle=0:0.01:2*pi; 
xp=r*cos(angle);
yp=r*sin(angle);
x_circle = x + xp; 
y_circle = y + yp; 
end

 
