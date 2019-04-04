% function Polynom4MinimizationJ;
% clc; clear all;

% ��������� � �������� ��������� �������
xData = [0,0,1,0]; % x = [x0, x'0, xt, x't]
yData = [0,0,1,0]; % y = [y0, y'0, yt, y't] 

% ������������ �����������
mechA = [0.3, 0.4, 0.08]; 
mechB = [0.6, 0.5, 0.08]; 
mechC = [0.7, 0.75, 0.08];

time = 3;
time_int = 0:0.01:time;  
maxVelocity = 0.8;
maxAcceleration = 2;
initApprox = [0,0];

maxVelocityList =  maxVelocity*ones(1, length(time_int));
maxAccelerationList =  maxAcceleration*ones(1, length(time_int));

% ����� ��� ����������� Px � Py 4-�� �������� a2=a2(a4)...b3=b3(b4)
syms a4 b4
LeftSideX = [time^2, time^3; 2*time, 3*time^2];
RightSideX = [xData(3)-xData(1)-xData(2)*time-a4*time^4; xData(4)-xData(2)-4*a4*time^3];
LeftSideY = [time^2, time^3; 2*time, 3*time^2];
RightSideY = [yData(3)-yData(1)-yData(2)*time-b4*time^4; yData(4)-yData(2)-4*b4*time^3];

coeffsPolynomX = LeftSideX \ RightSideX;
coeffsPolynomY = LeftSideY \ RightSideY;

% ���������� Px � Py
syms t
polynomX =  xData(1)+xData(2)*t+t.^2*coeffsPolynomX(1)+t.^3*coeffsPolynomX(2)+t.^4*a4;
polynomY = yData(1)+yData(2)*t+t.^2*coeffsPolynomY(1)+t.^3*coeffsPolynomY(2)+t.^4*b4;

% ���������� Px', Px", Py', Py"
dpolynomX = diff(polynomX, t);
ddpolynomX = diff(dpolynomX, t);
dpolynomY = diff(polynomY, t);
ddpolynomY = diff(dpolynomY, t);

% ����������� �� �������� � ���������
velConstr = sqrt(dpolynomX^2+dpolynomY^2);
accConstr = sqrt(ddpolynomX^2+ddpolynomY^2);
velConstrM = matlabFunction(velConstr);
accConstrM =  matlabFunction(accConstr);

% ������������ �����������
mechConstrA = (polynomX-mechA(1))^2+(polynomY-mechA(2))^2;
mechConstrB = (polynomX-mechB(1))^2+(polynomY-mechB(2))^2;
mechConstrC = (polynomX-mechC(1))^2+(polynomY-mechC(2))^2;

% ���������� ������
intExpression = polynomX^2+polynomY^2+dpolynomX^2+dpolynomY^2+ddpolynomX^2+ddpolynomY^2;
J = 0.5*int(intExpression,0,time);
global Jhandle; Jhandle = matlabFunction(J);

% ����������� �� �������� 
dvelConstr=diff(velConstr, t);
eqn = dvelConstr == 0;
solt = solve(eqn,t);
leftVel=sqrt(xData(2)^2+yData(2)^2);
rightVel=sqrt(xData(4)^2+yData(4)^2);

global vstar;
vstar = @(a,b) (-maxVelocity + vpa(max([leftVel,rightVel,... 
               subs(velConstr,[a4,b4,t],[a,b,subs(solt(1),[a4,b4],[a, b])]),...
               subs(velConstr,[a4,b4,t],[a,b,subs(solt(2),[a4,b4],[a, b])]),...
               subs(velConstr,[a4,b4,t],[a,b,subs(solt(3),[a4,b4],[a, b])])])));
           
% ����������� �� ���������                               
daccConstr = diff(accConstr, t);
eqn = daccConstr == 0;
solt = solve(eqn,t);
leftAcc = accConstrM(a4,b4,0);
rightAcc = accConstrM(a4,b4,time);

global astar;
astar = @(a,b) (-maxAcceleration + vpa(max([
               subs(leftAcc,[a4,b4],[a,b]),...
               subs(rightAcc,[a4,b4],[a,b]),...
               subs(accConstr,[a4,b4,t],[a,b,subs(solt(1),[a4,b4],[a, b])]),...
               subs(accConstr,[a4,b4,t],[a,b,subs(solt(2),[a4,b4],[a, b])]),...
               subs(accConstr,[a4,b4,t],[a,b,subs(solt(3),[a4,b4],[a, b])])])));
%% ������������ ����������� �   
 dmechConstrA = diff(mechConstrA, t);
 eqn = dmechConstrA == 0;
 eqnnd = mechConstrA == 0;
 solt = solve(eqn,t);
 soltnd] = solve(eqnnd,t);
 
 global mechAstar;
 mechAstar = @(a,b) (mechA(3)^2 - min([
                    isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(solt(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrA,[a4,b4,t],[a,b,subs(solt(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(solt(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(solt(2),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrA,[a4,b4,t],[a,b,subs(solt(2),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(solt(2),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(solt(3),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrA,[a4,b4,t],[a,b,subs(solt(3),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(solt(3),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(solt(4),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrA,[a4,b4,t],[a,b,subs(solt(4),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(solt(4),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(solt(5),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrA,[a4,b4,t],[a,b,subs(solt(5),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(solt(5),[a4,b4],[a, b])]))))*1,...
                    ...
                    isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrA,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1])^2);
%% ������������ ����������� B   
 % ������������ ����������� B
  dmechConstrB = diff(mechConstrB, t);
 eqn = dmechConstrB == 0;
 solt = solve(eqn,t);
 
 global mechBstar;
 mechBstar = @(a,b) (mechB(3)^2 - min([
                    isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(solt(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrB,[a4,b4,t],[a,b,subs(solt(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(solt(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(solt(2),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrB,[a4,b4,t],[a,b,subs(solt(2),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(solt(2),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(solt(3),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrB,[a4,b4,t],[a,b,subs(solt(3),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(solt(3),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(solt(4),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrB,[a4,b4,t],[a,b,subs(solt(4),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(solt(4),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(solt(5),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrB,[a4,b4,t],[a,b,subs(solt(5),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(solt(5),[a4,b4],[a, b])]))))*1,...
                    ...
                    isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrB,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1])^2);
%% ������������ ����������� C   
 % ������������ ����������� C

 dmechConstrC = diff(mechConstrC, t);
 eqn = dmechConstrC == 0;
 solt = solve(eqn,t);
 
 global mechCstar;
 mechCstar = @(a,b) (mechC(3)^2 - min([
                    isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(solt(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrC,[a4,b4,t],[a,b,subs(solt(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(solt(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(solt(2),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrC,[a4,b4,t],[a,b,subs(solt(2),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(solt(2),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(solt(3),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrC,[a4,b4,t],[a,b,subs(solt(3),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(solt(3),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(solt(4),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrC,[a4,b4,t],[a,b,subs(solt(4),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(solt(4),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(solt(5),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrC,[a4,b4,t],[a,b,subs(solt(5),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(solt(5),[a4,b4],[a, b])]))))*1,...
                    ...
                    isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1,...
                    isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))*...
                    vpa((subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])])))+...
                    (~isreal(vpa(subs(mechConstrC,[a4,b4,t],[a,b,subs(soltnd(1),[a4,b4],[a, b])]))))*1])^2);
%%
 optCoeffs = fmincon(@fun,[-0.01,-0.01],[],[],[],[],[],[],@constrglob);            
% gs = GlobalSearch;
% problem = createOptimProblem('fmincon', 'x0', [-0.01,-0.01], 'objective', @fun, 'nonlcon', @constrglob);
% minCoeffs = run(gs,problem);                                 

% ����������, ��������, ��������� ����� ���������� �����
polynomXopt=matlabFunction(polynomX);
polynomYopt=matlabFunction(polynomY);

velocity=matlabFunction(velConstr);
acceleration=matlabFunction(accConstr);

dotsPolynomX=polynomXopt(optCoeffs(1),time_int);
dotsPolynomY=polynomYopt(optCoeffs(2),time_int);

dotsVel=velocity(optCoeffs(1),optCoeffs(2),time_int);
dotsAcc=acceleration(optCoeffs(1),optCoeffs(2),time_int);

mechConstrAM = matlabFunction(mechConstrA);
mechConstrBM = matlabFunction(mechConstrB);

% ������������
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
plot(time_int,dotsVel,'Blue', 'LineWidth', 1.5);
plot(time_int,maxVelocityList,'Red', 'LineWidth', 0.5);
xlabel('t') 
ylabel('v') 
hold off

figure
hold on
plot(time_int,dotsAcc,'Blue', 'LineWidth', 1.5);
plot(time_int,maxAccelerationList,'Red', 'LineWidth', 0.5);
ylim([0 maxAcceleration+0.5])
xlabel('t') 
ylabel('a') 
hold off

%{ 
������������ ����������� + ��������
     Z=zeros(21,21); for X=-10:10
        for Y=-10:10
            if (not(X==0 && Y==0))
               Z(X+11,Y+11)=vstar(X,Y);
            end
        end
    end
X=-10:10; Y=-10:10; surf(-10:10,-10:10,Z) ; ind=1:length(X);

h=0.1; vstarp =@(a,b)(  (
Z(min(ind(X>=a)),min(ind(Y>=b)))*exp(-((X(min(ind(X>=a)))-a)/h)^2-((Y(min(ind(Y>=b)))-b)/h)^2)...
                     +Z(min(ind(X>=a)),max(ind(Y<b)))*exp(-((X(min(ind(X>=a)))-a)/h)^2-((Y(max(ind(Y<b)))-b)/h)^2)...
                     +Z(max(ind(X<a)),max(ind(Y<b)))*exp(-((X(max(ind(X<a)))-a)/h)^2-((Y(max(ind(Y<b)))-b)/h)^2)...
                     +Z(max(ind(X<a)),min(ind(Y>=b)))*exp(-((X(max(ind(X<a)))-a)/h)^2-((Y(min(ind(Y>=b)))-b)/h)^2))...
                     /(exp(-((X(min(ind(X>=a)))-a)/h)^2-((Y(min(ind(Y>=b)))-b)/h)^2)+exp(-((X(min(ind(X>=a)))-a)/h)^2-((Y(max(ind(Y<b)))-b)/h)^2)...
                       +exp(-((X(max(ind(X<a)))-a)/h)^2-((Y(max(ind(Y<b)))-b)/h)^2)+exp(-((X(max(ind(X<a)))-a)/h)^2-((Y(min(ind(Y>=b)))-b)/h)^2)));
vstarpx =@(x)([vstarp(x(1),x(2))-maxVelocity,0]);


Z2=zeros(19,19);
    for XX=-9:9
        for YY=-9:9
            if (not(XX==0 && YY==0))
               Z2(XX+10,YY+10)=vstarp(XX,YY);%-Z(XX+10,YY+10);
            end
        end
    end
figure; surf(-9:9,-9:9,Z2) ; figure; surf(-10:10,-10:10,Z) ;
%}


function [c,ceq] = constrglob(x)
global vstar;
global astar;
global mechAstar;
global mechBstar;
global mechCstar;
c =[double(vstar(x(1),x(2))); double(astar(x(1),x(2)));...
    double(mechAstar(x(1),x(2))); double(mechBstar(x(1),x(2))); double(mechCstar(x(1),x(2)))];
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

 