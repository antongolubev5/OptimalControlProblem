%% начальные данные
clc
clear all

% начальное и конечное положения системы (0,0)->(1,1)
x = [0,0,1,0]; % x = [x0, x'0, xt, x't]
y = [0,0,1,0]; % y = [y0, y'0, yt, y't] 

% вспомогательная точка 
z = [0.5, 0.5, 0.5, 0.5]; % z = [z_x, z_y, Vz_x, Vz_y] 

%% первая терминальная задача
% выбираем время и обнуляем старшую степень x_1(t)
time1 = 2*(z(1)-x(1))/(x(4)+z(3));
time_int1 = 0:0.001:time1; 

% накладываем ограничения на y'm и обнуляем старшую степень y_1(t)
z(4) = 2*(z(2)-y(1))/time1;

% находим коэффициенты для многочленов Px и Py 2-ых степеней
coeffPolynomX1 = (z(3)-x(2))/(2*time1);
coeffPolynomY1 = (z(4)-y(2))/(2*time1);

% вид многочленов 
syms t
polynomX1 =  x(1)+x(2)*t+t.^2*coeffPolynomX1;
polynomY1 =  y(1)+y(2)*t+t.^2*coeffPolynomY1;

polynomVX1 =  x(2)+2*t*coeffPolynomX1;
polynomVY1 =  y(2)+2*t*coeffPolynomY1;

% модуль скорости
polynomVAbs1 = sqrt(polynomVX1.^2+polynomVY1.^2);

%построение многочлена, который является траекторий системы из (y0 y'0) в (0 0)
polynomX1M = matlabFunction(polynomX1);
polynomY1M = matlabFunction(polynomY1);
polynomVX1M = matlabFunction(polynomVX1);
polynomVY1M = matlabFunction(polynomVY1);
polynomVAbs1M = matlabFunction(polynomVAbs1);

dotsPolynomX1  = polynomX1M(time_int1);
dotsPolynomY1  = polynomY1M(time_int1);
dotsPolynomVX1  = polynomVX1M(time_int1);
dotsPolynomVY1  = polynomVY1M(time_int1);
dotsPolynomVabs1  = polynomVAbs1M(time_int1);

%% вторая терминальная задача
% выбираем время и обнуляем старшую степень x2(t)
time2 = 2*(x(3)-z(1))/(x(4)+z(3));
time_int2 = 0:0.001:time2; 
time_int22 = 2:0.001:2+time2; 

% накладываем ограничения на x'm и обнуляем старшую степень y_2(t)
z(3)=z(4)*(x(3)-z(1))/(y(3)-z(2));

% находим коэффициенты для многочленов Px и Py 2-ых степеней
coeffPolynomX2 = (x(4)-z(3))/(2*time2);
coeffPolynomY2 = (y(4)-z(4))/(2*time2);

% вид многочленов
syms t
polynomX2 = z(1)+z(3)*t+t.^2*coeffPolynomX2;
polynomY2 = z(2)+z(4)*t+t.^2*coeffPolynomY2;

polynomVX2 = z(3)+2*t*coeffPolynomX2;
polynomVY2 = z(4)+2*t*coeffPolynomY2;

% модуль скорости
polynomVAbs2 = sqrt(polynomVX2.^2+polynomVY2.^2);

%построение многочлена, который является траекторий системы из (y0 y'0) в (0 0)
polynomX2M = matlabFunction(polynomX2);
polynomY2M = matlabFunction(polynomY2);
polynomVX2M = matlabFunction(polynomVX2);
polynomVY2M = matlabFunction(polynomVY2);
polynomVAbs2M = matlabFunction(polynomVAbs2);

dotsPolynomX2  = polynomX2M(time_int2);
dotsPolynomY2  = polynomY2M(time_int2);
dotsPolynomVX2  = polynomVX2M(time_int2);
dotsPolynomVY2  = polynomVY2M(time_int2);
dotsPolynomVabs2  = polynomVAbs2M(time_int2);

%% визуализация результатов
figure
hold on
grid on
legend('y(x)','Location','northwest')
plot(dotsPolynomX1,dotsPolynomY1,'Blue', 'LineWidth', 1.5, 'DisplayName','y(x)');
plot(dotsPolynomX2,dotsPolynomY2,'Blue', 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel('x') 
ylabel('y') 
hold off

figure
hold on
grid on
legend('x(t)','Location','northwest')
plot(time_int1,dotsPolynomX1,'Blue', 'LineWidth', 1.5, 'DisplayName','x(t)');
plot(time_int22,dotsPolynomX2,'Blue', 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel('t') 
ylabel('x') 
hold off

figure
hold on
grid on
legend('y(t)','Location','northwest')
plot(time_int1,dotsPolynomY1,'Blue', 'LineWidth', 1.5, 'DisplayName','y(t)');
plot(time_int22,dotsPolynomY2,'Blue', 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel('t') 
ylabel('y') 
hold off

figure
hold on
grid on
legend('v(t)','Location','northwest')
plot(time_int1,dotsPolynomVabs1,'Blue', 'LineWidth', 1.5, 'DisplayName','v(t)');
plot(time_int22,dotsPolynomVabs2,'Blue', 'LineWidth', 1.5, 'HandleVisibility','off');
xlabel('t') 
ylabel('v') 
hold off

