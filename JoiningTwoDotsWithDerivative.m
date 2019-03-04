function JoiningTwoDotsWithDerivative;
clc
clear all

% начальное и конечное положения системы (1,6)-->(3,8)
x = [1,2,3,4]; % x = [x0, x'0, xt, x't]
y = [6,7,8,9]; % y = [y0, y'0, yt, y't] 

% высчитанное время
time = 2*(x(3)-x(1))/(x(4)+x(2));
time_int = 0:0.001:time; 

% накладываем ограничения на y't (можно ли на другие?)
y(4) = 2*(y(3)-y(1))/time-y(2);

% находим коэффициенты для многочленов Px и Py 2-ых степеней
coeffPolynomX = (x(4)-x(2))/(2*time);
coeffPolynomY = (y(4)-y(2))/(2*time);

% вид многочленов
syms t
polynomX = @(t) x(1)+x(2)*t+t.^2*coeffPolynomX;
polynomY = @(t) y(1)+y(2)*t+t.^2*coeffPolynomY;

%построение многочлена, который является траекторий системы из (y0 y'0) в (0 0)
dotsPolynomX  = polynomX(time_int);
dotsPolynomY  = polynomY(time_int);

figure
hold on
plot(dotsPolynomX,dotsPolynomY,'Red', 'LineWidth', 1.5);
xlabel('x') 
ylabel('y') 
hold off

end