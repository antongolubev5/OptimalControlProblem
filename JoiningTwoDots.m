function JoiningTwoDots;
clc
clear all

% начальное положение системы
x0 = [1,1];
% конечное положение системы
xt = [2,2];

% время
Tdeg3 = 2*(xt(1)-x0(1))/(xt(2)+x0(2))-0.5; % заданное
Tdeg2 = 2*(xt(1)-x0(1))/(xt(2)+x0(2)); % высчитанное

t_int_deg3 = 0:0.001:Tdeg3; % временной интервал для многочлена deg=3
t_int_deg2 = 0:0.001:Tdeg2; % временной интервал для многочлена deg=2

% находим коэффициенты для многочлена 3-ей степени решив СЛАУ
LeftSidedeg3 = [Tdeg3^2, Tdeg3^3; 2*Tdeg3, 3*Tdeg3^2];
RightSidedeg3 = [xt(1)-x0(1)-x0(2)*Tdeg3; xt(2)-x0(2)];
coeffsdeg3 = LeftSidedeg3 \ RightSidedeg3;

% находим коэффициент для многочлена 2-ой степени 
coeffdeg2 = (xt(2)-x0(2))/(2*Tdeg2);

syms t
%многочлены второй и третьей степеней
pdeg3 = @(t) x0(1)+x0(2)*t+t.^2*coeffsdeg3(1)+(t.^3)*coeffsdeg3(2);
pdeg2 = @(t) x0(1)+x0(2)*t+t.^2*coeffdeg2;

% нахождение x' для построения фазового портрета
pdeg3diffs = diff(pdeg3, t);
pdeg2diffs = diff(pdeg2, t);

% преобразование 
pdeg3diff = matlabFunction(pdeg3diffs);
pdeg2diff = matlabFunction(pdeg2diffs);

%построение многочлена, который является траекторий системы из (y0 y'0) в (0 0)
xarray_deg2  = pdeg2(t_int_deg2);
yarray_deg2  = pdeg2diff(t_int_deg2);

xarray_deg3  = pdeg3(t_int_deg3);
yarray_deg3  = pdeg3diff(t_int_deg3);

figure
hold on
plot(xarray_deg2,yarray_deg2,'Red', 'LineWidth', 1.5);
plot(xarray_deg3,yarray_deg3, 'Blue', 'LineWidth', 1.5);
legend('deg2','deg3')
xlabel('x1') 
ylabel('x2') 
hold off

end
