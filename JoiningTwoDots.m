function JoiningTwoDots;
clc
clear all

% ��������� ��������� �������
x0 = [1,1];
% �������� ��������� �������
xt = [2,2];

% �����
Tdeg3 = 2*(xt(1)-x0(1))/(xt(2)+x0(2))-0.5; % ��������
Tdeg2 = 2*(xt(1)-x0(1))/(xt(2)+x0(2)); % �����������

t_int_deg3 = 0:0.001:Tdeg3; % ��������� �������� ��� ���������� deg=3
t_int_deg2 = 0:0.001:Tdeg2; % ��������� �������� ��� ���������� deg=2

% ������� ������������ ��� ���������� 3-�� ������� ����� ����
LeftSidedeg3 = [Tdeg3^2, Tdeg3^3; 2*Tdeg3, 3*Tdeg3^2];
RightSidedeg3 = [xt(1)-x0(1)-x0(2)*Tdeg3; xt(2)-x0(2)];
coeffsdeg3 = LeftSidedeg3 \ RightSidedeg3;

% ������� ����������� ��� ���������� 2-�� ������� 
coeffdeg2 = (xt(2)-x0(2))/(2*Tdeg2);

syms t
%���������� ������ � ������� ��������
pdeg3 = @(t) x0(1)+x0(2)*t+t.^2*coeffsdeg3(1)+(t.^3)*coeffsdeg3(2);
pdeg2 = @(t) x0(1)+x0(2)*t+t.^2*coeffdeg2;

% ���������� x' ��� ���������� �������� ��������
pdeg3diffs = diff(pdeg3, t);
pdeg2diffs = diff(pdeg2, t);

% �������������� 
pdeg3diff = matlabFunction(pdeg3diffs);
pdeg2diff = matlabFunction(pdeg2diffs);

%���������� ����������, ������� �������� ���������� ������� �� (y0 y'0) � (0 0)
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
