%% ������� ����������� ��������
% ����
del_alph = input('��������� ���� �����, ����: '); % ��������� ���� �����
del_gamma_a = input('��������� ���� �����, ����: '); % ��������� ���� �����
del_c_xa = input('��������� ������������ �������� �������������: '); % ��������� ������������ �������� �������������
del_c_ya = input('��������� ������������ ��������� ����: '); % ��������� ������������ ��������� ����
% ���������
Sigma = input('������� �������� ������������ �������� �����, �-1: '); % ������� �������� ������������ �������� �����, �-1
R_e = input('�������������� ������, �: '); % �������������� ������, �
R_r = input('�������� ������, �: '); % �������� ������, �
gamma_z = input('�������������� ���������� �����, ��^3/c^2: '); % �������������� ���������� �����, ��^3/c^2
a = input('�������� �����, �/�: '); % �������� �����, �/�
P_spec = input('�������� ������ ����, ��*�/�^2: '); % �������� ������ ����, ��*�/�^2
% ����������� �����������
%alph = alph+del_alph; % ���� �����, ����
gamma_a = gamma_a+del_gamma_a; % ���� �����, ����
% ����������� ������������� �� ���� ����� alph
%c_xa = c_xa + del_c_xa;
%c_ya = c_ya + del_c_ya;
% ����������
beth = input('�������� ������ �������, ��/�: '); % �������� ������ �������, ��/�
b = input('������� ������������, �: '); % ������� ������������, �
M = ones(1,length(alph));
M(:) = input('����� ����: ');
dV0 = a*M(i); % ������� ��������, �/�
dtheta0 = X(end,2); % ���� ������� ����������, ����
dhi0 = X(end,3); % ���� ����, ����
dR0 = X(end,4); % ������-������, ������ ������� - ����� �����, ����
dphi0 = X(end,5); % �������������� ������, ����
dlambda0 = X(end,6); % �������������� �������, ����
dm0 = X(end,7); % �����, ��
% time_interval_for_plot = (time_interval(end):0.1:time_interval(end)+10); % �
% time_interval_for_plot = (time_interval_for_plot(end):0.1:time_interval_for_plot(end)+10); % �
% time_interval_for_plot = (10:0.1:20); % �
% time_interval_for_plot = (0:0.1:10); % �
% time_interval_for_plot = (time_interval_for_plot(end)+0.1:0.1:time_interval_for_plot(end)+20); % �
start_cond = [dV0; dtheta0; dhi0; dR0; dphi0; dlambda0; dm0];
g_r = (gamma_z./dR0.^2)*(1+0.00162*(R_e./dR0).^2*(1-3*(sind(dphi0)).^2)); % ���������� ������������, �/�^2
g_z = -0.00162*gamma_z*R_e.^4.*sind(2*dphi0)./dR0.^4; % ��������������� ������������, �/�^2
