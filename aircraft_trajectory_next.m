%% ����������  ���������� ���� � ���������������� ���������� alph
% aircraft_trajectory - ��������� ����� ������� start_conditions � 
% ������� parameters_find, ��������� ���������� ���������� � ������������ 
% � ��������� alpha;
% start_conditions - ������, � ������� �������� ��������� �������;
% parameters_find - �������, ����������� ������� ���.
% � ���������� �� ������ ��������� ��������� ���������, ������ ����������
% �����, ����� ���������� �������� ������ ����
clc
close all
clearvars
start_conditions_next;
num = length(alph);
str_alpha = '';
str = {'��������, �/�';'���� ������� ����������, ����'; '���� ����, ����'; '������-������, �'};
% ������� ������� ���������������� ���������
for i = 1:num
    [t,X] = ode45(@(t,X) parameters_find(t,X,Sigma,g_r,g_z,ro,P_spec,c_xa(i),c_ya(i),b,alph(i),gamma_a,beth), time_interval_for_plot, start_cond);
    R = X(:,4).';
    phi = X(:,5).';
    lmb = X(:,6).';
    H = R-R_r./sqrt(1-0.0066934*(cosd(phi)).^2);
    % ���������� ��������
    figure(1)
    plot(time_interval_for_plot,H,time_interval_for_plot(1),H(1),'*'); grid on; hold on;
    %text(time_interval_for_plot(1),H(1)-1,['\leftarrow alpha= ', num2str(alph(i))]);
    title('���������� ����')
    xlabel('�����, �'); ylabel('������, �');
    % ������� � ����
    wgs84 = wgs84Ellipsoid('meters');
    [x,y,z] = geodetic2ecef(wgs84,phi,lmb,H);
    if i == 1
        next_conditions;
    end

%     figure(2)
% 	for k = 1:4
%         subplot(2,2,k)
%         plot(time_interval_for_plot,X(:,k)); grid on;hold on; ylabel(str(k)); xlabel('�����, �');
% 	end
end

% str = {'��������, �/�';'���� ������� ����������, ����'; '���� ����, ����'; '������-������, �'};
% figure
% for i = 1:4
%     subplot(2,2,i)
%     plot(X(:,i)); grid on; ylabel(str(i)); xlabel('�����, �');
% end