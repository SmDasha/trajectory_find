%% Построение  траектории ГЗЛА с последовательным изменением alph
% aircraft_trajectory - выполняет вызов скрипта start_conditions и 
% функции parameters_find, выполняет построение траекторий в соответствии 
% с заданными alpha;
% start_conditions - скрипт, в котором задаются начальные условия;
% parameters_find - функция, описывающая систему ОДУ.
% В параметрах не учтено изменение плотности атмосферы, задана постоянная
% масса, задан постоянный удельный вектор тяги
clc
close all
clearvars
start_conditions_next;
num = length(alph);
str_alpha = '';
str = {'Скорость, м/с';'Угол наклона траектории, град'; 'Угол пути, град'; 'Радиус-вектор, м'};
% Решение системы дифференциальных уравнений
for i = 1:num
    [t,X] = ode45(@(t,X) parameters_find(t,X,Sigma,g_r,g_z,ro,P_spec,c_xa(i),c_ya(i),b,alph(i),gamma_a,beth), time_interval_for_plot, start_cond);
    R = X(:,4).';
    phi = X(:,5).';
    lmb = X(:,6).';
    H = R-R_r./sqrt(1-0.0066934*(cosd(phi)).^2);
    % Построение графиков
    figure(1)
    plot(time_interval_for_plot,H,time_interval_for_plot(1),H(1),'*'); grid on; hold on;
    %text(time_interval_for_plot(1),H(1)-1,['\leftarrow alpha= ', num2str(alph(i))]);
    title('Траектория ГЗЛА')
    xlabel('Время, с'); ylabel('Высота, м');
    % Перевод в ПГСК
    wgs84 = wgs84Ellipsoid('meters');
    [x,y,z] = geodetic2ecef(wgs84,phi,lmb,H);
    if i == 1
        next_conditions;
    end

%     figure(2)
% 	for k = 1:4
%         subplot(2,2,k)
%         plot(time_interval_for_plot,X(:,k)); grid on;hold on; ylabel(str(k)); xlabel('Время, с');
% 	end
end

% str = {'Скорость, м/с';'Угол наклона траектории, град'; 'Угол пути, град'; 'Радиус-вектор, м'};
% figure
% for i = 1:4
%     subplot(2,2,i)
%     plot(X(:,i)); grid on; ylabel(str(i)); xlabel('Время, с');
% end