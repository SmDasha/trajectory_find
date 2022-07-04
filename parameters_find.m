function tX=parameters_find(t,X,Sigma,g_r,g_z,ro,P_spec,c_xa,c_ya,b,alph,gamma_a,beth)
    S = sqrt(3)*b^2/4; % характерная площадь ЛА, м^2
    % Коэффициенты sigma
    sigma_x = c_xa*S/(2*X(7));
    sigma_y = c_ya*S/(2*X(7));
    % Определение проекций вектора тяги
    P = P_spec*beth;
    P_x = P*cosd(alph);
    P_y = P*sind(alph)*cosd(gamma_a);
    P_z = P*sind(alph)*sind(gamma_a);
    % Проекции перегрузок
    % g0 = 9.81; % м/с^2
    % n_x = P./(g0.*X(7))+S*ro*X(1).^2.*(c_ya.*sin(alph)-c_xa.*cos(alph))./(g0*X(7)*2);
    % n_y = S*ro*X(1).^2.*(c_ya.*cos(alph)+c_xa.*sin(alph))./(2*g0*X(7));
    
    dV = -sigma_x*ro*X(1)^2-g_r*sind(X(2))+g_z*sind(X(3))*cosd(X(2))+P_x/X(7)+...
        X(4)*Sigma^2*cosd(X(5))*(sind(X(2))*cosd(X(5))-cosd(X(2))*sind(X(5))*sind(X(3)));
    dtheta = sigma_y*ro*X(1)*cosd(gamma_a)+(X(1)/X(4)-g_r/X(1))*cosd(X(2))-...
        (g_z/X(1))*sind(X(3))*sind(X(2))+P_y/(X(1)*X(7))+2*Sigma*cosd(X(5))*cosd(X(3))...
        +(X(4)*Sigma.^2/X(1))*cosd(X(5))*(cosd(X(2))*cosd(X(5))+sind(X(2))*sind(X(5))*sind(X(3)));
    dhi = -sigma_y*ro*X(1)*sind(gamma_a)/cosd(X(2))-X(1)*cosd(X(2))*tand(X(5))*cosd(X(3))/X(4)...
        +g_z*cosd(X(3))/(X(1)*cosd(X(2)))-P_z/(X(7)*X(1)*cosd(X(2)))...
        -2*Sigma*(sind(X(5))-cosd(X(5))*sind(X(3))*tand(X(2)))...
        -X(4)*Sigma.^2*sind(X(5))*cosd(X(5))*cosd(X(3))/(X(1)*cosd(X(2)));
    dR = X(1)*sind(X(2));
    dphi = X(1)*cosd(X(2))*sind(X(3))/X(4);
    dlambda = X(1)*cosd(X(2))*cosd(X(3))/(X(4)*cosd(X(5)));
    dm = -beth; 
    tX = [dV; dtheta; dhi; dR; dphi; dlambda; dm];    
end
