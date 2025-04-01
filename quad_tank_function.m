function x_dot = quad_tank_function(t,x,u,A_sez,a_for,g,k_1,k_2,gamma_1,gamma_2)
% Modello nello spazio degli stati quad tank

%   Equazioni di stato
x_dot = zeros(4,1);
x_dot(1) = -a_for(1)/A_sez(1)*sqrt(2*g*x(1)) + a_for(3)/A_sez(1)*sqrt(2*g*x(3)) + gamma_1*k_1*u(1)/A_sez(1);
x_dot(2) = -a_for(2)/A_sez(2)*sqrt(2*g*x(2)) + a_for(4)/A_sez(2)*sqrt(2*g*x(4)) + gamma_2*k_2*u(2)/A_sez(2);
x_dot(3) = -a_for(3)/A_sez(3)*sqrt(2*g*x(3)) + (1-gamma_2)*k_2*u(2)/A_sez(3);
x_dot(4) = -a_for(4)/A_sez(4)*sqrt(2*g*x(4)) + (1-gamma_1)*k_1*u(1)/A_sez(4);

end

