function dx = controlledTanks(t,x,u,A_sez,a_for,g,k_1,k_2,gamma_1,gamma_2)

% Stati
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

% Equazioni non lineari
dx1 = -a_for(1)/A_sez(1)*sqrt(2*g*x1) + a_for(3)/A_sez(1)*sqrt(2*g*x3) + gamma_1*k_1*u(1)/A_sez(1);
dx2 = -a_for(2)/A_sez(2)*sqrt(2*g*x2) + a_for(4)/A_sez(2)*sqrt(2*g*x4) + gamma_2*k_2*u(2)/A_sez(2);
dx3 = -a_for(3)/A_sez(3)*sqrt(2*g*x3) + (1-gamma_2)*k_2*u(2)/A_sez(3);
dx4 = -a_for(4)/A_sez(4)*sqrt(2*g*x4) + (1-gamma_1)*k_1*u(1)/A_sez(4);

dx = [dx1; dx2;dx3;dx4];

end