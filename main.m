clear;
clc;
close all
set(0,'DefaultLineLineWidth', 1.5);
set(0,'defaultAxesFontSize', 20)
set(0,'DefaultFigureWindowStyle', 'docked') 
% set(0,'DefaultFigureWindowStyle', 'normal')
set(0,'defaulttextInterpreter','latex')
rng('default');

%%  1. Inizializzazione parametri


T_sample = 40;%[s] Tempo di campionamento
A_sez = [28 , 32 , 28 , 32]; %[cm^2] area sezione
a_foro = [0.071, 0.057, 0.071, 0.057]; %[cm^2] area foro
gg = 981; %[cm/s^2] gravità
k_1 = 2.7; %[cm^3/sV] Flusso acqua generato da pompa 1
k_2 = 3.2; %[cm^3/sV] Flusso acqua generato da pompa 2
gamma_1 = 0.3; %definisce come il flusso generato dalla pompa 1 viene diviso.
gamma_2 = 0.4; %definisce come il flusso generato dalla pompa 2 viene diviso.

% Test con ingressi
u = @(t) [3.75;3]; % Ingressi costanti (calcolati a mano ponendo a 0 x_dot)

% ODE del sistema
dxdt = @(t,x) quad_tank_function(t,x,u(t),A_sez,a_foro,gg,k_1,k_2,gamma_1,gamma_2);

%%  2.Simulazione del sistema

% Valori iniziali dei livelli dei serbatoi
x0 = [1.3767 2.2772 0.8386 0.5604]; 

% Tempo di simulazione
tspan = [0 1000]; % Intervallo di tempo per la simulazione

% Simulazione con ode45
[t, x] = ode45(dxdt, tspan, x0);

% Visualizzazione dei risultati
figure(1)
plot(t, x);
title('Simulazione del sistema quad tank con ingressi costanti');
xlabel('Tempo [s]');
ylabel('Livelli dei serbatoi [cm]');
legend('x1', 'x2', 'x3', 'x4');
grid on;

%% 3.Linearizzazione del sistema nell'equilibrio

% Punto di equilibrio
x_e = [7.8253; 18.7323; 3.3545; 7.8801]; % Punto di equilibrio per x (dato dal testo)
u_e = [3.75; 3]; % Punto di equilibrio per u

% Sistema di equazioni
syms x1 x2 x3 x4 u1 u2
x = [x1; x2; x3; x4];
u = [u1; u2];

% Spazio degli stati
x_dot = [-a_foro(1)/A_sez(1)*sqrt(2*gg*x1) + a_foro(3)/A_sez(1)*sqrt(2*gg*x3) + gamma_1*k_1*u1/A_sez(1);
         -a_foro(2)/A_sez(2)*sqrt(2*gg*x2) + a_foro(4)/A_sez(2)*sqrt(2*gg*x4) + gamma_2*k_2*u2/A_sez(2);
         -a_foro(3)/A_sez(3)*sqrt(2*gg*x3) + (1-gamma_2)*k_2*u2/A_sez(3);
         -a_foro(4)/A_sez(4)*sqrt(2*gg*x4) + (1-gamma_1)*k_1*u1/A_sez(4)];

% Derivate parziali 
A_sym = jacobian(x_dot, x);
B_sym = jacobian(x_dot, u);

% Derivate parziali valutate nel punto di equilibrio
A = double(subs(A_sym, [x1, x2, x3, x4, u1, u2], [x_e', u_e']));
B = double(subs(B_sym, [x1, x2, x3, x4, u1, u2], [x_e', u_e']));

% Definizione delle matrici C e D
C = eye(4);  %
D = zeros(4, 2);

% Sistema linearizzato
continuous_time_ss = ss(A, B, C, D);

% Discretizzazione del sistema (Eulero in avanti)
discrete_time_ss = c2d(continuous_time_ss, T_sample);
A_disc = discrete_time_ss.A;
B_disc = discrete_time_ss.B;

%% Stabilità, poli e zeri (sistema continuo e discretizzato)

% Autovalori (sistema continuo)
lambda0 = eig(A);
figure(2);
pzmap(continuous_time_ss);
title('Mappa Poli-Zeri (Sistema tempo continuo)');

% Matrice di raggiungibilità (sistema discreto)
Mr = ctrb(discrete_time_ss);
Mr_rank = rank(Mr);

% Autovalori (Sistema discreto)
lambda1 = eig(A_disc);
figure(3);
pzmap(discrete_time_ss);
title('Mappa Poli-Zeri (Sistema tempo discreto)');

%% Definizione variabili per il CIS (Control Invariant Set)

% Definizione delle matrici Q ed R
Q = 1*eye(size(A_disc,1));
[n,m] = size(B_disc);
R = 1*eye(m); 

% Soluzione eq. di Riccati
[K,P] = dlqr(A_disc,B_disc,Q,R); 

% Definizioni matrici dei vincoli sugli stati
Hx = [eye(size(A_disc));-eye(size(A_disc))];
Xmax = [20;20;20;20];
Xmin = [0.5;0.5;0.5;0.5];
hx = [Xmax - x_e ; -Xmin + x_e];

% Definizioni matrici dei vincoli sugli ingressi
Hu = [eye(m);-eye(m)];
Umax = 4.5 * ones(m, 1);
Umin = zeros(m, 1);
hu = [Umax - u_e ; -Umin + u_e];

%% Computazione CIS

[G,g] = cis(A_disc,B_disc,[0;0;0;0],[0;0],Hx,hx,Hu,hu,Q,R);

% Rappresentazione del CIS
CIS = Polyhedron(G,g); 
CIS_projected_12 = projection(CIS , 1:2);
CIS_projected_34 = projection(CIS , 3:4);

figure(4);
CIS_projected_12.plot();
title('CIS del sistema linearizzato proiettato ($x_1$ - $x_2$)');
xlabel('$x_1$');
ylabel('$x_2$');
xlim([-20,20]);
ylim([-20,20]);
grid on;

figure(5);
CIS_projected_34.plot();
title('CIS del sistema linearizzato proiettato ($x_3$ - $x_4$)');
xlabel('$x_3$');
ylabel('$x_4$');
xlim([-20,20]);
ylim([-20,20]);
grid on;

%% N-Steps Controllable Set (o dominio di attrazione)

% Stato iniziale traslato rispetto all'equilibrio
x0_traslato = x0 - x_e;

% Specifica delle proiezioni degli stati
dim1 = [1, 2]; % Per x1 e x2
dim2 = [3, 4]; % Per x3 e x4

% Computazione N-Steps
[Np_steps_H, Np_steps_h, Np_steps_H_proj1, Np_steps_h_proj1, Np_steps_H_proj2, Np_steps_h_proj2, Np] = controllable_set(Hx, hx, Hu, hu, G, g, A_disc, B_disc, dim1, dim2, x0_traslato);
Np_step_set_proj1 = Polyhedron('A', Np_steps_H_proj1, 'b', Np_steps_h_proj1);
Np_step_set_proj2 = Polyhedron('A', Np_steps_H_proj2, 'b', Np_steps_h_proj2);

% Visualizzazion N-Steps Controllable Set + CIS
figure(6)
h1 = Np_step_set_proj1.plot('Alpha', 0);
set(h1, 'LineWidth', 2); 
title('N-step-controllable set: $x_1$ vs $x_2$');
xlabel('$x_1$');
ylabel('$x_2$');
hold on;
CIS_projected_12.plot();
xlim([-20, 20]);
ylim([-20, 20]);
legend({'N-step set', 'Control-invariant set'})

figure(7)
h2 = Np_step_set_proj2.plot('Alpha', 0);
set(h2, 'LineWidth', 2);
title('N-step-controllable set: $x_3$ vs $x_4$');
xlabel('$x_3$');
ylabel('$x_4$');
hold on;
CIS_projected_34.plot();
xlim([-20, 20]);
ylim([-20, 20]);
legend({'N-step set', 'Control-invariant set'})

% Se Np è basso e il vincolo di uguaglianza risulta "Infeasible" forzare Np
% ad un valore più alto rispetto a quanto ottenuto con la funzione
% "controllable_set"
Np = 30; 

%% Trasformazione dei vincoli 

% Vincoli X ed U in forma matriciale
Htilde_x = getMatrixVincoli(Hx, Np);
Htilde_u = getMatrixVincoli(Hu, Np-1);

htilde_x = getVecVincoli(hx,Np);
htilde_u = getVecVincoli(hu,Np-1);

%% Definizione della dinamica del sistema

[A_call,A_call_last] = getA_calligrafico(A_disc,Np);
[B_call, B_call_last] = getB_calligrafico(A_disc,B_disc,Np-1);

%% Trasformazione del costo

Q_tilde = getQ_tilde(Q,P,Np);
R_tilde = getR_tilde(R,Np);

% Vincolo di disuguaglianza
A_dis = [Htilde_u; 
        Htilde_x * B_call;
        G*B_call_last];

% B_dis è nel for della simulazione dell'MPC in quanto è necessario 
% utilizzare x(0) che ad ogni iterazione varia

% Definizione H,f
H = 2 * (B_call' * Q_tilde * B_call + R_tilde);
f_base = 2 * A_call' * Q_tilde * B_call;


%% Simulazione Vincolo terminale disuguaglianza

% Tempo di simulazione
T_sim = 100;

% Vettori di stato ed ingresso per salvataggio finale
x = zeros(n,T_sim+1);
u = zeros(m,T_sim);

% Stato iniziale
x(:,1) = x0';

% Simmetria forzata (evita il Warning)
H = (H+H')/2; 

t_tot = [];
x_tot = [];


for k = 1:T_sim

    % Aggiornamento stato x0
    x_current = x(:,k) - x_e;

    % Aggiornamento di f
    f = x_current' * f_base;

    % Vincolo di disuguaglianza
    b_dis = [htilde_u; 
            htilde_x - Htilde_x * A_call * x_current; 
            g - G*A_call_last * x_current];

    options = optimset('Display', 'off');
    [U, fval, exitflag, output] = quadprog(H,f,A_dis,b_dis,[],[],[],[],[],options);

    % receding horizon + equilibrio
    u(:,k) = U(1:m) + u_e; 

    % Simulazione del sistema con u0 trovata per il tempo di campionamento
    dxdt = @(t,x) controlledTanks(t,x,u(:,k),A_sez,a_foro,gg,k_1,k_2,gamma_1,gamma_2);
    [t, xx] = ode45(dxdt, [0 T_sample], x(:,k));

    % Memorizza dell'ultimo stato della simulazione
    x(:,k+1) = xx(end,:)';
    
    % Memorizzazione dei risultati ottenuti
    t_tot = [t_tot; t + (k-1)*T_sample];
    x_tot = [x_tot; xx];
end

% Visualizzazione degli stati 
figure;
plot(t_tot, x_tot(:, 1),'DisplayName', 'x1');
hold on;
plot(t_tot, x_tot(:, 2),'DisplayName', 'x2');
plot(t_tot, x_tot(:, 3),'DisplayName', 'x3');
plot(t_tot, x_tot(:, 4),'DisplayName', 'x4');
xlabel('Tempo (s)');
ylabel('Stati');
yline(x_e(1))
yline(x_e(2))
yline(x_e(3))
yline(x_e(4))
legend(["$X_1$" , "$X_2$" , "$X_3$" , "$X_4$","Obiettivo"], Interpreter="latex")
grid on;
title('Andamento degli stati nel tempo');

% Visualizzazione delle azioni di controllo
figure;
stairs(u(1, :), 'LineWidth', 1.5); 
hold on; 
stairs(u(2, :), 'LineWidth', 1.5); 
hold off; 
yline(Umax(1))
title('Andamento delle azioni di controllo');
xlabel('Tempo [s]');
ylabel('Controllo');
legend(["$u_1$", "$u_1$","Umax"],Interpreter="latex"); 


%% Simulazione Vincolo terminale uguaglianza

x = zeros(n,T_sim+1);
u = zeros(m,T_sim);

x(:,1) = x0';

t_tot = [];
x_tot = [];

% Vincolo su stati e controllo senza vincolo x(N) appartenente a CIS
A_dis = [Htilde_u; 
        Htilde_x * B_call];

% Vincolo di uguaglianza per x(N)
G = [eye(4);
    -eye(4)];
A_eq = G * B_call_last;

for k = 1:T_sim

    % Aggiornamento stato x0
    x_current = x(:,k) - x_e;

    % Aggiornamento di f
    f = x_current' * f_base;

    % Vincolo su stati e controllo
    b_dis = [htilde_u; 
            htilde_x - Htilde_x * A_call * x_current];

    % Vincolo di uguaglianza per x(N)
    b_eq = [zeros(8,1)] - G * A_call_last * x_current;

    options = optimset('Display', 'off');
    [U, fval, exitflag, output] = quadprog(H, f, A_dis, b_dis, A_eq, b_eq,[],[],[],options);

    % receding horizon + equilibrio
    u(:,k) = U(1:m) + u_e; 

    % Simulazione del sistema con u0 trovata per il tempo di campionamento
    dxdt = @(t,x) controlledTanks(t,x,u(:,k),A_sez,a_foro,gg,k_1,k_2,gamma_1,gamma_2);
    [t, xx] = ode45(dxdt, [0 T_sample], x(:,k));

    % Memorizza dell'ultimo stato della simulazione
    x(:,k+1) = xx(end,:)';
    
    % Memorizzazione dei risultati ottenuti
    t_tot = [t_tot; t + (k-1)*T_sample];
    x_tot = [x_tot; xx];

end

% Visualizzazione degli stati 
figure;
plot(t_tot, x_tot(:, 1),'DisplayName', 'x1');
hold on;
plot(t_tot, x_tot(:, 2),'DisplayName', 'x2');
plot(t_tot, x_tot(:, 3),'DisplayName', 'x3');
plot(t_tot, x_tot(:, 4),'DisplayName', 'x4');
xlabel('Tempo (s)');
ylabel('Stati');
yline(x_e(1))
yline(x_e(2))
yline(x_e(3))
yline(x_e(4))
legend(["$X_1$" , "$X_2$" , "$X_3$" , "$X_4$","Obiettivo"], Interpreter="latex")
grid on;
title('Andamento degli stati nel tempo');

% Visualizzazione delle azioni di controllo
figure;
stairs(u(1, :), 'LineWidth', 1.5); 
hold on; 
stairs(u(2, :), 'LineWidth', 1.5); 
hold off; 
yline(Umax(1))
title('Andamento delle azioni di controllo');
xlabel('Tempo [s]');
ylabel('Controllo');
legend(["$u_1$", "$u_1$","Umax"],Interpreter="latex"); 