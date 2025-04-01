function [G,g] = cis(A,B,x_ref,u_ref,Fx,fx,Fu,fu,Q,R)
%CIS Calcolo del control invariant set (CIS) di un sistema lineare
%   Questo metodo assume che un controllore LQR venga utilizzato
%   all'interno del CIS
%   Input:
%       - A, B: matrici del sistema
%       - x_ref: equilibrio attorno al quale costruire il CIS
%       - Fx*x<=fx: vincoli sullo stato
%       - Fu*u<=fu: vincoli sull'ingresso
%       - Q,R: matrici per LQR

%   Controllore LQR
% LQR è un oggetto astratto per dimostrare che in un apiccola regione
% esiste qualcosa che mi stabilizza il sistema. Lo calcolo SOLO per trovare
% i vincoli, poi uso solo MPC
K = -dlqr(A,B,Q,R); 

%   Matrice A del sistema controllato con LQR
A_lqr = A+B*K;

%   Vincoli sullo stato e sull'ingresso (F*x <= f)
F = [Fx; Fu*K];
f = [fx; fu + Fu*(K*x_ref - u_ref)];

%   Calcolo del CIS (G*x<=g)
%CIS_poly_prev = Polyhedron(); %inizialmente vuoto
CIS_poly_curr = Polyhedron(F,f);

primaIterazione= 1;
counter = 0;
% se nel progetto dà problemi, mettiamo un numero finito di iterazioni
while primaIterazione || CIS_poly_prev ~= CIS_poly_curr
    
    primaIterazione = 0;
    %   Memorizza vecchio candidato
    CIS_poly_prev = CIS_poly_curr;

    %   Calcola nuovo candidato (G_hat * x <= g_hat)
    G_hat = [CIS_poly_curr.A * A_lqr; F]; %metto i vincoli iniziali F in quanto il nuovo poly potrebbe uscire
    g_hat = [CIS_poly_curr.b + CIS_poly_curr.A*B*(K*x_ref - u_ref);f];

    CIS_poly_curr = Polyhedron(G_hat,g_hat);
    % potenzialmente il risultato finale potrebbe avere matrici o vettori
    % con moltissime righe e alcune potrebbero essere ridondanti. usando
    % minHRep e in oggetto il poliedro, si ottiene un risultato identico
    % geometricamente ma una rappresentazione con il minor numero di righe
    % possibili.
    CIS_poly_curr = CIS_poly_curr.minHRep();
    
    counter = counter+1;

end


%   Disequazioni che descrivono il CIS
G = CIS_poly_curr.A;
g = CIS_poly_curr.b;