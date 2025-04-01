function [H_nsteps, h_nsteps, H_proj1, h_proj1, H_proj2, h_proj2, Np] = controllable_set(Hx, hx, Hu, hu, H_target, h_target, A, B, dims1, dims2, x0_traslato)
%CONTROLLABLE_SET Calcolo dell'N-step-controllable set ad un set target
%descritto da H_target * x <= h_target

n = size(A,2);
m = size(B,2);

%   Candidato iniziale
H_ii_steps= H_target; 
h_ii_steps = h_target;

for ii=1:9999
    %   Computazione in R^(n+m) di stati e ingressi assieme
temp = Polyhedron('A',[H_ii_steps*A,H_ii_steps*B; ...
    zeros(size(Hu,1),n),Hu],'b',[h_ii_steps;hu]);
    
    %   Proiezione in R^n
    temp = projection(temp,1:n); %proiezione di temp da 1 ad n (stati)
    
    %Per ottenere una rappresentazione "minimizzata" che considera solo
    %vincoli non ridondanti
    temp.minHRep();

    if temp.contains(x0_traslato)
            % se contiene il punto usciamo dal ciclo
            break
    end


    %   Intersezione con X := {x | Hx*x <= hx} (basta aggiungere righe di vincoli)
    H_ii_steps = [temp.A; Hx];
    h_ii_steps = [temp.b; hx];
    
end

H_nsteps = H_ii_steps;
h_nsteps = h_ii_steps;

Np = ii;

% Proiettare il set controllabile sui due sottoinsiemi di dimensioni specificati per il plotting
[H_proj1, h_proj1] = projection_dim(H_nsteps, h_nsteps, dims1);
[H_proj2, h_proj2] = projection_dim(H_nsteps, h_nsteps, dims2);
end

function [H_proj, h_proj] = projection_dim(H, h, dims)
%Proiezione del set su un sottoinsieme di dimensioni specificate
%   H: Matrice dei vincoli del set originale
%   h: Vettore dei vincoli del set originale
%   dims: Indici delle dimensioni su cui proiettare

% Seleziona le colonne corrispondenti alle dimensioni specificate
H_proj = H(:, dims);
h_proj = h;
end
