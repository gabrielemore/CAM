function R_tilde = getR_tilde(R,N)

% Inizializza H_tilde come una matrice vuota
    R_tilde = [];
    
    % Aggiungi H alla diagonale N volte
    for i = 1:N
        R_tilde = blkdiag(R_tilde, R);
    end

end