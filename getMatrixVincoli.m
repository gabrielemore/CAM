function H_tilde = getMatrixVincoli(H, N)

% Inizializza H_tilde come una matrice vuota
    H_tilde = [];
    
    % Aggiungi H alla diagonale N volte
    for i = 1:N+1
        H_tilde = blkdiag(H_tilde, H);
    end

end