function Q_tilde = getQ_tilde(Q,P,N)

% Inizializza H_tilde come una matrice vuota
    Q_tilde = [];
    
%     % Aggiungi H alla diagonale N volte
%     for i = 1:N
%         Q_tilde = blkdiag(Q_tilde, Q);
%     end
%     Q_tilde = blkdiag(Q_tilde, P);

    % Crea una cella array con N-1 copie di Q seguite da P
    blocks = [repmat({Q}, 1, N), {P}];
    
    % Usa blkdiag per creare la matrice a blocchi diagonali
    Q_tilde = blkdiag(blocks{:});

end