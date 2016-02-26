function [F, G, ind] = OraclePG(qc, ind)
    
    // qc : vecteur des debits rduits
    // F : valeur du critère evalue au point qc
    // G : vecteur des derivees du critere par rapport à qc
    // ind : indicateur du type de calcul à effectuer
    //      = 2 : calcul de F uniquement
    //      = 3 : calcul de G uniquement
    //      = 4 : calcul de F et G
    
    // n : nombre total d'arcs, N
    // m : nombre total de noeuds, N
    // md : nombre de noeuds de demande, N
    // mr : nombre de noeuds reservoir, N
    // fd : demandes aux noeuds de demande, M(md,1)
    // pr : pressions aux noeuds reservoir, M(mr,1)
    // r : resistances des arcs, M(n,1)
    // q0 : vecteur initial des debits, M(n,1)
    // A : matrice d'incidence noeuds-arcs, M(m,n)
    // Ad : sous-matrice "demande" de A, M(md,n)
    // Ar : sous-matrice "reservoir" de A, M(mr,n)
    // AdT : sous-matrice "arbre" de Ad, M(md,md)
    // AdC : sous--matrice "coarbre" de Ad, M(md,n-md)
    // AdI : matrice inverse de AdT, M(md,md)
    // B : matrice d'incidence arcs-cycles, M(n,n-md)
    // q : debits des arcs, M(n,1)
    // qc : debits reduits des arcs, M(n-md,1)
    
    F = 0;
    G = zeros(n-md,1);
    ind = ind;
    
    A = ones(m,n); // donne
    p = ones(m,1); // donne
    
    AdT = A(m-md+1:m, 1:md);
    // AdT = A(mr+1:m, 1:md); // autre definition
    AdC = A(m-md+1:m, md+1:n);
    AdI = inv(AdT);
    Ar = A(1:mr,:);
    B = cat(1, -AdI*AdC, eye(n-md,n-md));
    q0 = cat(1, AdI*fd, zeros(n-md,1));
    pr = p(1:m-md,1);
    // pr = p(1:mr,1); // autre definition
    
    if ind == 2 then
        F = ((q0+B*qc)*(r.*(q0+B*qc).*(abs(q0+B*qc)).)/3 + pr*(Ar*(q0+B*qc))';
    elseif ind == 3 then
        G = B'*Ar'*pr + B'*(r.*(q0+B*qc).*(abs(q0+B*qc)).);
    elseif ind == 4 then
        F = ((q0+B*qc)*(r.*(q0+B*qc).*(abs(q0+B*qc)).)/3 + pr*(Ar*(q0+B*qc))';
        G = B'*Ar'*pr + B'*(r.*(q0+B*qc).*(abs(q0+B*qc)).);
    end
    
endfunction
