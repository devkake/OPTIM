function [F, G, H, ind] = OracleDH(lambda, ind)

// ---------------------------
// Arguments and Return values
// ---------------------------
    
// lambda : vecteur M(md, 1)
// F : valeur du critere evalue au point qc
// G : vecteur des derivees du critere par rapport a qc
// H : matrice des derivees secondes du critere par rapport a qc
// ind : indicateur du type de calcul à effectuer
//      = 2 : calcul de F uniquement
//      = 3 : calcul de G uniquement
//      = 4 : calcul de F et G
//      = 5 : calcul de H uniquement
//      = 6 : calcul de G et H
//      = 7 : calcul de F, G et H

// ---------------------------
// Variables
// ---------------------------

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

// ---------------------------
// Initialisation
// ---------------------------
    
    F = 0;
    G = zeros(n-md,1);
    H = zeros(n-md,n-md);
    ind = ind;
    
// A est donne
// pr est donne
// p n'est pas donne
// AdT = A(m-md+1:m, 1:md);
// AdT = A(mr+1:m, 1:md); // autre definition
// AdC = A(m-md+1:m, md+1:n);
// AdI = inv(AdT); // prob singularity
// Ar = A(1:mr,:);
// B = cat(1, -AdI*AdC, eye(n-md,n-md));
// q0 = cat(1, AdI*fd, zeros(n-md,1));

// ---------------------------
// Calculs de F et G
// ---------------------------

    z = -Ar'*pr-Ad'*lambda;
    qdiese = (z./abs(z)).*sqrt(abs(z)./r);
    if ind == 2 then
        // cas : seulement F
        F = -1 * (1./3.*qdiese'*(r.*qdiese.*abs(qdiese)) + pr'*(Ar*qdiese) + lambda'*(Ad*qdiese-fd));
    elseif ind == 3 then
        // cas : seulement G
        G = -1 * (Ad*qdiese - fd);
    elseif ind == 4 then
        // cas : F et G
        F = -1 * (1./3.*qdiese'*(r.*qdiese.*abs(qdiese)) + pr'*(Ar*qdiese) + lambda'*(Ad*qdiese-fd));
        G = -1 * (Ad*qdiese - fd);
    elseif ind == 5 then
        // cas : seulement H
        H = -1 * ((-1./2.)*Ad*diag((1. ./ sqrt(r.*abs(z))))*Ad');
    elseif ind == 6 then
        // cas : G et H
        G = -1 * (Ad*qdiese - fd);

        H = -1 * ((-1./2.)*Ad*diag((1. ./ sqrt(r.*abs(z))))*Ad');
    elseif ind == 7 then
        // cas : F, G et H
        F = -1 * (1./3.*qdiese'*(r.*qdiese.*abs(qdiese)) + pr'*(Ar*qdiese) + lambda'*(Ad*qdiese-fd));
        G = -1 * (Ad*qdiese - fd);
        H = -1 * ((-1./2.)*Ad*diag((1. ./ sqrt(r.*abs(z))))*Ad');
    end
    
endfunction
