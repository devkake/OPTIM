// ---------------------------
// OraclePH
// ---------------------------

function [F, G, ind] = OraclePG(qc, ind)

// ---------------------------
// Arguments and Return values
// ---------------------------
    
// qc : vecteur des debits reduits
// F : valeur du critere evalue au point qc
// G : vecteur des derivees du critere par rapport a qc
// ind : indicateur du type de calcul a effectuer
//      = 2 : calcul de F uniquement
//      = 3 : calcul de G uniquement
//      = 4 : calcul de F et G

// ---------------------------
// Variables
// ---------------------------

// ---------------------------
// Initialisation
// ---------------------------
    
    F = 0;
    G = zeros(n-md,1);
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
    
    if ind == 2 then
        // cas : seulement F
        F = ((q0+B*qc)'*(r.*(q0+B*qc).*(abs(q0+B*qc))))/3 + pr'*(Ar*(q0+B*qc));
    elseif ind == 3 then
        // cas : seulement G
        G = B'*Ar'*pr + B'*(r.*(q0+B*qc).*(abs(q0+B*qc)));
    elseif ind == 4 then
        // cas : F et G
        F = ((q0+B*qc)'*(r.*(q0+B*qc).*(abs(q0+B*qc))))/3 + pr'*(Ar*(q0+B*qc));
        G = B'*Ar'*pr + B'*(r.*(q0+B*qc).*(abs(q0+B*qc)));
    
endfunction


// ---------------------------
// OraclePH
// ---------------------------

function [F, G, H, ind] = OraclePH(qc, ind)

// ---------------------------
// Arguments and Return values
// ---------------------------
    
// qc : vecteur des debits reduits
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
// Initialisation
// ---------------------------
    
    F = 0;
    G = zeros(n-md,1);
    H = zeros(n-md,n-md);
    ind = ind;
    
    if ind == 2 then
        // cas : seulement F
        F = ((q0+B*qc)'*(r.*(q0+B*qc).*(abs(q0+B*qc))))/3 + pr'*(Ar*(q0+B*qc));
    elseif ind == 3 then
        // cas : seulement G
        G = B'*Ar'*pr + B'*(r.*(q0+B*qc).*(abs(q0+B*qc)));
    elseif ind == 4 then
        // cas : F et G
        F = ((q0+B*qc)'*(r.*(q0+B*qc).*(abs(q0+B*qc))))/3 + pr'*(Ar*(q0+B*qc));
        G = B'*Ar'*pr + B'*(r.*(q0+B*qc).*(abs(q0+B*qc)));
    elseif ind == 5 then
        // cas : seulement H
        H = 2*B'*diag(r.*(abs(q0+B*qc)))*B
    elseif ind == 6 then
        // cas : G et H
        G = B'*Ar'*pr + B'*(r.*(q0+B*qc).*(abs(q0+B*qc)));
        H = 2*B'*diag(r.*(abs(q0+B*qc)))*B
    elseif ind == 7 then
        // cas : F, G et H
        F = ((q0+B*qc)'*(r.*(q0+B*qc).*(abs(q0+B*qc))))/3 + pr'*(Ar*(q0+B*qc));
        G = B'*Ar'*pr + B'*(r.*(q0+B*qc).*(abs(q0+B*qc)));
        H = 2*B'*diag(r.*(abs(q0+B*qc)))*B
    end
    
endfunction

