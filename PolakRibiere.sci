function [fopt,xopt,gopt]=PolakRibiere(Oracle,xini)


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
//                                                                           //
//         Methode de gradient conjugue a pas non-fixe                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


// ------------------------
// Parametres de la methode
// ------------------------

   titre = "Parametres du gradient conjugue a pas non-fixe";
   labels = ["Nombre maximal d''iterations";...
             "Valeur du pas de gradient";...
             "Seuil de convergence sur ||G||"];
   typ = list("vec",1,"vec",1,"vec",1);
   default = ["5000";"1.0";"0.000001"];
   [ok,iter,alphai,tol] = getvalue(titre,labels,typ,default);

// ----------------------------
// Initialisation des variables
// ----------------------------
   
   logG = [];
   logP = [];
   Cout = [];

   timer();

// -------------------------
// Boucle sur les iterations
// -------------------------

   x = xini;
   [F, G_prec] = Oracle(xini,3);
   disp(G_prec)
   D = 0; // Premier D est le gradient (D = -G + beta*D;)
   
   kstar = iter;
   for k = 1:iter

//    - valeur du critere et du gradient

      ind = 4;
      [F,G] = Oracle(x,ind); // Oracle est l'argument

//    ind = 7;
//    [F,G,H] = OraclePH(x,ind);
//    disp(H);

//    - calcul de la direction de descente
        
      beta = (G - G_prec)'*G / (G_prec'*G_prec);

      D = -G + beta*D;

//    - calcul de la longueur du pas de gradient

      [alpha,ok]=Wolfe(alphai,x,D,Oracle); // les conditions de Wolfe
      if ok == 2 then
          disp('Un pas converge. Peut-etre, le pas ne satisfait pas les conditions de Wolfe')
      end
      // alpha = alphai: // pas fixe

//    - mise a jour des variables

      x = x + (alpha*D);
      
      G_prec = G;

//    - evolution du gradient, du pas et du critere

      logG = [ logG ; log10(norm(G)) ];
      logP = [ logP ; log10(alpha) ];
      Cout = [ Cout ; F ];
      
//    - test de convergence

      if norm(G) <= tol then
         kstar = k;
         break
      end
      
   end

   disp(x);
// ---------------------------
// Resultats de l'optimisation
// ---------------------------

   fopt = F;
   xopt = x;
   gopt = G;

   tcpu = timer();

   cvge = ['Iteration         : ' + string(kstar);...
           'Temps CPU         : ' + string(tcpu);...
           'Critere optimal   : ' + string(fopt);...
           'Norme du gradient : ' + string(norm(gopt))];

   disp('Fin de la methode de gradient a pas fixe')
   disp(cvge)

// - visualisation de la convergence

   Visualg(logG,logP,Cout);

endfunction
