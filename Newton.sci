function [fopt,xopt,gopt]=Newton(Oracle,xini)


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
//                                                                           //
//         Methode de Newton a pas non-fixe                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


// ------------------------
// Parametres de la methode
// ------------------------

   titre = "Parametres de Newton a pas non-fixe";
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
   tol = 1e-10 // tolerance for inverse of hessian

   kstar = iter;
   for k = 1:iter

//    - valeur du critere et du gradient

      ind = 7;
      [F,G,H] = Oracle(x,ind); // Oracle est l'argument


//    - calcul de la direction de descente
      
      regularity = rcond(H) > tol; // H is regular (inversible) or not
      if regularity then
          D = -H\G;
          // D = -inv(H)*G
      else
          // H is not regular (inversible) or not
          D = -G; // Descent du gradient
          
//          // @https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization
//          // Levenberg–Marquardt algorithm
//          mu = 0.1
//          while ~regularity
//              H = H + mu*eye(H);
//              regularity = rcond(H) > tol
//          end
//          D = -H\G;
          
      end

//    - calcul de la longueur du pas de gradient

      [alpha,ok]=Wolfe(alphai,x,D,Oracle); // les conditions de Wolfe
      if ok == 2 then
          disp('Un pas converge. Peut-etre, le pas ne satisfait pas les conditions de Wolfe')
      end
      // alpha = alphai: // pas fixe

//    - mise a jour des variables

      x = x + (alpha*D);

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
