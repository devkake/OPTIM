function [alphan,ok]=Wolfe(alpha,x,D,Oracle)


//////////////////////////////////////////////////////////////
//                                                          //
//   RECHERCHE LINEAIRE SUIVANT LES CONDITIONS DE WOLFE     //
//                                                          //
//                                                          //
//  Arguments en entree                                     //
//  -------------------                                     //
//    alpha  : valeur initiale du pas                       //
//    x      : valeur initiale des variables                //
//    D      : direction de descente                        //
//    Oracle : nom de la fonction Oracle                    //
//                                                          //
//  Arguments en sortie                                     //
//  -------------------                                     //
//    alphan : valeur du pas apres recherche lineaire       //
//    ok     : indicateur de reussite de la recherche       //
//             = 1 : conditions de Wolfe verifiees          //
//             = 2 : indistinguabilite des iteres           //
//                                                          //
//                                                          //
//    omega1 : coefficient pour la 1-ere condition de Wolfe //
//    omega2 : coefficient pour la 2-eme condition de Wolfe //
//                                                          //
//////////////////////////////////////////////////////////////


// -------------------------------------
// Coefficients de la recherche lineaire
// -------------------------------------

   omega1 = 0.1; // coefficient pour la 1-ere condition de Wolfe
   omega2 = 0.9; // coefficient pour la 2-eme condition de Wolfe

   alphamin = 0.0; // initial alpha inf
   alphamax = %inf; // initial alpha sup

   ok = 0; // pour le boucle de while
   dltx = 0.00000001; // tolerance

// ---------------------------------
// Algorithme de Fletcher-Lemarechal
// ---------------------------------

   // Appel de l'oracle au point initial
   
   ind = 4;
   [F,G] = Oracle(x,ind);
   beta = G'*D; // <(nabla)J(x), D>

   // Initialisation de l'algorithme

   alphan = alpha; // initialisation de alpha
   // alphan = 2*alpha
   // alphan = -2*(alpha)/beta // pas de Fletcher
   xn     = x; // initialisation de xn

   // Boucle de calcul du pas
   //
   // xn represente le point pour la valeur courante du pas,
   // xp represente le point pour la valeur precedente du pas.

   while ok == 0
      
      xp = xn; // le point pour la valeur precedente du pas
      xn = x + (alphan*D); // mettre xn a jour

      // Calcul des conditions de Wolfe

      // -----> A completer...
      // -----> A completer...
      [Fn,Gn] = Oracle(xn, 4);
      cond1 = Fn-F <= omega1*alphan*beta;
      cond2 = Gn'*D >= omega2*beta;
      

      // Test de la valeur de alphan :
      // - si les deux conditions de Wolfe sont verifiees,
      //   faire ok = 1 : on sort alors de la boucle while
      // - sinon, modifier la valeur de alphan : on reboucle.

      // -----> A completer...
      // -----> A completer...
      if ~cond1 then
          alphamax = alphan;
          alphan = (alphamin+alphamax)/2.0;
      else
          if ~cond2 then
              alphamin = alphan;
              if isinf(alphamax)  then
                  alphan = 2*alphamin;
              else
                  alphan = (alphamin+alphamax)/2.0;
              end
          else
              ok = 1;
          end
      end

      // Test d'indistinguabilite

      if norm(xn-xp) < dltx then
        ok = 2;
      end

   end

endfunction
