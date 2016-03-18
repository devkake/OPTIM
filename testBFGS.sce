
disp('This test is for OraclePG, merged');
exec('Probleme_R.sce');
exec('Structures_R.sce');
exec('Visualg.sci');
exec('HydrauliqueP.sci'); 
exec('Verification.sci');
// exec('Oracle.sci'); 
exec('OraclePG.sci'); 
exec('OraclePH.sci'); 
exec('Wolfe_Skel.sci');
exec('BFGS.sci');
qini = 0.1 * rand(n-md,1);
[Fopt,qopt,Gopt] = BFGS(OraclePG,qini);
[q,z,f,p] = HydrauliqueP(qopt); 
Verification(q,z,f,p);
