exec('Probleme_R.sce');
exec('Structures_R.sce');
exec('Visualg.sci');
exec('HydrauliqueP.sci');
exec('HydrauliqueD.sci');  
exec('Verification.sci');
// exec('Oracle.sci'); 
exec('OraclePG.sci'); 
exec('OraclePH.sci'); 
exec('OracleDG.sci'); 
exec('OracleDH.sci'); 
exec('Wolfe_Skel.sci');
exec('Gradient.sci');
// qini = 0.1 * rand(n-md,1);
// qini = 0.2 * ones(n-md, 1);
// qini = 0.1 * rand(md,1);
qini = 0.1 * ones(md, 1);
for i = 1:md
    qini(i) = 0.01*i;
end
// [Fopt,qopt,Gopt] = Gradient(OraclePG,qini);
[Fopt,qopt,Gopt] = Gradient(OracleDG,qini);
// [q,z,f,p] = HydrauliqueP(qopt); 
[q,z,f,p] = HydrauliqueD(qopt); 
Verification(q,z,f,p);
