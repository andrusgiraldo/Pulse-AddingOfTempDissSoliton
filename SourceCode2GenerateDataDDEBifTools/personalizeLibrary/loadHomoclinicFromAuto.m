function [hcli_br] = loadHomoclinicFromAuto(nameFile,funcs,indPar1,indPar2,dS,epsAux)
%LOADHOMOCLINICFROMAUTO Summary of this function goes here
%   Detailed explanation goes here

% loading data from auto and the corresponding parameters.
solAuto      =  load([nameFile '_solution']);
parAuto      =  load([nameFile '_PARAMS']);

% create a point structure to load after 
pt.period    =   parAuto(10);
pt.mesh      =   solAuto(:, 1);
pt.mesh      =   linspace(0,1,900);
pt.profile   =   interp1(solAuto(:, 1),solAuto(:, [2, 3, 4]),pt.mesh,'pchip').';
pt.kind      =   'psol';
pt.degree    =   1; 
pt.parameter =   parAuto'; 
pt           =   p_remesh(pt,4,301);
method       =   df_mthod(funcs,'psol');
[pt,suc]     =   p_correc(funcs,pt,[],[],method.point);

% create the homoclinic object.
hcli          =  p_tohcli(funcs,pt);
mhcli         =  df_mthod(funcs,'hcli');
hcli.epsilon  =  epsAux;
[hcli, s]     =  p_correc(funcs,hcli,indPar2,[],mhcli.point);   % correct 
hcli          =  p_remesh(hcli,4,301);                          % remesh it and
[hcli, s]     =  p_correc(funcs,hcli,indPar2,[],mhcli.point);   % correct it again

% create the homoclinic branch object.
hcli_br          =  df_brnch(funcs,[indPar1, indPar2],'hcli');
hcli_br.point    =  hcli;
hcli.parameter(indPar2) = hcli.parameter(indPar2)+dS;
[hcli,s]                = p_correc(funcs,hcli,indPar2,[],mhcli.point);
hcli_br.point(2) = hcli;
end

