function [SNPfuncs, SNPbranch] = loadSaddleNodePerFromAuto(nameFile,funcs,indPar1,indPar2,dS)
%LOADHOMOCLINICFROMAUTO Summary of this function goes here
%   Detailed explanation goes here

% loading data from auto and the corresponding parameters.
solAuto      =  load([nameFile '_solution']);
parAuto      =  load([nameFile '_PARAMS']);

% create a point structure to load after 
pt.period    =   parAuto(10);
pt.mesh      =   linspace(0,1,900);
pt.profile   =   interp1(solAuto(:, 1),solAuto(:, [2, 3, 4]),pt.mesh,'pchip').';
pt.kind      =   'psol';
pt.degree    =   1; 
pt.parameter =   [parAuto' 0]; 
pt           =   p_remesh(pt,3,121);
method       =   df_mthod(funcs,'psol');

% create the period
br           =   df_brnch(funcs,indPar2,'psol');
br.point     =   pt;   

% create the SNP Object
[SNPfuncs,SNPbranch]  = SetupPOfold(funcs,br,1,'contpar',[indPar1 indPar2],'dir', indPar2 ,'step',dS);
end

