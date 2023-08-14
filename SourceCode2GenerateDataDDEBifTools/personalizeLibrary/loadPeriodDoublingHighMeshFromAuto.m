function [PDfuncs, PDbranch] = loadPeriodDoublingHighMeshFromAuto(nameFile,funcs,indPar1,indPar2,dS)
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
pt.parameter =   [parAuto' 0]; 
pt           =   p_remesh(pt,4,200);
method       =   df_mthod(funcs,'psol');
[pt,suc]     =   p_correc(funcs,pt,[],[],method.point);

% create the period
br           =   df_brnch(funcs,indPar2,'psol');
br.point     =   pt;   
pt.parameter(indPar2) =  1e-6;
[pt,suc]     =   p_correc(funcs,pt,[],[],method.point);
br.point(2)  =   pt;


% create the PD Object
br                 = br_stabl(funcs,br,0,1);
[PDfuncs,PDbranch] = SetupPeriodDoubling(funcs,br,2,'contpar',[indPar1 indPar2],'dir',indPar1,'step',dS);
end

