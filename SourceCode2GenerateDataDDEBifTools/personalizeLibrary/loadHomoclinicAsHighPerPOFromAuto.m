function [per_br] = loadHomoclinicAsHighPerPOFromAuto(nameFile,funcs_per,indPar1,indPar2,indper,dS)
%LOADHOMOCLINICASHIGHPERPOFROMAUTO Summary of this function goes here
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
pt.parameter(indper) = parAuto(10); % for constant period continuation
pt           =   p_remesh(pt,3,121);
method       =   df_mthod(funcs_per,'psol');
[pt,suc]     =   p_correc(funcs_per,pt,[indPar1, indPar2],[],method.point);

% create the homoclinic branch object.
per_br          =  df_brnch(funcs_per,[indPar1, indPar2],'psol');
per_br.point    =  pt;
pt.parameter(indPar2) = pt.parameter(indPar2)+dS;
[pt,suc]        = p_correc(funcs_per,pt,[indPar1, indPar2],[],method.point);
per_br.point(2) = pt;
end

