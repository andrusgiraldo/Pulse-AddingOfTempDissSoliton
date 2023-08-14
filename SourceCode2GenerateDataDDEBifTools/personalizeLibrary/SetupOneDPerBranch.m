function [per_br] = SetupOneDPerBranch(pt, funcs, indPar, indParFix, setParFix, dS)
%SETUPONEDPERBRANCH break homoclinic in one parameter
%   Detailed explanation goes here

% create the psol branch object.
per_br          =   df_brnch(funcs,indPar,'psol');
method          =   df_mthod(funcs,'psol');

if isequal(pt.kind,'hcli')
    pt           =   p_topsol(funcs,pt);
else
    pt.degree    =   3; 
end
pt.parameter(indParFix) = setParFix;
[pt,suc]        =   p_correc(funcs,pt, [],[],method.point);
per_br.point    =   pt;
pt.parameter(indPar) = pt.parameter(indPar) + dS;
[pt,suc]        =   p_correc(funcs,pt, [],[],method.point);
per_br.point(2) =   pt;
end

