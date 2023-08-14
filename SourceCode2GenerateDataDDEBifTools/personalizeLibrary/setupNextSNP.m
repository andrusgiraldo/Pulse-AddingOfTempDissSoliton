function [PDfuncs,PDbranch] = setupNextSNP(pt,funcs,indPar1,indPar2,dS,suiton)
%SETUPNEXTPD Summary of this function goes here
%   Detailed explanation goes here
method       =   df_mthod(funcs,'psol');
pt.parameter     =   pt.parameter(1:11);
pt.parameter(10) =   pt.period;
pt.profile   =   pt.profile(1:3,:);

if suiton>0
    pt.parameter(9) = pt.parameter(9) + suiton*pt.period;
end

auxMesh      =   pt.mesh;
pt.mesh      =   linspace(0,1,900);
pt.profile   =   interp1(auxMesh',pt.profile',pt.mesh','pchip')';
pt.degree    =   1; 

pt           =   p_remesh(pt,3,121);
method       =   df_mthod(funcs,'psol');
[pt,suc]     =   p_correc(funcs,pt,[],[],method.point);

br           =   df_brnch(funcs,indPar2,'psol');
br.point     =   pt;   

% create the PD Object
[PDfuncs,PDbranch] = SetupPOfold(funcs,br,1,'contpar',[indPar1 indPar2],'dir',indPar1,'step',dS);
end

