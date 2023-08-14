function [hcli_br] = loadHomoclinicFromPerDDETool(pt,funcs,indPar1,indPar2,dS)
%LOADHOMOCLINICFROMAUTO Summary of this function goes here
%   Detailed explanation goes here

% create the homoclinic object.
hcli          =  p_tohcli(funcs,pt);
mhcli         =  df_mthod(funcs,'hcli');
[hcli, s]     =  p_correc(funcs,hcli,indPar2,[],mhcli.point);   % correct 
hcli          =  p_remesh(hcli,3,121);                          % remesh it and
[hcli, s]     =  p_correc(funcs,hcli,indPar2,[],mhcli.point);   % correct it again

% create the homoclinic branch object.
hcli_br          =  df_brnch(funcs,[indPar1, indPar2],'hcli');
hcli_br.point    = hcli;
hcli.parameter(indPar2) = hcli.parameter(indPar2)+dS;
[hcli,s]                = p_correc(funcs,hcli,indPar2,[],mhcli.point);
hcli_br.point(2) = hcli;
end

