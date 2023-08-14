%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%  - Script that generates the data for the bifurcation diagram -  %%%
% %%%%%%%%%%    - in parameter plane used in Figures 5, 8 12  %%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all;
close all;

% Adding DDE Biftools in the search path: 
% CHANGE this line to the path where your DDE Biftools installation is
addpath('../../../../../MATLAB/dde_biftool_v3.1.1/ddebiftool/',...
        '../../../../../MATLAB/dde_biftool_v3.1.1/demos/phase_oscillator',...
        '../../../../../MATLAB/dde_biftool_v3.1.1/ddebiftool_extra_psol/',...
        '../../../../../MATLAB/dde_biftool_v3.1.1/ddebiftool_utilities/',...
        '../../../../../MATLAB/dde_biftool_v3.1.1/ddebiftool_extra_nmfm/',...
        '../../../../../MATLAB/dde_biftool_v3.1.1/ddebiftool_extra_rotsym'); 


% Adding a collection of functions that are used to initialize different
% variables for DDE BifTools
addpath('./personalizeLibrary')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%      Definition of the Delay differential model     %%%%%%%%%%%
% %%%%%%      used in the manuscript: Sandstede model with Delay     %%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cleaning the screen and variables
clc
clear all

% Vectors with the names of the variables, and their corresponding index
parnames    =   {'a','b','c','alpha','beta', 'gamma', 'mu', 'mu_tilde', 'tau', 'T','kappa'};
c_ind       =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(c_ind{:});

% Definition of Sandstede's model with Delay
f           =   @(zz,p) [...
                % first component
                p(ind.a).*zz(1,1,:) + p(ind.b).*zz(2,1,:) - p(ind.a).*zz(1,1,:).^2 + zz(1,1,:).*(p(ind.mu_tilde)-p(ind.alpha).*zz(3,1,:)).*(2-3.*zz(1,1,:)); ...
                % second component
                p(ind.b).*zz(1,1,:) + p(ind.a).*zz(2,1,:) - 3/2.*p(ind.b).*zz(1,1,:).^2 - 3/2.*p(ind.a).*zz(1,1,:).*zz(2,1,:) + p(ind.kappa).*zz(1,2,:).*zz(2,2,:) - 2*zz(2,1,:).*(p(ind.mu_tilde)-p(ind.alpha).*zz(3,1,:)); ...
                % third component
                p(ind.c).*zz(3,1,:) + p(ind.mu).*zz(1,1,:) + p(ind.gamma).*zz(1,1,:).*zz(3,1,:) + p(ind.alpha).*p(ind.beta).*( zz(1,1,:).^2.*(1-zz(1,1,:))-zz(2,1,:).^2 )];

% Definition of the DDE function for equilibria for DDE Biftools
funcs       =   set_funcs('sys_rhs', f, 'sys_tau', @()ind.tau,...
                'x_vectorized',true); 
% Definition of the DDE function for periodic solutions for DDE Biftools
funcs_per   =   set_funcs('sys_rhs', f, 'sys_tau', @()ind.tau,...
                'x_vectorized',true ,'sys_cond', @(p)copy_period(p,ind.T));    

prefixFolderName  =   '../AutoResonantBifurcationSolution/Matlab/BifPlaneSolution_';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%             Bifurcation diagram for the non-orientable         %%%%%%
%%%%%%%%%%%%%%       resonant case near Orbit bifurcation          %%%%%%%%
%%%%%%%%%%%%%%%%        Left resonant point in Figure 8            %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variable to store the curves. 
auxCurves         =  cell(2,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the homoclinic solution from Auto07p using the fixed period 
% formulation and run the continuation scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sufixFolderName   =  'OrbitFlipNoNOriPaper/HOM_1_-2.0620_0.0550';
auxPerHom_br      =   loadHomoclinicAsHighPerPOFromAuto([prefixFolderName, sufixFolderName],funcs_per,ind.mu,ind.kappa,ind.T, -1e-01);
auxPerHom_br.parameter.max_step  = [ind.kappa  0.05;ind.mu 0.05];
auxPerHom_br.parameter.max_bound = [ind.mu   0.2];
auxPerHom_br.parameter.min_bound = [ind.mu  -0.2];
figure(1); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,200);

% Now Move kappa
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.kappa,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.tau) = pt.parameter(ind.tau)-0.1;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.kappa,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.kappa  0.1;ind.tau 0.1];
auxPerHom_br.parameter.max_bound = [ind.kappa  2];
auxPerHom_br.parameter.min_bound = [ind.kappa -1];
figure(2); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,200);

% Now Move a
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.a,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.a) = pt.parameter(ind.a)+0.1;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.a,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.a  0.025;ind.tau 0.025];
auxPerHom_br.parameter.max_bound = [];
auxPerHom_br.parameter.min_bound = [];
figure(3); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,300);
auxPerHom_br          =   br_rvers(auxPerHom_br);
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,300);
auxCurves{1}      = auxPerHom_br;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the two-homoclinic solution using the fixed period formulation 
% and run the continuation scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sufixFolderName   =  'OrbitFlipNoNOriPaper/HOM2_1_-0.7490_0.0550';
auxPerHom_br      =   loadHomoclinicAsHighPerPOFromAuto([prefixFolderName, sufixFolderName],funcs_per,ind.mu,ind.kappa,ind.T, +1e-01);
auxPerHom_br.parameter.max_step  = [ind.kappa  0.05;ind.mu 0.05];
auxPerHom_br.parameter.max_bound = [ind.mu   0.2];
auxPerHom_br.parameter.min_bound = [ind.mu  -0.2];
figure(1); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,200);

% Now Move kappa
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.kappa,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.tau) = pt.parameter(ind.tau)-0.1;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.kappa,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.kappa  0.1;ind.tau 0.1];
auxPerHom_br.parameter.max_bound = [ind.kappa  1];
auxPerHom_br.parameter.min_bound = [ind.kappa -1];
figure(2); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,200);

% Now Move a
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.a,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.a) = pt.parameter(ind.a)+0.01;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.a,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.a  0.03;ind.tau 0.03];
auxPerHom_br.parameter.max_bound = [];
auxPerHom_br.parameter.min_bound = [];
figure(3); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,300);
auxPerHom_br          =   br_rvers(auxPerHom_br);
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,400);
auxCurves{2}   = auxPerHom_br;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load of the Period-doubling bifuraction run the continuation scheme 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sufixFolderName           =  'OrbitFlipNoNOriPaper/PD_1_-0.6414_0.0550';
[auxPDfuncs, auxPDbranch] = loadPeriodDoublingFromAuto([prefixFolderName, sufixFolderName],funcs,ind.mu,ind.kappa, 1e-03);
auxPDbranch.parameter.max_step  = [ind.kappa  0.05;ind.mu 0.05];
auxPDbranch.parameter.max_bound = [ind.mu   0.2];
auxPDbranch.parameter.min_bound = [ind.mu  -0.2];
figure(1);hold on;
auxPDbranch    = br_contn(auxPDfuncs,auxPDbranch,200);

% Here move kappa
pt                         =  auxPDbranch.point(end);
[auxPDfuncs,auxPDbranch]   = setupNextPD(pt,funcs,ind.kappa,ind.tau,-1e-03,0);
auxPDbranch.parameter.max_step  = [ind.kappa  0.1;ind.tau 0.1];
auxPDbranch.parameter.max_bound = [ind.kappa  1];
auxPDbranch.parameter.min_bound = [ind.kappa -1];
figure(2);hold on;
auxPDbranch    = br_contn(auxPDfuncs,auxPDbranch,200);

% Here move a, but we need to map tau to positive, so we can continue
pt                           =  auxPDbranch.point(end);
[auxPDfuncs1,auxPDbranch1]   = setupNextPD(pt,funcs,ind.a,ind.tau,1e-03,2); % you need the period double to reappear
auxPDbranch1.parameter.max_step  = [ind.a  0.025;ind.tau 0.025];
auxPDbranch1.parameter.max_bound = [];
auxPDbranch1.parameter.min_bound = [];
figure(3);hold on;
auxPDbranch1    = br_contn(auxPDfuncs1,auxPDbranch1,400);
auxPDbranch1          =   br_rvers(auxPDbranch1);
[auxPDbranch1,s,f,r]  =   br_contn(auxPDfuncs1,auxPDbranch1,400);

% Remapping everything back to negative delay
auxPDbranch = auxPDbranch1;
for i =1:length(auxPDbranch1.point)
    auxPDbranch.point(i).parameter(ind.tau) = auxPDbranch1.point(i).parameter(ind.tau) - 2*auxPDbranch1.point(i).period;
end
auxCurves{3}    = auxPDbranch;


aux      = auxCurves{1};  % This is the curve!!!
figure(10)
hold on
plot3(aux.point(1).profile(1,:),aux.point(1).profile(2,:),aux.point(1).profile(3,:))
plot3(aux.point(end).profile(1,:),aux.point(end).profile(2,:),aux.point(end).profile(3,:),'linewidth',2.0)
scatter3(0,0,0,100,'filled')

% Undocument the next line to save the curves
%save('./DDEResults/BifCurvesNonOriResPaper_a_tau_closeOrbitFlip','auxCurves') 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%             Bifurcation diagram for the Orientable             %%%%%%
%%%%%%%%%%%%%%       resonant case near Orbit bifurcation          %%%%%%%%
%%%%%%%%%%%%%%%%                  Figure 5                         %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

auxCurves         =  cell(2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load of the Homoclinic bifuraction run the continuation scheme 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sufixFolderName       =  'OrbitFlipOriPaper/HOM_1_-2.1220_-0.0550';
auxPerHom_br          =   loadHomoclinicAsHighPerPOFromAuto([prefixFolderName, sufixFolderName],funcs_per,ind.a,ind.kappa,ind.T,-1e-01);
auxPerHom_br.parameter.max_step  = [ind.kappa  0.05;ind.a 0.05];
auxPerHom_br.parameter.max_bound = [ind.a  -0.5];
auxPerHom_br.parameter.min_bound = [ind.a  -2.5];
figure(1); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,200);

% Now Move mu, I have to stop at mu_tilde
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.mu_tilde,ind.kappa],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.kappa) = pt.parameter(ind.kappa)-0.1;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.mu_tilde,ind.kappa],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.mu_tilde  0.05;ind.kappa 0.05];
auxPerHom_br.parameter.max_bound = [ind.kappa   1.0];
auxPerHom_br.parameter.min_bound = [ind.kappa  -1.0];
figure(2); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,200);

% Now Move mu_tilde back
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.mu_tilde,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.tau) = pt.parameter(ind.tau)-0.01;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.mu_tilde,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.mu_tilde  0.1;ind.tau 0.1];
auxPerHom_br.parameter.max_bound = [ind.mu_tilde   -0.055];
auxPerHom_br.parameter.min_bound = [ind.mu_tilde  -1];
figure(3); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,200);


% Now Move mu
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.mu,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.tau) = pt.parameter(ind.tau)-0.1;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.mu,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.mu   0.1;  ind.tau 0.1];
auxPerHom_br.parameter.max_bound = [ind.mu   0.2];
auxPerHom_br.parameter.min_bound = [ind.mu  -0.2];
figure(4); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,200);


% Now Move a
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.a,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.a) = pt.parameter(ind.a)+0.1;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.a,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.a  0.01;ind.tau 0.01];
auxPerHom_br.parameter.max_bound = [];
auxPerHom_br.parameter.min_bound = [];
figure(5); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,400);

auxCurves{1}          = auxPerHom_br;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load of the Saddle-node run the continuation scheme 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sufixFolderName             =  'OrbitFlipOriPaper/SNP_1_-0.3852_-0.0550';
[auxSNPfuncs, auxSNPbranch] = loadSaddleNodePerFromAuto([prefixFolderName, sufixFolderName],funcs,ind.a,ind.kappa,-1e-03);
auxSNPbranch.parameter.max_step  = [ind.a   0.05;ind.kappa 0.05];
auxSNPbranch.parameter.max_bound = [ind.a  -0.01];
auxSNPbranch.parameter.min_bound = [ind.a  -0.8];
figure(1);hold on;
[auxSNPbranch,s,f,r]   =   br_contn(auxSNPfuncs,auxSNPbranch,200);

% Here move kappa
pt                           = auxSNPbranch.point(end);
[auxSNPfuncs,auxSNPbranch]   = setupNextSNP(pt,funcs,ind.mu_tilde,ind.kappa,-1e-03,0);
auxSNPbranch.parameter.max_step  = [ind.mu_tilde  0.05;ind.tau 0.1];
auxSNPbranch.parameter.max_bound = [ind.kappa  1];
auxSNPbranch.parameter.min_bound = [ind.kappa -1];
figure(2);hold on;
auxSNPbranch    = br_contn(auxSNPfuncs,auxSNPbranch,200);

% Here move mu_tilde
pt                           = auxSNPbranch.point(end);
[auxSNPfuncs,auxSNPbranch1]   = setupNextSNP(pt,funcs,ind.mu_tilde,ind.tau,1e-03,0);
auxSNPbranch1.parameter.max_step  = [ind.mu_tilde  0.1;ind.tau 0.1];
auxSNPbranch1.parameter.max_bound = [ind.mu_tilde   -0.055];
auxSNPbranch1.parameter.min_bound = [ind.mu_tilde   -1];
figure(3);hold on;
auxSNPbranch1    = br_contn(auxSNPfuncs,auxSNPbranch1,200);
auxSNPbranch     = auxSNPbranch1;

% Here move mu
pt                           = auxSNPbranch.point(end);
[auxSNPfuncs,auxSNPbranch1]   = setupNextSNP(pt,funcs,ind.mu,ind.tau,-1e-03,0);
auxSNPbranch1.parameter.max_step  = [ind.mu  0.1;ind.tau 0.1];
auxSNPbranch1.parameter.max_bound = [ind.mu   0.2];
auxSNPbranch1.parameter.min_bound = [ind.mu  -0.2];
figure(4);hold on;
auxSNPbranch1    = br_contn(auxSNPfuncs,auxSNPbranch1,200);
auxSNPbranch     = auxSNPbranch1;


% Here move a
pt                           = auxSNPbranch.point(end);
[auxSNPfuncs,auxSNPbranch1]   = setupNextSNP(pt,funcs,ind.a,ind.tau,-1e-03,0);
auxSNPbranch1.parameter.max_step  = [ind.a  0.1;ind.tau 0.1];
auxSNPbranch1.parameter.max_bound = [];
auxSNPbranch1.parameter.min_bound = [];
figure(5);hold on;
auxSNPbranch1         = br_contn(auxSNPfuncs,auxSNPbranch1,100);
auxSNPbranch1         = br_rvers(auxSNPbranch1);
auxSNPbranch1.parameter.max_step  = [ind.a  0.01;ind.tau 0.01];
[auxSNPbranch1,s,f,r]  =   br_contn(auxSNPfuncs,auxSNPbranch1,300);
auxSNPbranch     = auxSNPbranch1;

% Here move mu_tilde

auxCurves{2}           =   auxSNPbranch;

aux      = auxSNPbranch;  % This is the curve!!!
figure(10)
hold on
plot3(aux.point(1).profile(1,:),aux.point(1).profile(2,:),aux.point(1).profile(3,:))
plot3(aux.point(end).profile(1,:),aux.point(end).profile(2,:),aux.point(end).profile(3,:),'linewidth',2.0)
scatter3(0,0,0,100,'filled')

% Undocument this line to save 
% save('./DDEResults/BifCurvesOriResPaper_a_tau_closeOrbitFlip','auxCurves') 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%        Bifurcation diagram for the Orbit bifurcation         %%%%%%%%
%%%%%%%%%%%%%%%%                Figure 12                          %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
auxCurves         =  cell(3,1);

% load the homoclinic solution using the fixed period formulation and run the continuation scheme
sufixFolderName   =  'OrbitFlipPaper/HOM_1_0.0595_0.0300';%'OrbitFlipNoNOri/HOM_1_-2.2976_0.0700';
auxPerHom_br      =   loadHomoclinicAsHighPerPOFromAuto([prefixFolderName, sufixFolderName],funcs_per,ind.mu,ind.kappa,ind.T,-1e-02);
auxPerHom_br.parameter.max_step  = [ind.kappa  0.05;ind.mu 0.05];
auxPerHom_br.parameter.max_bound = [ind.mu   0.2];
auxPerHom_br.parameter.min_bound = [ind.mu  -0.2];
figure(1); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,200);

% Now Move kappa
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.kappa,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.tau) = pt.parameter(ind.tau)-0.1;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.kappa,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.kappa  0.05;ind.tau 0.05];
auxPerHom_br.parameter.max_bound = [ind.kappa  2];
auxPerHom_br.parameter.min_bound = [ind.kappa -1];
figure(2); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,200);

% Moving mu_tilde to zero
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.mu_tilde,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.tau) = pt.parameter(ind.tau)+0.01;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.mu_tilde,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.mu_tilde  0.05;ind.tau 0.05];
auxPerHom_br.parameter.max_bound = [];
auxPerHom_br.parameter.min_bound = [ind.mu_tilde  0];
figure(3); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,100);

% Here is the orbitFlip point
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.mu,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.tau) = pt.parameter(ind.tau)-0.01;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.mu,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.mu  0.01;ind.tau 0.01];
auxPerHom_br.parameter.max_bound = [ind.tau  6];
auxPerHom_br.parameter.min_bound = [ind.tau -6];
figure(4); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,500);
auxPerHom_br          =   br_rvers(auxPerHom_br);
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,500);

auxCurves{1}      = auxPerHom_br;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the two-homoclinic solution using the fixed period formulation 
% and run the continuation scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sufixFolderName   =  'OrbitFlipNoNOriPaper/HOM2_1_-0.7490_0.0550';
auxPerHom_br      =   loadHomoclinicAsHighPerPOFromAuto([prefixFolderName, sufixFolderName],funcs_per,ind.mu,ind.kappa,ind.T, +1e-01);
auxPerHom_br.parameter.max_step  = [ind.kappa  0.05;ind.mu 0.05];
auxPerHom_br.parameter.max_bound = [ind.mu   0.2];
auxPerHom_br.parameter.min_bound = [ind.mu  -0.2];
figure(1); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,200);

% Now Move kappa
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.kappa,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.tau) = pt.parameter(ind.tau)-0.1;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.kappa,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.kappa  0.1;ind.tau 0.1];
auxPerHom_br.parameter.max_bound = [ind.kappa  1];
auxPerHom_br.parameter.min_bound = [ind.kappa -1];
figure(2); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,200);

% Now Move a
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.a,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.a) = pt.parameter(ind.a)+0.01;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.a,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.a  0.05;ind.tau 0.05];
auxPerHom_br.parameter.max_bound = [ind.a -0.5];
auxPerHom_br.parameter.min_bound = [];
figure(10); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,200);

% Moving mu_tilde to zero
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.mu_tilde,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.tau) = pt.parameter(ind.tau)+0.01;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.mu_tilde,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.mu_tilde  0.05;ind.tau 0.05];
auxPerHom_br.parameter.max_bound = [];
auxPerHom_br.parameter.min_bound = [ind.mu_tilde  0];
figure(3); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,100);

% Here is the orbitFlip point
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.mu,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.tau) = pt.parameter(ind.tau)-0.01;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.mu,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.mu  0.01;ind.tau 0.01];
auxPerHom_br.parameter.max_bound = [ind.tau  6];
auxPerHom_br.parameter.min_bound = [ind.tau -6];
figure(4); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,200);
auxPerHom_br          =   br_rvers(auxPerHom_br);
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,200);
auxCurves{2}      = auxPerHom_br;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load of the Period-doubling bifuraction run the continuation scheme 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sufixFolderName           =  'OrbitFlipNoNOriPaper/PD_1_-0.6414_0.0550';
[auxPDfuncs, auxPDbranch] = loadPeriodDoublingFromAuto([prefixFolderName, sufixFolderName],funcs,ind.mu,ind.kappa, 1e-03);
auxPDbranch.parameter.max_step  = [ind.kappa  0.05;ind.mu 0.05];
auxPDbranch.parameter.max_bound = [ind.mu   0.2];
auxPDbranch.parameter.min_bound = [ind.mu  -0.2];
figure(1);hold on;
auxPDbranch    = br_contn(auxPDfuncs,auxPDbranch,200);

% Here move kappa
pt                         =  auxPDbranch.point(end);
[auxPDfuncs,auxPDbranch]   =  setupNextPD(pt,funcs,ind.kappa,ind.tau,-1e-03,0);
auxPDbranch.parameter.max_step  = [ind.kappa  0.1;ind.tau 0.1];
auxPDbranch.parameter.max_bound = [ind.kappa  1];
auxPDbranch.parameter.min_bound = [ind.kappa -1];
figure(2);hold on;
auxPDbranch    = br_contn(auxPDfuncs,auxPDbranch,200);

% Here move a, but I need to map tau to positive, so I can continue
pt                           =  auxPDbranch.point(end);
[auxPDfuncs1,auxPDbranch1]   =  setupNextPD(pt,funcs,ind.a,ind.tau,1e-03,2); % you need the period double to reappear
auxPDbranch1.parameter.max_step  = [ind.a  0.05;ind.tau 0.05];
auxPDbranch1.parameter.max_bound = [ind.a -0.5];
auxPDbranch1.parameter.min_bound = [];
figure(10);hold on;
auxPDbranch1    = br_contn(auxPDfuncs1,auxPDbranch1,200);
auxPDbranch     = auxPDbranch1;

% Here move mu_tilde, I am postive tau
pt                           =  auxPDbranch.point(end);
[auxPDfuncs1,auxPDbranch1]   =  setupNextPD(pt,funcs,ind.mu_tilde,ind.tau,-1e-03,0); 
auxPDbranch1.parameter.max_step  = [ind.mu_tilde  0.05;ind.tau 0.05];
auxPDbranch1.parameter.max_bound = [];
auxPDbranch1.parameter.min_bound = [ind.mu_tilde  0];
figure(3);hold on;
auxPDbranch1    = br_contn(auxPDfuncs1,auxPDbranch1,200);
auxPDbranch     = auxPDbranch1;

% Connecting Orbit flip at postivive (Stefan would love this)
pt                           =  auxPDbranch.point(end);
[auxPDfuncs1,auxPDbranch1]   =  setupNextPD(pt,funcs,ind.mu,ind.tau,-1e-03,0); 
auxPDbranch1.parameter.max_step  = [ind.mu  0.025;ind.tau 0.025];
auxPDbranch1.parameter.max_bound = [];
auxPDbranch1.parameter.min_bound = [];
figure(4);hold on;
auxPDbranch1    = br_contn(auxPDfuncs1,auxPDbranch1,300);
auxPDbranch1          =   br_rvers(auxPDbranch1);
[auxPDbranch1,s,f,r]  =   br_contn(auxPDfuncs1,auxPDbranch1,300);

% Remapping everything back to negative delay
auxPDbranch = auxPDbranch1;
for i =1:length(auxPDbranch1.point)
    auxPDbranch.point(i).parameter(ind.tau) = auxPDbranch1.point(i).parameter(ind.tau) - 2*auxPDbranch1.point(i).period;
end
auxCurves{3}    = auxPDbranch;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load of the Saddle-node run the continuation scheme 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sufixFolderName             =  'OrbitFlipOriPaper/SNP_1_-0.3852_-0.0550';
[auxSNPfuncs, auxSNPbranch] = loadSaddleNodePerFromAuto([prefixFolderName, sufixFolderName],funcs,ind.a,ind.kappa,-1e-03);
auxSNPbranch.parameter.max_step  = [ind.a   0.05;ind.kappa 0.05];
auxSNPbranch.parameter.max_bound = [ind.a   0];
auxSNPbranch.parameter.min_bound = [ind.a  -0.5];
figure(10);hold on;
[auxSNPbranch,s,f,r]   =   br_contn(auxSNPfuncs,auxSNPbranch,200);

% Here move kappa
pt                           = auxSNPbranch.point(end);
[auxSNPfuncs,auxSNPbranch]   = setupNextSNP(pt,funcs,ind.mu_tilde,ind.kappa,-1e-03,0);
auxSNPbranch.parameter.max_step  = [ind.mu_tilde  0.05;ind.tau 0.1];
auxSNPbranch.parameter.max_bound = [ind.kappa  1];
auxSNPbranch.parameter.min_bound = [ind.kappa -1];
figure(2);hold on;
auxSNPbranch    = br_contn(auxSNPfuncs,auxSNPbranch,200);

% Here move mu_tilde
pt                           = auxSNPbranch.point(end);
[auxSNPfuncs,auxSNPbranch1]   = setupNextSNP(pt,funcs,ind.mu_tilde,ind.tau,1e-03,0);
auxSNPbranch1.parameter.max_step  = [ind.mu_tilde  0.1;ind.tau 0.1];
auxSNPbranch1.parameter.max_bound = [ind.mu_tilde   0];
auxSNPbranch1.parameter.min_bound = [];
figure(3);hold on;
auxSNPbranch1    = br_contn(auxSNPfuncs,auxSNPbranch1,200);
auxSNPbranch     = auxSNPbranch1;

% Here move mu
pt                           = auxSNPbranch.point(end);
[auxSNPfuncs,auxSNPbranch1]   = setupNextSNP(pt,funcs,ind.mu,ind.tau,-1e-03,0);
auxSNPbranch1.parameter.max_step  = [ind.mu  0.01;ind.tau 0.01];
auxSNPbranch1.parameter.max_bound = [];
auxSNPbranch1.parameter.min_bound = [ind.tau -6];
figure(4);hold on;
auxSNPbranch1    = br_contn(auxSNPfuncs,auxSNPbranch1,300);
auxSNPbranch1          =   br_rvers(auxSNPbranch1);
[auxSNPbranch1,s,f,r]  =   br_contn(auxSNPfuncs,auxSNPbranch1,300);

auxSNPbranch     = auxSNPbranch1;

auxCurves{4}           =   auxSNPbranch;

% Undocument this line to save 
%save('./DDEResults/BifCurvesOrbitFlipPaper_mu_tau_closeOrbitFlip','auxCurves') 
