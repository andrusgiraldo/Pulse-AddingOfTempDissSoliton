%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%  --- Script that generates the data for Figures 11 - 13 ---  %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all;
close all;

% Adding DDE Biftools in the search path: 
% CHANGE this line to the path where your DDE Biftools installation is
addpath('../../../../../../MATLAB/dde_biftool_v3.1.1/ddebiftool/',...
        '../../../../../../MATLAB/dde_biftool_v3.1.1/demos/phase_oscillator',...
        '../../../../../../MATLAB/dde_biftool_v3.1.1/ddebiftool_extra_psol/',...
        '../../../../../../MATLAB/dde_biftool_v3.1.1/ddebiftool_utilities/',...
        '../../../../../../MATLAB/dde_biftool_v3.1.1/ddebiftool_extra_nmfm/',...
        '../../../../../../MATLAB/dde_biftool_v3.1.1/ddebiftool_extra_rotsym'); 


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
%%%%%     Computing the homoclinic solutions used in Figure 13        %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

auxCurves         =  cell(3,1);

% load the homoclinic solution using the fixed period formulation and run the continuation scheme
sufixFolderName   =  'OrbitFlipPaper/HOM_1_0.0595_0.0300';
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
auxPerHom_br.parameter.max_bound = [ind.mu  1];
auxPerHom_br.parameter.min_bound = [ind.mu  0.1];
figure(4); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,500);

% Sample the solution branch
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.mu,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.tau) = pt.parameter(ind.tau)-0.0001;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.mu,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.mu  0.01;ind.tau 0.01];
auxPerHom_br.parameter.max_bound = [ind.mu  0.1];
auxPerHom_br.parameter.min_bound = [ind.mu -0.1];
figure(4); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,500);

homOrbitBranch      = auxPerHom_br;

% Get the orbit
pt                    =  auxPerHom_br.point(end);
auxPerHom_br          =  df_brnch(funcs_per,[ind.mu,ind.tau],'psol');
auxPerHom_br.point    =  pt;
pt.parameter(ind.tau) = pt.parameter(ind.tau)+0.01;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.mu,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.mu  0.01;ind.tau 0.01];
auxPerHom_br.parameter.max_bound = [ind.mu  0.0];
auxPerHom_br.parameter.min_bound = [ind.mu -0.1];
figure(4); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,500);

homOrbitFlipBranch      = auxPerHom_br;

% Undocument the next line to save the variable
%save('./DDEResults/homOrbitBranch','homOrbitBranch','homOrbitFlipBranch')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Compute solutions as they approach tau infinity.             %%%%%
%%%%%%%%%%%%% These solutions are used in Figure 11       %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


auxCurves         =  cell(2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the homoclinic solution using the fixed period formulation and run 
% the continuation scheme
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
pt.parameter(ind.a) = pt.parameter(ind.a)-0.1;
method                     = df_mthod(funcs_per,'psol');
[pt,suc]              = p_correc(funcs_per,pt,[ind.a,ind.tau],[],method.point);
auxPerHom_br.point(2) = pt;
auxPerHom_br.parameter.max_step  = [ind.a  0.025;ind.tau 0.025];
auxPerHom_br.parameter.max_bound = [];
auxPerHom_br.parameter.min_bound = [];
figure(3); hold on;
[auxPerHom_br,s,f,r]  =   br_contn(funcs_per,auxPerHom_br,800); %600 is good


approxInfCycle = auxPerHom_br;

% Undocument the next line to save the variable
%save('./DDEResults/approxInfCycle','approxInfCycle')

