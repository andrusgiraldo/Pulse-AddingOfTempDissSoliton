%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%  --- Script that generates the data for Figures   2 ---  %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  Compute TDS with tau=100 near NonOriRes %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computation of the profiles of two temporal dissipative solitons at parameter
% values tau=100 and a= -1.56 corresponding to POs close to a 1-hom and
% 2-hom in Sandstedes model with delay

% Load data precomputed in mainSandBifDiagPaper.m 
FileName = 'BifCurvesNonOriResPaper_a_tau_closeOrbitFlip';
load(['./DDEResults/' FileName]);
slice = -1.56;

auxPoints = cell(1,1);

setDataFigure.tag         = 'Points tau=100 from NonOrientable Resonant Bifurcation';
setDataFigure.paramPoints = slice;
setDataFigure.Points      = cell(1,1);
setDataFigure.Curves      = cell(1,1);
setDataFigure.PD          = cell(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find intersections of curves with slices and correct points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Continuation of the periodic orbit in tau for Figure 2.

% find point along continuation branch auxCurves{1} close to a=-1.56. 
% auxCurves{1} has the solution branch corresponding to the periodic orbit
% that goes homoclinic.
[~,auxLoc]= find(diff(sign(arrayfun(@(x) x.parameter(ind.a),auxCurves{1}.point)-slice))~=0);
% correct point and set up auxauxCurves
pt =  auxCurves{1}.point(auxLoc(1));
pt.parameter(ind.tau)=pt.parameter(ind.tau)+pt.period;
auxPer_br = SetupOneDPerBranch(pt, funcs, ind.tau, ind.a, slice, -1e-3);
auxPer_br.parameter.max_step  = [0 1.; ind.tau  1.];
auxPer_br.parameter.max_bound = [ind.tau 100.];
auxPer_br.parameter.min_bound = [];
figure(10); hold on;
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,1000);
auxPer_br = br_rvers(auxPer_br);
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,250);  
setDataFigure.Points{1} = auxPer_br.point(end);
setDataFigure.Curves{1} = auxPer_br;



% find point along continuation branch auxCurves{1} close to a=-1.56. 
% auxCurves{1} has the solution branch corresponding to the periodic orbit
% that goes homoclinic.
[~,auxLoc]= find(diff(sign(arrayfun(@(x) x.parameter(ind.a),auxCurves{1}.point)-slice))~=0);
% Here, the computation of the solution with 2-pulse-per-delay-interval orbit  
% with tau=100 is computed.
pt =  auxCurves{1}.point(auxLoc(1));
pt.parameter(ind.tau)=pt.parameter(ind.tau)+2*pt.period;
auxPer_br = SetupOneDPerBranch(pt, funcs, ind.tau, ind.a, slice, -1e-3);
auxPer_br.parameter.max_step  = [0 1.; ind.tau  1.];
auxPer_br.parameter.max_bound = [];
auxPer_br.parameter.min_bound = [ind.tau 100.];
figure(10); hold on;
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,1000);
setDataFigure.Points{4} = auxPer_br.point(end);


% find points along continuation branch auxCurves{2} close to a=-1.56. 
% auxCurves{2} has the solution branch corresponding to the periodic orbit
% that goes two-homoclinic and comes from the period-doubling bifurcation
[~,auxLoc]= find(diff(sign(arrayfun(@(x) x.parameter(ind.a),auxCurves{2}.point)-slice))~=0);

% correct points and set up auxauxCurves for the first PD bifurcation
% branch
pt =  auxCurves{2}.point(auxLoc(1));
pt.parameter(ind.tau)=pt.parameter(ind.tau)+pt.period;
auxPer_br = SetupOneDPerBranch(pt, funcs, ind.tau, ind.a, slice, -1e-3);
auxPer_br.parameter.max_step  = [0 1.; ind.tau  1.];
auxPer_br.parameter.max_bound = [ind.tau 200.];
auxPer_br.parameter.min_bound = [];
figure(10); hold on;
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,1000);
auxPer_br = br_rvers(auxPer_br);
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,1000);  
setDataFigure.Points{2} = auxPer_br.point(end);
setDataFigure.Curves{2} = auxPer_br;

% correct points and set up auxauxCurves for the second PD bifurcation
% branch
pt =  auxCurves{2}.point(auxLoc(2));
pt.parameter(ind.tau)=pt.parameter(ind.tau)+pt.period;
auxPer_br = SetupOneDPerBranch(pt, funcs, ind.tau, ind.a, slice, -1e-3);
auxPer_br.parameter.max_step  = [0 1.; ind.tau  1.];
auxPer_br.parameter.max_bound = [ind.tau 200.];
auxPer_br.parameter.min_bound = [];
figure(10); hold on;
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,1000);
auxPer_br = br_rvers(auxPer_br);
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,1000);  
setDataFigure.Points{3} = auxPer_br.point(end);
setDataFigure.Curves{3} = auxPer_br;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Period-doubling bifurcation curve in parameter plane %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find point close a=-1.56;
[~,auxLoc]= find(diff(sign(arrayfun(@(x) x.parameter(ind.a),auxCurves{3}.point)-slice))~=0);

% correct points and set up auxauxCurves
pt =  auxCurves{3}.point(auxLoc(1));
[auxPDfuncs1,auxPDbranch1]   = setupNextPD(pt,funcs,ind.tau,ind.a,1e-01,2); % you need the period double to reappear
auxPDbranch1.parameter.max_step  = [0 1., ind.a  0.025, ind.tau 10.];
auxPDbranch1.parameter.max_bound = [ind.tau, 200, ind.a, -1];
auxPDbranch1.parameter.min_bound = [ind.a, -2];
figure(3);hold on;
auxPDbranch1 = br_contn(auxPDfuncs1,auxPDbranch1,100);
auxPDbranch1 = br_rvers(auxPDbranch1);
auxPDbranch1 = br_contn(auxPDfuncs1,auxPDbranch1,100);
setDataFigure.PD = auxPDbranch1;

% Undocument the next line to save the variable
%save('./DDEResults/PointsResNonOri','setDataFigure')

