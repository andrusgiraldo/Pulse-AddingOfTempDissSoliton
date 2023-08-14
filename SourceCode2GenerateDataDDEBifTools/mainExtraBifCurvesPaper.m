%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%  --- Script that generates the data for Hopf bifurcation ---  %%%%%
% %%%%  ---  curve for the orientable resonant and orbit flip   ---  %%%%%
% %%%%                  ---  bifurcation curve ---                   %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

prefixFolderName  =   './AutoResonantBifurcationSolution/Matlab/BifPlaneSolution_';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Hopf curve can be picked up from slice
% Orientable Resonant case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('./DDEResults/SlicesResOri','setDataFigure')

% pick nontrivial equilibrium from PO br
pt_hopf  = setDataFigure.curvesSlice{6}{2}.point(1).profile(:,1); % Check that the four positio correspond to slice -1.7
par_hopf = setDataFigure.curvesSlice{6}{2}.point(1).parameter;
par_hopf(ind.tau) = 20;
% set up equilibrium branch
[branch1,suc]=SetupStst(funcs,'parameter',par_hopf,'x', pt_hopf,...
    'contpar',ind.tau,'step', -1e-2, 'max_step', [0, .3]); %,'min_bound',[],'max_bound',1

% continue
figure(1); clf; hold on;
branch1.method.continuation.plot=1; % don't plot prgress
[branch1,s,f,r]=br_contn(funcs,branch1,1000);
hold off;

% GetStability
branch1.method.stability.minimal_real_part=-10;
[nunst,~,~,branch1.point] =GetStability(branch1,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2);

% set up Hopf curve
[branch2,suc]=SetupHopf(funcs,branch1,indhopf(end-2),'contpar',[ind.tau, ind.a],...
'dir',ind.tau,'step',1e-2, 'max_step', [0 1. ind.a 0.01], 'min_bound', [], 'plot',1);

% continue
figure(11); hold on; 
[branch2,s,f,r]=br_contn(funcs,branch2,100);
branch2=br_rvers(branch2);
branch2.parameter.max_step = [0 0.01 ind.a 0.01];
[branch2,s,f,r]=br_contn(funcs,branch2,600);
branch2.parameter.max_step = [0 1. ind.a 0.01];
[branch2,s,f,r]=br_contn(funcs,branch2,100); 
xlabel('\tau');ylabel('\a');
hold off;

% map once to negative delay
for i=1:length(branch2.point)
    branch2.point(i).parameter(ind.tau)=branch2.point(i).parameter(ind.tau)...
                                        -2*2*pi./branch2.point(i).omega;
end

% plot Hopf curves
figure(4); clf; hold on;
for j=1:20
    tauAux = arrayfun( @(x) x.parameter(ind.tau), branch2.point);
    aAux = arrayfun( @(x) x.parameter(ind.a), branch2.point);
    omAux = arrayfun( @(x) x.omega, branch2.point);
    plot(tauAux+(j-10)*2*pi./omAux, aAux,'k','linewidth',2)
    axis([-5  5  -3  1])
end
plot(tauAux, aAux,'r','linewidth',2)
hold off;

auxauxCurves = branch2;

% Undocument the next line to save the variable
%save('./DDEResults/BifCurvesOriResPaper_a_tau_closeOrbitFlip_OneDBifD_Hopf.mat','auxauxCurves');
%save(['./DDEResults/' FileName '_OneDBifD_Hopf.mat'],'auxauxCurves');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Hopf curve can be picked up from slice 
% Orbit Flip case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('../DDEPaper/DDEResults/SlicesOrbitFlip')
% I need to find a first approximation of the Hopf point
auxSliceCurves = setDataFigure.curvesSlice{1};
dampAux = arrayfun( @(x) norm(diff(x.profile.')), auxSliceCurves{1}.point);
auxInd   = find(dampAux<0.01);

pt_hopf  = auxSliceCurves{1}.point(auxInd(4)).profile(:,1); 
par_hopf = auxSliceCurves{1}.point(auxInd(4)).parameter;


[branch1,suc]=SetupStst(funcs,'parameter',par_hopf,'x', pt_hopf,...
    'contpar',ind.tau,'step', -1e-2, 'max_step', [0, .3]);

% continue
figure(1); clf; hold on;
branch1.method.continuation.plot=1; % don't plot prgress
[branch1,s,f,r]=br_contn(funcs,branch1,100);
branch1=br_rvers(branch1);
[branch1,s,f,r]=br_contn(funcs,branch1,20);
hold off;

% GetStability
branch1.method.stability.minimal_real_part=-10;
[nunst,~,~,branch1.point] =GetStability(branch1,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2);

% set up Hopf curve
[branch2,suc]=SetupHopf(funcs,branch1,indhopf(end-2),'contpar',[ind.tau, ind.mu],...
'dir',ind.tau,'step',1e-2, 'max_step', [0 0.01 ind.mu 0.01], 'min_bound', [], 'plot',1);

% continue
figure(11); hold on; 
[branch2,s,f,r]=br_contn(funcs,branch2,200);
branch2=br_rvers(branch2);
branch2.parameter.max_step = [0 0.01 ind.mu 0.01];
[branch2,s,f,r]=br_contn(funcs,branch2,300);
%branch2.parameter.max_step = [0 1. ind.mu 0.01];
%[branch2,s,f,r]=br_contn(funcs,branch2,100); 
xlabel('\tau');ylabel('\a');
hold off;

% map once to negative delay
for i=1:length(branch2.point)
    branch2.point(i).parameter(ind.tau)=branch2.point(i).parameter(ind.tau)...
                                        -2*2*pi./branch2.point(i).omega;
end

% plot Hopf curves
figure(4); clf; hold on;
for j=1:20
    tauAux = arrayfun( @(x) x.parameter(ind.tau), branch2.point);
    muAux = arrayfun( @(x) x.parameter(ind.mu), branch2.point);
    omAux = arrayfun( @(x) x.omega, branch2.point);
    plot(tauAux+(j-10)*2*pi./omAux, muAux,'k','linewidth',2)
    axis([-5  5  -3  1])
end
plot(tauAux, muAux,'r','linewidth',2)
hold off;

auxauxCurves = branch2;

% Undocument the next line to save the variable
%save('./DDEResults/BifCurvesOrbitFlipPaper_mu_tau_closeOrbitFlip_OneDBifD_Hopf.mat','auxauxCurves');

