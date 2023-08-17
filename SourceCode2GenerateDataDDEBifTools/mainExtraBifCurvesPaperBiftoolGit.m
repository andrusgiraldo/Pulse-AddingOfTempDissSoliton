%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%  --- Script that generates the data for Figures   6 - 7 ---  %%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all;
close all;

% For this set of script to run, please download and use the latest commit
% of DDE BifTools from GitHub

% addpath('../../../../MATLAB/dde_biftool-git/ddebiftool/',...
%         '../../../../MATLAB/dde_biftool-git/demos/phase_oscillator',...
%         '../../../../MATLAB/dde_biftool-git/ddebiftool_extra_psol/',...
%         '../../../../MATLAB/dde_biftool-git/ddebiftool_utilities/',...
%         '../../../../MATLAB/dde_biftool-git/ddebiftool_extra_nmfm/',...
%         '../../../../MATLAB/dde_biftool-git/ddebiftool_extra_rotsym'); 

addpath('../../../Dropbox/Work/Research/MATLAB/ddebiftool-git/ddebiftool/',...
        '../../../Dropbox/Work/Research/MATLAB/ddebiftool-git/demos/phase_oscillator',...
        '../../../Dropbox/Work/Research/MATLAB/ddebiftool-git/ddebiftool_extra_psol/',...
        '../../../Dropbox/Work/Research/MATLAB/ddebiftool-git/ddebiftool_utilities/',...
        '../../../Dropbox/Work/Research/MATLAB/ddebiftool-git/ddebiftool_extra_nmfm/',...
        '../../../Dropbox/Work/Research/MATLAB/ddebiftool-git/ddebiftool_extra_rotsym');

% Adding a collection of functions that are used to initialize different
% variables for DDE-Biftool
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


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Extra PD curve for Fig 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('./DDEResults/SlicesResOri','setDataFigure')
% pick nontrivial equilibrium from Hopf branch in slice (c)
pt_hopf  = setDataFigure.curvesSlice{6}{2}.point(1).profile(:,1); % Check that the four position correspond to slice -1.7
par_hopf = setDataFigure.curvesSlice{6}{2}.point(1).parameter;

% set up equilibrium branch in tau
[branch1,suc]=SetupStst(funcs,'parameter',par_hopf,'x', pt_hopf,...
    'contpar',ind.tau,'step', 1e-3, 'max_step', [0, 10.],...
    'min_bound', [], 'max_bound', [ind.tau, 0]);

% continue equilibrium to tau=0
figure(1); clf; hold on;
branch1.method.continuation.plot=1; % don't plot prgress
[branch1,s,f,r]=br_contn(funcs,branch1,1000);
hold off;

% set up equilibrium branch in a 
[branch2,suc]=SetupStst(funcs,'parameter',branch1.point(end).parameter,...
    'x', branch1.point(end).x,'contpar',ind.a,'step', 1e-3,...
    'max_step', [0, 0.1],'min_bound', [], 'max_bound', [ind.a, -1.575]);
    
% continue equilibrium to a=-1.575
figure(2); clf; hold on;
branch2.method.continuation.plot=1; % don't plot prgress
[branch2,s,f,r]=br_contn(funcs,branch2,1000);
hold off;
    
% set up equilibrium branch in tau (locally) and branch off from first
% reappearance
[branch3,suc]=SetupStst(funcs,'parameter',branch2.point(end).parameter,...
    'x', branch2.point(end).x,'contpar',ind.tau,'step', 1e-3,...
    'max_step', [0, 0.1],'min_bound', [], 'max_bound', [ind.tau, 10]);
    
% continue equilibrium
figure(3); clf; hold on;
branch3.method.continuation.plot=1; % don't plot prgress
[branch3,s,f,r]=br_contn(funcs,branch3,1000);
hold off;    

% GetStability
branch3.method.stability.minimal_real_part=-10;
[nunst,~,~,branch3.point] =GetStability(branch3,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2);% set up Hopf curve
indhopf=indhopf(1);

% Branch off from Hopf and continue to tau=50
[perbranch1,suc]=SetupPsolFrom_stst(funcs,branch3,indhopf(1),...
    'degree', 11, 'intervals', 67, 'matrix','sparse','eigmatrix','sparse',...
    'contpar',ind.tau,'step', 1e-3,...
    'max_step', [0, 10.0], 'min_bound', [], 'max_bound', [ind.tau, 50]);

% continue Psol
figure(4); clf; hold on;
perbranch1.method.continuation.plot=1; % don't plot prgress
[perbranch1,s,f,r]=br_contn(funcs,perbranch1,20);
hold off;

perbranch2 = perbranch1;   

%set up Psol branch in mu (but map once forward)
perbranch2.point(end).parameter(ind.tau)=perbranch2.point(end).parameter(ind.tau)+perbranch2.point(end).period;
[perbranch3,suc]=SetupPsolFrom_psol(funcs,perbranch2,length(perbranch2.point),...
    'contpar',ind.mu,'step', 1e-3,'max_step', [0, 1.0],'min_bound', [], 'max_bound', [ind.mu, 0.2]);
    
% continue Psol to mu=0.2
figure(6); clf; hold on;
perbranch3.method.continuation.plot=1; % don't plot prgress
[perbranch3,s,f,r]=br_contn(funcs,perbranch3,1000);
hold off;    
    
%set up Psol branch in mu_tilde
[perbranch4,~]=SetupPsolFrom_psol(funcs,perbranch3,length(perbranch3.point),...
    'contpar',ind.mu_tilde,'step', 1e-3,'max_step', [0, 1.0],'min_bound', [], 'max_bound', [ind.mu_tilde, 0.03]);
    
% continue Psol to mu=0.2
figure(7); clf; hold on;
perbranch4.method.continuation.plot=1; % don't plot prgress
[perbranch4,s,f,r]=br_contn(funcs,perbranch4,1000);
hold off;  
    
% Now we locally look like Fig.10(h) in the first reappearance. Let us now apply
% the reappearance map once more and continue in tau, then compute
% stability and from there we can continue 
perbranch4.point(end).parameter(ind.tau)=perbranch4.point(end).parameter(ind.tau)+perbranch4.point(end).period;
[perbranch5,suc]=SetupPsolFrom_psol(funcs,perbranch4,length(perbranch4.point),...
    'contpar',ind.tau,'step', -1e-1,'compute_eigenfuncs', true,...
    'degree', 11, 'intervals', 331, 'matrix','sparse','eigmatrix','sparse',...
    'max_step', [0, .1],'min_bound', [], 'max_bound', [ind.tau, 20.]);
    
% continue Psol to mu=0.2
figure(8); clf; hold on;
perbranch5.method.continuation.plot=1; % don't plot prgress
[perbranch5,s,f,r]=br_contn(funcs,perbranch5,1000);
perbranch5 = br_rvers(perbranch5);
[perbranch5,s,f,r]=br_contn(funcs,perbranch5,1000);
hold off;  
    
%GetStability
[pernunst,~,~,perbranch5.point] =GetStability(perbranch5,'funcs',funcs,...
    'exclude_trivial',true, 'locate_trivial', @(p)1);
indpd=find(abs(diff(pernunst))==1);
indpd=indpd(end); % many points founds pick first
    
       
% Now set up the PD curve 
[pdfuncs1,pdbranch1,suc]=SetupPeriodDoubling(funcs,perbranch5,indpd,...
    'contpar', [ind.tau ind.a], 'dir', ind.a, 'step', 1e-2,...
    'max_step', [0, 0.1; ind.a, 0.1; ind.tau, 1.],'min_bound', [],...
    'max_bound', [ind.a -1.2; ind.tau 35.], 'locate_trivial', @(p)[-1,1],...
    'newton_nmon_iterations',5,'newton_max_iterations',10,...
    'matrix','sparse','eigmatrix','sparse','degree',11,'intervals',34);

% continue PD
figure(9); clf; hold on;
pdbranch1.method.continuation.plot=1; 
[pdbranch1,s,f,r]=br_contn(pdfuncs1,pdbranch1,1000);
pdbranch1 = br_rvers(pdbranch1);
[pdbranch1,s,f,r]=br_contn(pdfuncs1,pdbranch1,1000);
hold off;
    
auxPDbranch = pdbranch1;   

% Undocument the next line to save the variable
save('./DDEResults/SpecialPDBranchFig3','auxPDbranch')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Extra FOLD Curve in Fig 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('./DDEResults/SlicesResOri','setDataFigure')
% pick up equilibrium from Hopf branch in slice (b)
pt_hopf  = setDataFigure.curvesSlice{6}{2}.point(1).profile(:,1); % Check that the four position correspond to slice -1.7
par_hopf = setDataFigure.curvesSlice{6}{2}.point(1).parameter;

% set up equilibrium branch in tau
[branch1,suc]=SetupStst(funcs,'parameter',par_hopf,'x', pt_hopf,...
    'contpar',ind.tau,'step', 1e-3, 'max_step', [0, 10.],...
    'min_bound', [], 'max_bound', [ind.tau, 0.]);

% continue equilibrium to tau=0
figure(1); clf; hold on;
branch1.method.continuation.plot=1; % don't plot prgress
[branch1,s,f,r]=br_contn(funcs,branch1,1000);
hold off;   

% set up equilibrium branch in tau (locally) and compute (1st RA of) Hopf points
[branch3,suc]=SetupStst(funcs,'parameter',branch1.point(end).parameter,...
    'x', branch1.point(end).x,'contpar',ind.tau,'step', 1e-3,...
    'max_step', [0, 0.1],'min_bound', [], 'max_bound', [ind.tau, 10]);
    
% continue equilibrium
figure(3); clf; hold on;
branch3.method.continuation.plot=1; % don't plot prgress
[branch3,s,f,r]=br_contn(funcs,branch3,1000);
hold off;    

% GetStability
branch3.method.stability.minimal_real_part=-10;
[nunst,~,~,branch3.point] =GetStability(branch3,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2);% set up Hopf curve
indhopf=indhopf(2);

% Branch off from Hopf, map back once and continue
[perbranch1,suc]=SetupPsolFrom_stst(funcs,branch3,indhopf,...
    'degree', 11, 'intervals', 67, 'matrix','sparse','eigmatrix','sparse',...
    'contpar',ind.tau,'step', 1e-3,...
    'max_step', [0, 0.1], 'min_bound', [], 'max_bound', [ind.tau, 50]);

perbranch1.point(1).parameter(ind.tau)=perbranch1.point(1).parameter(ind.tau)-perbranch1.point(1).period;
perbranch1.point(2).parameter(ind.tau)=perbranch1.point(2).parameter(ind.tau)-perbranch1.point(2).period;

% continue Psol
figure(4); clf; hold on;
perbranch1.method.continuation.plot=1; % don't plot prgress
[perbranch1,s,f,r]=br_contn(funcs,perbranch1,50);
hold off;
    
% find fold and setup Fold branch
[~,indfold] = max(arrayfun( @(x) x.parameter(ind.tau), perbranch1.point));
[foldfuncs1,foldbranch1,suc]=SetupPOfold(funcs,perbranch1,indfold,...
    'contpar', [ind.tau ind.a], 'dir', ind.a, 'step', 1e-2,...
    'max_step', [0, 0.1; ind.a, 0.1; ind.tau, 0.01],...
    'min_bound',[ind.tau, -2.; ind.a, -3.],...
    'max_bound', [ind.tau, 0; ind.a 1.], 'locate_trivial', @(p)[1,1],...
    'newton_nmon_iterations',5,'newton_max_iterations',10,...
    'matrix','sparse','eigmatrix','sparse','degree',11,'intervals',34);

% continue fold
figure(9); clf; hold on;
foldbranch1.method.continuation.plot=1; % don't plot prgress
[foldbranch1,s,f,r]=br_contn(foldfuncs1,foldbranch1,100);
foldbranch1 = br_rvers(foldbranch1);
[foldbranch1,s,f,r]=br_contn(foldfuncs1,foldbranch1,300);
hold off;
    
auxFOLDbranch = foldbranch1;   

% Undocument the next line to save the variable
save('./DDEResults/SpecialFOLDBranchFig5','auxFOLDbranch')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Extra FOLD Curve in Fig 12 (emanating from B2)
% strategy is similar to before pick up equilibrium in Oriressclice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('./DDEResults/SlicesResOri','setDataFigure')
% pick nontrivial equilibrium from Hopf branch in slice (c)
pt_hopf  = setDataFigure.curvesSlice{6}{2}.point(1).profile(:,1); % Check that the four position correspond to slice -1.7
par_hopf = setDataFigure.curvesSlice{6}{2}.point(1).parameter;

% set up equilibrium branch in tau
[branch1,suc]=SetupStst(funcs,'parameter',par_hopf,'x', pt_hopf,...
    'contpar',ind.tau,'step', 1e-3, 'max_step', [0, 10.],...
    'min_bound', [], 'max_bound', [ind.tau, 0]);

% continue equilibrium to tau=0
figure(1); clf; hold on;
branch1.method.continuation.plot=1; % don't plot prgress
[branch1,s,f,r]=br_contn(funcs,branch1,1000);
hold off;

% set up equilibrium branch in a 
[branch2,suc]=SetupStst(funcs,'parameter',branch1.point(end).parameter,...
    'x', branch1.point(end).x,'contpar',ind.a,'step', 1e-3,...
    'max_step', [0, 0.1],'min_bound', [], 'max_bound', [ind.a, -0.5]);
    
% continue equilibrium to a=-0.5
figure(2); clf; hold on;
branch2.method.continuation.plot=1; % don't plot prgress
[branch2,s,f,r]=br_contn(funcs,branch2,1000);
hold off;

% set up equilibrium branch in mu
[branch3,suc]=SetupStst(funcs,'parameter',branch2.point(end).parameter,...
    'x', branch2.point(end).x,'contpar',ind.mu,'step', 1e-3,...
    'max_step', [0, 0.1],'min_bound', [], 'max_bound', [ind.mu, -0.11]);
    
% continue equilibrium to mu=-0.11
figure(3); clf; hold on;
branch3.method.continuation.plot=1; % don't plot prgress
[branch3,s,f,r]=br_contn(funcs,branch3,1000);
hold off;    

% set up equilibrium branch in mu
[branch4,suc]=SetupStst(funcs,'parameter',branch3.point(end).parameter,...
    'x', branch3.point(end).x,'contpar',ind.mu_tilde,'step', 1e-3,...
    'max_step', [0, 0.1],'min_bound', [], 'max_bound', [ind.mu_tilde, -0.0]);
    
% continue equilibrium to mu_tilde=0.0
figure(4); clf; hold on;
branch4.method.continuation.plot=1; % don't plot prgress
[branch4,s,f,r]=br_contn(funcs,branch4,1000);
hold off; 
    
% set up equilibrium branch in tau (locally) and branch off from first
% reappearance
[branch5,suc]=SetupStst(funcs,'parameter',branch4.point(end).parameter,...
    'x', branch4.point(end).x,'contpar',ind.tau,'step', 1e-3,...
    'max_step', [0, 0.1],'min_bound', [], 'max_bound', [ind.tau, 10]);
    
% continue equilibrium
figure(5); clf; hold on;
branch5.method.continuation.plot=1; % don't plot prgress
[branch5,s,f,r]=br_contn(funcs,branch5,1000);
hold off;
    
% GetStability
branch5.method.stability.minimal_real_part=-10;
[nunst,~,~,branch5.point] =GetStability(branch5,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2);% set up Hopf curve
indhopf=indhopf(1);
    
    
% Branch off from Hopf, map back once and continue
[perbranch1,suc]=SetupPsolFrom_stst(funcs,branch5,indhopf,...
    'degree', 11, 'intervals', 67, 'matrix','sparse','eigmatrix','sparse',...
    'contpar',ind.tau,'step', 1e-3,...
    'max_step', [0, 0.1], 'min_bound', [], 'max_bound', [ind.tau, 10]);

    
% continue Psol
figure(6); clf; hold on;
perbranch1.method.continuation.plot=1; % don't plot prgress
[perbranch1,s,f,r]=br_contn(funcs,perbranch1,50);
hold off;
    
% find fold and setup Fold branch
[~,indfold] = max(arrayfun( @(x) x.parameter(ind.tau), perbranch1.point));
[foldfuncs1,foldbranch1,suc]=SetupPOfold(funcs,perbranch1,indfold,...
    'contpar', [ind.tau ind.mu], 'dir', ind.a, 'step', 1e-2,...
    'max_step', [0, 0.1; ind.a, 0.1; ind.tau, 0.1],...
    'min_bound',[ind.tau, -4.; ind.mu, -1.],...
    'max_bound', [ind.tau, 4.; ind.mu 1.], 'locate_trivial', @(p)[1,1],...
    'newton_nmon_iterations',5,'newton_max_iterations',10,...
    'matrix','sparse','eigmatrix','sparse','degree',11,'intervals',34);

% continue PD
figure(7); clf; hold on;
foldbranch1.method.continuation.plot=1;
foldbranch1 = br_rvers(foldbranch1);
[foldbranch1,s,f,r]=br_contn(foldfuncs1,foldbranch1,100);
hold off;
    
auxFOLDbranch = foldbranch1; 

% Undocument the next line to save the variable
save('./DDEResults/SpecialFOLDBranchFig12','auxFOLDbranch')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Extra period-doubled FOLD Curve in Fig 10(h) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('./DDEResults/SlicesResOri','setDataFigure')
% pick nontrivial equilibrium from Hopf branch in slice (c)
pt_hopf  = setDataFigure.curvesSlice{6}{2}.point(1).profile(:,1); % Check that the four position correspond to slice -1.7
par_hopf = setDataFigure.curvesSlice{6}{2}.point(1).parameter;

% set up equilibrium branch in tau
[branch1,suc]=SetupStst(funcs,'parameter',par_hopf,'x', pt_hopf,...
    'contpar',ind.tau,'step', 1e-3, 'max_step', [0, 10.],...
    'min_bound', [], 'max_bound', [ind.tau, 0]);

% continue equilibrium to tau=0
figure(1); clf; hold on;
branch1.method.continuation.plot=1; % don't plot prgress
[branch1,s,f,r]=br_contn(funcs,branch1,1000);
hold off;

% set up equilibrium branch in a 
[branch2,suc]=SetupStst(funcs,'parameter',branch1.point(end).parameter,...
    'x', branch1.point(end).x,'contpar',ind.a,'step', 1e-3,...
    'max_step', [0, 0.1],'min_bound', [], 'max_bound', [ind.a, -1.48]);
    
% continue equilibrium to a=-1.48
figure(2); clf; hold on;
branch2.method.continuation.plot=1; % don't plot prgress
[branch2,s,f,r]=br_contn(funcs,branch2,1000);
hold off;
    
% set up equilibrium branch in tau (locally) and branch off from first
% reappearance
[branch3,suc]=SetupStst(funcs,'parameter',branch2.point(end).parameter,...
    'x', branch2.point(end).x,'contpar',ind.tau,'step', 1e-3,...
    'max_step', [0, 0.1],'min_bound', [], 'max_bound', [ind.tau, 10]);
    
% continue equilibrium
figure(3); clf; hold on;
branch3.method.continuation.plot=1; % don't plot prgress
[branch3,s,f,r]=br_contn(funcs,branch3,1000);
hold off;    

% GetStability
branch3.method.stability.minimal_real_part=-10;
[nunst,~,~,branch3.point] =GetStability(branch3,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2);% set up Hopf curve
indhopf=indhopf(1);

% Branch off from Hopf and continue to tau=50
[perbranch1,suc]=SetupPsolFrom_stst(funcs,branch3,indhopf(1),...
    'degree', 11, 'intervals', 67, 'matrix','sparse','eigmatrix','sparse',...
    'contpar',ind.tau,'step', 1e-3,...
    'max_step', [0, 10.0], 'min_bound', [], 'max_bound', [ind.tau, 50]);

% continue Psol
figure(4); clf; hold on;
perbranch1.method.continuation.plot=1; % don't plot prgress
[perbranch1,s,f,r]=br_contn(funcs,perbranch1,20);
hold off;

perbranch2 = perbranch1;   

%set up Psol branch in mu (but map once forward)
perbranch2.point(end).parameter(ind.tau)=perbranch2.point(end).parameter(ind.tau)+perbranch2.point(end).period;
[perbranch3,suc]=SetupPsolFrom_psol(funcs,perbranch2,length(perbranch2.point),...
    'contpar',ind.mu,'step', 1e-3,'max_step', [0, 1.0],'min_bound', [], 'max_bound', [ind.mu, 0.2]);
    
% continue Psol to mu=0.2
figure(6); clf; hold on;
perbranch3.method.continuation.plot=1; % don't plot prgress
[perbranch3,s,f,r]=br_contn(funcs,perbranch3,1000);
hold off;    
    
%set up Psol branch in mu_tilde
[perbranch4,suc]=SetupPsolFrom_psol(funcs,perbranch3,length(perbranch3.point),...
    'contpar',ind.mu_tilde,'step', 1e-3,'max_step', [0, 1.0],'min_bound', [], 'max_bound', [ind.mu_tilde, 0.03]);
    
% continue Psol to mu=0.2
figure(7); clf; hold on;
perbranch4.method.continuation.plot=1; % don't plot prgress
[perbranch4,s,f,r]=br_contn(funcs,perbranch4,1000);
hold off;  
    
% Now we locally look like Fig.10(h) in the first reappearance. Let us now apply
% the reappearance map once more and continue in tau, then compute
% stability and from there we can continue 
perbranch4.point(end).parameter(ind.tau)=perbranch4.point(end).parameter(ind.tau)+perbranch4.point(end).period;
[perbranch5,suc]=SetupPsolFrom_psol(funcs,perbranch4,length(perbranch4.point),...
    'contpar',ind.tau,'step', -1e-1,'compute_eigenfuncs', true,...
    'degree', 11, 'intervals', 331, 'matrix','sparse','eigmatrix','sparse',...
    'max_step', [0, 1.1],'min_bound', [], 'max_bound', [ind.tau, 30.]);
    
% continue Psol to mu=0.2
figure(8); clf; hold on;
perbranch5.method.continuation.plot=1; % don't plot prgress
[perbranch5,s,f,r]=br_contn(funcs,perbranch5,1000);
perbranch5 = br_rvers(perbranch5);
[perbranch5,s,f,r]=br_contn(funcs,perbranch5,1000);
hold off;  
    
%GetStability
[pernunst,~,~,perbranch5.point] =GetStability(perbranch5,'funcs',funcs,...
    'exclude_trivial',true, 'locate_trivial', @(p)1);
indpd=find(abs(diff(pernunst))==1);
indpd=indpd(1); % many points founds pick first
    

% Now set up period doubled branch
% it is possible to switch back to negative delay here
[dperbranch1,suc]=DoublePsol(funcs,perbranch5,indpd,...
    'contpar', ind.tau, 'dir', ind.tau, 'step', 1e-2,...
    'max_step', [0, 1.1; ind.tau, 1.],'min_bound', [],...
    'max_bound', [ind.tau 100], 'locate_trivial', @(p)[-1,1],...
    'newton_nmon_iterations',5,'newton_max_iterations',10,...
    'matrix','sparse','eigmatrix','sparse','degree',11,'intervals',34);

% continue PD
figure(9); clf; hold on;
dperbranch1.method.continuation.plot=1; % don't plot prgress
[dperbranch1,s,f,r]=br_contn(funcs,dperbranch1,1000);
hold off;
    
auxDoublebranch = dperbranch1;   

% Undocument the next line to save the variable
save('./DDEResults/SpecialDoubleBranchFig10(h)','auxDoublebranch')
  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Extra period-doubled FOLD Curve in Fig 8 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('./DDEResults/SlicesResOri','setDataFigure')
% pick nontrivial equilibrium from Hopf branch in slice (c)
pt_hopf  = setDataFigure.curvesSlice{6}{2}.point(1).profile(:,1); % Check that the four position correspond to slice -1.7
par_hopf = setDataFigure.curvesSlice{6}{2}.point(1).parameter;

% set up equilibrium branch in tau
[branch1,suc]=SetupStst(funcs,'parameter',par_hopf,'x', pt_hopf,...
    'contpar',ind.tau,'step', 1e-3, 'max_step', [0, 10.],...
    'min_bound', [], 'max_bound', [ind.tau, 0]);

% continue equilibrium to tau=0
figure(1); clf; hold on;
branch1.method.continuation.plot=1; % don't plot prgress
[branch1,s,f,r]=br_contn(funcs,branch1,1000);
hold off;

% set up equilibrium branch in a 
[branch2,suc]=SetupStst(funcs,'parameter',branch1.point(end).parameter,...
    'x', branch1.point(end).x,'contpar',ind.a,'step', 1e-3,...
    'max_step', [0, 0.1],'min_bound', [], 'max_bound', [ind.a, -1.48]);
    
% continue equilibrium to a=-1.48
figure(2); clf; hold on;
branch2.method.continuation.plot=1; % don't plot prgress
[branch2,s,f,r]=br_contn(funcs,branch2,1000);
hold off;
    
% set up equilibrium branch in tau (locally) and branch off from first
% reappearance
[branch3,suc]=SetupStst(funcs,'parameter',branch2.point(end).parameter,...
    'x', branch2.point(end).x,'contpar',ind.tau,'step', 1e-3,...
    'max_step', [0, 0.1],'min_bound', [], 'max_bound', [ind.tau, 10]);
    
% continue equilibrium
figure(3); clf; hold on;
branch3.method.continuation.plot=1; % don't plot prgress
[branch3,s,f,r]=br_contn(funcs,branch3,1000);
hold off;    

% GetStability
branch3.method.stability.minimal_real_part=-10;
[nunst,~,~,branch3.point] =GetStability(branch3,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2);% set up Hopf curve
indhopf=indhopf(1);

% Branch off from Hopf and continue to tau=50
[perbranch1,suc]=SetupPsolFrom_stst(funcs,branch3,indhopf(1),...
    'degree', 11, 'intervals', 67, 'matrix','sparse','eigmatrix','sparse',...
    'contpar',ind.tau,'step', 1e-3,...
    'max_step', [0, 10.0], 'min_bound', [], 'max_bound', [ind.tau, 50]);

% continue Psol
figure(4); clf; hold on;
perbranch1.method.continuation.plot=1; % don't plot prgress
[perbranch1,s,f,r]=br_contn(funcs,perbranch1,20);
hold off;

perbranch2 = perbranch1;   

%set up Psol branch in mu (but map once forward)
perbranch2.point(end).parameter(ind.tau)=perbranch2.point(end).parameter(ind.tau)+perbranch2.point(end).period;
[perbranch3,suc]=SetupPsolFrom_psol(funcs,perbranch2,length(perbranch2.point),...
    'contpar',ind.mu,'step', 1e-3,'max_step', [0, 1.0],'min_bound', [], 'max_bound', [ind.mu, 0.2]);
    
% continue Psol to mu=0.2
figure(6); clf; hold on;
perbranch3.method.continuation.plot=1; % don't plot prgress
[perbranch3,s,f,r]=br_contn(funcs,perbranch3,1000);
hold off;    
    
%set up Psol branch in mu_tilde
[perbranch4,suc]=SetupPsolFrom_psol(funcs,perbranch3,length(perbranch3.point),...
    'contpar',ind.mu_tilde,'step', 1e-3,'max_step', [0, 1.0],'min_bound', [], 'max_bound', [ind.mu_tilde, 0.03]);
    
% continue Psol to mu=0.2
figure(7); clf; hold on;
perbranch4.method.continuation.plot=1; % don't plot prgress
[perbranch4,s,f,r]=br_contn(funcs,perbranch4,1000);
hold off;  
    
% Now we locally look like Fig.10(h) in the first reappearance. Let us now apply
% the reappearance map once more and continue in tau, then compute
% stability and from there we can continue 
perbranch4.point(end).parameter(ind.tau)=perbranch4.point(end).parameter(ind.tau)+perbranch4.point(end).period;
[perbranch5,suc]=SetupPsolFrom_psol(funcs,perbranch4,length(perbranch4.point),...
    'contpar',ind.tau,'step', -1e-1,'compute_eigenfuncs', true,...
    'degree', 11, 'intervals', 331, 'matrix','sparse','eigmatrix','sparse',...
    'max_step', [0, 1.1],'min_bound', [], 'max_bound', [ind.tau, 20.]);
    
% continue Psol to mu=0.2
figure(8); clf; hold on;
perbranch5.method.continuation.plot=1; % don't plot prgress
[perbranch5,s,f,r]=br_contn(funcs,perbranch5,1000);
perbranch5 = br_rvers(perbranch5);
[perbranch5,s,f,r]=br_contn(funcs,perbranch5,1000);
hold off;  
    
%GetStability
[pernunst,~,~,perbranch5.point] =GetStability(perbranch5,'funcs',funcs,...
    'exclude_trivial',true, 'locate_trivial', @(p)1);
indpd=find(abs(diff(pernunst))==1);
indpd=indpd(end); % many points founds pick first

% Now set up period doubled branch
% it is possible to switch back to negative delay here
perbranch5.point(indpd).parameter(ind.tau)=perbranch5.point(indpd).parameter(ind.tau)-2*perbranch5.point(indpd).period;
[dperbranch1,suc]=DoublePsol(funcs,perbranch5,indpd,...
    'contpar', ind.tau, 'dir', ind.tau, 'step', -1e-2,...
    'max_step', [0, 1.; ind.tau, 1.],'min_bound', [],...
    'max_bound', [ind.tau 100], 'locate_trivial', @(p)[-1,1],...
    'newton_nmon_iterations',5,'newton_max_iterations',10,...
    'matrix','sparse','eigmatrix','sparse','degree',11,'intervals',34);

% continue PD
figure(9); clf; hold on;
dperbranch1.method.continuation.plot=1; % don't plot prgress
[dperbranch1,s,f,r]=br_contn(funcs,dperbranch1,100);
hold off;
    
% find fold and setup Fold branch
[~,indfold] = min(arrayfun( @(x) x.parameter(ind.tau), dperbranch1.point));
[foldfuncs1,foldbranch1,suc]=SetupPOfold(funcs,dperbranch1,indfold,...
    'contpar', [ind.tau ind.a], 'dir', ind.a, 'step', 1e-2,...
    'max_step', [0, 0.1; ind.a, 0.1; ind.tau, 0.1],...
    'min_bound',[ind.tau, -2.; ind.a, -1.7],...
    'max_bound', [ind.tau, 0.; ind.a -1.], 'locate_trivial', @(p)[1,1],...
    'newton_nmon_iterations',5,'newton_max_iterations',10,...
    'matrix','sparse','eigmatrix','sparse','degree',11,'intervals',34);

% continue PD
figure(7); clf; hold on;
foldbranch1.method.continuation.plot=1; % don't plot prgress
[foldbranch1,s,f,r]=br_contn(foldfuncs1,foldbranch1,100);
foldbranch1 = br_rvers(foldbranch1);
[foldbranch1,s,f,r]=br_contn(foldfuncs1,foldbranch1,100);
hold off;
    
auxFOLDbranch = foldbranch1;   

% Undocument the next line to save the variable
save('./DDEResults/SpecialFOLDBranchFig8','auxFOLDbranch')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Figure 7 panel (h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('./DDEResults/SlicesResOri','setDataFigure')
% pick nontrivial equilibrium from Hopf branch in slice (c)
pt_hopf  = setDataFigure.curvesSlice{6}{2}.point(1).profile(:,1); % Check that the four position correspond to slice -1.7
par_hopf = setDataFigure.curvesSlice{6}{2}.point(1).parameter;

% set up equilibrium branch in tau
[branch1,suc]=SetupStst(funcs,'parameter',par_hopf,'x', pt_hopf,...
    'contpar',ind.tau,'step', 1e-3, 'max_step', [0, 10.],...
    'min_bound', [], 'max_bound', [ind.tau, 0]);

% continue equilibrium to tau=0
figure(1); clf; hold on;
branch1.method.continuation.plot=1; % don't plot prgress
[branch1,s,f,r]=br_contn(funcs,branch1,1000);
hold off;

% set up equilibrium branch in a 
[branch2,suc]=SetupStst(funcs,'parameter',branch1.point(end).parameter,...
    'x', branch1.point(end).x,'contpar',ind.a,'step', 1e-3,...
    'max_step', [0, 0.1],'min_bound', [], 'max_bound', [ind.a, 0.1]);
    
% continue equilibrium to a=-1.48
figure(2); clf; hold on;
branch2.method.continuation.plot=1; % don't plot prgress
[branch2,s,f,r]=br_contn(funcs,branch2,1000);
hold off;
    
% set up equilibrium branch in tau (locally) and branch off from first
% reappearance
[branch3,suc]=SetupStst(funcs,'parameter',branch2.point(end).parameter,...
    'x', branch2.point(end).x,'contpar',ind.tau,'step', 1e-3,...
    'max_step', [0, 0.1],'min_bound', [], 'max_bound', [ind.tau, 10]);
    
% continue equilibrium
figure(3); clf; hold on;
branch3.method.continuation.plot=1; % don't plot prgress
[branch3,s,f,r]=br_contn(funcs,branch3,1000);
hold off;    

% GetStability
branch3.method.stability.minimal_real_part=-10;
[nunst,~,~,branch3.point] =GetStability(branch3,'funcs',funcs);
indhopf=find(abs(diff(nunst))==2);% set up Hopf curve
indhopf=indhopf(1);

% Branch off from Hopf and continue to tau=50
[perbranch1,suc]=SetupPsolFrom_stst(funcs,branch3,indhopf(1),...
    'degree', 11, 'intervals', 67, 'matrix','sparse','eigmatrix','sparse',...
    'contpar',ind.tau,'step', 1e-3,...
    'max_step', [0, 10.0], 'min_bound', [], 'max_bound', [ind.tau, 50]);

% continue Psol
figure(4); clf; hold on;
perbranch1.method.continuation.plot=1; % don't plot prgress
[perbranch1,s,f,r]=br_contn(funcs,perbranch1,20);
hold off;

perbranch2 = perbranch1;   

%set up Psol branch in a
[perbranch3,suc]=SetupPsolFrom_psol(funcs,perbranch2,length(perbranch2.point),...
    'contpar',ind.a,'step', 1e-3,'max_step', [0, 1.0],'min_bound', [], 'max_bound', [ind.a, 0.25]);
    
% continue Psol to mu=0.2
figure(6); clf; hold on;
perbranch3.method.continuation.plot=1; % don't plot prgress
[perbranch3,s,f,r]=br_contn(funcs,perbranch3,100);  
    
%set up Psol branch in tau
[perbranch4,suc]=SetupPsolFrom_psol(funcs,perbranch3,length(perbranch3.point),...
    'contpar',ind.tau,'step', 1e-3,'max_step', [0, 0.1],'min_bound', [], 'max_bound', []);
    
% continue Psol to mu=0.2
figure(7); clf; hold on;
perbranch4.method.continuation.plot=1; 
[perbranch4,s,f,r]=br_contn(funcs,perbranch4,1000);
perbranch4 = br_rvers(perbranch4);
[perbranch4,s,f,r]=br_contn(funcs,perbranch4,1000);
hold off; 
     
auxPerbranch = perbranch4;
% Undocument the next line to save the variable
save('./DDEResults/SpecialPerBranchFig7(h)','auxPerbranch')