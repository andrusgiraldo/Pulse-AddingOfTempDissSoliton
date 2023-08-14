%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%  --- Script that generates the data for Figures   9 -10 ---  %%%%%
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
%%%%%  Slices for the orientable resonant case near Orbit bifurcation %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data for the orientable resonant bifurcation diagram, choose slices
FileName = 'BifCurvesNoNOriResPaper_a_tau_closeOrbitFlip';
load(['./DDEResults/' FileName]);
slices = [-2.4 -2.3 -2.1 -2.0 -1.9 -1.7 -1.584 -1.56 -1.48];

setDataFigure.tag         = 'NonOrientable Resonant Bifurcation Slices';
setDataFigure.paramSlices = slices;
setDataFigure.curvesSlice = cell(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the corresponding tau value of the slice and correct points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pick slice
for k=2:length(slices)
    slice = slices(k);
    auxSliceCurves = cell(1,1);

    % find the corresponding tau value of the slice
    [~,auxLoc]= find(diff(sign(arrayfun(@(x) x.parameter(ind.a),auxCurves{1}.point)-slice))~=0);
    % correct points and set up auxauxCurves
    for i=1:length(auxLoc)
        pt =  auxCurves{1}.point(auxLoc(i));
        pt.parameter(ind.tau)=pt.parameter(ind.tau)+pt.period;
        auxPer_br = SetupOneDPerBranch(pt, funcs, ind.tau, ind.a, slice, -1e-3);
        auxPer_br.parameter.max_step  = [0 1.; ind.tau  1.];
        auxPer_br.parameter.max_bound = [ind.tau 100.];
        auxPer_br.parameter.min_bound = [];
        figure(k); hold on;
        [auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,250);
        auxPer_br.point(1:5) = [];
        auxPer_br         =  br_rvers(auxPer_br);
        [auxPer_br,s,f,r] =  br_contn(funcs,auxPer_br,20);
        auxSliceCurves{i} =  auxPer_br;
    end
    setDataFigure.curvesSlice{k} = auxSliceCurves;
    disp(k)
end

% Undocument the next line to save the variable
%save('./DDEResults/SlicesResNonOri','setDataFigure')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing the period-two branch from the homoclinic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setDataFigure.Hom2curvesSlice = cell(1,length(slices));
for k=[7,8,9]%length(slices)
    slice = slices(k);
    auxSliceCurves = cell(1,1);

    % find the corresponding tau value of the slice
    [~,auxLoc]= find(diff(sign(arrayfun(@(x) x.parameter(ind.a),auxCurves{2}.point)-slice))~=0);
    % correct points and set up auxauxCurves
    for i=1:length(auxLoc)
        pt =  auxCurves{2}.point(auxLoc(i));
        pt.parameter(ind.tau)=pt.parameter(ind.tau)+pt.period;
        auxPer_br = SetupOneDPerBranch(pt, funcs, ind.tau, ind.a, slice, -1e-3);
        auxPer_br.parameter.max_step  = [0 1.; ind.tau  1.];
        auxPer_br.parameter.max_bound = [ind.tau 100.];
        auxPer_br.parameter.min_bound = [];
        figure(k); hold on;
        [auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,200);
        auxPer_br.point(1:5) = [];
        auxPer_br         =  br_rvers(auxPer_br);
        [auxPer_br,s,f,r] =  br_contn(funcs,auxPer_br,150);
        auxSliceCurves{i} =  auxPer_br;
    end
    setDataFigure.Hom2curvesSlice{k} = auxSliceCurves;
    disp(k)
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting a single slice related to the PD when there does not exist 
% a homoclinic double orbit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=7;
auxA   = slices(k);
auxSliceCurves = cell(1,1);


% find the corresponding tau value of the slice
auxTau    = -0.6825;
[~,auxLoc]= find(diff(sign(arrayfun(@(x) x.parameter(ind.tau),auxCurves{2}.point)-(auxTau)))~=0);

% correct points and set up auxauxCurves  
pt =  auxCurves{2}.point(auxLoc(1));
pt.parameter(ind.tau)=pt.parameter(ind.tau) + pt.period;
auxPer_br = SetupOneDPerBranch(pt, funcs, ind.tau, ind.a, pt.parameter(ind.a), -1e-3);
auxPer_br.parameter.max_step  = [0 1.; ind.tau   1];
auxPer_br.parameter.max_bound = [];
auxPer_br.parameter.min_bound = [ind.tau 13];
figure(2); hold on;
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,200);

pt =  auxPer_br.point(end);
auxPer_br = SetupOneDPerBranch(pt, funcs, ind.a, ind.tau, 13, -1e-04);
auxPer_br.parameter.max_step  = [0 1.; ind.a   0.05];
auxPer_br.parameter.max_bound = [];
auxPer_br.parameter.min_bound = [ind.a auxA];
figure(3); hold on;
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,30);

pt =  auxPer_br.point(end);
auxPer_br = SetupOneDPerBranch(pt, funcs, ind.tau, ind.a, pt.parameter(ind.a), -1e-04);
auxPer_br.parameter.max_step  = [0 1.; ind.tau   0.5];
auxPer_br.parameter.max_bound = [];
auxPer_br.parameter.min_bound = [];
figure(4); hold on;
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,120);
auxPer_br.point(1:5) = [];
auxPer_br         =  br_rvers(auxPer_br);
[auxPer_br,s,f,r] =  br_contn(funcs,auxPer_br,120);
auxSliceCurves{1} =  auxPer_br;

% Move to the other side of the branch to find intersections with curves
auxTau    = -0.696;
[~,auxLoc]= find(diff(sign(arrayfun(@(x) x.parameter(ind.tau),auxCurves{2}.point)-(auxTau)))~=0);

% Correct points and set up auxauxCurves  
pt =  auxCurves{2}.point(auxLoc(1));
auxPer_br = SetupOneDPerBranch(pt, funcs, ind.a, ind.tau, auxTau, -1e-4);
pt.parameter(ind.tau)=pt.parameter(ind.tau) + pt.period;
auxPer_br = SetupOneDPerBranch(pt, funcs, ind.tau, ind.a, pt.parameter(ind.a), -1e-3);
auxPer_br.parameter.max_step  = [0 1.; ind.tau   1];
auxPer_br.parameter.max_bound = [];
auxPer_br.parameter.min_bound = [ind.tau 15.2];
figure(2); hold on;
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,200);

pt =  auxPer_br.point(end);
auxPer_br = SetupOneDPerBranch(pt, funcs, ind.a, ind.tau, 15.2, -1e-04);
auxPer_br.parameter.max_step  = [0 1.; ind.a   0.05];
auxPer_br.parameter.max_bound = [];
auxPer_br.parameter.min_bound = [ind.a auxA];
figure(3); hold on;
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,30);

pt =  auxPer_br.point(end);
auxPer_br = SetupOneDPerBranch(pt, funcs, ind.tau, ind.a, pt.parameter(ind.a), -1e-04);
auxPer_br.parameter.max_step  = [0 1.; ind.tau   0.5];
auxPer_br.parameter.max_bound = [];
auxPer_br.parameter.min_bound = [];
figure(4); hold on;
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,120);
auxPer_br.point(1:5) = [];
auxPer_br         =  br_rvers(auxPer_br);
[auxPer_br,s,f,r] =  br_contn(funcs,auxPer_br,120);
auxSliceCurves{2} =  auxPer_br;

setDataFigure.Hom2curvesSlice{k} = auxSliceCurves;
disp(k)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finding the solution branches of periodic solutions for the slice of 
% a=-2.4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

auxA   = slices(1); %-2.4
auxSliceCurves = cell(1,1);

% find intersections with curves
auxTau    = -1.5;
[~,auxLoc]= find(diff(sign(arrayfun(@(x) x.parameter(ind.tau),auxCurves{1}.point)-(auxTau)))~=0);

% correct points and set up auxauxCurves  
pt =  auxCurves{1}.point(auxLoc(2));
pt.parameter(ind.tau)=pt.parameter(ind.tau) + pt.period;
auxPer_br = SetupOneDPerBranch(pt, funcs, ind.tau, ind.a, pt.parameter(ind.a), -1e-3);
auxPer_br.parameter.max_step  = [0 1.; ind.tau   1];
auxPer_br.parameter.max_bound = [];
auxPer_br.parameter.min_bound = [ind.tau 13];
figure(2); hold on;
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,200);

pt =  auxPer_br.point(end);
auxPer_br = SetupOneDPerBranch(pt, funcs, ind.a, ind.tau, 13, -1e-04);
auxPer_br.parameter.max_step  = [0 1.; ind.a   0.05];
auxPer_br.parameter.max_bound = [];
auxPer_br.parameter.min_bound = [ind.a auxA];
figure(3); hold on;
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,30);

pt =  auxPer_br.point(end);
auxPer_br = SetupOneDPerBranch(pt, funcs, ind.tau, ind.a, pt.parameter(ind.a), -1e-04);
auxPer_br.parameter.max_step  = [0 1.; ind.tau   0.5];
auxPer_br.parameter.max_bound = [ind.tau 80];
auxPer_br.parameter.min_bound = [];
figure(4); hold on;
[auxPer_br,s,f,r]=   br_contn(funcs,auxPer_br,200);
auxPer_br.point(1:5) = [];
auxPer_br         =  br_rvers(auxPer_br);
[auxPer_br,s,f,r] =  br_contn(funcs,auxPer_br,350);
auxSliceCurves{1} =  auxPer_br;
setDataFigure.curvesSlice{1} = auxSliceCurves;

% Undocument the next line to save the variable
%save('./DDEResults/SlicesResNonOri','setDataFigure')
