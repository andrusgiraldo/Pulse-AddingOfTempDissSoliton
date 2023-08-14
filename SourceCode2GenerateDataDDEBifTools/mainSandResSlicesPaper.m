%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%  --- Script that generates the data for Figures   6 - 7 ---  %%%%%
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

prefixFolderName  =   '../AutoResonantBifurcationSolution/Matlab/BifPlaneSolution_';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  Slices for the orientable resonant case near Orbit bifurcation %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data for the orientable resonant bifurcation diagram, choose slices
FileName = 'BifCurvesOriResPaper_a_tau_closeOrbitFlip';
load(['./DDEResults/' FileName]);
load(['./DDEResults/' FileName '_OneDBifD_Hopf']);
hopfCurves = cell(1,1); % Hopf curve in second file called auxauxCurves
hopfCurves{1} = auxauxCurves; % Hopf curve in second file called auxauxCurves

slices = [0.1 -0.08 -0.35 -0.38 -1.1 -1.7 -2.0];

setDataFigure.tag         = 'Orientable Resonant Bifurcation Slices';
setDataFigure.paramSlices = slices;
setDataFigure.curvesSlice = cell(1,length(slices));
setDataFigure.hopfCurvesSlice = cell(1,length(slices));
load('./DDEResults/SlicesResOri','setDataFigure');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find intersections of curves with slices and correct points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=[2,3,4,5,6]%length(slices)
    % pick slice
    slice = slices(k);
    auxSliceCurves = cell(1,1);
    
    % find the corresponding tau value of the slice
    [~,auxLoc]= find(diff(sign(arrayfun(@(x) x.parameter(ind.a),auxCurves{1}.point)-slice))~=0);
    % correct points and set up auxauxCurves
    for i=1:2 % aucLoc may have repetition. Nevertheless, we only need that the first two points which are distinct. More is repetition (the homoclinic curves is a circle)
        pt =  auxCurves{1}.point(auxLoc(i));
        pt.parameter(ind.a)= slice;
        pt =  p_correc(funcs, pt, ind.tau, [], auxCurves{1}.method.point);
        pt.parameter(ind.tau)=pt.parameter(ind.tau)+pt.period;
        auxPer_br = SetupOneDPerBranch(pt, funcs, ind.tau, ind.a, slice, -1e-3);
        auxPer_br.parameter.max_step  = [0 1.1; ind.tau  10.];
        auxPer_br.parameter.max_bound = [ind.tau 100.];
        auxPer_br.parameter.min_bound = [ind.tau -100.];
        figure(1); hold on;
        [auxPer_br,s,f,r] =   br_contn(funcs,auxPer_br,200);
        auxPer_br.point(1:5)  = [];
        auxPer_br         =   br_rvers(auxPer_br);
        [auxPer_br,s,f,r] =   br_contn(funcs,auxPer_br,50);
        % Save as auxauxCurves{i}
        auxSliceCurves{i} = auxPer_br;
    end
    setDataFigure.curvesSlice{k} = auxSliceCurves;
    disp(k)
end

% Undocument the next line to save the variable
%save('./DDEResults/SlicesResOri','setDataFigure')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the solution branches of periodic solutions that emanating from
% the Hopf bifurcation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for k=1:length(slices)
    % pick slice
    slice = slices(k);
    auxSliceHopf = cell(1,1);

    % find the corresponding tau value of the slice
    [~,auxLoc]= find(diff(sign(arrayfun(@(x) x.parameter(ind.a),hopfCurves{1}.point)-slice))~=0);
    %correct points and set up auxauxCurves
    for i=1:length(auxLoc) % aucLoc may has repetitions 
        pt =  hopfCurves{1}.point(auxLoc(i));
        pt.parameter(ind.a)= slice;
        pt =  p_correc(funcs, pt, ind.tau, [], hopfCurves{1}.method.point);
        hopfCurves{1}.point(auxLoc(i)) = pt;
        auxPer_br = SetupPsol(funcs, hopfCurves{1}, auxLoc(i), ...
                    'print_residual_info',1,'intervals',31,'degree',3,...
                    'free_par', ind.tau);
        auxPer_br.parameter.free = ind.tau;
        auxPer_br.parameter.max_bound = [ind.tau 100.];
        auxPer_br.parameter.min_bound = [ind.tau -100.];
        figure(1); clf; hold on;
        auxPer_br.parameter.max_step  = [0 1.; ind.tau  10.];
        [auxPer_br,s,f,r] =   br_contn(funcs,auxPer_br,100);
        %Save as auxauxCurves{i}
        auxSliceHopf{i} = auxPer_br;
    end
    setDataFigure.hopfCurvesSlice{k} = auxSliceHopf;
end

% Undocument the next line to save the variable
%save('./DDEResults/SlicesResOri','setDataFigure')


%% Plotting 
p = 7;

figure(5+p); clf; hold on;
auxauxCurves = setDataFigure.curvesSlice{p};
if ~isempty(auxauxCurves)
    for i=1:length(auxauxCurves)
        if ~isempty(auxauxCurves{i})
            tauAux = arrayfun( @(x) x.parameter(ind.tau), auxauxCurves{i}.point);
            perAux = arrayfun( @(x) x.period, auxauxCurves{i}.point);
            dampAux = arrayfun( @(x) norm(diff(x.profile.')), auxauxCurves{i}.point);
            for j=1:20
                plot(tauAux(dampAux>0.01)+(j-3)*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2)
            end
        end
    end
end
auxauxCurves = setDataFigure.hopfCurvesSlice{p};
for i=1:length(auxauxCurves)
    tauAux = arrayfun( @(x) x.parameter(ind.tau), auxauxCurves{i}.point);
    perAux = arrayfun( @(x) x.period, auxauxCurves{i}.point);
    dampAux = arrayfun( @(x) norm(diff(x.profile.')), auxauxCurves{i}.point);
    for j=1:20
        plot(tauAux(dampAux>0.01)+(j-3)*perAux(dampAux>0.01),perAux(dampAux>0.01),'-b','linewidth',2)
        
    end
end
axis([-4.5  20  0  30])
hold off;

%save(['./DDEResults/' FileName '_OneDBifD_' num2str(setParFix) '.mat'],'auxauxCurves');
% % manipulate P_2 branch if it turns around at the PD point
% %[~,loc] = min(arrayfun( @(x) x.parameter(indPar), auxauxCurves{2}.point));
% %auxauxCurves{2}.point(loc+1:end)=[];