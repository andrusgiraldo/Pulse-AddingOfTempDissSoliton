%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%  --- Script to create the differentbifurcation slices   ---  %%%%%%%
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

prefixFolderName  =   './AutoResonantBifurcationSolution/Matlab/BifPlaneSolution_';

            
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of the bifurcation slices for Figures 6 and 7  
%  Orientable resonant delay case 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
load('../SourceCode2GenerateDataDDEBifTools/DDEResults/SlicesResOri')

ticklengthUn = 0.15;
labelSize    = 27;
ratioFigure  = 6.9/6;          % ratio of this figure
colorAux     =[0.5*0.5,0.5*0.5,0.8;...
               1 0 1; ...
               1 0 1;...
               0.6    0.5    0.75]; 
colorAuxRep  =[0.64*0.8,0.808*0.8,0.922*0.8]; 

colorAuxHopf =[1 0 1;...
               0.5*0.5,0.5*0.5,0.8;...
               0.5*0.5,0.5*0.5,0.8;...
               0.6    0.5    0.75]; 

% Definition of the colors used for the solution branches
colorAuxCur    = cell(1,7);
colorAuxCur{1} = [0.5*0.5,0.5*0.5,0.8; 0.6    0.5    0.75]; 
colorAuxCur{2} = [0.6    0.5    0.75; 1 0 1]; 
colorAuxCur{3} = [0.6    0.5    0.75; 1 0 1];
colorAuxCur{4} = [0.6    0.5    0.75; 0.1    0.6    0.1]; 
colorAuxCur{5} = [0.1    0.6    0.1;0.1    0.6    0.1];
colorAuxCur{6} = [0.5*0.5,0.5*0.5 0.8;0 0.8 0.8]; 
colorAuxCur{7} = []; 

colorAuxCurHB    = cell(1,7);
colorAuxCurHB{1} = [0.5*0.5,0.5*0.5,0.8; 0.6    0.5    0.75]; 
colorAuxCurHB{2} = [0.5*0.5,0.5*0.5,0.8; 0.6    0.5    0.75; 1 0 1]; 
colorAuxCurHB{3} = [0.5*0.5,0.5*0.5 0.8; 0.5*0.5,0.5*0.5 0.8; 0.6    0.5    0.75; 0.6    0.5    0.75;1 0 1]; 
colorAuxCurHB{4} = [0.5*0.5,0.5*0.5 0.8; 0.5*0.5,0.5*0.5 0.8; 0 0.8 0.8; 0 0.8 0.8]; 
colorAuxCurHB{5} = [0.5*0.5,0.5*0.5 0.8; 0 0.8 0.8; 0 0.8 0.8]; 
colorAuxCurHB{6} = []; 
colorAuxCurHB{7} = [0.5*0.5,0.5*0.5 0.8; 0.5*0.5,0.5*0.5 0.8]; 

for k=1:length(setDataFigure.paramSlices)
    auxSliceCurves = setDataFigure.curvesSlice{k};
    if ~(isempty(auxSliceCurves))
        figure(10+k); clf; hold on;
        for i=1:length(auxSliceCurves)
            for j=1:15
                tauAux = arrayfun( @(x) x.parameter(ind.tau), auxSliceCurves{i}.point);
                perAux = arrayfun( @(x) x.period, auxSliceCurves{i}.point);
                dampAux = arrayfun( @(x) norm(diff(x.profile.')), auxSliceCurves{i}.point);
                auxInd = dampAux>0.01;
                if ((j-5)==-1)
                    auxColor = colorAuxCur{k};
                    if ~(isempty(colorAuxCur{k}))
                        plot(tauAux(auxInd)+(j-5)*perAux(auxInd),perAux(auxInd),'k','linewidth',2,'Color',auxColor(i,:))
                    else
                        plot(tauAux(auxInd)+(j-5)*perAux(auxInd),perAux(auxInd),'k','linewidth',2,'Color',colorAux(i,:))
                    end
                else
                    plot(tauAux(auxInd)+(j-5)*perAux(auxInd),perAux(auxInd),'k','linewidth',2,'Color',colorAuxRep)
                end
                scatter(tauAux(auxInd(1))+(j-5)*perAux(auxInd(1)),perAux(auxInd(1)),30,'s','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',1.5)
                scatter(tauAux(auxInd(end))+(j-5)*perAux(auxInd(end)),perAux(auxInd(end)),30,'s','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',1.5)
                if (k == 6) 
                    auxPunto1 = -2.11972;
                    auxPunto2 =  2.94881;
                    scatter(auxPunto1+(j-5)*auxPunto2,auxPunto2,30,'s','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',1.5) 
                    auxPunto1 = -1.28545;
                    auxPunto2 =  3.78136;
                    scatter(auxPunto1+(j-5)*auxPunto2,auxPunto2,30,'s','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',1.5) 
                end
            end
        end
    end

    auxSliceCurves = setDataFigure.hopfCurvesSlice{k};
    if ~(isempty(auxSliceCurves))
        figure(10+k); hold on;
        for i=1:length(auxSliceCurves)
            tauAux = arrayfun( @(x) x.parameter(ind.tau), auxSliceCurves{i}.point);
            perAux = arrayfun( @(x) x.period, auxSliceCurves{i}.point);
            dampAux = arrayfun( @(x) norm(diff(x.profile.')), auxSliceCurves{i}.point);
            for j=1:15
                auxInd = dampAux>0.01;
                if ((j-5)==0)
                    auxColor = colorAuxCurHB{k};
                    if ~(isempty(colorAuxCurHB{k}))
                        plot(tauAux(auxInd)+(j-5)*perAux(auxInd),perAux(auxInd),'k','linewidth',2,'Color',auxColor(i,:))
                    else
                        plot(tauAux(auxInd)+(j-5)*perAux(auxInd),perAux(auxInd),'k','linewidth',2,'Color',colorAuxHopf(i,:))
                    end
                else
                    plot(tauAux(auxInd)+(j-5)*perAux(auxInd),perAux(auxInd),'k','linewidth',2,'Color',colorAuxRep)
                end
                scatter(tauAux(auxInd(1))+(j-5)*perAux(auxInd(1)),perAux(auxInd(1)),30,'s','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',1.5)
                scatter(tauAux(auxInd(end))+(j-5)*perAux(auxInd(end)),perAux(auxInd(end)),30,'s','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',1.5) 
                if (k==7)
                    auxPunto1 = -2.0998;
                    auxPunto2 =  2.89632;
                    scatter(auxPunto1+(j-5)*auxPunto2,auxPunto2,30,'s','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',1.5) % Generalised Hopf bifurcation
                end
            end
        end
    end 
    if (k==3 || k==4 || k==5)
        patch([-1.8336 -0.8463 -0.8463 -1.8336],[2.7502 2.7502 7.7726 7.7726],[0.2,0.2,0.2],'FaceAlpha',0.15,'LineStyle','none')
    end
    box on
    axis([-9.5  25  0  30])
    xAuxTicks = [-10, 0 10 20];
    yAuxTicks = linspace( 0,   30,4);
    set(gcf,'Color',[1 1 1]);
    set(gcf,'units','centimeters','pos', [5,20,6.9*3,6*3])
    set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
    xticklabels({'','','',''})
    yticklabels({'','','',''})
    %hgexport(gcf, ['./Figures/SlicesResOri_' sprintf('%.4f',setDataFigure.paramSlices(k)) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of the bifurcation slices for 7(h)
%  Orientable resonant delay case 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('../SourceCode2GenerateDataDDEBifTools/DDEResults/SpecialPerBranchFig7(h)')
ticklengthUn = 0.15;
labelSize    = 27;
ratioFigure  = 6.9/6; % ration of this figure
colorAux     =[0.5*0.5,0.5*0.5,0.8;...
               1 0 1; ...
               1 0 1;...
               0.6    0.5    0.75]; 
colorAuxRep  =[0.64*0.8,0.808*0.8,0.922*0.8]; 

colorAuxHopf =[1 0 1;...
               0.5*0.5,0.5*0.5,0.8;...
               0.5*0.5,0.5*0.5,0.8;...
               0.6    0.5    0.75]; 


colorAuxCurHB    = cell(1,7);
colorAuxCurHB{1} = [0.5*0.5,0.5*0.5,0.8; 0.6    0.5    0.75]; % Colors of the first family

figure(10); clf; hold on;

auxSliceCurves = auxPerbranch;
if ~(isempty(auxSliceCurves))
    figure(10); hold on;
    tauAux = arrayfun( @(x) x.parameter(ind.tau), auxSliceCurves.point);
    perAux = arrayfun( @(x) x.period, auxSliceCurves.point);
    dampAux = arrayfun( @(x) norm(diff(x.profile.')), auxSliceCurves.point);
    for j=1:15
        auxInd = dampAux>0.01;
        if ((j-4)==0)
            auxColor = colorAuxCurHB{1};
            if ~(isempty(colorAuxCurHB{1}))
                plot(tauAux(auxInd)+(j-5)*perAux(auxInd),perAux(auxInd),'k','linewidth',2,'Color',auxColor(1,:))
            else
                plot(tauAux(auxInd)+(j-5)*perAux(auxInd),perAux(auxInd),'k','linewidth',2,'Color',colorAuxHopf(1,:))
            end
        else
            plot(tauAux(auxInd)+(j-5)*perAux(auxInd),perAux(auxInd),'k','linewidth',2,'Color',colorAuxRep)
        end
        scatter(tauAux(auxInd(1))+(j-5)*perAux(auxInd(1)),perAux(auxInd(1)),30,'s','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',1.5)
        scatter(tauAux(auxInd(end))+(j-5)*perAux(auxInd(end)),perAux(auxInd(end)),30,'s','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',1.5) 
    end
end 

box on
axis([-9.5  25  0  30])
xAuxTicks = [-10, 0 10 20];
yAuxTicks = linspace( 0,   30,4);
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,6.9*3,6*3])
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/SlicesResOri_' sprintf('%.4f',auxSliceCurves.point(1).parameter(ind.a)) '.eps'], hgexport('factorystyle'), 'Format', 'eps');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of the bifurcation slices for Figures 9 and 10  
%  NonOrientable resonant delay case 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
load('../SourceCode2GenerateDataDDEBifTools/DDEResults/SlicesResNonOri')

ticklengthUn = 0.15;
labelSize    = 27;
ratioFigure  = 6.9/6; 
colorAux     =[0.5*0.5,0.5*0.5,0.8;...
               1 0 1; ...
               1 0 1;...
               0.6    0.5    0.75]; 
colorAuxRep  =[0.64*0.8,0.808*0.8,0.922*0.8]; 
for k=1:length(setDataFigure.paramSlices)
    auxSliceCurves = setDataFigure.curvesSlice{k};
    if ~(isempty(auxSliceCurves))
        figure(10+k); clf; hold on;
        for i=1:length(auxSliceCurves)
            for j=1:15
                tauAux = arrayfun( @(x) x.parameter(ind.tau), auxSliceCurves{i}.point);
                perAux = arrayfun( @(x) x.period, auxSliceCurves{i}.point);
                dampAux = arrayfun( @(x) norm(diff(x.profile.')), auxSliceCurves{i}.point);
                if ((j-5)==-1)
                    plot(tauAux(dampAux>0.01)+(j-5)*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2,'Color',colorAux(i,:))
                else
                    plot(tauAux(dampAux>0.01)+(j-5)*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2,'Color',colorAuxRep)
                end
            end
        end
    end
    
    auxSliceCurves = setDataFigure.Hom2curvesSlice{k};
    if ~(isempty(auxSliceCurves))
        figure(10+k); hold on;
        for i=1:length(auxSliceCurves)
            for j=1:15
                tauAux = arrayfun( @(x) x.parameter(ind.tau), auxSliceCurves{i}.point);
                perAux = arrayfun( @(x) x.period, auxSliceCurves{i}.point);
                dampAux = arrayfun( @(x) norm( interp1(x.mesh,x.profile',0:0.05:0.95) - interp1(x.mesh,x.profile',[0.5:0.05:0.95,0:0.05:0.45]) ), auxSliceCurves{i}.point);
                plot(tauAux(dampAux>0.01)+(j-5)*perAux(dampAux>0.01),0.5*perAux(dampAux>0.01),'k','linewidth',2,'Color',[0.7,0,0])
            end
        end
    end

    if(k==9)
        load('../SourceCode2GenerateDataDDEBifTools/DDEResults/SpecialDoubleBranchFig10(h)')
        tauAux = arrayfun( @(x) x.parameter(ind.tau), auxDoublebranch.point);
        perAux = arrayfun( @(x) x.period, auxDoublebranch.point);
        dampAux = arrayfun( @(x) norm( interp1(x.mesh,x.profile',0:0.05:0.95) - interp1(x.mesh,x.profile',[0.5:0.05:0.95,0:0.05:0.45]) ), auxDoublebranch.point);
        for j=1:15
            plot(tauAux(dampAux>0.01)+(j-5)*perAux(dampAux>0.01),0.5*perAux(dampAux>0.01),'k','linewidth',2,'Color',[0.7,0,0])
        end
    end
    if(k==7 || k==8 || k==9)
        patch([-0.9 -0.25 -0.25 -0.9],[4.2236 4.2236 12.0645 12.0645],[0.2,0.2,0.2],'FaceAlpha',0.2,'LineStyle','none')
    end
    box on
    axis([-16  30  0  40])
    xAuxTicks = [-10, 0 10 20 30];
    yAuxTicks = linspace( 0,   40,5);
    set(gcf,'Color',[1 1 1]);
    set(gcf,'units','centimeters','pos', [5,20,6.9*3,6*3])
    set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
    xticklabels({'','','',''})
    yticklabels({'','','',''})
    hgexport(gcf, ['./Figures/SlicesResNonOri_' sprintf('%.4f',setDataFigure.paramSlices(k)) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of the bifurcation slices for Figures 14 and 15  
%  Orbit flip delay case 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
load('../DDEPaper/DDEResults/SlicesOrbitFlip')

% Matrix to store the Hopf bifurcation for slices k=1,...,6
auxPuntoHB = [-1.57303, 2.32522, -1.03604, 2.9113;...
            -1.58199, 2.30682, -1.00718, 2.92659;...
            -1.5903,  2.30365, -0.992089, 2.9301;...
            -1.65182, 2.34723, -0.964932, 3.07916;...
            -1.71788, 2.30529, -0.811336, 3.0831;...
            -0.497408, 3.02666, 0.387704, 2.29777];

% Matrix to store the patches!
auxPatchesX = zeros(10,4);
auxPatchesX(2,:) = [1.0915 1.5774 1.5774 1.0915];
auxPatchesX(3,:) = [-1.4341 -1.05 -1.05 -1.4341];
auxPatchesX(4,:) = [-3.2613 3.2882 3.2882 -3.2613];
auxPatchesX(6,:) = [0.40 0.80 0.80 0.40]; 
auxPatchesX(7,:) = [-0.6 0.6 0.6 -0.6];
auxPatchesX(8,:) = [-0.3 0.35 0.35 -0.3];

auxPatchesY = zeros(10,4);
auxPatchesY(2,:) = [4.7711 4.7711 11.5756 11.5756];
auxPatchesY(3,:) = [4.7711 4.7711 11.5756 11.5756];
auxPatchesY(4,:) = [0.9125 0.9125 14.9732 14.9732];
auxPatchesY(6,:) = [2.5 2.5 8 8];
auxPatchesY(7,:) = [2.5 2.5 8 8];
auxPatchesY(8,:) = [3 3 8 8];

colorAuxCur    = cell(1,9);
colorAuxCur{1} = [1 0 1]; 
colorAuxCur{2} = [1 0 1;0.6    0.5    0.75]; 
colorAuxCur{3} = [0.0    0.8    0.8; 0.0    0.8    0.8;1 0 1;0.6    0.5    0.75]; 
colorAuxCur{4} = [0.0    0.8    0.8; 0.0    0.8    0.8;1 0 1;0.6    0.5    0.75];
colorAuxCur{5} = [0.0    0.8    0.8; 0.0    0.8    0.8;1 0 1;0.6    0.5    0.75]; 
colorAuxCur{6} = [0.0    0.8    0.8;1 0 1]; 
colorAuxCur{7} = [0.5*0.5,0.5*0.5,0.8];
colorAuxCur{8} = [0.5*0.5,0.5*0.5,0.8]; 
colorAuxCur{9} = [0.5*0.5,0.5*0.5,0.8];

ticklengthUn = 0.15;
labelSize    = 27;
ratioFigure  = 6.9/6; % ration of this figure
colorAux     =[0.5*0.5,0.5*0.5,0.8;...
               0.5*0.5,0.5*0.5,0.8; ...
               1 0 1;...
               0.6    0.5    0.75;...
               0.0    0.8    0.8]; 
colorAuxRep  =[0.64*0.8,0.808*0.8,0.922*0.8]; 

for k=1:length(setDataFigure.paramSlices)
    auxSliceCurves = setDataFigure.curvesSlice{k};
    if ~(isempty(auxSliceCurves))
        figure(k); clf; hold on;
        for i=1:length(auxSliceCurves)
            for j=1:20
                tauAux    = arrayfun( @(x) x.parameter(ind.tau), auxSliceCurves{i}.point);
                perAux    = arrayfun( @(x) x.period, auxSliceCurves{i}.point);
                dampAux   = arrayfun( @(x) norm(diff(x.profile.')), auxSliceCurves{i}.point);
                splitAux1 = find(diff(tauAux)<=0);
                splitAux2 = find(diff(tauAux)>=0);
                if ((j-8)==-1 && i~=1)
                    switch k
                        case 2
                            plot(tauAux(dampAux>0.01)+(j-8)*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2,'Color',colorAuxCur{k}(i,:))
                        case 3
                            plot(tauAux(dampAux>0.01)+(j-8)*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2,'Color',colorAuxCur{k}(i,:))
                        case 4
                            plot(tauAux(dampAux>0.01)+(j-8)*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2,'Color',colorAuxCur{k}(i,:))
                        case 5
                            plot(tauAux(dampAux>0.01)+(j-8)*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2,'Color',colorAuxCur{k}(i,:))
                        case 6
                            plot(tauAux(dampAux>0.01)+(j-8)*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2,'Color',colorAuxCur{k}(i,:))
                        case 7
                            disp('I am not allowed here')
                        case 8
                            disp('I am not allowed here')
                        case 9
                            disp('I am not allowed here')
                        otherwise
                            plot(tauAux(dampAux>0.01)+(j-8)*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2,'Color',colorAux(3,:))
                    end
                elseif ((j-8)==-1)
                    switch k
                        case 1
                            plot(tauAux(dampAux>0.01)+(j-8)*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2,'Color',colorAuxCur{k}(i,:))
                        case 2 
                            plot(tauAux(dampAux>0.01)+(j-8)*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2,'Color',colorAuxCur{k}(i,:))
                        case 6
                            plot(tauAux(dampAux>0.01)+(j-8)*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2,'Color',colorAuxCur{k}(i,:))
                        case 7
                            plot(tauAux(dampAux>0.01)+(j-8)*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2,'Color',colorAuxCur{k}(1,:))
                        case 8
                            plot(tauAux(dampAux>0.01)+(j-8)*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2,'Color',colorAuxCur{k}(1,:))
                        case 9
                            plot(tauAux(dampAux>0.01)+(j-8)*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2,'Color',colorAuxCur{k}(1,:))
                        otherwise
                            auxInd = dampAux>0.01;
                            plot(tauAux(auxInd)+(j-8)*perAux(auxInd),perAux(auxInd),'k','linewidth',2,'Color',colorAuxRep)
                    end
                else
                    auxInd = dampAux>0.01;
                    plot(tauAux(auxInd)+(j-8)*perAux(auxInd),perAux(auxInd),'k','linewidth',2,'Color',colorAuxRep)
                end
            end
        end
    end

    % This is for figure 4
    if (k==4)
        tauAux    = arrayfun( @(x) x.parameter(ind.tau), auxSliceCurves{4}.point);
        perAux    = arrayfun( @(x) x.period, auxSliceCurves{4}.point);
        dampAux   = arrayfun( @(x) norm(diff(x.profile.')), auxSliceCurves{4}.point);
        splitAux1 = find(diff(tauAux)<=0);
        splitAux2 = find(diff(tauAux)>=0);
        plot(tauAux(dampAux>0.01)-2*perAux(dampAux>0.01),perAux(dampAux>0.01),'k','linewidth',2,'Color',colorAuxCur{k}(i,:))
        plot(tauAux(auxInd)-1*perAux(auxInd),perAux(auxInd),'k','linewidth',2,'Color',colorAuxRep)
    end

    for j=1:20
        if(k<=6)
            auxPunto1 = auxPuntoHB(k,1);
            auxPunto2 = auxPuntoHB(k,2);
            scatter(auxPunto1+(j-5)*auxPunto2,auxPunto2,30,'s','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',1.5)
            auxPunto1 = auxPuntoHB(k,3);
            auxPunto2 = auxPuntoHB(k,4);
            scatter(auxPunto1+(j-5)*auxPunto2,auxPunto2,30,'s','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',1.5) % Generalised Hopf bifurcation
        end
    end
    
    patch(auxPatchesX(k,:),auxPatchesY(k,:),[0.2,0.2,0.2],'FaceAlpha',0.2,'LineStyle','none')

    auxSliceCurves = setDataFigure.Hom2curvesSlice{k};
    if ~(isempty(auxSliceCurves))
        figure(k); hold on;
        for i=1:length(auxSliceCurves)
            for j=1:20
                tauAux = arrayfun( @(x) x.parameter(ind.tau), auxSliceCurves{i}.point);
                perAux = arrayfun( @(x) x.period, auxSliceCurves{i}.point);
                dampAux = arrayfun( @(x) norm( interp1(x.mesh,x.profile',0:0.05:0.95) - interp1(x.mesh,x.profile',[0.5:0.05:0.95,0:0.05:0.45]) ), auxSliceCurves{i}.point);
                plot(tauAux(dampAux>0.01)+(j-8)*perAux(dampAux>0.01),0.5*perAux(dampAux>0.01),'k','linewidth',2,'Color',[0.7,0,0])
                
            end
        end
    end
    box on
    axis([-4.5  30  0  30])
    xAuxTicks = [-10, 0 10 20];
    yAuxTicks = linspace( 0,   30,4);
    set(gcf,'Color',[1 1 1]);
    set(gcf,'units','centimeters','pos', [5,20,6.9*3,6*3])
    set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
    xticklabels({'','','',''})
    yticklabels({'','','',''})
    hgexport(gcf, ['./Figures/SlicesOrbitFlip_' sprintf('%.4f',setDataFigure.paramSlices(k)) '.eps'], hgexport('factorystyle'), 'Format', 'eps');
end


