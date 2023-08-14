%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%       --- Script that generates Figure 5, 8 and 12 --        %%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of one of the insets of Figure 5 
%  Bifurcation diagram of orientable resonant delay case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

ticklengthUn = 0.15;
labelSize = 27;
% Parameter indices
parnames    =   {'a','b','c','alpha','beta', 'gamma', 'mu', 'mu_tilde', 'tau', 'T'};
c_ind       =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(c_ind{:});

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesOriResPaper_a_tau_closeOrbitFlip')

colorAux =[0.8    0.4    0.1; 0    1  0]
figure(5); clf; hold on;

auxManyPoints = 10;
plot(linspace(-2.0,-0.0,auxManyPoints), -1.5006*ones(1,auxManyPoints), 'linewidth',3, 'LineStyle','--','Color',[255, 191, 0]/255)

for i=length(auxCurves):-1:1
    aAux   = arrayfun( @(x) x.parameter(ind.a), auxCurves{i}.point);
    tauAux = arrayfun( @(x) x.parameter(ind.tau), auxCurves{i}.point);
    plot(tauAux,aAux,'Color',colorAux(i,:),'linewidth',5)
end

scatter([-1.56044],[-1.50018],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)

box on
axis([-1.6   -1.5   -2.0   -0.5])

xAuxTicks = [];
yAuxTicks = [];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,3*3,2.5*3])
set(gca,'position',[0.002,0.002,0.996,0.996],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/FigureResDelayOriZoom2.eps'], hgexport('factorystyle'), 'Format', 'eps');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of the second inset of Figure 5 
%  Bifurcation diagram in parameter plane of orientable resonant delay case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

ticklengthUn = 0.15;
labelSize = 27;
% Parameter indices
parnames    =   {'a','b','c','alpha','beta', 'gamma', 'mu', 'mu_tilde', 'tau', 'T'};
c_ind       =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(c_ind{:});

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesOriResPaper_a_tau_closeOrbitFlip')
load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesOriResPaper_a_tau_closeOrbitFlip_OneDBifD_Hopf.mat')

colorAux =[0.8    0.4    0.1; 0    1  0];
figure(5); clf; hold on;
auxManyPoints = 10;
plot(linspace(-2.0,-0.0,auxManyPoints), -1.5006*ones(1,auxManyPoints), 'linewidth',3, 'LineStyle','--','Color',[255, 191, 0]/255)


% Hopf curve of notrivial ss
aAux = arrayfun( @(x) x.parameter(ind.a), auxauxCurves.point);
tauAux = arrayfun( @(x) x.parameter(ind.tau), auxauxCurves.point);
plot(tauAux, aAux, '-', 'Color', [0.3922         0    0.3922],'linewidth',5)
% Curves
for i=1:length(auxCurves)%:-1:1
    aAux   = arrayfun( @(x) x.parameter(ind.a), auxCurves{i}.point);
    tauAux = arrayfun( @(x) x.parameter(ind.tau), auxCurves{i}.point);
    plot(tauAux,aAux,'Color',colorAux(i,:),'linewidth',5)
end
box on
axis([-1.21   -1.08   -1.65   -1.15])

scatter([-1.10849],[-1.30961],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0) % Generalised Hopf bifurcation
scatter([-1.16877],[-1.5006],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0) % Second Resonant


xAuxTicks = [];
yAuxTicks = [];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,3*3,2.5*3])
set(gca,'position',[0.002,0.002,0.996,0.996],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/FigureResDelayOriZoom3.eps'], hgexport('factorystyle'), 'Format', 'eps');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of the main diagram of Figure 5 
%  Bifurcation diagram in parameter plane of orientable resonant delay case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

ticklengthUn = 0.15;
labelSize = 27;
% Parameter indices
parnames    =   {'a','b','c','alpha','beta', 'gamma', 'mu', 'mu_tilde', 'tau', 'T'};
c_ind       =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(c_ind{:});

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesOriResPaper_a_tau_closeOrbitFlip')
load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesOriResPaper_a_tau_closeOrbitFlip_OneDBifD_Hopf.mat')

colorAux =[0.8    0.4    0.1; 0    1  0];
figure(5); clf; hold on;

auxManyPoints = 10;
plot(linspace(-2.0,-0.0,auxManyPoints), -1.5006*ones(1,auxManyPoints), 'linewidth',3, 'LineStyle','--','Color',[255, 191, 0]/255)


% Extra Fold Branch coming from the resonant point
load('../SourceCode2GenerateDataDDEBifTools/DDEResults/SpecialFOLDBranchFig5.mat')
aAux = arrayfun( @(x) x.parameter(ind.a), auxFOLDbranch.point);
tauAux = arrayfun( @(x) x.parameter(ind.tau), auxFOLDbranch.point);
plot(tauAux, aAux, '-', 'Color', [0         1    0],'linewidth',5) % Hopf


% Hopf curve of notrivial ss
aAux = arrayfun( @(x) x.parameter(ind.a), auxauxCurves.point);
tauAux = arrayfun( @(x) x.parameter(ind.tau), auxauxCurves.point);
plot(tauAux, aAux, '-', 'Color', [0.3922         0    0.3922],'linewidth',5) % Hopf
% Curves
for i=length(auxCurves):-1:1
    aAux   = arrayfun( @(x) x.parameter(ind.a), auxCurves{i}.point);
    tauAux = arrayfun( @(x) x.parameter(ind.tau), auxCurves{i}.point);
    plot(tauAux,aAux,'Color',colorAux(i,:),'linewidth',5)
end


scatter([-1.10849],[-1.30961],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0) % Generalised Hopf bifurcation

scatter([-1.56044],[-1.50018],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)% First resonant
scatter([-1.16877],[-1.5006],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0) % Second Resonant


plot(linspace(-2.0,-0.0,auxManyPoints), 0.25*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-2.0,-0.0,auxManyPoints), 0.1*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-2.0,-0.0,auxManyPoints),-0.08*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-2.0,-0.0,auxManyPoints),-0.35*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-2.0,-0.0,auxManyPoints),-0.38*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-2.0,-0.0,auxManyPoints),-1.1*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-2.0,-0.0,auxManyPoints),-1.7*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-2.0,-0.0,auxManyPoints),-2.0*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])

% Maxima and minima for curves
scatter([-1.28206],[0.00170113],150,'^','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0) % Second Resonant
scatter([-1.31017],[-0.368925],150,'^','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0) % Second Resonant
scatter([-1.31017],[-0.368925],150,'^','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0) % Second Resonant
scatter([-1.28141],[0.211409],150,'^','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0) % Second Resonant
scatter([-1.40104],[-1.81686],200,'sq','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0) % Second Resonant


box on
axis([-1.75   -0.75   -2.2000    0.3])

xAuxTicks = linspace(-1.75,   -0.75,5);
yAuxTicks = linspace(-2.2000,    0.3,5);
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,14*3,6*3])
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/12,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/FigureResDelayOri.eps'], hgexport('factorystyle'), 'Format', 'eps');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of the first zoom of Figure 8 
%  Bifurcation diagram in parameter plane of 
%  nonorientable resonant delay case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

ticklengthUn = 0.15;
labelSize = 27;
% Parameter indices
parnames    =   {'a','b','c','alpha','beta', 'gamma', 'mu', 'mu_tilde', 'tau', 'T'};
c_ind       =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(c_ind{:});

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesNoNOriResPaper_a_tau_closeOrbitFlip')

colorAux =[0.8    0.4    0.1; 0    0.7843    0.7843; 1, 0, 0];
figure(1); clf; hold on;

auxManyPoints = 10;
plot(linspace(-2.0,-0.0,auxManyPoints), -1.5007*ones(1,auxManyPoints), 'linewidth',3, 'LineStyle','--','Color',[255, 191, 0]/255)


load('../SourceCode2GenerateDataDDEBifTools/DDEResults/SpecialPDBranchFig3')
aAux   = arrayfun( @(x) x.parameter(ind.a), auxPDbranch.point);
tauAux = arrayfun( @(x) x.parameter(ind.tau), auxPDbranch.point);
plot(tauAux,aAux,'Color',colorAux(3,:),'linewidth',5)

for i=length(auxCurves):-1:1
    aAux   = arrayfun( @(x) x.parameter(ind.a), auxCurves{i}.point);
    tauAux = arrayfun( @(x) x.parameter(ind.tau), auxCurves{i}.point);
    plot(tauAux,aAux,'Color',colorAux(i,:),'linewidth',5)
end

% Extra Fold Branch coming from the resonant point
load('../SourceCode2GenerateDataDDEBifTools/DDEResults/SpecialFOLDBranchFig8.mat')
aAux = arrayfun( @(x) x.parameter(ind.a), auxFOLDbranch.point);
tauAux = arrayfun( @(x) x.parameter(ind.tau), auxFOLDbranch.point);
plot(tauAux, aAux, '-', 'Color', [0         0.5    0],'linewidth',5) 

scatter([-0.6868],[-1.5007],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)
scatter([-0.6268],[-1.5007],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)

box on
axis([-0.7496   -0.5786   -1.6546   -1.4000])%[-0.8   -0.55   -1.7   -1.4])

xAuxTicks = [];
yAuxTicks = [];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,3*3,2.5*3])
set(gca,'position',[0.002,0.002,0.996,0.996],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/FigureResDelayNoNOriZoom.eps'], hgexport('factorystyle'), 'Format', 'eps');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of the second zoom of Figure 8 
%  Bifurcation diagram in parameter plane of 
%  nonorientable resonant delay case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

ticklengthUn = 0.15;
labelSize = 27;
% Parameter indices
parnames    =   {'a','b','c','alpha','beta', 'gamma', 'mu', 'mu_tilde', 'tau', 'T'};
c_ind       =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(c_ind{:});

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesNoNOriResPaper_a_tau_closeOrbitFlip')

colorAux =[0.8    0.4    0.1; 0    0.7843    0.7843; 1, 0, 0];
figure(1); clf; hold on;

auxManyPoints = 20;
plot(linspace(-2.0,2.0,auxManyPoints), -1.5007*ones(1,auxManyPoints), 'linewidth',3, 'LineStyle','--','Color',[255, 191, 0]/255)


for i=length(auxCurves):-1:1
    aAux   = arrayfun( @(x) x.parameter(ind.a), auxCurves{i}.point);
    tauAux = arrayfun( @(x) x.parameter(ind.tau), auxCurves{i}.point);
    plot(tauAux,aAux,'Color',colorAux(i,:),'linewidth',5)
end

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesNoNOriResMirrorPaper_a_tau_closeOrbitFlip')
for i=length(auxCurves):-1:1
    aAux   = arrayfun( @(x) x.parameter(ind.a), auxCurves{i}.point);
    tauAux = arrayfun( @(x) x.parameter(ind.tau), auxCurves{i}.point);
    plot(tauAux,aAux,'Color',colorAux(i,:),'linewidth',5)
end
scatter([0.6685],[-1.5007],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)


auxManyPoints = 10;

box on
axis([0.5240    0.7679   -1.6411   -1.2788])

xAuxTicks = [];
yAuxTicks = [];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,3*3,2.5*3])
set(gca,'position',[0.002,0.002,0.996,0.996],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/FigureResDelayNoNOriZoom2.eps'], hgexport('factorystyle'), 'Format', 'eps');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of the main diagram of Figure 8 
%  Bifurcation diagram in parameter plane of 
%  nonorientable resonant delay case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

ticklengthUn = 0.15;
labelSize = 27;

% Parameter indices
parnames    =   {'a','b','c','alpha','beta', 'gamma', 'mu', 'mu_tilde', 'tau', 'T'};
c_ind       =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(c_ind{:});

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesNonOriResPaper_a_tau_closeOrbitFlip')

colorAux =[0.8    0.4    0.1; 0    0.7843    0.7843; 1, 0, 0];
figure(1); clf; hold on;

auxManyPoints = 20;
plot(linspace(-8.0,3.0,auxManyPoints), -1.5007*ones(1,auxManyPoints), 'linewidth',3, 'LineStyle','--','Color',[255, 191, 0]/255)


load('../SourceCode2GenerateDataDDEBifTools/DDEResults/SpecialPDBranchFig3')
aAux   = arrayfun( @(x) x.parameter(ind.a), auxPDbranch.point);
tauAux = arrayfun( @(x) x.parameter(ind.tau), auxPDbranch.point);
plot(tauAux,aAux,'Color',colorAux(3,:),'linewidth',5)

for i=length(auxCurves):-1:1
    aAux   = arrayfun( @(x) x.parameter(ind.a), auxCurves{i}.point);
    tauAux = arrayfun( @(x) x.parameter(ind.tau), auxCurves{i}.point);
    plot(tauAux,aAux,'Color',colorAux(i,:),'linewidth',5)
end

auxManyPoints = 10;
scatter([-0.6868],[-1.5007],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesNoNOriResMirrorPaper_a_tau_closeOrbitFlip')
for i=length(auxCurves):-1:1
    aAux   = arrayfun( @(x) x.parameter(ind.a), auxCurves{i}.point);
    tauAux = arrayfun( @(x) x.parameter(ind.tau), auxCurves{i}.point);
    plot(tauAux,aAux,'Color',colorAux(i,:),'linewidth',5)
end
scatter([0.6685],[-1.5007],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)


% Dashes for the slices
plot(linspace(-8.0,3.0,auxManyPoints),-1.48*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-8.0,3.0,auxManyPoints),-1.56*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-8.0,3.0,auxManyPoints),-1.584*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-8.0,3.0,auxManyPoints),-1.7*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])

plot(linspace(-8.0,3.0,auxManyPoints),-2.0*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-8.0,3.0,auxManyPoints),-2.1*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-8.0,3.0,auxManyPoints),-2.3*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-8.0,3.0,auxManyPoints),-2.4*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])

% Maxima and minima for curves
scatter(-1.5391,-2.38438,200,'sq','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0) % Second Resonant
scatter(-2.71476,-1.84282,150,'^','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0) % Second Resonant
scatter(-4.16497,-2.12259,200,'sq','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0) % Second Resonant
scatter(-5.43999,-2.06996,150,'^','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0) % Second Resonant

box on
axis([-6.5   2.0   -2.45   -1.0])

yAuxTicks = linspace(-2.5,   -1,5);
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,14*3,6*3])
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',[-5 -2.5 0],'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/12,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/FigureResDelayNoNOri.eps'], hgexport('factorystyle'), 'Format', 'eps');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of the first zoom of Figure 12 
%  Bifurcation diagram in parameter plane of 
%  orbit flip delay case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

ticklengthUn = 0.15;
labelSize = 27;

% Parameter indices
parnames    =   {'a','b','c','alpha','beta', 'gamma', 'mu', 'mu_tilde', 'tau', 'T'};
c_ind       =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(c_ind{:});

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesOrbitFlipPaper_mu_tau_closeOrbitFlip')
load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesOrbitFlipPaper_mu_tau_closeOrbitFlip_OneDBifD_Hopf.mat')

figure(5); clf; hold on;
% Hopf curve of notrivial ss
muAux = arrayfun( @(x) x.parameter(ind.mu), auxauxCurves.point);
tauAux = arrayfun( @(x) x.parameter(ind.tau), auxauxCurves.point);
omAux = arrayfun( @(x) x.omega, auxauxCurves.point);
for j=1:20
    plot(tauAux+(j-10)*2*pi./omAux, muAux, '-', 'Color', [.7 .7 1.],'linewidth',5)
end


colorAux =[0.8    0.4    0.1; 0    0.7843    0.7843;1, 0, 0; 0    1  0];
for i=length(auxCurves):-1:1
    aAux   = arrayfun( @(x) x.parameter(ind.mu), auxCurves{i}.point);
    tauAux = arrayfun( @(x) x.parameter(ind.tau), auxCurves{i}.point);
    plot(tauAux,aAux,'Color',colorAux(i,:),'linewidth',4)
end

auxManyPoints = 10;
scatter([-0.7828],[0],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)

box on
axis([-1  -0.55   -0.1    0.1])

xAuxTicks = linspace(-1,  -0.55,5);
yAuxTicks = linspace(-0.1,    0.1,5);
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,6*3,6*3])
hgexport(gcf, ['./Figures/FigureDelayOrbitZoom.eps'], hgexport('factorystyle'), 'Format', 'eps');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of the second inset of Figure 12 
%  Bifurcation diagram in parameter plane of 
%  orbit flip delay case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all;
close all;

figure(1)
hold on

x1 = -0.7483;
y1 =  0.01858;

x2 = -0.8491;
y2 = -0.03199;

m  = (y1-y2)/(x1-x2);


ticklengthUn = 0.15;
labelSize = 27;
fPlane         = @(planeCurve) planeCurve(:,1);
gPlane         = @(planeCurve) planeCurve(:,2)- (m*(planeCurve(:,1)-x1)+y1);

% Parameter indices
parnames    =   {'a','b','c','alpha','beta', 'gamma', 'mu', 'mu_tilde', 'tau', 'T'};
c_ind       =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(c_ind{:});

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesOrbitFlipPaper_mu_tau_closeOrbitFlip')

figure(1); clf; hold on;

colorAux   = [0.8    0.4    0.1; 0    0.7843    0.7843;1, 0, 0; 0    1  0];
interMes   = linspace(-1.0,-0.1021,1000);

tauAuxRef1 = arrayfun( @(x) x.parameter(ind.tau), auxCurves{1}.point);
muAuxRef1  = arrayfun( @(x) x.parameter(ind.mu), auxCurves{1}.point);
auxIndRef  = find((tauAuxRef1'>-1.0).*(tauAuxRef1'<-0.1021));
cs         = spline(tauAuxRef1(auxIndRef),muAuxRef1(auxIndRef));
for i=length(auxCurves):-1:1
    tauAux = arrayfun( @(x) x.parameter(ind.tau), auxCurves{i}.point);
    muAux  = arrayfun( @(x) x.parameter(ind.mu), auxCurves{i}.point);
    
    auxInd    = find((tauAux'>-1.0).*(tauAux'<-0.1021));
    plot(tauAux(auxInd),muAux(auxInd)-ppval(cs,tauAux(auxInd)),'Color',colorAux(i,:),'linewidth',4)
end

auxManyPoints = 10;
scatter(-0.7842,0,300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)

box on
axis([-0.86   -0.71   -2e-3    2e-3])

xAuxTicks = linspace(-1,  -0.55,5);
yAuxTicks = linspace(-0.1,    0.1,5);
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,4*3,2*3])
set(gca,'position',[0.002,0.002,0.996,0.996],'XTick',[],'YTick',[],'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/FigureDelayOrbitZoomSlanted.eps'], hgexport('factorystyle'), 'Format', 'eps');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of the main diagram of Figure 12 
%  Bifurcation diagram in parameter plane of 
%  orbit flip delay case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

ticklengthUn = 0.15;
labelSize = 27;

% Parameter indices
parnames    =   {'a','b','c','alpha','beta', 'gamma', 'mu', 'mu_tilde', 'tau', 'T'};
c_ind       =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(c_ind{:});

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesOrbitFlipPaper_mu_tau_closeOrbitFlip')
load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesOrbitFlipPaper_mu_tau_closeOrbitFlip_OneDBifD_Hopf.mat')

figure(5); clf; hold on;

% Saddle-node periodic bifurcation curve coming from the Orbit Flip
load('../SourceCode2GenerateDataDDEBifTools/DDEResults/SpecialFOLDBranchFig12')
aAux = arrayfun( @(x) x.parameter(ind.mu), auxFOLDbranch.point);
tauAux = arrayfun( @(x) x.parameter(ind.tau), auxFOLDbranch.point);
plot(tauAux, aAux, '-', 'Color', [0         1    0],'linewidth',5) 

% Hopf curve of notrivial ss
muAux = arrayfun( @(x) x.parameter(ind.mu), auxauxCurves.point);
tauAux = arrayfun( @(x) x.parameter(ind.tau), auxauxCurves.point);
omAux = arrayfun( @(x) x.omega, auxauxCurves.point);
for j=1:20
    plot(tauAux+(j-10)*2*pi./omAux, muAux, '-', 'Color', [0.3922         0    0.3922],'linewidth',5) % Hopf
end

colorAux =[0.8    0.4    0.1; 0    0.7843    0.7843;1, 0, 0; 0    1  0];
plot(tauAux, muAux, '-', 'Color', [0.3922         0    0.3922],'linewidth',5) % Hopf

% Plotting the other curves
for i=length(auxCurves):-1:1
    aAux   = arrayfun( @(x) x.parameter(ind.mu), auxCurves{i}.point);
    tauAux = arrayfun( @(x) x.parameter(ind.tau), auxCurves{i}.point);
    plot(tauAux,aAux,'Color',colorAux(i,:),'linewidth',4)
end

auxManyPoints = 10;
scatter([-0.7828],[0],200,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)
scatter([ 0.7305],[0],200,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)

% Dashes for the slices
plot(linspace(-8.0,5.0,auxManyPoints),0.27*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-8.0,5.0,auxManyPoints),0.26*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-8.0,5.0,auxManyPoints),0.2*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-8.0,5.0,auxManyPoints),0.05*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])

plot(linspace(-8.0,5.0,auxManyPoints),-0.10*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-8.0,5.0,auxManyPoints),-0.115*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-8.0,5.0,auxManyPoints),-0.12*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
plot(linspace(-8.0,5.0,auxManyPoints),-0.13*ones(1,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])

% Maxima and minima for curves
scatter(-1.24375,-0.116717,200,'sq','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0) % Second Resonant
scatter( 1.32845,-0.123776,200,'sq','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0) % Second Resonant

scatter(0.0565943,0.264527,120,'^','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0) % Second Resonant
scatter(0.0683301,0.345798,120,'^','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0) % Second Resonant

patch([-1.4387 -1.0516 -1.0516 -1.4387],[-0.1351 -0.1351 -0.0920 -0.0920],[0.2,0.2,0.2],'FaceAlpha',0.2,'LineStyle','none')

box on
axis([-3.75    3.75   -0.22    0.36])

xAuxTicks = linspace(-3.75,   3.75,5);
yAuxTicks = linspace(-0.22,   0.36,5);
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,14*3,6*3])
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/12,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/FigureDelayOrbit.eps'], hgexport('factorystyle'), 'Format', 'eps');

set(gcf,'units','centimeters','pos', [5,20,4*3,2*3])
axis([-1.4387   -1.0516   -0.1351   -0.0920])
hgexport(gcf, ['./Figures/FigureDelayOrbitZoomMinimum.eps'], hgexport('factorystyle'), 'Format', 'eps');