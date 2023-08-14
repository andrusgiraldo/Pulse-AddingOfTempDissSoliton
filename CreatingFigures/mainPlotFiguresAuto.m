%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%          --- Script that generates Figure 4 ---              %%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Inset inside the Orbit Flip Panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all;
close all;

figure(1)
hold on

ticklengthUn = 0.15;
labelSize = 27;
fPlane         = @(planeCurve) planeCurve(:,1);
gPlane         = @(planeCurve) planeCurve(:,2);

load('PaperCurvesAuto')

auxTag    = 'OrbitFlipPaper';
auxPlanes = bifurcationDiagram.planes(auxTag).curves;
auxKeys   = {'SNP','PD','HOM2','HOM'};
for i=1:length(auxKeys)
    auxCurve = auxPlanes(auxKeys{i});
    x = fPlane(auxCurve{1});
    y = gPlane(auxCurve{1});
    plot(x,y,'Color',bifurcationDiagram.colorMap(auxKeys{i}),'LineWidth',3);
end
scatter([0],[0],200,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',2.0)

box on
axis([-0.075 0.075 -0.04 0.04])

set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,4*3,2*3])
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',[],'YTick',[],'FontSize',labelSize,'ticklength',[ticklengthUn/12,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/FigureOriOrbitFlip.eps'], hgexport('factorystyle'), 'Format', 'eps');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Orbit Flip Panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

figure(1)
hold on

x1 = -0.04325;
y1 = -0.02183;

x2 =  0.02965;
y2 =  0.01496;

m  = (y1-y2)/(x1-x2);


ticklengthUn = 0.15;
labelSize = 27;
fPlane         = @(planeCurve) planeCurve(:,1);
gPlane         = @(planeCurve) planeCurve(:,2)- (m*(planeCurve(:,1)-x1)+y1);


load('PaperCurvesAuto')

auxTag    = 'OrbitFlipPaper';
auxPlanes = bifurcationDiagram.planes(auxTag).curves;
auxKeys   = {'SNP','PD','HOM2','HOM'};
for i=1:length(auxKeys)
    auxCurve = auxPlanes(auxKeys{i});
    x = fPlane(auxCurve{1});
    y = gPlane(auxCurve{1});
    plot(x,y,'Color',bifurcationDiagram.colorMap(auxKeys{i}),'LineWidth',5);
end
scatter([0],[0],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)

box on
axis([-0.075 0.075 -0.005 0.005])

xAuxTicks = linspace(-0.075,0.075,5);
yAuxTicks = linspace(-0.005,0.005,5);
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,14*3,6*3])
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/12,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/FigureSlantedOrbitFlip.eps'], hgexport('factorystyle'), 'Format', 'eps');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Orientable Resonant bifurcation Panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all;
close all;

figure(1)
hold on

ticklengthUn = 0.15;
labelSize = 27;
fPlane         = @(planeCurve) planeCurve(:,1);
gPlane         = @(planeCurve) planeCurve(:,2);

load('PaperCurvesAuto')

auxTag    = 'OrbitFlipOriPaper';
auxPlanes = bifurcationDiagram.planes(auxTag).curves;
%auxKeys   = keys(auxPlanes);
auxKeys   = {'SNP','HOM'};
for i=1:length(auxKeys)
    auxCurve = auxPlanes(auxKeys{i});
    x = fPlane(auxCurve{1});
    y = gPlane(auxCurve{1});
    plot(x,y,'Color',bifurcationDiagram.colorMap(auxKeys{i}),'LineWidth',5);
end
scatter([-1.50072],[-5.92615E-02],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)

box on
axis([-2.0   -0.5   -0.061   -0.056])

xAuxTicks = linspace(-2.0,   -0.5,5);
yAuxTicks = linspace(-0.061, -0.056,5);
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,6*3,4*3])
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/FigureResOriOrbitFlip.eps'], hgexport('factorystyle'), 'Format', 'eps');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Nonorientable Resonant bifurcation Panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s
clc
clear all;
close all;

figure(1)
hold on

ticklengthUn = 0.15;
labelSize = 27;
fPlane         = @(planeCurve) planeCurve(:,1);
gPlane         = @(planeCurve) planeCurve(:,2);

load('PaperCurvesAuto')

auxTag    = 'OrbitFlipNoNOriPaper';
auxPlanes = bifurcationDiagram.planes(auxTag).curves;
auxKeys   = {'PD','HOM2','HOM'};
for i=1:length(auxKeys)
    auxCurve = auxPlanes(auxKeys{i});
    x = fPlane(auxCurve{1});
    y = gPlane(auxCurve{1});
    plot(x,y,'Color',bifurcationDiagram.colorMap(auxKeys{i}),'LineWidth',5);
end
scatter([-1.50072],[5.89599E-02],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)

box on
axis([-2.0   -0.5   0.0526    0.0626])

xAuxTicks = linspace(-2.0,   -0.5,5);
yAuxTicks = linspace(0.0526,    0.0626,5);
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,6*3,4*3])
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/FigureResNoNOriOrbitFlip.eps'], hgexport('factorystyle'), 'Format', 'eps');
