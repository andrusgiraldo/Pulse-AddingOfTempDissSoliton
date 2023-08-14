%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%         --- Script to create the Figures 11 and 13   ---     %%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of Figure 13(a) : Homoclinic solution before orbit flip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/homOrbitBranch')

x_axis        = @(x) x(:,1);
y_axis        = @(x) x(:,2);
z_axis        = @(x) x(:,3);

auxManyPoints = 20;
[X,Y] = meshgrid([-1.2,1.2],[-1.2,1.2]);

figure(1)
hold on
plot3(x_axis(homOrbitBranch.point(1).profile'),y_axis(homOrbitBranch.point(1).profile'), z_axis(homOrbitBranch.point(1).profile'),'Color',[1,0,0],'linewidth',1.0*5)
scatter3(0,0,0,150,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',2.0)
plot3(zeros(1,auxManyPoints),zeros(1,auxManyPoints),linspace(-0.5,0.5,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
surf(X,Y,0*X,'lines','none','facecolor',[0,0.7,0],'facealpha',0.3);


axis([-1.2000    1.2000   -0.5000    0.5000   -0.2000    0.2000])
xAuxTicks = [];
yAuxTicks = [];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,4*3,4*3])
set(gca,'position',[-0.0,-0.6,1.5,2.0],'XTick',xAuxTicks,'YTick',yAuxTicks,'linewidth',3) 
xticklabels({'','','',''})
yticklabels({'','','',''})
axis off
view(167,5)

hgexport(gcf, ['./Figures/FigureHomPreOrbit.eps'], hgexport('factorystyle'), 'Format', 'eps');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of Figure 13(c) : Homoclinic solution after orbit flip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
load('../SourceCode2GenerateDataDDEBifTools/DDEResults/homOrbitBranch')

x_axis        = @(x) x(:,1);
y_axis        = @(x) x(:,2);
z_axis        = @(x) x(:,3);

auxManyPoints = 20;
[X,Y] = meshgrid([-1.2,1.2],[-1.2,1.2]);

figure(2)
hold on
plot3(x_axis(homOrbitBranch.point(end).profile'),y_axis(homOrbitBranch.point(end).profile'), z_axis(homOrbitBranch.point(end).profile'),'Color',[1,0,0],'linewidth',1.0*5)
scatter3(0,0,0,150,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',2.0)
plot3(zeros(1,auxManyPoints),zeros(1,auxManyPoints),linspace(-0.5,0.5,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
surf(X,Y,0*X,'lines','none','facecolor',[0,0.7,0],'facealpha',0.3);

axis([-1.2000    1.2000   -0.5000    0.5000   -0.2000    0.2000])
xAuxTicks = [];
yAuxTicks = [];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,4*3,4*3])
set(gca,'position',[-0.0,-0.6,1.5,2.0],'XTick',xAuxTicks,'YTick',yAuxTicks,'linewidth',3) 
xticklabels({'','','',''})
yticklabels({'','','',''})
axis off
view(167,5)

hgexport(gcf, ['./Figures/FigureHomPosOrbit.eps'], hgexport('factorystyle'), 'Format', 'eps');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of Figure 13(b) : Homoclinic solution at orbit flip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/homOrbitBranch')

x_axis        = @(x) x(:,1);
y_axis        = @(x) x(:,2);
z_axis        = @(x) x(:,3);

auxManyPoints = 20;
[X,Y] = meshgrid([-1.2,1.2],[-1.2,1.2]);

figure(3)
hold on
plot3(x_axis(homOrbitFlipBranch.point(end).profile'),y_axis(homOrbitFlipBranch.point(end).profile'), z_axis(homOrbitFlipBranch.point(end).profile'),'Color',[1,0,0],'linewidth',1.0*5)
scatter3(0,0,0,150,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',2.0)
plot3(zeros(1,auxManyPoints),zeros(1,auxManyPoints),linspace(-0.5,0.5,auxManyPoints), 'linewidth',1, 'LineStyle','--','Color',[0,0,0])
surf(X,Y,0*X,'lines','none','facecolor',[0,0.7,0],'facealpha',0.3);

axis([-1.2000    1.2000   -0.5000    0.5000   -0.2000    0.2000])
xAuxTicks = [];
yAuxTicks = [];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,4*3,4*3])
set(gca,'position',[-0.0,-0.6,1.5,2.0],'XTick',xAuxTicks,'YTick',yAuxTicks,'linewidth',3) 
xticklabels({'','','',''})
yticklabels({'','','',''})
axis off
view(167,5)

hgexport(gcf, ['./Figures/FigureHomAtOrbit.eps'], hgexport('factorystyle'), 'Format', 'eps');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of Figure 11(b) : Homoclinic as tau goes to infty
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

ticklengthUn = 0.15;
labelSize = 27;

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/approxInfCycle')

x_axis        = @(x) x(:,1);
y_axis        = @(x) x(:,2);
z_axis        = @(x) x(:,3);

auxInd        = 600;
auxSolution   = approxInfCycle.point(auxInd).profile';
[X,Y] = meshgrid([-1.2,1.2],[-1.2,1.2]);

figure(4)
hold on
plot3(x_axis(auxSolution),y_axis(auxSolution), z_axis(auxSolution),'Color',[1,0,0],'linewidth',1.0*3)
scatter3(0,0,0,150,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',2.0)
scatter3(x_axis(auxSolution(end,:)),y_axis(auxSolution(end,:)), z_axis(auxSolution(end,:)),300,'*','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0)
scatter3(x_axis(auxSolution(100,:)),y_axis(auxSolution(100,:)), z_axis(auxSolution(100,:)),300,'*','MarkerEdgeColor',[0,0.7,0.7],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0)

auxManyPoints = 20;
[X,Y] = meshgrid([-1.2,1.2],[-1.2,1.2]);
plot3(zeros(1,auxManyPoints),zeros(1,auxManyPoints),linspace(-0.5,0.5,auxManyPoints), 'linewidth',2, 'LineStyle','--','Color',[0,0,0])
surf(X,Y,0*X,'lines','none','facecolor',[0,0.7,0],'facealpha',0.3);

axis([-1.2000    1.2000   -0.5000    0.5000   -0.3000    0.3000])
xAuxTicks = [];
yAuxTicks = [];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,6*3,6.4*3])
set(gca,'position',[-0.0,-0.85,1.7,2.0],'XTick',xAuxTicks,'YTick',yAuxTicks,'linewidth',3) 
xticklabels({'','','',''})
yticklabels({'','','',''})
axis off
view(167,5)

hgexport(gcf, ['./Figures/FigureAtInfty.eps'], hgexport('factorystyle'), 'Format', 'eps');


% 258
auxPos          = 264;
auxMesh         = approxInfCycle.point(auxInd).mesh;
auxMesh(1:(auxPos-1))  = auxMesh(1:(auxPos-1))+1.0;
auxMesh         = [auxMesh(auxPos:end) auxMesh(1:(auxPos-1))];

auxMeasure      = x_axis(auxSolution);
auxMeasure      = [auxMeasure(auxPos:end);auxMeasure(1:(auxPos-1))];
figure(5)
hold on
set(gcf,'units','centimeters','Color',[1 1 1],'pos', [20,20,6*3,2*3]);
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',[],'YTick',[],'ticklength',[ticklengthUn/6,0.50],'FontSize',24,'linewidth',3)%[0.16,0.15,0.80,0.84]
plot([-3 3],0.696967*[1 1],'linewidth',4.0, 'LineStyle','--','Color',[0,0.7,0.7])
plot([-3 3],0.154054*[1 1],'linewidth',4.0, 'LineStyle','--','Color',[0.5,0.5,0.5])
plot(auxMesh , auxMeasure,'Color',[1,0,0],'linewidth',1.0*3)
plot([1.52254 10],[0,0],'Color',[1,0,0],'linewidth',1.0*3)
plot([-10 0.493461],[0,0],'Color',[1,0,0],'linewidth',1.0*3)
plot([-3 3],0.0*[1 1],'linewidth',4.0, 'LineStyle','--','Color',0*[0.5,0.5,0.5])
axis off
box off 
axis([0.54    1.52         -0.1    1.0000])

hgexport(gcf, ['./Figures/FigureAtInftyXComp.eps'], hgexport('factorystyle'), 'Format', 'eps');


auxMeasure      = y_axis(auxSolution);
auxMeasure      = [auxMeasure(auxPos:end);auxMeasure(1:(auxPos-1))];
figure(6)
hold on
set(gcf,'units','centimeters','Color',[1 1 1],'pos', [20,20,6*3,2*3]);
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',[],'YTick',[],'ticklength',[ticklengthUn/6,0.50],'FontSize',24,'linewidth',3)%[0.16,0.15,0.80,0.84]
plot([-3 3],0.172323*[1 1],'linewidth',4.0, 'LineStyle','--','Color',[0,0.7,0.7])
plot([-3 3],0.109541*[1 1],'linewidth',4.0, 'LineStyle','--','Color',[0.5,0.5,0.5])
plot(auxMesh , auxMeasure,'Color',[1,0,0],'linewidth',1.0*3)
plot([1.52254 10],[0,0],'Color',[1,0,0],'linewidth',1.0*3)
plot([-10 0.493461],[0,0],'Color',[1,0,0],'linewidth',1.0*3)
plot([-3 3],0.0*[1 1],'linewidth',4.0, 'LineStyle','--','Color',0*[0.5,0.5,0.5])
axis off
box off 
axis([0.54    1.52         -0.4    0.4])

hgexport(gcf, ['./Figures/FigureAtInftyYComp.eps'], hgexport('factorystyle'), 'Format', 'eps');


auxMeasure      = z_axis(auxSolution);
auxMeasure      = [auxMeasure(auxPos:end);auxMeasure(1:(auxPos-1))];
figure(7)
hold on
set(gcf,'units','centimeters','Color',[1 1 1],'pos', [20,20,6*3,2*3]);
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',[],'YTick',[],'ticklength',[ticklengthUn/6,0.50],'FontSize',24,'linewidth',3)%[0.16,0.15,0.80,0.84]
plot([-3 3],0.213668*[1 1],'linewidth',4.0, 'LineStyle','--','Color',[0,0.7,0.7])
plot([-3 3],0.0333472*[1 1],'linewidth',4.0, 'LineStyle','--','Color',[0.5,0.5,0.5])
plot(auxMesh , auxMeasure,'Color',[1,0,0],'linewidth',1.0*3)
plot([1.52254 10],[0,0],'Color',[1,0,0],'linewidth',1.0*3)
plot([-10 0.493461],[0,0],'Color',[1,0,0],'linewidth',1.0*3)
plot([-3 3],0.0*[1 1],'linewidth',4.0, 'LineStyle','--','Color',0*[0.5,0.5,0.5])
axis off
box off 
axis([0.54    1.52         -0.025    0.28])

hgexport(gcf, ['./Figures/FigureAtInftyZComp.eps'], hgexport('factorystyle'), 'Format', 'eps');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of Figure 11(a) : Periodic orbit converging to period-two solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

ticklengthUn = 0.15;
labelSize = 27;

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/HetBranch')

x_axis        = @(x) x(:,1);
y_axis        = @(x) x(:,2);
z_axis        = @(x) x(:,3);

auxInd        = 1;
auxSolution   = setDataFigure.Points{1}.profile';
[X,Y] = meshgrid([-1.2,1.2],[-1.2,1.2]);

equil1 = [0.677338, 0.208242, 0.203802];
equil2 = [0.274986, 0.195758, 0.0637645];

figure(4)
hold on
plot3(x_axis(auxSolution),y_axis(auxSolution), z_axis(auxSolution),'Color',[1,0,0],'linewidth',1.0*3)
scatter3(0,0,0,150,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',2.0)
scatter3(x_axis(equil1),y_axis(equil1), z_axis(equil1),300,'*','MarkerEdgeColor',[0,0.7,0.7],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0)
scatter3(x_axis(equil2),y_axis(equil2), z_axis(equil2),300,'*','MarkerEdgeColor',[0.5,0.5,0.5],'MarkerFaceColor',0*[1,1,1],'LineWidth',2.0)

auxManyPoints = 20;
[X,Y] = meshgrid([-1.2,1.2],[-1.2,1.2]);
plot3(zeros(1,auxManyPoints),zeros(1,auxManyPoints),linspace(-0.5,0.5,auxManyPoints), 'linewidth',2, 'LineStyle','--','Color',[0,0,0])
surf(X,Y,0*X,'lines','none','facecolor',[0,0.7,0],'facealpha',0.3);

axis([-1.2000    1.2000   -0.5000    0.5000   -0.3000    0.3000])
xAuxTicks = [];
yAuxTicks = [];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,6*3,6.4*3])
set(gca,'position',[-0.0,-0.85,1.7,2.0],'XTick',xAuxTicks,'YTick',yAuxTicks,'linewidth',3) 
xticklabels({'','','',''})
yticklabels({'','','',''})
axis off
view(167,5)
hgexport(gcf, ['./Figures/FigureAtInftyPeriodicTwoPeriod.eps'], hgexport('factorystyle'), 'Format', 'eps');


% This is for the temporal profile
pt1     =  setDataFigure.Points{1};
ind.tau =  9;
figure(5)
hold on
set(gcf,'units','centimeters','Color',[1 1 1],'pos', [20,20,6*3,2*3]);
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',[],'YTick',[],'ticklength',[ticklengthUn/6,0.50],'FontSize',24,'linewidth',3)%[0.16,0.15,0.80,0.84]

for i=-1:3
    plot(-0.3+(pt1.mesh+i)*pt1.period/abs(pt1.parameter(ind.tau)),pt1.profile(1,:),'-', 'Color',[1,0,0],'linewidth',3);
end
plot([-3 3],x_axis(equil2)*[1 1],'linewidth',4.0, 'LineStyle','--','Color',[0.5,0.5,0.5])
plot([-3 3],x_axis(equil1)*[1 1],'linewidth',4.0, 'LineStyle','--','Color',[0,0.7,0.7])
plot([-3 3],0.0*[1 1],'linewidth',4.0, 'LineStyle','--','Color',0*[0.5,0.5,0.5])

axis off
box on 
axis([-1    3         -0.1    1.0000])

hgexport(gcf, ['./Figures/FigureAtInftyXCompPeriodicTwoPeriod.eps'], hgexport('factorystyle'), 'Format', 'eps');


figure(6)
hold on
set(gcf,'units','centimeters','Color',[1 1 1],'pos', [20,20,6*3,2*3]);
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',[],'YTick',[],'ticklength',[ticklengthUn/6,0.50],'FontSize',24,'linewidth',3)%[0.16,0.15,0.80,0.84]
for i=-1:3
    plot(-0.3+(pt1.mesh+i)*pt1.period/abs(pt1.parameter(ind.tau)),pt1.profile(2,:),'-', 'Color',[1,0,0],'linewidth',3);
end
plot([-3 3],y_axis(equil2)*[1 1],'linewidth',4.0, 'LineStyle','--','Color',[0.5,0.5,0.5])
plot([-3 3],y_axis(equil1)*[1 1],'linewidth',4.0, 'LineStyle','--','Color',[0,0.7,0.7])
plot([-3 3],0.0*[1 1],'linewidth',4.0, 'LineStyle','--','Color',0*[0.5,0.5,0.5])
axis off
box on
axis([-1    3          -0.05    0.32])
hgexport(gcf, ['./Figures/FigureAtInftyYCompPeriodicTwoPeriod.eps'], hgexport('factorystyle'), 'Format', 'eps');


figure(7)
hold on
set(gcf,'units','centimeters','Color',[1 1 1],'pos', [20,20,6*3,2*3]);
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',[],'YTick',[],'ticklength',[ticklengthUn/6,0.50],'FontSize',24,'linewidth',3)%[0.16,0.15,0.80,0.84]
for i=-1:3
    plot(-0.3+(pt1.mesh+i)*pt1.period/abs(pt1.parameter(ind.tau)),pt1.profile(3,:),'-', 'Color',[1,0,0],'linewidth',3);
end
plot([-3 3],z_axis(equil2)*[1 1],'linewidth',4.0, 'LineStyle','--','Color',[0.5,0.5,0.5])
plot([-3 3],z_axis(equil1)*[1 1],'linewidth',4.0, 'LineStyle','--','Color',[0,0.7,0.7])
plot([-3 3],0.0*[1 1],'linewidth',4.0, 'LineStyle','--','Color',0*[0.5,0.5,0.5])
axis off
box on
axis([-1    3         -0.025    0.28])
hgexport(gcf, ['./Figures/FigureAtInftyZCompPeriodicTwoPeriod.eps'], hgexport('factorystyle'), 'Format', 'eps');



