%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%        --- Script that generates Figure 1 to 3 ---           %%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of Figure 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
load('../SourceCode2GenerateDataDDEBifTools/DDEResults/PointsResNonOri','setDataFigure')

parnames    =   {'a','b','c','alpha','beta', 'gamma', 'mu', 'mu_tilde', 'tau', 'T','kappa'};
c_ind       =   [parnames;num2cell(1:length(parnames))];
ind         =   struct(c_ind{:});

ticklengthUn = 0.15;
labelSize    = 27;
Colors = [[0, 0, 1]; [0.4940, 0.1840, 0.5560]; [0.75, 0, 0.75]];
ColorsGradient = [ [0.5, 1.447/2, 1.741/2]; [1.494/2, 1.184/2, 1.556/2]; [1.75/2, 0.5, 1.75/2]];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting first profile
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pt=setDataFigure.Points{1};
figure(11); clf; hold on;
for i=-1:3
    plot(0.15+(pt.mesh+i)*pt.period/pt.parameter(ind.tau),pt.profile(1,:),'-', 'Color', Colors(1,:),'linewidth',4);
end
hold off;
box on
axis([0 3 -0.025 1])
xAuxTicks = [0, 1, 2];
yAuxTicks = [];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,9*3,3*3])
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','',''})
yticklabels({'',''})
hgexport(gcf, ['./Figures/Figure1a1.eps'], hgexport('factorystyle'), 'Format', 'eps');

Colors = [[0, 0.4470, 0.7410]; [0.4940, 0.1840, 0.5560]; [0.75, 0, 0.75]];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting second profile
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pt=setDataFigure.Points{2};
figure(12); clf; hold on;
for i=-1:3
    plot(-0.45+(pt.mesh+i)*pt.period/pt.parameter(ind.tau),pt.profile(1,:),'-', 'Color', Colors(2,:),'linewidth',4);
end
hold off;
xlim([0,3])
box on
axis([0 3 -0.025 1])
xAuxTicks = [0, 1, 2];
yAuxTicks = [];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,9*3,3*3])
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','',''})
yticklabels({'',''})
hgexport(gcf, ['./Figures/Figure1b1.eps'], hgexport('factorystyle'), 'Format', 'eps');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting third profile
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pt=setDataFigure.Points{4};
figure(13); clf; hold on;
for i=-1:10
    plot(-0.45+(pt.mesh+i)*pt.period/pt.parameter(ind.tau),pt.profile(1,:),'-', 'Color', ColorsGradient(1,:),'linewidth',4);
end
hold off;
box on
axis([0 3 -0.025 1])
xAuxTicks = [0, 1, 2];
yAuxTicks = [];
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,9*3,3*3])
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','',''})
yticklabels({'',''})
hgexport(gcf, ['./Figures/Figure1c1.eps'], hgexport('factorystyle'), 'Format', 'eps');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Space-Time Representations for the one-pulse TDS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Colors = [[0, 0, 1]; [0.4940, 0.1840, 0.5560]; [0.75, 0, 0.75]];
ColorsGradient = [ [0.8, 0.8, 1]; [1.494/2, 1.184/2, 1.556/2]; [1.75/2, 0.5, 1.75/2]];
cMap = makeColorMap([1 1 1], ColorsGradient(1,:), Colors(1,:), 256);

pt       = setDataFigure.Points{1};
tau      = pt.parameter(ind.tau);
per      = pt.period;

resy     = 100;                 % number of delay intervals to be plotted
aspct    = 1.;                  % aspect ratio of figure 
resx     = floor(resy/aspct);   % resolution of horizontal component (one delay interval)
resT     = ceil(per/tau*resx);  % resolution of period
nprd     = ceil(per/tau)+1;     % number of periods appended to the profile

drft     = resT-resx;           % drift with respect to delay interval
t        = linspace(0,1,resT);  % interpolate solution along equidistant mesh
Z1       = nan*zeros(resy,resx);% grid

I        = interp1(pt.mesh,pt.profile(1,:),t);
I        = repmat(I,1,nprd);
I        = circshift(I, 14);    % set initial drift manually
for i=1:1:resy
    Z1(i,:) = I(1:resx);
end
figure(14); clf; hold on;
s= surf(Z1);
view(0,90);
    xl = xlim();                % Get existing range.
    yl = ylim();
    % Make 2 tick marks.  (Existing # of tick marks may not = 5)
    xTickLocations = linspace(xl(1)+1, xl(2)-2, 2);
    yTickLocations = linspace(yl(1)+1, yl(2)-2, 3);
    set(gca,'XTick', xTickLocations);
    set(gca,'YTick', yTickLocations);
    % Make new labels for the new tick marks
    set(gca,'XTickLabel', {'0' ,'\tau+\delta'});
    set(gca,'YTickLabel', {'0' ,resy/2, resy});
    xlabel('$\theta$', 'interpreter', 'latex')
    ylabel('$n$', 'interpreter', 'latex')
set(gca,'FontSize', 18);
s.EdgeColor = 'none';
s.FaceColor = 'interp';
axis tight;
axis equal;
set(gca,'DataAspectRatio',[aspct*resx/resy 1 1])
colormap(cMap);
box on
axis on
hcb = colorbar();
hcb.Title.String = 'I';
set(gca,'linewidth',3)
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,3*3,3*3])
set(gca,'position',[0.005,0.005,0.99,0.99],'linewidth',3,'YTick',[50],'ticklength',[ticklengthUn/6,0.50])
hgexport(gcf, ['./Figures/Figure1a2.eps'], hgexport('factorystyle'), 'Format', 'eps');

Colors = [[0, 0.4470, 0.7410]; [0.4940, 0.1840, 0.5560]; [0.75, 0, 0.75]];
ColorsGradient = [ [0.5, 1.447/2, 1.741/2]; [1.494/2, 1.184/2, 1.556/2]; [1.75/2, 0.5, 1.75/2]];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Space-Time Representations for the bound two-pulse TDS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cMap = makeColorMap([1 1 1], ColorsGradient(2,:), Colors(2,:), 256);

pt       = setDataFigure.Points{2};
tau      = pt.parameter(ind.tau);
per      = pt.period;

resy     = 100;                     % number of delay intervals to be plotted
aspct    = 1.;                      % aspect ratio of figure 
resx     = floor(resy/aspct);       % resolution of horizontal component (one delay interval)
resT     = ceil(per/tau*resx);      % resolution of period
nprd     = ceil(per/tau)+1;         % number of periods appended to the profile

drft     = resT-resx;               % drift with respect to delay interval
t        = linspace(0,1,resT);      % interpolate solution along equidistant mesh
Z1       = nan*zeros(resy,resx);    % grid

I        = interp1(pt.mesh,pt.profile(1,:),t);
I        = repmat(I,1,nprd);
I        = circshift(I, 60);        % set initial drift manually
for i=1:1:resy
    Z1(i,:) = I(1:resx);
end
figure(15); clf; hold on;
s= surf(Z1);
view(0,90);
    xl = xlim();                    % Get existing range.
    yl = ylim();
    % Make 2 tick marks.  (Existing # of tick marks may not = 5)
    xTickLocations = linspace(xl(1)+1, xl(2)-2, 2);
    yTickLocations = linspace(yl(1)+1, yl(2)-2, 3);
    set(gca,'XTick', xTickLocations);
    set(gca,'YTick', yTickLocations);
    % Make new labels for the new tick marks
    set(gca,'XTickLabel', {'0' ,'\tau+\delta'});
    set(gca,'YTickLabel', {'0' ,resy/2, resy});
    xlabel('$\theta$', 'interpreter', 'latex')
    ylabel('$n$', 'interpreter', 'latex')
set(gca,'FontSize', 18);
s.EdgeColor = 'none';
s.FaceColor = 'interp';
axis tight;
axis equal;
set(gca,'DataAspectRatio',[aspct*resx/resy 1 1])
colormap(cMap);
box on
hcb = colorbar();
hcb.Title.String = 'I';
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,3*3,3*3])
set(gca,'position',[0.005,0.005,0.99,0.99],'linewidth',3,'YTick',[50],'ticklength',[ticklengthUn/6,0.50])
hgexport(gcf, ['./Figures/Figure1b2.eps'], hgexport('factorystyle'), 'Format', 'eps');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Space-Time Representations for the equidistant two-pulse TDS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cMap = makeColorMap([1 1 1], ColorsGradient(1,:), Colors(1,:), 256);

pt       = setDataFigure.Points{4};
tau      = pt.parameter(ind.tau);
per      = pt.period;

resy     = 100;               % number of delay intervals to be plotted
aspct    = 1.;%1.6/1;          % aspect ratio of figure 
resx     = floor(resy/aspct);         % resolution of horizontal component (one delay interval)
resT     = ceil(per/tau*resx);% resolution of period
nprd     = ceil(per/tau)+1;   % number of periods appended to the profile

drft     = resT-resx;         % drift with respect to delay interval
t        = linspace(0,1,resT);% interpolate solution along equidistant mesh
Z1       = nan*zeros(resy,resx); % grid

I        = interp1(pt.mesh,pt.profile(1,:),t);
I        = repmat(I,1,nprd);
I        = circshift(I, 60);    % set initial drift manually
for i=1:1:resy
    %I       = circshift(I, 1);
    Z1(i,:) = I(1:resx);
end
figure(16); clf; hold on;
s= surf(Z1);%surf(Z1(:,1:idtau));
view(0,90);
    xl = xlim(); % Get existing range.
    yl = ylim();
    % Make 2 tick marks.  (Existing # of tick marks may not = 5)
    xTickLocations = linspace(xl(1)+1, xl(2)-2, 2);
    yTickLocations = linspace(yl(1)+1, yl(2)-2, 3);
    set(gca,'XTick', xTickLocations);
    set(gca,'YTick', yTickLocations);
    % Make new labels for the new tick marks
    set(gca,'XTickLabel', {'0' ,'\tau+\delta'});
    set(gca,'YTickLabel', {'0' ,resy/2, resy});
    xlabel('$\theta$', 'interpreter', 'latex')
    ylabel('$n$', 'interpreter', 'latex')
set(gca,'FontSize', 18);
s.EdgeColor = 'none';
s.FaceColor = 'interp';
axis tight;
axis equal;
set(gca,'DataAspectRatio',[aspct*resx/resy 1 1])
colormap(cMap);
box on
hcb = colorbar();
hcb.Title.String = 'I';
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,3*3,3*3])
set(gca,'position',[0.005,0.005,0.99,0.99],'linewidth',3,'YTick',[50],'ticklength',[ticklengthUn/6,0.50])
hgexport(gcf, ['./Figures/Figure1c2.eps'], hgexport('factorystyle'), 'Format', 'eps');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of Figure 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Colors = [[0, 0, 1]; 0.9*[156, 0, 151]/255; [0.75, 0, 0.75]];

per1 = arrayfun( @(x) x.period, setDataFigure.Curves{1}.point);
tau1 = arrayfun( @(x) x.parameter(ind.tau), setDataFigure.Curves{1}.point);
per2 = arrayfun( @(x) x.period, setDataFigure.Curves{2}.point);
tau2 = arrayfun( @(x) x.parameter(ind.tau), setDataFigure.Curves{2}.point);
per3 = arrayfun( @(x) x.period, setDataFigure.Curves{3}.point);
tau3 = arrayfun( @(x) x.parameter(ind.tau), setDataFigure.Curves{3}.point);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting panel (a)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(21); clf; hold on;
for i=1:100
    plot(tau1+i*per1,per1,'linewidth', 2, 'Color', [0, 255, 200]/255)%ColorsGradient(1,:))
end
for i=-1:100
    plot(tau2+i*per2,per2/2, 'linewidth',2, 'Color', ColorsGradient(2,:)*0.7)
end
for i=-1:100
    plot(tau3+i*per3,per3/2, 'linewidth',2, 'Color', ColorsGradient(2,:)*0.7)
end 
plot(tau1,per1,'linewidth',2, 'Color', Colors(1,:))
plot(tau2,per2/2,'linewidth',2, 'Color', Colors(2,:))
plot(tau3,per3/2,'linewidth',2, 'Color', Colors(2,:))
hold off;
box on
set(gca,'linewidth',1);
axis([0  100   0  100])
xlabel('$\tau$','interpreter','latex')
ylabel('$T$','interpreter','latex')
xAuxTicks = linspace(0, 100, 5);
yAuxTicks = linspace(0, 100, 5);
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,6.5*3,6.5*3])
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
xticklabels({'','','',''})
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/Figure2a.eps'], hgexport('factorystyle'), 'Format', 'eps');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting panel (b)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(22); clf; hold on;
plot(tau2-per2,per2/2,'linewidth',3, 'Color', Colors(2,:))
plot(tau3-per3,per3/2,'linewidth',3, 'Color', Colors(2,:))
plot(tau1-per1,per1,'linewidth',3, 'Color', Colors(1,:))
patch([-0.72 -0.64 -0.64 -0.72],[0 0 50 50],[0.2,0.2,0.2],'FaceAlpha',0.2,'LineStyle','none')
box on
xAuxTicks = linspace(-0.8, 0.8, 5);
axis([-0.8  0.8   0    100])
xlabel('$\tau-T$','interpreter','latex')
ylabel('$T$','interpreter','latex')
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,6.5*3,6.5*3])
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',xAuxTicks,'YTick',yAuxTicks,'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/Figure2b.eps'], hgexport('factorystyle'), 'Format', 'eps');

scatter([-0.699891 -0.64904],[8.71691 5.72097],200,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)
hold off;
axis([-0.72  -0.64   0    50])
set(gcf,'units','centimeters','pos', [5,20,3*3,3*3])
set(gca,'XTick',[],'YTick',[]) %[0.07,0.10,0.92,0.88]
hgexport(gcf, ['./Figures/Figure2bZoom.eps'], hgexport('factorystyle'), 'Format', 'eps');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Creation of Figure 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('../SourceCode2GenerateDataDDEBifTools/DDEResults/SpecialPDBranchFig3')
aPD   = arrayfun( @(x) x.parameter(ind.a), auxPDbranch.point);
tauPD = arrayfun( @(x) x.parameter(ind.tau), auxPDbranch.point);
perPD = arrayfun( @(x) x.period, auxPDbranch.point);

x = linspace(-5,100,100);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting panel (a)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(31); clf; hold on;
plot(tauPD,aPD,'linewidth',5, 'Color', [1. 0. 0.]);
plot(x, -1.5007*x./x, '--','linewidth',3, 'Color', [0.5 0.5 0.5]);
hold off;
xlabel('$\tau$','interpreter','latex')
ylabel('$a$','interpreter','latex')
axis([5 35 -1.65 -1.35])
box on
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,6.5*3,6.5*3])
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',[15,25],'YTick',[-1.575 -1.5 -1.425],'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
yticklabels({'','','',''})
hgexport(gcf, ['./Figures/FigurePDComputation_a.eps'], hgexport('factorystyle'), 'Format', 'eps');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting panel (b)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(32); clf; hold on;
plot(tauPD-2*perPD,aPD,'linewidth',5, 'Color', [1. 0.0 0.0]);
plot(x, -1.5*x./x, ' --','linewidth',3, 'Color', [0.5 0.5 0.5]);

% Ploting the homoclinic here!
load('../SourceCode2GenerateDataDDEBifTools/DDEResults/BifCurvesNonOriResPaper_a_tau_closeOrbitFlip')
colorAux =[0.8    0.4    0.1; 0    0.7843    0.7843; 1, 0, 0];
for i=1
    aAux   = arrayfun( @(x) x.parameter(ind.a), auxCurves{i}.point);
    tauAux = arrayfun( @(x) x.parameter(ind.tau), auxCurves{i}.point);
    plot(tauAux,aAux,'Color',colorAux(i,:)*0.75,'linewidth',5)
end

scatter([-0.687],[-1.5001],300,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1],'LineWidth',3.0)

hold off;

xlabel('$\tau-T$','interpreter','latex')
ylabel('$a$','interpreter','latex')
axis([-0.75 -0.55 -1.65 -1.35])
box on
set(gcf,'Color',[1 1 1]);
set(gcf,'units','centimeters','pos', [5,20,6.5*3,6.5*3])
set(gca,'position',[0.001,0.001,0.998,0.998],'XTick',[-0.7 -0.65 -0.6],'YTick',[-1.575 -1.5 -1.425],'FontSize',labelSize,'ticklength',[ticklengthUn/6,0.50],'linewidth',3) %[0.07,0.10,0.92,0.88]
yticklabels({'','','',''})

hgexport(gcf, ['./Figures/FigurePDComputation_b.eps'], hgexport('factorystyle'), 'Format', 'eps');