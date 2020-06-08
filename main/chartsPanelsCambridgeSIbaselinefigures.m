%% plot CAMBRIDGE  results

%% load toolboxes
path(pathdef)

addpath ../matlabbox/emtools/
addpath ../matlabbox/emtexbox/
addpath ../matlabbox/emgibbsbox/
addpath ../matlabbox/emeconometrics/
addpath ../matlabbox/emstatespace/

%% init
initscript


MODELS = {'SIthetaCONSTlambdaTVP', ...
    'SIthetaCONSTlambdaCONST', ...
    'SIthetaTVPlambdaTVP', ...
    'SIthetaTVPlambdaCONST', ...
    };

m      = 3;
alt_m  = 1;
alt2_m = 2;

vintageChoices   = {'GDPD', 'CPI'};
v                = 1;
doInflationNoise = true;


datadir = 'datGDPD';
% datadir = pwd;
% datadir = '~/jam/datGDPD';
picdir  = [];

%% setting
doWrap           = false;

fontsize = 16;

%% loop over models, vintageChoices and doInflationNoise



modellabel  = MODELS{m};
alt_modellabel = MODELS{alt_m};
alt2_modellabel = MODELS{alt2_m};



vintagechoice = vintageChoices{v};


%% process settings
datasetlabel = sprintf('cambridge2018%s', vintagechoice);


datalabel = sprintf('%s.%s', datasetlabel, modellabel);
if ~doInflationNoise
    datalabel = strcat(datalabel, '.nonoise');
end

alt_datalabel = sprintf('%s.%s', datasetlabel, alt_modellabel);
if ~doInflationNoise
    alt_datalabel = strcat(alt_datalabel, '.nonoise');
end
alt2_datalabel = sprintf('%s.%s', datasetlabel, alt2_modellabel);
if ~doInflationNoise
    alt2_datalabel = strcat(alt2_datalabel, '.nonoise');
end

Nsurveys = 5;

SurveyLabels = cell(Nsurveys,1);
for n = 1 : Nsurveys
    SurveyLabels{n} = sprintf('h = %d', n);
end

timestamp   = [];

titlename = strrep(datalabel, '.', '-');
if doWrap
    initwrap %#ok<*UNRCH>
end

tvpLambda = ~(strcmpi(modellabel, 'SIconstlambda') || strcmpi(modellabel, 'SIconstthetalambda') ...
    || strcmpi(modellabel, 'REdriftAR1'));

%% get parameters
filext = sprintf('particles.%s.dat', datalabel);
alt_filext = sprintf('particles.%s.dat', alt_datalabel);
alt2_filext = sprintf('particles.%s.dat', alt2_datalabel);


%% get data

y       = importdata(fullfile(datadir, sprintf('YDATA.%s', filext)))';
yNaN    = logical(importdata(fullfile(datadir, sprintf('YNAN.%s', filext))))';
y(yNaN) = NaN;
dates = importdata(sprintf('%s.dates.txt', datasetlabel));



%% read results
T   = size(y,1);
Ny  = size(y,2);
Nstates = 4;
Nsv     = 2;


type(fullfile(datadir, strcat('settings.', filext)));
if ~isempty(wrap)
    copyfile(fullfile(datadir, strcat('settings.', filext)), fullfile(wrap.dir, strcat('settings.', filext)))
    latexwrapper(wrap, 'add', 'listing', strcat('settings.', filext))
end

%% collect results
% linear states
TAUhatRE = importdata(fullfile(datadir, sprintf('TAUHATRE.%s', filext)));
TAUhatSI = importdata(fullfile(datadir, sprintf('TAUHATSI.%s', filext)));
GAPhatRE = importdata(fullfile(datadir, sprintf('GAPHATRE.%s', filext)));
GAPhatSI = importdata(fullfile(datadir, sprintf('GAPHATSI.%s', filext)));

GAPNOISEhatRE = loaddat(fullfile(datadir, sprintf('GAPNOISEHATRE.%s', filext)));
GAPNOISEhatSI = loaddat(fullfile(datadir, sprintf('GAPNOISEHATSI.%s', filext)));

TAURE = importdata(fullfile(datadir, sprintf('TAURE.%s', filext)));
TAUSI = importdata(fullfile(datadir, sprintf('TAUSI.%s', filext)));
GAPRE = importdata(fullfile(datadir, sprintf('GAPRE.%s', filext)));
GAPSI = importdata(fullfile(datadir, sprintf('GAPSI.%s', filext)));

INFLATIONRE = loaddat(fullfile(datadir, sprintf('INFLATIONRE.%s', filext)));
INFLATIONSI = loaddat(fullfile(datadir, sprintf('INFLATIONSI.%s', filext)));



% SV
SV = NaN(T,12,Nsv);
for s = 1 : Nsv
    SV(:,:,s)   = importdata(fullfile(datadir, sprintf('SV%d.%s', s, filext)));
end

SVhat = NaN(T,1,Nsv);
for s = 1 : Nsv
    SVhat(:,:,s)   = importdata(fullfile(datadir, sprintf('SVHAT%d.%s', s, filext)));
end

% lambda
LAMBDA     = importdata(fullfile(datadir, sprintf('LAMBDA.%s', filext)));
LAMBDAhat  = importdata(fullfile(datadir, sprintf('LAMBDAHAT.%s', filext)));
XLAMBDA    = loaddat(fullfile(datadir, sprintf('XLAMBDA.%s', filext)));
XLAMBDAhat = loaddat(fullfile(datadir, sprintf('XLAMBDAHAT.%s', filext)));

if ~isempty(XLAMBDA)
    checkdiff(normcdf(XLAMBDA(:,2:end)), LAMBDA(:,2:end)); % check that normcdf-relationship holds for every quantile (though not mean)
end

alt_LAMBDA     = importdata(fullfile(datadir, sprintf('LAMBDA.%s', alt_filext)));
alt_LAMBDAhat  = importdata(fullfile(datadir, sprintf('LAMBDAHAT.%s', alt_filext)));

alt2_LAMBDA     = importdata(fullfile(datadir, sprintf('LAMBDA.%s', alt2_filext)));
alt2_THETA     = importdata(fullfile(datadir, sprintf('A.%s', alt2_filext)));

% theta
THETA    = loaddat(fullfile(datadir, sprintf('A.%s', filext)));
THETAhat = loaddat(fullfile(datadir, sprintf('AHAT.%s', filext)));

%% load SCALE PARAMETERS

% siga
SIGA    = loaddat(fullfile(datadir, sprintf('SIGA.%s', filext)));
SIGAhat = loaddat(fullfile(datadir, sprintf('SIGAHAT.%s', filext)));

% siglambda
SIGLAMBDA    = loaddat(fullfile(datadir, sprintf('SIGLAMBDA.%s', filext)));
SIGLAMBDAhat = loaddat(fullfile(datadir, sprintf('SIGLAMBDAHAT.%s', filext)));

% hInno
HINNO = NaN(T,12,Nsv);
for s = 1 : Nsv
    HINNO(:,:,s) = importdata(fullfile(datadir, sprintf('HINNO%d.%s', s, filext)));
end
HINNOhat = NaN(T,1,Nsv);
for s = 1 : Nsv
    HINNOhat(:,:,s) = importdata(fullfile(datadir, sprintf('HINNOHAT%d.%s', s, filext)));
end


% measurement error
SIGMA = NaN(T,12,Ny);
for s = 1 : Ny
    SIGMA(:,:,s) = importdata(fullfile(datadir, sprintf('SIGMA%d.%s', s, filext)));
end
SIGMAhat = NaN(T,1,Ny);
for s = 1 : Ny
    SIGMAhat(:,:,s) = importdata(fullfile(datadir, sprintf('SIGMAHAT%d.%s', s, filext)));
end

%% load state uncertainty measure
XSIGall         = loaddat(fullfile(datadir, sprintf('XSIGall.%s', filext)));
XSIGsrv         = loaddat(fullfile(datadir, sprintf('XSIGsrv.%s', filext)));
XSIGinf         = loaddat(fullfile(datadir, sprintf('XSIGinf.%s', filext)));
XSIGinfspflong  = loaddat(fullfile(datadir, sprintf('XSIGinfspflong.%s', filext)));
XSIGspflong     = loaddat(fullfile(datadir, sprintf('XSIGspflong.%s', filext)));

STATEa       = loaddat(fullfile(datadir, sprintf('STATEa.%s', filext)));
STATElambda  = loaddat(fullfile(datadir, sprintf('STATElambda.%s', filext)));
STATEsvol    = loaddat(fullfile(datadir, sprintf('STATEsvol.%s', filext)));
STATEsvol    = STATEsvol';

%% smoother data

% linear states
smootherTAURE = importdata(fullfile(datadir, sprintf('smootherTAURE.%s', filext)));
smootherTAUSI = importdata(fullfile(datadir, sprintf('smootherTAUSI.%s', filext)));
smootherGAPRE = importdata(fullfile(datadir, sprintf('smootherGAPRE.%s', filext)));
smootherGAPSI = importdata(fullfile(datadir, sprintf('smootherGAPSI.%s', filext)));

smootherTAUREdelta  = importdata(fullfile(datadir, sprintf('smootherTAUREdelta.%s', filext)));

if tvpLambda
    LAMBDAsmoother      = loaddat(fullfile(datadir, sprintf('smootherLAMBDA.%s', filext)));
    deltaLAMBDAsmoother = loaddat(fullfile(datadir, sprintf('smootherDELTALAMBDA.%s', filext)));

    alt_LAMBDAsmoother      = loaddat(fullfile(datadir, sprintf('smootherLAMBDA.%s', alt_filext)));
    deltaalt_LAMBDAsmoother = loaddat(fullfile(datadir, sprintf('smootherDELTALAMBDA.%s', alt_filext)));
else
    LAMBDAsmoother = [];
    deltaLAMBDAsmoother = [];
    alt_LAMBDAsmoother = [];
    deltaalt_LAMBDAsmoother = [];
end

THETAsmoother           = loaddat(fullfile(datadir, sprintf('smootherA.%s', filext)));

SVsmoother      = NaN(T,12,Nsv);
deltaSVsmoother = NaN(T,12,Nsv);
for s = 1 : Nsv
    SVsmoother(:,:,s)        = importdata(fullfile(datadir, sprintf('smootherSV%d.%s', s, filext)));
    deltaSVsmoother(:,:,s)   = importdata(fullfile(datadir, sprintf('smootherDELTASV%d.%s', s, filext)));
end


%% settings
ndxmean     = 1;
ndxmedian   = 2;
ndxtails    = 2 + [3 4 7 8];
ndxtail68   = 2 + [4 7];
ndxtail90   = 2 + [3 8];
ndxmid      = 2;

% fractiles = [0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995] * 100;

xtickdates = datenum(1968:5:2018,1,1);


%% DATA (Figure 1)
close all
figname = 'figDATA';
% panellabels = {'a', 'b', 'c', 'd'};
for n = [2 3 5 6]


    figure
    set(gcf, 'defaultLegendAutoUpdate','off')
    set(gca, 'fontsize', fontsize)
    hold on

    h1 = plot(dates, y(:,1), 'r--', 'linewidth', 2);
    h2 = plot(dates, y(:,n), 'b-', 'linewidth', 3);
    ylim([-1 15])
    set(gca, 'ytick', -2 : 2 : 16)
    set(gca, 'xtick', xtickdates)
    if n == 2
        legend([h1 h2], 'Realized Inflation \pi_t', 'SPF Nowcast \pi_{t,1}^{SPF}', 'location', 'northeast')
        % title(sprintf('(%s) Realized Inflation and SPF Inflation Nowcast', panellabels{panel}))
    else
        legend([h1 h2], 'Realized Inflation \pi_t', sprintf('%d-Q Ahead SPF: \\pi_{t,%d}^{SPF}', n-2,n-1), 'location', 'northeast')
        % title(sprintf('(%s) Realized Inflation and %d-Quarter Ahead SPF Prediction', panellabels{panel}, n-2))
    end
    nbershades(xtickdates)
    orient landscape
    panelname = sprintf('%spanelh%d', figname, n-1);
    saveas(gcf, panelname, 'epsc')

end


%% FIGURE 2: static parameters


figname = 'figSTATICPARAMETER';

% panel a: vol of trend SV
n = 1;
panel = 1;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plotCI(HINNO(:,ndxmid,n), HINNO(:,ndxtails,n), dates, 0)
ylim([0 .35])
set(gca, 'ytick', 0 : .05 : 2)
set(gca, 'xtick', xtickdates)
% title('(a) PLE Path of Shock Variance Parameter of ln \varsigma^2_{\eta,t}, \sigma^2_\eta')
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')

% panel b: vol of trend SV
n = 2;
panel = 2;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plotCI(HINNO(:,ndxmid,n), HINNO(:,ndxtails,n), dates, 0)
ylim([0 .2])
set(gca, 'ytick', 0 : .05 : 2)
set(gca, 'xtick', xtickdates)
% title('(b) PLE Path of Shock Variance Parameter of ln \varsigma^2_{\nu,t}, \sigma^2_\nu')
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')

% panel c: gap AR1 coeffcient
panel = 3;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plotCI(SIGA(:,ndxmid), SIGA(:,ndxtails), dates, -1)
ylim([0 .05])
set(gca, 'ytick', 0 : .01 : 1)
set(gca, 'xtick', xtickdates)
% title('(c) PLE Path of Shock Variance Parameter of Gap AR(1) coefficient \theta_{t}')
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')


% panel d: vol of lambda
panel = 4;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plotCI(SIGLAMBDA(:,ndxmid), SIGLAMBDA(:,ndxtails), dates, 0)
ylim([0 .06])
set(gca, 'ytick', 0 : .02 : 1)
set(gca, 'xtick', xtickdates)
% title('(d) PLE Path of Shock Variance Parameter of \lambda_{t}, \sigma^2_\kappa')
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')


%% FIG3


figname = 'figTRENDandGAP';

% panel a: SI Trend and SPF Nowcast
panel = 1;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plotCI(TAUSI(:,ndxmid), TAUSI(:,ndxtail68), dates, 0, 'k-.', 'linewidth', 3)
h2 = plot(dates, y(:,2), 'r-', 'linewidth', 2);
h1 = plot(dates, TAUSI(:,ndxmid), 'k-.', 'linewidth', 3);
ylim([0 10])
set(gca, 'ytick', 0 : 2 : 20)
set(gca, 'xtick', xtickdates)
% title('(a) Filtered SI Trend and SPF Inflation Nowcast')
legend([h1 h2], 'Filtered SI Trend: F_{t|t} \tau_t', 'SPF Nowcast: \pi_{t,t+1}^{SPF}', 'location', 'northeast')
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')


%% panel b: RE Trend and SPF-4Q
panel = 2;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plotCI(TAURE(:,ndxmid), TAURE(:,ndxtail68), dates, 0, 'k-.', 'linewidth', 3)
h2 = plot(dates, y(:,6), 'r-.', 'linewidth', 3);
h1 = plot(dates, TAURE(:,ndxmid), 'k:', 'linewidth', 3);
h3 = plot(dates, y(:,1), 'b-', 'linewidth', 1);
ylim([-2 14])
set(gca, 'ytick', -4 : 2 : 20)
set(gca, 'xtick', xtickdates)
% title('(b) Filtered SI Trend and 4-Quarter Ahead SPF Inflation Prediction')
legend([h3 h1 h2], 'Realized Inflation', 'Filtered Trend: \tau_{t|t}', '4-Q Ahead SPF: \pi_{t,t+5}^{SPF}', 'location', 'northeast')
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')


%% panel c: Inflation, SI and RE Trend
panel = 3;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plotCI(TAUSI(:,ndxmid), TAUSI(:,ndxtail68), dates, 0, 'k-', 'linewidth', 3)
h1 = plot(dates, y(:,1), 'b-.', 'linewidth', 2);
h2 = plot(dates, TAUSI(:,ndxmid), 'k-', 'linewidth', 3);
h3 = plot(dates, TAURE(:,ndxmid), 'r--', 'linewidth', 3);
ylim([-2 14])
set(gca, 'ytick', -2 : 2 : 20)
set(gca, 'xtick', xtickdates)
% title('(c) Realized Inflation, Filtereed SI Trend and Filtered RE Trend')
legend([h1 h2 h3], 'Realized Inflation: \pi_t', 'Filtered SI Trend: F_{t|t} \tau_t', 'Filtered RE Trend: \tau_{t|t}')
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')


%% panel d: Gap
panel = 4;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
% plotCI(GAPSI(:,ndxmid), GAPSI(:,ndxtail68), dates, 0, 'k-', 'linewidth', 2)
h1 = plot(dates, GAPSI(:,ndxmid), 'k-', 'linewidth', 2);
h2 = plot(dates, GAPRE(:,ndxmid), 'r--', 'linewidth', 2);
plotOrigin('k-')
ylim([-5 10])
set(gca, 'ytick', -10 : 2.5 : 20)
set(gca, 'xtick', xtickdates)
% title('(d) Filtered SI Gap Inflation and Filtered RE Gap Inflation')
legend([h1 h2], 'Filtered SI Gap: F_{t|t} \epsilon_t', 'Filtered RE Gap: \epsilon_{t|t}')
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')

%% panel e: Gap (incl noise)
panel = 5;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
h2 = plot(dates, GAPNOISEhatRE, 'r-', 'linewidth', 3);
h1 = plot(dates, GAPNOISEhatSI, 'k:', 'linewidth', 4);
plotOrigin('k-')
ylim([-5 10])
set(gca, 'ytick', -10 : 2.5 : 20)
set(gca, 'xtick', xtickdates)
% title('(d) Filtered SI Gap Inflation and Filtered RE Gap Inflation')
legend([h1 h2], 'SI', 'RE')
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')


%% FIG4

figname = 'figSV';

% panel a: Filtered Trend SV
n = 1;
panel = 1;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plot(dates, SV(:,ndxmid,n), 'r:', 'linewidth', 2);
plot(dates, SV(:,ndxtail90,n), 'k-', 'linewidth', 1);
ylim([0 0.8])
set(gca, 'ytick', 0 : .2 : 1.2)
set(gca, 'xtick', xtickdates)
% title('(a) Filtered Trend SV, \varsigma_{\eta,t|t}, and 90% Uncertainty Bands')
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')


% panel b: Smoothed Trend SV
n = 1;
panel = 2;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plot(dates, SVsmoother(:,ndxmid,n), 'r:', 'linewidth', 2);
plot(dates, SVsmoother(:,ndxtail90,n), 'k-', 'linewidth', 1);
ylim([0 0.8])
set(gca, 'ytick', 0 : .2 : 1.2)
set(gca, 'xtick', xtickdates)
% title('(b) Smoothed Trend SV, \varsigma_{\eta,t|T}, and 90% Uncertainty Bands')
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')

% panel c: Filtered Gap SV
n = 2;
panel = 3;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plot(dates, SV(:,ndxmid,n), 'r:', 'linewidth', 2);
plot(dates, SV(:,ndxtail90,n), 'k-', 'linewidth', 1);
ylim([0 3.5])
set(gca, 'ytick', 0 : .5 : 20)
set(gca, 'xtick', xtickdates)
% title('(c) Filtered Gap SV, \varsigma_{\nu,t|t}, and 90% Uncertainty Bands')
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')

% panel d: Smoothed Gap SV
n = 2;
panel = 4;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plot(dates, SVsmoother(:,ndxmid,n), 'r:', 'linewidth', 2);
plot(dates, SVsmoother(:,ndxtail90,n), 'k-', 'linewidth', 1);
ylim([0 3.5])
set(gca, 'ytick', 0 : .5 : 20)
set(gca, 'xtick', xtickdates)
% title('(b) Smoothed Gap SV, \varsigma_{\nu,t|T}, and 90% Uncertainty Bands')
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')



%% FIG 5 -- lambda

figname = 'figLAMBDA';

% panel a: baseline lambda
panel = 1;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plotCI(LAMBDA(:,ndxmid), LAMBDA(:,ndxtail90), dates, 0, 'k-.', 'linewidth', 2)
plot(dates, LAMBDAsmoother(:,ndxmid), 'r-.', 'linewidth', 2);
plot(dates, LAMBDAsmoother(:,ndxtail90), 'r-', 'linewidth', 1);
ylim([0 1])
set(gca, 'ytick', 0 : .1 : 1)
set(gca, 'xtick', xtickdates)
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')

% panel b: baseline delta-lambda
panel = 2;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plotCI(deltaLAMBDAsmoother(:,ndxmid), deltaLAMBDAsmoother(:,ndxtails), dates, 0, 'w-.', 'linewidth', 2)
plot(xtickdates([1 end]), [0 0], 'k-')
ylim([-.3 1])
set(gca, 'ytick', -1 : .1 : 1)
set(gca, 'xtick', xtickdates)
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')

% panel c: alt lambda
panel = 3;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plotCI(alt_LAMBDA(:,ndxmid), alt_LAMBDA(:,ndxtail90), dates, 0, 'k-.', 'linewidth', 2)
plot(dates, alt_LAMBDAsmoother(:,ndxmid), 'r-.', 'linewidth', 2);
plot(dates, alt_LAMBDAsmoother(:,ndxtail90), 'r-', 'linewidth', 1);
ylim([0 1])
set(gca, 'ytick', 0 : .1 : 1)
set(gca, 'xtick', xtickdates)
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')

% panel b: alt delta-lambda
panel = 4;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plotCI(deltaalt_LAMBDAsmoother(:,ndxmid), deltaalt_LAMBDAsmoother(:,ndxtails), dates, 0, 'w-.', 'linewidth', 2)
plot(xtickdates([1 end]), [0 0], 'k-')
ylim([-.3 1])
set(gca, 'ytick', -1 : .1 : 1)
set(gca, 'xtick', xtickdates)
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')

%% FIG 5b -- lambda and theta
figname = 'figLAMBDATHETA';
figure
set(gcf, 'defaultLegendAutoUpdate','off')

% panel b: baseline lambda
panel = 2;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
plotCI(LAMBDA(:,ndxmid), LAMBDA(:,ndxtail90), dates, 0, 'k-.', 'linewidth', 2)
plot(dates, alt2_LAMBDA(:,ndxmid), 'r-', 'linewidth', 2);
plot(dates, alt2_LAMBDA(:,ndxtail90), 'r-.', 'linewidth', 1);
ylim([0 1])
set(gca, 'ytick', 0 : .1 : 1)
set(gca, 'xtick', xtickdates)
% title({'(b) Models M_2 and M_1:', '\lambda_{t|t} vs PLE path of \lambda'})
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')

% panel a: baseline delta-lambda
panel = 1;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on
h1 = plotCI(THETA(:,ndxmid), THETA(:,ndxtail90), dates, 0, 'k-.', 'linewidth', 2);
h2 = plot(dates, alt2_THETA(:,ndxmid), 'r-', 'linewidth', 2);
plot(dates, alt2_THETA(:,ndxtail90), 'r-.', 'linewidth', 1);
ylim([-.6 1])
plotOrigin
set(gca, 'ytick', -1 : .2 : 1)
set(gca, 'xtick', xtickdates)
% title({'(a) Models M_2 and M_1:', '\theta_{t|t} vs PLE path of \theta'})
legend([h1 h2], 'M_2', 'M_1', 'location', 'southeast')
nberlines(xtickdates)
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')

%% FIG 6 -- uncertainty

figname = 'figTRENDUNCERTAINTY';

hanni = NaN(4,1);

XSIGall(1,:) = NaN;
XSIGinf(1,:) = NaN;
XSIGsrv(1,:) = NaN;
XSIGinfspflong(1,:) = NaN;
XSIGspflong(1,:) = NaN;

% panel a: Uncertainty about RE Trend
n = 1;
panel = 1;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on

hanni(1) = plot(dates, XSIGall(:,n), 'k-', 'linewidth', 3);
hanni(2) = plot(dates, XSIGinf(:,n), 'k:', 'linewidth', 2);
hanni(3) = plot(dates, XSIGspflong(:,n), 'b-.', 'linewidth', 2);
hanni(4) = plot(dates, XSIGinfspflong(:,n), 'r--', 'linewidth', 2);

ylim([0 2])
set(gca, 'ytick', 0 : .5 : 20)
set(gca, 'xtick', xtickdates)
% title('(a) Volatility of Posterior for RE Trend Inflation, \tau_t, Conditional on Different Information Sets')
nberlines(xtickdates)
legend(hanni, 'Full Information Set', 'Realized Inflation: \pi_t', '4-Q Ahead SPF: \pi_{t,t+5}^{SPF}', 'Combination: \pi_t and \pi_{t,t+5}^{SPF}', 'location', 'northeast')
pbaspect([2 1 1])
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')

% panel b: Uncertainty about SI Trend
n = 1;
panel = 2;
figure
set(gcf, 'defaultLegendAutoUpdate','off')
set(gca, 'fontsize', fontsize)
hold on

hanni(1) = plot(dates, XSIGall(:,n), 'k-', 'linewidth', 3);
hanni(2) = plot(dates, XSIGinf(:,n), 'k:', 'linewidth', 2);
hanni(3) = plot(dates, XSIGspflong(:,n), 'b-.', 'linewidth', 2);
hanni(4) = plot(dates, XSIGinfspflong(:,n), 'r--', 'linewidth', 2);

ylim([0 2])
set(gca, 'ytick', 0 : .5 : 20)

set(gca, 'xtick', xtickdates)
legend(hanni, 'Full Information Set', 'Realized Inflation: \pi_t', '4-Q Ahead SPF: \pi_{t,t+5}^{SPF}', 'Combination: \pi_t and \pi_{t,t+5}^{SPF}', 'location', 'northeast')

% title('(b) Volatility of Posterior for SI Trend Inflation, F_t \tau_t, Conditional on Different Information Sets')
nberlines(xtickdates)
pbaspect([2 1 1])
orient landscape
panelname = sprintf('%spanel%d', figname, panel);
saveas(gcf, panelname, 'epsc')




%% finish
finishwrap
dockAllFigures
finishscript
