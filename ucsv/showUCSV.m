%% plot ucsv fortran results

%% load toolboxes
path(pathdef)

addpath ../matlabbox/emtools/
addpath ../matlabbox/emtexbox/
addpath ../matlabbox/emgibbsbox/
addpath ../matlabbox/emeconometrics/
addpath ../matlabbox/emstatespace/

%% settings
initscript
showGains    = true;

datasetlabel = 'cambridge2018GDPD';
datalabel    = sprintf('%s.UCSV', datasetlabel);

titlename = datalabel;
initwrap

timestamp   = [];

datadir = pwd;

%% get parameters
if isempty(timestamp)
    filext = sprintf('particles.%s.dat', datalabel);
else
    filext = sprintf('%s.particles.%s.dat', timestamp, datalabel);
end


%% get data

y     = importdata(fullfile(datadir, sprintf('YDATA.%s', filext)))';
yNaN  = logical(importdata(fullfile(datadir, sprintf('YNAN.%s', filext))))';
y(yNaN) = NaN;
dates = importdata(sprintf('%s.dates.txt', datasetlabel));

y = y(:,1);


%% read results
T   = size(y,1);
Ny  = size(y,2);
Nstates = 2;
Nsv     = 2;


type(fullfile(datadir, strcat('settings.', filext)));

% linear states
TAUhatRE = importdata(fullfile(datadir, sprintf('TAUHATRE.%s', filext)));
GAPhatRE = importdata(fullfile(datadir, sprintf('GAPHATRE.%s', filext)));

TAURE = importdata(fullfile(datadir, sprintf('TAURE.%s', filext)));
GAPRE = importdata(fullfile(datadir, sprintf('GAPRE.%s', filext)));

if showGains
    GAIN = NaN(T,Ny,Nstates);
    for s = 1 : Nstates
        GAIN(:,:,s) = importdata(fullfile(datadir, sprintf('GAIN%d.%s', s, filext)));
    end
end

% loglike and ESS
LOGLIKE         = importdata(fullfile(datadir, sprintf('LOGLIKE.%s', filext)));
ESS             = importdata(fullfile(datadir, sprintf('ESS.%s', filext)));

% SV
SV = NaN(T,12,Nsv);
for s = 1 : Nsv
    SV(:,:,s)   = importdata(fullfile(datadir, sprintf('SV%d.%s', s, filext)));
end

SVhat = NaN(T,1,Nsv);
for s = 1 : Nsv
    SVhat(:,:,s)   = importdata(fullfile(datadir, sprintf('SVHAT%d.%s', s, filext)));
end

%% load SCALE PARAMETERS

% hInno
HINNO = NaN(T,12,Nsv);
for s = 1 : Nsv
    HINNO(:,:,s) = importdata(fullfile(datadir, sprintf('HINNO%d.%s', s, filext)));
end
HINNOhat = NaN(T,1,Nsv);
for s = 1 : Nsv
    HINNOhat(:,:,s) = importdata(fullfile(datadir, sprintf('HINNOHAT%d.%s', s, filext)));
end




%% settings
ndxmean     = 1;
ndxmedian   = 2;
ndxtails    = 2 + [3 8]; % 90 percent
% fractiles = [0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995] * 100;

fontsize = 12;

%% show some data
close all
figure
hold on
set(gca, 'fontsize', fontsize)
plot(dates, y(:,1), 'k-.', 'linewidth', 2)
xtickdates(dates)
box off
wrapcf('YDATA', wrap)

%% ESS
figure
hold on
set(gca, 'fontsize', fontsize)
plot(dates, ESS, 'k-', 'linewidth', 2)
ylim([0 1])
xtickdates(dates)
wrapcf('ESS', wrap)

%% LLF

figure
plot(dates, cumsum(LOGLIKE) ./ (1:T)', 'k-')
xtickdates(dates)
title(sum(LOGLIKE))
wrapcf('llf', wrap);

%% filter: TAU and data
hanni = NaN(3,1);
newfigure('TAU')
set(gcf, 'Renderer', 'painters')
hold on
set(gca, 'fontsize', fontsize)
hanni(1) = plotCI(TAUhatRE, TAURE(:,ndxtails), dates, -2, 'k--', 'linewidth', 2);
% plot(dates, TAURE(:,ndxmean), 'r--')
hanni(2) = plot(dates, y(:,1), 'r:', 'linewidth', 2);
xtickdates(dates)
legend(hanni(1:2), 'trend inflation', 'realized inflation')
plotOrigin
box off
wrapcf('TAU', wrap);

newfigure('GAP')
hold on
set(gcf, 'Renderer', 'painters')
set(gca, 'fontsize', fontsize)
plot(dates, GAPRE(:,ndxmean), 'k-', 'linewidth', 2);
xtickdates(dates)
wrapcf('GAP', wrap);


%% GAIN
if showGains
    nanny       = repmat(yNaN, [1 1 Nstates]); %#ok<*UNRCH>
    GAIN(nanny) = NaN;
    
    hanni = NaN(2,1);
    for n = 1 : Ny
        figure
        hold on
        hanni(1) = plot(dates, squeeze(GAIN(:,n,1)), 'b-', 'linewidth', 2);
        hanni(2) = plot(dates, squeeze(GAIN(:,n,2)),'b--', 'linewidth', 2);
        xtickdates(dates)
        legend(hanni, 'Trend', 'Gap')
        %     title(sprintf('Kalman Gains from %s', Ynames{n}))
        wrapcf(sprintf('KalmanGain%d', n), wrap)
    end
end

%% SV
hanni = NaN(2,1);
for n = 1 : Nsv
    figure
    hold on;
    set(gca, 'fontsize', fontsize)
    hanni(1) = plot(dates, SVhat(:,1,n), 'k', 'linewidth', 3);
    plot(dates, SV(:,ndxtails,n), 'k', 'linewidth', 1);
    xtickdates(dates)
    wrapcf(sprintf('SVfilter%d', n), wrap);
end


%% REPORT SCALE PARAMETERS
% hInno
for n = 1 : Nsv
    figure
    hold on;
    set(gca, 'fontsize', fontsize)
    % plot(dates, HINNOhat(:,1,n), 'b', 'linewidth', 3);
    plot(dates, HINNO(:,ndxmedian,n), 'k', 'linewidth', 3);
    plot(dates, HINNO(:,ndxtails,n), 'k', 'linewidth', 1);
    ylim([0 max(ylim)])
    xtickdates(dates)
    wrapcf(sprintf('HINNO%d', n), wrap);
end


%% finish
finishwrap
dockAllFigures
finishscript
