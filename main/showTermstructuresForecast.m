%% plot CAMBRIDGE  results

%% load toolboxes
path(pathdef)

addpath ../matlabbox/emtools/
addpath ../matlabbox/emtexbox/
addpath ../matlabbox/emgibbsbox/
addpath ../matlabbox/emeconometrics/
addpath ../matlabbox/emstatespace/

%% setup
initscript

MODELS = {'SIthetaCONSTlambdaTVP', ...
    'SIthetaCONSTlambdaCONST', ...
    'SIthetaTVPlambdaTVP', ...
    'SIthetaTVPlambdaCONST', ...
    'SIthetaTVPlambdaTVPnormcdf', ...
    };

m = 3;
vintageChoices = {'GDPD', 'CPI'};
v = 1;

doInflationNoise = true;
dofAPFfattail    = 0;

lambdaN = 1;

switch v
    case 1
        datadir   = 'datGDPD';
    case 2
        datadir   = 'datCPI';
    otherwise
        datadir = pwd;
end

%% setting

doWrap   = true;
fontsize = 10;

%% quantile settings
ndxmean     = 1;
ndxmedian   = 2;

ndxmid      = ndxmedian;

ndxtailIQR    = 2 + [5 6];
% ndxtails    = 2 + [3 8];
% fractiles = [0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995] * 100;


ndxtail90   = 2 + [3 8];
ndxtails    = ndxtail90;
ndxtails9068  = 2 + [3 4 7 8];


%% loop over models, vintageChoices and doInflationNoise



modellabel = MODELS{m};



vintagechoice = vintageChoices{v};


%% process settings
tvpA      = contains(modellabel, 'thetaTVP');
tvpLambda = contains(modellabel, 'lambdaTVP');

datasetlabel = sprintf('cambridge2018%s', vintagechoice);


datalabel = sprintf('%s.%s', datasetlabel, modellabel);
if ~doInflationNoise
    datalabel = strcat(datalabel, '.nonoise');
end

Nsurveys = 5;

SurveyLabels = cell(Nsurveys,1);
for n = 1 : Nsurveys
    SurveyLabels{n} = sprintf('h = %d', n);
end

timestamp   = [];

titlename = strrep(datalabel, '.', '-');
titlename = strcat('termstructuresForecasts-',titlename);
if dofAPFfattail > 0
    titlename = sprintf('%s-dofAPFfattail%d', titlename, dofAPFfattail);
end

if doWrap
    initwrap
end


%% get parameters
filext = sprintf('particles.%s', datalabel);
if dofAPFfattail > 0
    filext = sprintf('%s.dofAPFfattail%d', filext, dofAPFfattail);
end
if ~tvpLambda && lambdaN > 1
    filext = sprintf('%s.lambdaN%d', filext, lambdaN);
end
filext = strcat(filext, '.dat');
if ~isempty(timestamp)
    filext = strcat(timestamp, '.', filext);
end


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

% inflation forecasts
PIhatRE = loaddat(fullfile(datadir, sprintf('PIHAT.%s', filext)));
PIhatSI = loaddat(fullfile(datadir, sprintf('PIFCST.%s', filext)));

% smoother forecasts
% counterfactual forecasts
smootherPIhatSI = NaN(T, 12, Nsurveys);
for h = 1 : Nsurveys
    filename = sprintf('smootherSIFORECASTh%d.%s', h, filext);
    smootherPIhatSI(:,:,h) = importdata(fullfile(datadir, filename));
end
smootherPIhatRE = NaN(T, 12, Nsurveys);
for h = 1 : Nsurveys
    filename = sprintf('smootherREFORECASTh%d.%s', h, filext);
    smootherPIhatRE(:,:,h) = importdata(fullfile(datadir, filename));
end

% counterfactual forecasts
ndxCounterfactuals = [2 4 6 8];
Ncounterfactuals   = length(ndxCounterfactuals);
smootherPIhatSIcounterfact = NaN(T, 12, Nsurveys, Ncounterfactuals);
for n = 1 : Ncounterfactuals
    for h = 1 : Nsurveys
        filename = sprintf('counterfactualSIFORECASTh%dlambda%d.%s', h, ndxCounterfactuals(n), filext);
        smootherPIhatSIcounterfact(:,:,h,n) = importdata(fullfile(datadir, filename));
    end
end

% lambda
LAMBDA            = importdata(fullfile(datadir, sprintf('LAMBDA.%s', filext)));
LAMBDAhat         = importdata(fullfile(datadir, sprintf('LAMBDAHAT.%s', filext)));
LAMBDAsmoother    = loaddat(fullfile(datadir, sprintf('smootherLAMBDA.%s', filext)));

%% plot model-implied forecasts: RE vs SI (filtered)
for t =  [26 114 198] % 1  : T
    
    thisLambda = LAMBDAsmoother(t,ndxmid);
    
    hanni = NaN(3,1);
    thisfig = figure;
    set(gca, 'fontsize', fontsize)
    hold on
    hanni(1) = plot(1:Nsurveys, squeeze(smootherPIhatRE(t,ndxmid,:)), 'k-', 'linewidth', 2);
    hanni(2) = plot(1:Nsurveys, squeeze(smootherPIhatSI(t,ndxmid,:)), 'r--', 'linewidth', 2);
    %     hanni(3) = plot(1:Nsurveys, squeeze(smootherPIhatSIcounterfact(t,ndxmid,:,3)), 'b--', 'linewidth', 2);
    hanni(3) = plot(1:Nsurveys, y(t,2:end), 'r:d', 'linewidth', 1);
    %     plot(1:Nsurveys, y(t+(0:Nsurveys-1),1), 'bx', 'linewidth', 1);
    lambs = ylim;
    lambs(2) = floor(lambs(2) * 1 + 1);
    %     ylim([min(0, lambs(1)), lambs(2)])
    %     if lambs(2) > 4
    %         set(gca, 'ytick', 0 : 2 : 10)
    %     else
    %         set(gca, 'ytick', 0 : 10)
    %     end
    set(gca, 'xtick', 0 : 10)
    %     xlabel('h')
    pbaspect([4 1 1])
    orient(thisfig, 'landscape')
    saveas(thisfig, sprintf('forecastCrossectionT%s', datestr(dates(t), 'yyyyqq')), 'epsc');
    legend(hanni, 'RE (model)', 'SI (model)', 'SPF data', 'location', 'southwest')
    saveas(thisfig, sprintf('forecastCrossectionT%s-withlegend', datestr(dates(t), 'yyyyqq')), 'epsc');
    
    fprintf('lambda(%s) = %5.2f\n', datestr(dates(t), 'yyyy:qq'), thisLambda);
end

%% plot all
if doWrap
    
    close all
    %% catalogue of all term structures
    for t =  1  : T
        
        thisLambda = LAMBDAsmoother(t,ndxmid);
        
        hanni = NaN(3,1);
        thisfig = figure;
        set(gca, 'fontsize', fontsize)
        hold on
        hanni(1) = plot(1:Nsurveys, squeeze(smootherPIhatRE(t,ndxmid,:)), 'k-', 'linewidth', 2);
        hanni(2) = plot(1:Nsurveys, squeeze(smootherPIhatSI(t,ndxmid,:)), 'r--', 'linewidth', 2);
        %     hanni(3) = plot(1:Nsurveys, squeeze(smootherPIhatSIcounterfact(t,ndxmid,:,3)), 'b--', 'linewidth', 2);
        hanni(3) = plot(1:Nsurveys, y(t,2:end), 'r:d', 'linewidth', 1);
        %     plot(1:Nsurveys, y(t+(0:Nsurveys-1),1), 'bx', 'linewidth', 1);
        lambs = ylim;
        lambs(2) = floor(lambs(2) * 10 + 5) / 10;
        %         ylim([min(0, lambs(1)), lambs(2)])
        % legend(hanni, 'RE', 'SI (actual)', 'SI with \lambda=0.6', 'SPF', 'location', 'best')
        legend(hanni, 'RE (model)', 'SI (model)', 'SPF data', 'location', 'northeast')
        title(sprintf('Forecasts for t=%s when \\lambda=%4.2f', datestr(dates(t), 'yyyy-qq'), thisLambda))
        wrapthisfigure(thisfig,sprintf('forecastCrossectionT%s-withtitle', datestr(dates(t), 'yyyyqq')), wrap)
        close(thisfig)
        
    end
    
    %% catalogue of all term structures w/counterfactual
    for t =  1  : T
        
        thisLambda = LAMBDAsmoother(t,ndxmid);
        
        hanni = NaN(3,1);
        thisfig = figure;
        set(gca, 'fontsize', fontsize)
        hold on
        hanni(1) = plot(1:Nsurveys, squeeze(smootherPIhatRE(t,ndxmid,:)), 'k-', 'linewidth', 2);
        hanni(2) = plot(1:Nsurveys, squeeze(smootherPIhatSI(t,ndxmid,:)), 'r--', 'linewidth', 2);
        hanni(3) = plot(1:Nsurveys, squeeze(smootherPIhatSIcounterfact(t,ndxmid,:,3)), 'b--', 'linewidth', 2);
        hanni(4) = plot(1:Nsurveys, y(t,2:end), 'r:d', 'linewidth', 1);
        %     plot(1:Nsurveys, y(t+(0:Nsurveys-1),1), 'bx', 'linewidth', 1);
        lambs = ylim;
        lambs(2) = floor(lambs(2) * 10 + 5) / 10;
        %         ylim([min(0, lambs(1)), lambs(2)])
        % legend(hanni, 'RE', 'SI (actual)', 'SI with \lambda=0.6', 'SPF', 'location', 'best')
        legend(hanni, 'RE (model)', 'SI (model)', 'SI (\lambda=0.6)', 'SPF data', 'location', 'northeast')
        title(sprintf('Forecasts for t=%s when \\lambda=%4.2f', datestr(dates(t), 'yyyy-qq'), thisLambda))
        wrapthisfigure(thisfig,sprintf('forecastCrossectionT%s-withcounterfactual', datestr(dates(t), 'yyyyqq')), wrap)
%         close(thisfig)
        
    end
    
    %% memo: lambda
    thisfig = figure;
    hold on;
    %         plot(dates, LAMBDA(:,ndxmid), 'b', 'linewidth', 3);
    %         plot(dates, LAMBDA(:,ndxtails), 'b', 'linewidth', 1);
    
    
    plot(dates, LAMBDAsmoother(:,ndxmid), 'r-', 'linewidth', 3);
    plot(dates, LAMBDAsmoother(:,ndxtails), 'r-', 'linewidth', 1);
    
    ylim([0 1])
    
    nbershades(dates)
    wrapthisfigure(thisfig,'LAMBDAsmoothed', wrap);
    
    %% Compare time series of forecasts
    for n = 1 : Nsurveys
        if n == 1
            limerick = [0 12];
        else
            limerick = [0 10];
        end
        hanni = NaN(2,1);
        thisfig = figure;
        %           subplot(2,1,1)
        hold on
        hanni(1) = plot(dates, PIhatRE(:,n), 'k-', 'linewidth', 1);
        hanni(2) = plot(dates, PIhatSI(:,n), 'r--', 'linewidth', 1);
        plot(dates, y(:,n+1), 'b:', 'linewidth', 2)
        ylim(limerick)
        nbershades(dates)
        legend(hanni, 'RE', 'SI')
        title(sprintf('Forecasts for h=%d, RE vs SI', n))
        wrapthisfigure(thisfig,sprintf('forecastREvsSI%d', n), wrap)
    end
    
    
end
%% finish
finishwrap
dockAllFigures
finishscript
