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

%% set up models
MODELS = {'SIthetaCONSTlambdaTVP', ...
    'SIthetaCONSTlambdaCONST', ...
    'SIthetaTVPlambdaTVP', ...
    'SIthetaTVPlambdaCONST', ...
    };

m      = 3;

vintageChoices   = {'GDPD', 'CPI'};
v                = 1;
doInflationNoise = true;


switch v
    case 1
        datadir   = 'datGDPD';
    case 2
        datadir   = 'datCPI';
    otherwise
        datadir = pwd;
end


Nlags = 16;

fontsize = 16;

%% quantile settings
ndxmean     = 1; %#ok<*NASGU>
ndxmedian   = 2;

ndxmid      = ndxmean;

ndxtailIQR    = 2 + [3 6];
% ndxtails    = 2 + [3 8];
% fractiles = [0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995] * 100;


ndxtails90   = 2 + [3 8];
ndxtails    = ndxtails90;
ndxtails9068  = 2 + [3 4 7 8];
ndxtails68  = 2 + [4 7];


%% loop over models, vintageChoices and doInflationNoise
modellabel = MODELS{m};
vintagechoice = vintageChoices{v};

%% process settings
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


%% set up latexwrapper
titlename = strrep(datalabel, '.', '-');
titlename = strcat('SILAMBDA-', titlename);

%% get parameters
if isempty(timestamp)
    filext = sprintf('particles.%s.dat', datalabel);
else
    filext = sprintf('%s.particles.%s.dat', timestamp, datalabel);
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

%% collect results

% RE and SI forecasts
REFORECASTsmoother      = NaN(T,12,Nsv);
for n = 1 : Nsurveys
    REFORECASTsmoother(:,:,n)  = loaddat(fullfile(datadir, sprintf('smootherREFORECASTh%d.%s', n, filext)));
end
SIFORECASTsmoother      = NaN(T,12,Nsv);
for n = 1 : Nsurveys
    SIFORECASTsmoother(:,:,n)  = loaddat(fullfile(datadir, sprintf('smootherSIFORECASTh%d.%s', n, filext)));
end

smootherPIpersistence = loaddat(fullfile(datadir, sprintf('smootherPIpersistence.%s', filext)));
smootherTRENDpersistence = loaddat(fullfile(datadir, sprintf('smootherTRENDpersistence.%s', filext)));
smootherGAPpersistence = loaddat(fullfile(datadir, sprintf('smootherGAPpersistence.%s', filext)));

LAMBDAsmoother      = loaddat(fullfile(datadir, sprintf('smootherLAMBDA.%s', filext)));
BETAsmoother        = loaddat(fullfile(datadir, sprintf('smootherBETA.%s', filext)));
LOGBETAsmoother     = loaddat(fullfile(datadir, sprintf('smootherLOGBETA.%s', filext)));


FEVsmoother = NaN(T,12,Nsurveys);
for n = 1 : Nsurveys
    FEVsmoother(:,:,n)  = loaddat(fullfile(datadir, sprintf('smootherFEVh%d.%s', n, filext)));
end

FEVTRENDsmoother = NaN(T,12,Nsurveys);
for n = 1 : Nsurveys
    FEVTRENDsmoother(:,:,n) = loaddat(fullfile(datadir, sprintf('smootherFEVTRENDh%d.%s', n, filext)));
end
FEVGAPsmoother = NaN(T,12,Nsurveys);
for n = 1 : Nsurveys
    FEVGAPsmoother(:,:,n) = loaddat(fullfile(datadir, sprintf('smootherFEVGAPh%d.%s', n, filext)));
end
FEVNOISEsmoother = NaN(T,12,Nsurveys);
for n = 1 : Nsurveys
    FEVNOISEsmoother(:,:,n) = loaddat(fullfile(datadir, sprintf('smootherFEVNOISEh%d.%s', n, filext)));
end

FEVTRENDSHAREsmoother = NaN(T,12,Nsurveys);
for n = 1 : Nsurveys
    FEVTRENDSHAREsmoother(:,:,n) = loaddat(fullfile(datadir, sprintf('smootherFEVTRENDSHAREh%d.%s', n, filext)));
end
FEVGAPSHAREsmoother = NaN(T,12,Nsurveys);
for n = 1 : Nsurveys
    FEVGAPSHAREsmoother(:,:,n) = loaddat(fullfile(datadir, sprintf('smootherFEVGAPSHAREh%d.%s', n, filext)));
end
FEVNOISESHAREsmoother = NaN(T,12,Nsurveys);
for n = 1 : Nsurveys
    FEVNOISESHAREsmoother(:,:,n) = loaddat(fullfile(datadir, sprintf('smootherFEVNOISESHAREh%d.%s', n, filext)));
end

checkdiff(FEVNOISESHAREsmoother(:,ndxmean,:)+FEVTRENDSHAREsmoother(:,ndxmean,:)+FEVGAPSHAREsmoother(:,ndxmean,:) - 1);

VEsmoother = NaN(T,12,Nsurveys);
for n = 1 : Nsurveys
    VEsmoother(:,:,n) = loaddat(fullfile(datadir, sprintf('smootherVEh%d.%s', n, filext)));
end
R2smoother = NaN(T,12,Nsurveys);
for n = 1 : Nsurveys
    R2smoother(:,:,n) = loaddat(fullfile(datadir, sprintf('smootherR2h%d.%s', n, filext)));
end


SIBIASsmoother = NaN(T,12,Nsurveys);
for n = 1 : Nsurveys
    SIBIASsmoother(:,:,n) = loaddat(fullfile(datadir, sprintf('smootherSIBIASh%d.%s', n, filext)));
end

altlambda = 2:2:8;
altNlambda = length(altlambda);
altSIBIASsmoother = NaN(T,12,Nsurveys, altNlambda);
for m = 1 : altNlambda
    for n = 1 : Nsurveys
        altSIBIASsmoother(:,:,n,m) = loaddat(fullfile(datadir, sprintf('counterfactualSIBIASh%dlambda%d.%s', n, altlambda(m), filext)));
    end
    altSIFORECASTsmoother = NaN(T,12,Nsurveys);
    for n = 1 : Nsurveys
        altSIFORECASTsmoother(:,:,n,m) = loaddat(fullfile(datadir, sprintf('counterfactualSIFORECASTh%dlambda%d.%s', n, altlambda(m), filext)));
    end
end

% Nsmoother = 1e2;
% SIBIASdraws = NaN(Nsmoother,T,Nsurveys);
% for n = 1 : Nsurveys
%     SIBIASdraws(:,:,n) = loaddat(fullfile(datadir, sprintf('smootherSIBIASh%d.draws.%s', n, filext)));
% end

SIMSEsmoother = NaN(T,12,Nsurveys);
for n = 1 : Nsurveys
    SIMSEsmoother(:,:,n) = loaddat(fullfile(datadir, sprintf('smootherSIMSEh%d.%s', n, filext)));
end

SIR2smoother = NaN(T,12,Nsurveys);
for n = 1 : Nsurveys
    SIR2smoother(:,:,n) = loaddat(fullfile(datadir, sprintf('smootherSIR2h%d.%s', n, filext)));
end

SIR2BIASsmoother = NaN(T,12,Nsurveys);
for n = 1 : Nsurveys
    SIR2BIASsmoother(:,:,n) = loaddat(fullfile(datadir, sprintf('smootherSIR2BIASh%d.%s', n, filext)));
end

SIR2PERSISTENCEsmoother = NaN(T,12,Nsurveys);
for n = 1 : Nsurveys
    SIR2PERSISTENCEsmoother(:,:,n) = loaddat(fullfile(datadir, sprintf('smootherSIR2PERSISTENCEh%d.%s', n, filext)));
end


%% FEV stuff

thisfig = figure;
set(gca, 'fontsize', fontsize)
hold on
h = plot(dates, squeeze(FEVNOISESHAREsmoother(:,ndxmid,:)), 'linewidth', 1);
hl = plot(dates, squeeze(LAMBDAsmoother(:,ndxmid)), 'k-.', 'linewidth', 3);
ylim([0 1])
nbershades(dates)
% title('FEVnoiseshare')
legend([hl; h], cat(1, {'\lambda'}, SurveyLabels{:}), 'location', 'Northwest')
figname = 'figFEVnoiseshare';
orient landscape
saveas(gcf, figname, 'epsc')


% thisfig = figure;
% set(gca, 'fontsize', fontsize)
% hold on
% h = plot(dates, squeeze(FEVNOISESHAREsmoother(:,ndxmid,:)), 'linewidth', 1);
% ylim([0 1])
% nbershades(dates)
% % title('FEVnoiseshare')
% legend(h, SurveyLabels{:}, 'location', 'Northwest')
% figname = 'figFEVnoiseshareStep1';
% orient landscape
% print(figname, '-loose', '-depsc')



%% stacked bar charts for SIMSE
% 
% for j = 1 : Nsurveys
%     
%     data = [squeeze(FEVsmoother(:,ndxmid,j)) squeeze(SIBIASsmoother(:,ndxmid,j)) ];
%     
%     thisfig = figure;
%     set(gca, 'fontsize', fontsize)
%     
%     hold on
%     h = bar(dates, data, 1, 'stacked');
%     % make colors b/w compatible
%     set(h(1), 'FaceColor', 0.7 * [1 1 1 ]);
%     set(h(1), 'EdgeColor', 0.7 * [1 1 1 ]);
%     set(h(2), 'FaceColor', 0 * [1 1 1 ]);
%     set(h(2), 'EdgeColor', 0 * [1 1 1 ]);
%     
%     
%     xtickdates(dates)
%     % title(sprintf('SI-MSE (h=%d)', j))
%     % legend(h, 'Bias share', 'Persistence share')
%     legend(h, 'FEV', 'Squared Bias')
%     figname = sprintf('figSIMSEbarsh%d', j);
%     orient landscape
%     print(figname, '-loose', '-depsc')
%     
% end

%%  compare squared bias across counterfactuals, bar plot
% n = 1;
% for thisALT = 2:4
%     thisfig = figure;
%     set(gca, 'fontsize', fontsize)
%     hold on
%     h1 = bar(dates, SIBIASsmoother(:,ndxmid,n), 1, 'facecolor', .2 * [1 1 1]);
%     h2 = plot(dates, squeeze(altSIBIASsmoother(:,ndxmid,n,thisALT)), 'r:',...
%         'linewidth', 2);
%     xtickdates(dates)
%     legend([h1 h2], 'actual SI forecasts', sprintf('counterfactual with \\lambda = %4.2f', altlambda(thisALT)/10))
%     figname = sprintf('figBIASaltSIvsSIforecastH%dLambda%d', n, altlambda(thisALT));
%     orient landscape
%     print(figname, '-loose', '-depsc')
% end


%% stacked bar charts for SIMSE with counterfactual bias
for thisALT = 3 % 2 : 4
    for n = 1 % : Nsurveys
        
        % data  = [squeeze(FEVsmoother(:,ndxmid,n)) squeeze(altSIBIASsmoother(:,ndxmid,n, thisALT)) ];
        data  = [squeeze(FEVsmoother(:,ndxmid,n)) squeeze(SIBIASsmoother(:,ndxmid,n)) ];
        
        REmse = squeeze(FEVsmoother(:,ndxmid,n)) + squeeze(SIBIASsmoother(:,ndxmid,n));
        SImse = squeeze(FEVsmoother(:,ndxmid,n)) + squeeze(altSIBIASsmoother(:,ndxmid,n, thisALT));
        
        thisfig = figure;
        set(gca, 'fontsize', fontsize)
        
        hold on
        h = bar(dates, data, 1, 'stacked');
        % make colors b/w compatible
        set(h(1), 'FaceColor', 0.7 * [1 1 1 ]);
        set(h(1), 'EdgeColor', 0.7 * [1 1 1 ]);
        set(h(2), 'FaceColor', 0 * [1 0 0 ]);
        set(h(2), 'EdgeColor', 0 * [1 0 0 ]);
        
        %     plot(dates, REmse, 'k-.', 'linewidth', 3)
        hl = plot(dates, SImse, 'r-.', 'linewidth', 3);
        
        xtickdates(dates)
        % title(sprintf('SI-MSE (h=%d)', n))
        % legend(h, 'Bias share', 'Persistence share')
        legend([h, hl], 'FEV', 'Squared Bias', sprintf('Counterfactual SI-MSE w/\\lambda=%4.2f' ,altlambda(thisALT)/10))
        figname = sprintf('figSIMSEbarsh%dALTlambda%d', n, altlambda(thisALT));
        orient landscape
        saveas(gcf, figname, 'epsc')
    end
end


%% R2 CPS measure (for appendix)
% for n = 1 : Nsurveys
%     thisfig = figure;
%     hold on
%     plot(dates, squeeze(R2smoother(:,ndxmid,n)), 'k-', 'linewidth', 2);
%     plot(dates, squeeze(R2smoother(:,ndxtails68,n)), 'k--', 'linewidth', 2);
%     plot(dates, squeeze(R2smoother(:,ndxtails90,n)), 'k--', 'linewidth', 1);
%     ylim([0 max(ylim)])
%     nbershades(dates)
%     %     title(sprintf('R2-PERSISTENCE h=%d', n))
%     figname = sprintf('R2PERSISTENCEh%d', n);
%     orient landscape
%     print(figname, '-loose', '-depsc')
% end


%% finish
dockAllFigures
finishscript
