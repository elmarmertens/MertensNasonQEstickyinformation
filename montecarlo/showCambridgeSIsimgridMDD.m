%% plot MDD detection rates

%% load toolboxes
path(pathdef)

addpath ../matlabbox/emtools/
addpath ../matlabbox/emtexbox/
addpath ../matlabbox/emgibbsbox/
addpath ../matlabbox/emeconometrics/
addpath ../matlabbox/emstatespace/

%% setting

initscript

modellabel = 'SIthetaTVPlambdaTVP';

datadir = pwd;
Nmodels = 4;

simlabel = ''; 
% simlabel = '.lownoise'; 

doInflationNoise = true;
T          = 200;
Ngrid      = 640;
MDDgrid    = 10;
Nparticles = 1e4;
simtype    = '';


%#ok<*UNRCH>

if doInflationNoise
    noiselabel = [];
else
    noiselabel = 'nonoise';
end

doDateAxis = false;
leadTime = 40;

%% comparisons

% Models

MODELS = {'SIthetaCONSTlambdaTVP', ...
    'SIthetaCONSTlambdaCONST', ...
    'SIthetaTVPlambdaTVP', ...
    'SIthetaTVPlambdaCONST', ...
    };

comparisons(1).baseline    = 3;
comparisons(1).alternative = 1;

comparisons(2).baseline    = 3;
comparisons(2).alternative = 2;

comparisons(3).baseline    = 3;
comparisons(3).alternative = 4;

comparisons(4).baseline    = 2;
comparisons(4).alternative = 1;

comparisons(5).baseline    = 3;
comparisons(5).alternative = 1;

comparisons(6).baseline    = 4;
comparisons(6).alternative = 1;


for c = 1 : length(comparisons)
    comparisons(c).label = sprintf('%s-vs-%s', ...
        MODELS{comparisons(c).baseline}, ...
        MODELS{comparisons(c).alternative});
end

%% settings
datasetlabel = sprintf('simdataT%d', T);
datalabel    = sprintf('%s.%s', datasetlabel, modellabel);
if ~doInflationNoise
    datalabel = strcat(datalabel, '.nonoise');
end
datalabel    = sprintf('%s.Ngrid%d.Nparticles%d', datalabel, Ngrid, Nparticles);

titlename = strcat('delta2MDD', simtype, datalabel);
titlename = strrep(titlename, '.', '_');
initwrap

timestamp   = [];

%% prep objects
filext = sprintf('particles.%sMDDgrid%d.%s%s.dat', simtype, MDDgrid, datalabel, simlabel);

if doDateAxis
    dates = genrQdates(1960, 1960 + ceil(T / 4));
    dates = dates(1:T);
else
    dates  = 1:T;
    xticks = 0 : 40 : T;
end

type(fullfile(datadir, sprintf('settings.%s', filext)))

%% load results
LOGMDD = NaN(T,Ngrid,Nmodels);

for n = 1 : Nmodels
    LOGMDD(:,:,n) = importdata(fullfile(datadir, sprintf('LOGMDD%d.%s', n-1, filext)));
end



%% plot pairwise comparisons
for c = 1 : length(comparisons)
    %% compare
    deltaLOGMDD = 2 * (LOGMDD(:,:,comparisons(c).baseline) - LOGMDD(:,:,comparisons(c).alternative));
    
    
    %% report strong evidence
    positive = sum(deltaLOGMDD > 2,2) / Ngrid * 100;
    negative = sum(deltaLOGMDD < -2,2) / Ngrid * 100;
    
    negative = round(negative);
    positive = round(positive);

    thisfig = figure;
    hold on
    set(gca, 'fontsize', 20)
    hp = plot(dates, positive, 'k-', 'linewidth', 2);
    hn = plot(dates, negative, 'r--', 'linewidth', 1);
    legend([hp, hn], sprintf('pro M_%d', comparisons(c).baseline-1), ...
        sprintf('against M_%d', comparisons(c).baseline-1)',...
        'autoupdate', 'off', 'location', 'best')
    if doDateAxis
        xtickdates(dates(leadTime:end))
    else
        xlim([leadTime T])
        set(gca, 'xtick', xticks);
    end
    ylim([-0.1 max(ylim)])
    
    title(sprintf(' M_%d vs  M_%d', comparisons(c).baseline-1, comparisons(c).alternative-1))
    wrapthisfigure(thisfig,sprintf('sharebiggerthan2-%s%s%s-withtitle', comparisons(c).label, noiselabel, simlabel), wrap)
    
end


%% finish
finishwrap
dockAllFigures
finishscript

