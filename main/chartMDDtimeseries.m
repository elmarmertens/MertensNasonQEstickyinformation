%% collect and plot MDD distributions

initscript
initwrap
wrap = diary2wrap(wrap);
%% settings
MODELS = {'SIthetaCONSTlambdaTVP', ...
    'SIthetaCONSTlambdaCONST', ...
    'SIthetaTVPlambdaTVP', ...
    'SIthetaTVPlambdaCONST', ...
    };

vintagechoice = 'GDPD';
Nparticles = [1e3 1e4 1e5];
labelNparticles = cell(length(Nparticles),1);
for p = 1 : length(Nparticles)
    labelNparticles{p} = sprintf('%d', Nparticles(p));
end

Ngrid      = 250;
datadir    = 'datGDPDstderr';

showESS = false;
xrange        = [-3.4 -2.8];

%% model labels
MODELlabels = cell(length(MODELS),1);
for m = 1 : length(MODELS)
    MODELlabels{m} = sprintf('M_%d', m-1);
end

%% load dates
datasetlabel = sprintf('cambridge2018%s', vintagechoice);
dates        = importdata(sprintf('%s.dates.txt', datasetlabel));
T            = length(dates);

%% allocate memory

dummy     = struct('LOGMDD', NaN(T, Ngrid), 'ESS', NaN(T, Ngrid));
MDDmodels = repmat(dummy, length(MODELS), length(Nparticles), 2); % w/ and w/o noise
clear dummy

%% collect data
for m = 1 : length(MODELS)

    modellabel = MODELS{m};

    for p = 1 : length(Nparticles)
        for i = 1 : 2
            doInflationNoise = i == 1;

            % build filext
            datasetlabel = sprintf('cambridge2018%s', vintagechoice);
            datalabel = sprintf('%s.%sSTDERR', datasetlabel, modellabel);
            if ~doInflationNoise
                datalabel = strcat(datalabel, '.nonoise');
            end

            filext = sprintf('particles.%s.Ngrid%d.Nparticles%d.dat', datalabel, Ngrid, Nparticles(p));

            % load objects
            settingsfile = fullfile(datadir, strcat('settings.', filext));
            if exist(settingsfile, 'file')
                type(settingsfile)
            end

            LOGMDD = loaddat(fullfile(datadir, sprintf('LOGMDD.%s', filext)));
            ESS    = loaddat(fullfile(datadir, sprintf('ESS.%s', filext)));

            if isempty(LOGMDD) || isempty(ESS)
                warning('missing data w/%s', filext)
            end

            % paste into structure
            MDDmodels(m,p,i).LOGMDD = LOGMDD;
            MDDmodels(m,p,i).ESS    = ESS;
        end
    end
end

% save MDD MDDmodels Nparticles MODELS Ngrid datadir vintagechoice

%% allocate memory for LLF
x    = xrange(1) : .001 : xrange(2);
Nx   = length(x);

dummy = struct('llf', NaN(T,Ngrid), 'pdf', NaN(T, Nx), 'mean', NaN(T,1), 'stderr', NaN(T,1), ...
    'meancontrib', NaN(T,1), 'stderrcontrib', NaN(T,1));
LLF   = repmat(dummy, length(MODELS), length(Nparticles), 2); % w/ and w/o noise
clear dummy

%% compute LLF
for m = 1 : length(MODELS)

    modellabel = MODELS{m};

    for p = 1 : length(Nparticles)
        for i = 1 : 2
            doInflationNoise = i == 1;

            % collect LOGMDD
            LOGMDD = MDDmodels(m,p,i).LOGMDD;

            if ~isempty(LOGMDD)

                if any(isnan(LOGMDD(:)))
                    hrulefill
                    if doInflationNoise
                        fprintf('WARNING \n \t %s, Nparticles %d: %d NaN in LOGMDD\n', modellabel, Nparticles(p), sum(isnan(LOGMDD(:))))
                    else
                        fprintf('WARNING \n \t %s (nonoise), Nparticles %d: %d NaN in LOGMDD\n', modellabel, Nparticles(p), sum(isnan(LOGMDD(:))))
                    end
                    hrulefill
                end
                % compute
                %                 logmdd   = bsxfun(@rdivide, cumsum(LOGMDD), (1:T)');
                logmdd   = cumsum(LOGMDD);
                contrib = bsxfun(@rdivide, logmdd, (1:T)');

                pdf   = NaN(T, Nx);
                for t = 1 : T
                    pdf(t,:)  = ksdensity(contrib(t,:), x);
                end

                % paste into structure
                LLF(m,p,i).logmdd    = logmdd;
                LLF(m,p,i).mean   = nanmean(logmdd,2);
                LLF(m,p,i).stderr = nanstd(logmdd,[],2);

                LLF(m,p,i).llfcontrib    = contrib;
                LLF(m,p,i).meancontrib   = nanmean(contrib,2);
                LLF(m,p,i).stderrcontrib = nanstd(contrib,[],2);

                LLF(m,p,i).pdf    = pdf;
            end
        end
    end
end


%% compare MDD over time

LINESTYLES = {'-', '--', '-.', ':'};

p = length(Nparticles);
for i = 1 : 2

    doInflationNoise = i == 1;
    if ~doInflationNoise
        inflationlabel = '-nonoise';
    else
        inflationlabel = '';
    end

    MDDmid = NaN(T,length(MODELS));
    for m= 1 : length(MODELS)
        MDDmid(:,m) = LLF(m,p,i).mean;
    end

    % scale by t
    MDDmid = bsxfun(@rdivide, MDDmid, (1:T)') * T;

    MDDmid(1:8,:) = NaN; % ignore first few obs
    thisfig = figure;
    set(gca, 'fontsize', 12)
    hold on
    for m = 1 : length(MODELS)

        linestyle = LINESTYLES{m};

        plot(dates,MDDmid(:,m), 'linestyle', linestyle, 'linewidth', 2)
    end
    xtickdates(dates)
    legend(MODELlabels, 'location', 'best')
    orient landscape
    wrapthisfigure(thisfig, sprintf('MDDthroughtime%s-particles%d', inflationlabel, Nparticles(p)), wrap)

end % i



%% finish
dockAllFigures
finishwrap
finishscript
