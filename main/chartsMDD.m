%% collect and plot MDD distributions

initscript
initwrap

%% settings
MODELS = {'SIthetaTVPlambdaTVP', ...
    'SIthetaCONSTlambdaTVP', 'SIthetaTVPlambdaCONST' , 'SIthetaCONSTlambdaCONST'};


vintagechoice = 'GDPD';
% vintagechoice = 'CPI';


Nparticles = [1e3 1e4 1e5];
labelNparticles = cell(length(Nparticles),1);
for p = 1 : length(Nparticles)
    labelNparticles{p} = sprintf('%3d,000', Nparticles(p)*1e-3);
end

noiselabel = {'','-nonoise'};

Ngrid      = 250;

yrange        = [0 0.1];

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
            
            
            switch upper(vintagechoice)
                case 'GDPD'
                    datadir    = 'datGDPDstderr';
                    
                    if i == 2
                        xrange = [-700 -525];
                    else
                        xrange = [-650 -500];
                    end
                    matfilename = 'cambridgeMDD';
                case 'CPI'
                    datadir    = 'datCPIstderr';
                    xrange     = [-525 -425];
                    matfilename = 'cambridgeMDD-CPI';
                otherwise
                    error houston
            end
            
            
            
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
            else
                error houston
            end
            
            LOGMDD = loaddat(fullfile(datadir, sprintf('LOGMDD.%s', filext)));
            ESS    = loaddat(fullfile(datadir, sprintf('ESS.%s', filext)));
            
            % paste into structure
            MDDmodels(m,p,i).LOGMDD = LOGMDD;
            MDDmodels(m,p,i).ESS    = ESS;
        end
    end
end

% save MDD MDDmodels Nparticles MODELS Ngrid datadir vintagechoice

%% allocate memory for MDD
x    = xrange(1) : .001 : xrange(2);
Nx   = length(x);

dummy = struct('logmdd', NaN(T,Ngrid), 'pdf', NaN(T, Nx), 'mean', NaN(T,1), 'stderr', NaN(T,1), ...
    'meancontrib', NaN(T,1), 'stderrcontrib', NaN(T,1));
MDD   = repmat(dummy, length(MODELS), length(Nparticles), 2); % w/ and w/o noise
clear dummy

%% compute MDD
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
                t = T;
                pdf(t,:)  = ksdensity(logmdd(t,:), x);
                
                % paste into structure
                MDD(m,p,i).logmdd    = logmdd;
                MDD(m,p,i).mean   = nanmean(logmdd,2);
                nobs = sum(~isnan(logmdd),2);
                MDD(m,p,i).stderr = nanstd(logmdd,[],2) ./ sqrt(nobs);
                
                MDD(m,p,i).pdf    = pdf;
            end
        end
    end
end

%% plot MDD per model, single panel


colors = {[0 0 1], [1 0 0], [0 0 0]};
close all
fontsize = 20;

for m = 1 : length(MODELS)
    
    modellabel = MODELS{m};
    switch upper(vintagechoice)
        case 'CPI'
            modellabel = strcat(modellabel, '-CPI');
        otherwise % do nothing
    end
    for i = 1 : 2
        
        thisfig = figure;
        set(gcf, 'defaultLegendAutoUpdate','off')
        set(gca, 'fontsize', fontsize)
        
        hold on
        % loop over Nparticles
        hanni = NaN(length(Nparticles),1);
        for p = 1 : length(Nparticles)
            switch p
                case 1
                    linestyle = ':';
                case 3
                    linestyle = '-';
                case 2
                    linestyle = '--';
            end
            hanni(p) = plot(x, MDD(m,p,i).pdf(end,:), ...
                'linewidth', 3, 'linestyle', linestyle, 'color', colors{p});
            
        end
        % mean
        for p = length(Nparticles)
            
            % values
            mddhat = MDD(m,p,i).mean(end);
            mddhi  = MDD(m,p,i).mean(end) + 2 * MDD(m,p,i).stderr(end);
            mddlo  = MDD(m,p,i).mean(end) - 2 * MDD(m,p,i).stderr(end);
            
            % mean
            hmid = plot(repmat(mddhat,2,1), ylim, ...
                'linewidth', 2, 'linestyle', ':', 'color', colors{p});
            ci = x >= mddlo & x < mddhi;
            if isequal(colors{p}, [0 0 0])
                shadecolor = .75 * [1 1 1];
            else
                shadecolor = .75 * colors{p};
            end
            cishades(1) = bar(x, max(ylim) * ci , 1, 'EdgeColor', shadecolor, 'FaceColor', shadecolor);
            cishades(2) = bar(x, min(ylim) * ci , 1, 'EdgeColor', shadecolor, 'FaceColor', shadecolor);
            uistack(cishades, 'bottom')
            set(gca, 'layer', 'top')
        end
        set(gca, 'xtick', -750 : 25 : -400)
        xlim(xrange)
        
        hl = legend([hanni; cishades(1)], labelNparticles{:}, 'mean +/- (2 x NSE) w/ 100,000', 'location', 'best');
        legend('boxoff')
        
        wrapthisfigure(thisfig,sprintf('figMDDsim-%s%s', modellabel, noiselabel{i}), wrap)
        delete(hl)
        wrapthisfigure(thisfig,sprintf('figMDDsim-%s%s-nolegend', modellabel, noiselabel{i}), wrap)
        
    end % noise choice
end


%% finish
varlist = {'vintagechoice', 'Nparticles', 'labelNparticles', 'MODELS', 'MDD', 'xrange', 'x'};
save(matfilename, varlist{:}, '-v7.3');

dockAllFigures
finishwrap
finishscript
