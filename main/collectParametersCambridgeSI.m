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

vintageChoices = {'GDPD', 'CPI'};
v = 1;

switch v
    case 1
        datadir   = 'datGDPD';
    case 2
        datadir   = 'datCPI';
    otherwise
        datadir = pwd;
end

%% setting

doWrap           = false;


%% quantile settings
ndxmean     = 1;
ndxmedian   = 2;

ndxmid      = ndxmedian;


% for tables

ndxtail90   = 2 + [3 8];
ndxtails = ndxtail90;

%% loop over models, vintageChoices and doInflationNoise

for m = 1 : 4 % length(MODELS)
    
    modellabel = MODELS{m};
    
    
    for doInflationNoise = [ true false ]
        
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
        
        titlename = strrep(datalabel, '.', '-');
        if doWrap
            initwrap %#ok<*UNRCH>
        end
        
        tvpA      = contains(modellabel, 'thetaTVP');
        tvpLambda = contains(modellabel, 'lambdaTVP');
        
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
        if ~isempty(wrap)
            copyfile(fullfile(datadir, strcat('settings.', filext)), fullfile(wrap.dir, strcat('settings.', filext)))
            latexwrapper(wrap, 'add', 'listing', strcat('settings.', filext))
        end
        
        %% collect results
        % lambda
        LAMBDA     = importdata(fullfile(datadir, sprintf('LAMBDA.%s', filext)));
        LAMBDAhat  = importdata(fullfile(datadir, sprintf('LAMBDAHAT.%s', filext)));
        
        % a
        A    = loaddat(fullfile(datadir, sprintf('A.%s', filext)));
        Ahat = loaddat(fullfile(datadir, sprintf('AHAT.%s', filext)));
        
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
        
        
        %% create structure to hold posterior moments of parameters
        params = struct('siga', [], 'siglambda', [], 'logSVvol', [], 'lambda', [], 'a', [], 'noisevar', [], 'ndxmid', []);
        
        params.ndxmid = ndxmid;
        
        %% report THETA / Lambda if need be
        if ~tvpA
            params.a.mid   = A(end,ndxmid);
            params.a.tails = A(end,ndxtail90);
        end
        if ~tvpLambda
            params.lambda.mid   = LAMBDA(end,ndxmid);
            params.lambda.tails = LAMBDA(end,ndxtail90);
        end
        
        %% REPORT SCALE PARAMETERS
        % siga
        if ~isempty(SIGA)
            figure
            hold on;
            plot(dates, SIGA(:,ndxmid), 'b', 'linewidth', 3);
            plot(dates, SIGA(:,ndxtails), 'b', 'linewidth', 1);
            ylim([0 max(ylim)])
            nbershades(dates)
            wrapcf(sprintf('SIGA'), wrap);
            
            params.siga.mid   = SIGA(end,ndxmid);
            params.siga.tails = SIGA(end,ndxtail90);
            
        end
        
        % siglambda
        if ~isempty(SIGLAMBDA)
            figure
            hold on;
            plot(dates, SIGLAMBDA(:,ndxmid), 'b', 'linewidth', 3);
            plot(dates, SIGLAMBDA(:,ndxtails), 'b', 'linewidth', 1);
            ylim([0 max(ylim)])
            nbershades(dates)
            wrapcf(sprintf('SIGLAMBDA'), wrap);
            
            params.siglambda.mid   = SIGLAMBDA(end,ndxmid);
            params.siglambda.tails = SIGLAMBDA(end,ndxtail90);
        end
        
        % hInno
        for n = 1 : Nsv
            figure
            hold on;
            plot(dates, HINNO(:,ndxmid,n), 'b', 'linewidth', 3);
            plot(dates, HINNO(:,ndxtails,n), 'b', 'linewidth', 1);
            ylim([0 max(ylim)])
            nbershades(dates)
            wrapcf(sprintf('HINNO%d', n), wrap);
            
            params.logSVvol.mid   = permute(HINNO(end,ndxmid,:), [2 3 1]);
            params.logSVvol.tails = permute(HINNO(end,ndxtail90,:), [2 3 1]);
            
        end
        
        % Sigma
        for n = 1 : Ny
            figure
            hold on;
            plot(dates, SIGMA(:,ndxmid,n), 'b', 'linewidth', 3);
            plot(dates, SIGMA(:,ndxtails,n), 'b', 'linewidth', 1);
            ylim([0 max(ylim)])
            nbershades(dates)
            wrapcf(sprintf('SIGMA%d', n), wrap);
            
            params.noisevar.mid   = permute(SIGMA(end,ndxmid,:), [2 3 1]);
            params.noisevar.tails = permute(SIGMA(end,ndxtail90,:), [2 3 1]);
            
        end
        
        %% finish
        
        finishwrap
        if ~isempty(wrap)
            close all
            matdir = wrap.dir;
        else
            dockAllFigures
            matdir = pwd;
        end
        
        matfilename = strcat('params-', strrep(datalabel, '.', '-'));
        save(fullfile(matdir, matfilename), '-struct', 'params');
        
    end
end

finishscript
