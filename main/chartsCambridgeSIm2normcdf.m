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


MODELS = {'SIthetaTVPlambdaTVPnormcdf'};


vintageChoices = {'GDPD', 'CPI'};
v = 2;



switch v
    case 1
        datadir   = 'datGDPD';
    case 2
        datadir   = 'datCPI';
    otherwise
        datadir = pwd;
end

picdir    = []; % fullfile(localtemp, 'foo');

%% setting
doWrap           = false;

fontsize = 10;

%% loop over models, vintageChoices and doInflationNoise

for m = 1
    for doInflationNoise = [true false]
        
        
        modellabel  = MODELS{m};
        
        
        
        vintagechoice = vintageChoices{v};
        
        
        %% process settings
        datasetlabel = sprintf('cambridge2018%s', vintagechoice);
        
        
        datalabel = sprintf('%s.%s', datasetlabel, modellabel);
        if ~doInflationNoise
            datalabel = strcat(datalabel, '.nonoise');
        end
        
        % append CPI to modellabel if needed
        if v == 2
            modellabel = strcat(modellabel, '-CPI');
        end
        
        Nsurveys = 5;
        
        SurveyLabels = cell(Nsurveys,1);
        for n = 1 : Nsurveys
            SurveyLabels{n} = sprintf('h = %d', n);
        end
        
        timestamp   = [];
        
        titlename = strrep(datalabel, '.', '-');
        
        tvpTheta  = contains(modellabel, 'thetaTVP');
        tvpLambda = contains(modellabel, 'lambdaTVP');
        
        
        %% get parameters
        filext = sprintf('particles.%s.dat', datalabel);
        
        
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
        
        TAURE = importdata(fullfile(datadir, sprintf('TAURE.%s', filext)));
        TAUSI = importdata(fullfile(datadir, sprintf('TAUSI.%s', filext)));
        GAPRE = importdata(fullfile(datadir, sprintf('GAPRE.%s', filext)));
        GAPSI = importdata(fullfile(datadir, sprintf('GAPSI.%s', filext)));
        
        INFLATIONRE = loaddat(fullfile(datadir, sprintf('INFLATIONRE.%s', filext)));
        INFLATIONSI = loaddat(fullfile(datadir, sprintf('INFLATIONSI.%s', filext)));
        
        % inflation forecasts
        PIhatRE = loaddat(fullfile(datadir, sprintf('PIHAT.%s', filext)));
        PIhatSI = loaddat(fullfile(datadir, sprintf('PIFCST.%s', filext)));
        PIpersistence = loaddat(fullfile(datadir, sprintf('PIPERSISTENCE.%s', filext)));
        TRENDpersistence = loaddat(fullfile(datadir, sprintf('TRENDPERSISTENCE.%s', filext)));
        GAPpersistence = loaddat(fullfile(datadir, sprintf('GAPPERSISTENCE.%s', filext)));
        
        % loglike and ESS
        % LOGMDD         = importdata(fullfile(datadir, sprintf('LOGMDD.%s', filext)));
        % ESS             = importdata(fullfile(datadir, sprintf('ESS.%s', filext)));
        % PARTICLEWEIGHTS = loaddat(fullfile(datadir, sprintf('PARTICLEWEIGHTS.%s', filext)));
        
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
        
        % theta
        THETA    = loaddat(fullfile(datadir, sprintf('A.%s', filext)));
        THETAhat = loaddat(fullfile(datadir, sprintf('AHAT.%s', filext)));
        
        %% load SCALE PARAMETERS
        
        % siga
        SIGA    = loaddat(fullfile(datadir, sprintf('SIGA.%s', filext)));
        SIGAhat = loaddat(fullfile(datadir, sprintf('SIGAHAT.%s', filext)));
        
        % siglambda
        % SIGLAMBDA    = loaddat(fullfile(datadir, sprintf('SIGLAMBDA.%s', filext)));
        % SIGLAMBDAhat = loaddat(fullfile(datadir, sprintf('SIGLAMBDAHAT.%s', filext)));
        XLAMBDA_SIGMA    = loaddat(fullfile(datadir, sprintf('XLAMBDA_SIGMA.%s', filext)));
        XLAMBDA_SIGMAhat = loaddat(fullfile(datadir, sprintf('XLAMBDA_SIGMAHAT.%s', filext)));

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
        % if showSmoother
        %     XSIGall         = loaddat(fullfile(datadir, sprintf('XSIGall.%s', filext)));
        %     XSIGsrv         = loaddat(fullfile(datadir, sprintf('XSIGsrv.%s', filext)));
        %     XSIGinf         = loaddat(fullfile(datadir, sprintf('XSIGinf.%s', filext)));
        %     XSIGinfspflong  = loaddat(fullfile(datadir, sprintf('XSIGinfspflong.%s', filext)));
        %     XSIGspflong     = loaddat(fullfile(datadir, sprintf('XSIGspflong.%s', filext)));
        %
        %     smootherPIpersistence = loaddat(fullfile(datadir, sprintf('smootherPIpersistence.%s', filext)));
        %     smootherTRENDpersistence = loaddat(fullfile(datadir, sprintf('smootherTRENDpersistence.%s', filext)));
        %     smootherGAPpersistence = loaddat(fullfile(datadir, sprintf('smootherGAPpersistence.%s', filext)));
        %
        %     STATEa       = loaddat(fullfile(datadir, sprintf('STATEa.%s', filext)));
        %     STATElambda  = loaddat(fullfile(datadir, sprintf('STATElambda.%s', filext)));
        %     STATEsvol    = loaddat(fullfile(datadir, sprintf('STATEsvol.%s', filext)));
        %     STATEsvol    = STATEsvol';
        % end
        
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
            
        else
            LAMBDAsmoother = [];
            deltaLAMBDAsmoother = [];
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
        
        if ~doInflationNoise
            modellabel = strcat(modellabel, '-nonoise');
        end
        
        %% FIGURE 2: static parameters
        
        % datalabels = {'Mid', '5', '15.86', '84.13', '95'};
        
        figname = strcat('figSTATICPARAMETER-', modellabel);
        figure
        set(gcf, 'defaultLegendAutoUpdate','off')
        
        % panel a: vol of trend SV
        n = 1;
        panel = 1;
        subplot(2,2,panel)
        set(gca, 'fontsize', fontsize)
        hold on
        plotCI(HINNO(:,ndxmid,n), HINNO(:,ndxtails,n), dates, 0)
        ylim([0 .2])
        set(gca, 'ytick', 0 : .05 : 2)
        set(gca, 'xtick', xtickdates)
        title('(a) PLE Path of Shock Variance Parameter of ln \varsigma^2_{\eta,t}, \sigma^2_\eta')
        nberlines(xtickdates)
        
        % panel b: vol of trend SV
        n = 2;
        panel = 2;
        subplot(2,2,panel)
        set(gca, 'fontsize', fontsize)
        hold on
        plotCI(HINNO(:,ndxmid,n), HINNO(:,ndxtails,n), dates, 0)
        ylim([0 .2])
        set(gca, 'ytick', 0 : .05 : 2)
        set(gca, 'xtick', xtickdates)
        title('(b) PLE Path of Shock Variance Parameter of ln \varsigma^2_{\nu,t}, \sigma^2_\nu')
        nberlines(xtickdates)
        
        if tvpTheta
            % panel c: gap AR1 coeffcient
            panel = 3;
            subplot(2,2,panel)
            set(gca, 'fontsize', fontsize)
            hold on
            plotCI(SIGA(:,ndxmid), SIGA(:,ndxtails), dates, -1)
            ylim([0 .05])
            set(gca, 'ytick', 0 : .01 : 1)
            set(gca, 'xtick', xtickdates)
            title('(c) PLE Path of Shock Variance Parameter of Gap AR(1) coefficient \theta_{t}')
            nberlines(xtickdates)
            
        else
            % panel c:  AR1 coeffcient
            panel = 3;
            subplot(2,2,panel)
            set(gca, 'fontsize', fontsize)
            hold on
            plotCI(THETA(:,ndxmid), THETA(:,ndxtails), dates, -1)
            ylim([0 1])
            set(gca, 'ytick', -1 : .2 : 1)
            set(gca, 'xtick', xtickdates)
            title('(c) PLE Path of  Gap AR(1) coefficient \theta')
            nberlines(xtickdates)
            
            
        end
        
        
        if tvpLambda % panel d: vol of lambda
            panel = 4;
            subplot(2,2,panel)
            set(gca, 'fontsize', fontsize)
            hold on
            plotCI(XLAMBDA_SIGMA(:,ndxmid), XLAMBDA_SIGMA(:,ndxtails), dates, 0)
            set(gca, 'xtick', xtickdates)
            title('(d) PLE Path of Shock Variance Parameter of \lambda_{t}, \sigma^2_\kappa')
            nberlines(xtickdates)
            
            
        else
            panel = 4;
            subplot(2,2,panel)
            set(gca, 'fontsize', fontsize)
            hold on
            plotCI(LAMBDA(:,ndxmid), LAMBDA(:,ndxtails), dates, 0)
            ylim([0 1])
            set(gca, 'ytick', 0 : .2 : 1)
            set(gca, 'xtick', xtickdates)
            title('(d) PLE Path of \lambda')
            nberlines(xtickdates)
            
        end
        
        % panels complete, now store the figure
        orient landscape
        print(figname, '-fillpage', '-dpdf')
        
        %% FIG3
        
        
        figname = strcat('figTRENDandGAP-', modellabel);
        
        figure
        set(gcf, 'defaultLegendAutoUpdate','off')
        
        % panel a: SI Trend and SPF Nowcast
        panel = 1;
        subplot(2,2,panel)
        set(gca, 'fontsize', fontsize)
        hold on
        plotCI(TAUSI(:,ndxmid), TAUSI(:,ndxtail68), dates, 0, 'k-.', 'linewidth', 2)
        h2 = plot(dates, y(:,2), 'r-', 'linewidth', 1);
        h1 = plot(dates, TAUSI(:,ndxmid), 'k-.', 'linewidth', 2);
        ylim([0 10])
        set(gca, 'ytick', 0 : 2 : 20)
        set(gca, 'xtick', xtickdates)
        title('(a) Filtered SI Trend and SPF Inflation Nowcast')
        legend([h1 h2], 'Filtered SI Trend: F_{t|t} \tau_t', 'SPF Nowcast: \pi_{t,t+1}^{SPF}', 'location', 'northeast')
        nberlines(xtickdates)

        
        % panel b: SI Trend and SPF-4Q
        panel = 2;
        subplot(2,2,panel)
        set(gca, 'fontsize', fontsize)
        hold on
        plotCI(TAUSI(:,ndxmid), TAUSI(:,ndxtail68), dates, 0, 'k-.', 'linewidth', 2)
        h2 = plot(dates, y(:,6), 'r-', 'linewidth', 1);
        h1 = plot(dates, TAUSI(:,ndxmid), 'k-.', 'linewidth', 2);
        ylim([0 10])
        set(gca, 'ytick', 0 : 2 : 20)
        set(gca, 'xtick', xtickdates)
        title('(b) Filtered SI Trend and 4-Quarter Ahead SPF Inflation Prediction')
        legend([h1 h2], 'Filtered SI Trend: F_{t|t} \tau_t', '4-Q Ahead SPF: \pi_{t,t+5}^{SPF}', 'location', 'northeast')
        nberlines(xtickdates)
        
        
        % panel c: Inflation, SI and RE Trend
        panel = 3;
        subplot(2,2,panel)
        set(gca, 'fontsize', fontsize)
        hold on
        h1 = plot(dates, y(:,1), 'b-.', 'linewidth', 1);
        h2 = plot(dates, TAUSI(:,ndxmid), 'k-', 'linewidth', 2);
        h3 = plot(dates, TAURE(:,ndxmid), 'r--', 'linewidth', 2);
        ylim([0 14])
        set(gca, 'ytick', 0 : 2 : 20)
        set(gca, 'xtick', xtickdates)
        title('(c) Realized Inflation, Filtereed SI Trend and Filtered RE Trend')
        legend([h1 h2 h3], 'Realized Inflation: \pi_t', 'Filtered SI Trend: F_{t|t} \tau_t', 'Filtered RE Trend: \tau_{t|t}')
        nberlines(xtickdates)
        
        
        % panel d: Gap
        panel = 4;
        subplot(2,2,panel)
        set(gca, 'fontsize', fontsize)
        hold on
        h1 = plot(dates, GAPSI(:,ndxmid), 'k-', 'linewidth', 1);
        h2 = plot(dates, GAPRE(:,ndxmid), 'r--', 'linewidth', 1);
        ylim([-4 14])
        set(gca, 'ytick', -20 : 2 : 20)
        set(gca, 'xtick', xtickdates)
        title('(d) Filtered SI Gap Inflation and Filtered RE Gap Inflation')
        legend([h1 h2], 'Filtered SI Gap: F_{t|t} \epsilon_t', 'Filtered RE Gap: \epsilon_{t|t}')
        nberlines(xtickdates)
        
       
        % panels complete, now store the figure
        orient landscape
        print(figname, '-fillpage', '-dpdf')
        
        
        %% FIG4
        
        figname = strcat('figSV-', modellabel);
        
        
        figure
        set(gcf, 'defaultLegendAutoUpdate','off')
        
        % panel a: Filtered Trend SV
        n = 1;
        panel = 1;
        subplot(2,2,panel)
        set(gca, 'fontsize', fontsize)
        hold on
        plot(dates, SV(:,ndxmid,n), 'r:', 'linewidth', 2);
        plot(dates, SV(:,ndxtail90,n), 'k-', 'linewidth', 1);
        ylim([0 max(ylim)])
        % set(gca, 'ytick', 0 : .2 : 5)
        set(gca, 'xtick', xtickdates)
        title('(a) Filtered Trend SV, \varsigma_{\eta,t|t}, and 90% Uncertainty Bands')
        nberlines(xtickdates)
        
       
        % panel b: Smoothed Trend SV
        n = 1;
        panel = 2;
        subplot(2,2,panel)
        set(gca, 'fontsize', fontsize)
        hold on
        plot(dates, SVsmoother(:,ndxmid,n), 'r:', 'linewidth', 2);
        plot(dates, SVsmoother(:,ndxtail90,n), 'k-', 'linewidth', 1);
        ylim([0 max(ylim)])
        % set(gca, 'ytick', 0 : .2 : 5)
        set(gca, 'xtick', xtickdates)
        title('(b) Smoothed Trend SV, \varsigma_{\eta,t|T}, and 90% Uncertainty Bands')
        nberlines(xtickdates)
                
        % panel c: Filtered Gap SV
        n = 2;
        panel = 3;
        subplot(2,2,panel)
        set(gca, 'fontsize', fontsize)
        hold on
        plot(dates, SV(:,ndxmid,n), 'r:', 'linewidth', 2);
        plot(dates, SV(:,ndxtail90,n), 'k-', 'linewidth', 1);
        % ylim([0 3.5])
        % set(gca, 'ytick', 0 : .5 : 20)
        set(gca, 'xtick', xtickdates)
        title('(c) Filtered Gap SV, \varsigma_{\nu,t|t}, and 90% Uncertainty Bands')
        nberlines(xtickdates)
                
        % panel d: Smoothed Gap SV
        n = 2;
        panel = 4;
        subplot(2,2,panel)
        set(gca, 'fontsize', fontsize)
        hold on
        plot(dates, SVsmoother(:,ndxmid,n), 'r:', 'linewidth', 2);
        plot(dates, SVsmoother(:,ndxtail90,n), 'k-', 'linewidth', 1);
        set(gca, 'xtick', xtickdates)
        title('(b) Smoothed Gap SV, \varsigma_{\nu,t|T}, and 90% Uncertainty Bands')
        nberlines(xtickdates)
        
        
        % panels complete, now store the figure
        orient landscape
        print(figname, '-fillpage', '-dpdf')
        
        %% FIG 5 -- lambda
        
        if tvpLambda
            figname = strcat('figLAMBDA-', modellabel);
            
            figure
            set(gcf, 'defaultLegendAutoUpdate','off')
            
            % panel a: baseline lambda
            panel = 1;
            subplot(2,1,panel)
            set(gca, 'fontsize', fontsize)
            hold on
            plotCI(LAMBDA(:,ndxmid), LAMBDA(:,ndxtail90), dates, 0, 'k-.', 'linewidth', 2)
            plot(dates, LAMBDAsmoother(:,ndxmid), 'r-.', 'linewidth', 2);
            plot(dates, LAMBDAsmoother(:,ndxtail90), 'r-', 'linewidth', 1);
            ylim([0 1])
            set(gca, 'ytick', 0 : .1 : 1)
            set(gca, 'xtick', xtickdates)
            title({'(a) Filtered and Smoothed SI Updating, \lambda_{t|t} and \lambda_{t|T}'})
            nberlines(xtickdates)
                        
            % panel b: baseline delta-lambda
            panel = 2;
            subplot(2,1,panel)
            set(gca, 'fontsize', fontsize)
            hold on
            plotCI(deltaLAMBDAsmoother(:,ndxmid), deltaLAMBDAsmoother(:,ndxtails), dates, -1, 'k-.', 'linewidth', 2)
            plot(xtickdates([1 end]), [0 0], 'k-')
            ylim([-1 1])
            set(gca, 'ytick', -1 : .2 : 1)
            set(gca, 'xtick', xtickdates)
            title({'(b) Accumulated Changes in Smoothed SI Updated, \lambda_{t|T} - \lambda_{1|T}'})
            nberlines(xtickdates)
            
            % panels complete, now store the figure
            orient landscape
            print(figname, '-fillpage', '-dpdf')
            
        end
        
        %% FIG 5 -- lambda
        
        if tvpTheta
            figname = strcat('figTHETA-', modellabel);
            
            figure
            set(gcf, 'defaultLegendAutoUpdate','off')
            
            % panel a: baseline lambda
            panel = 1;
            %     subplot(2,1,panel)
            set(gca, 'fontsize', fontsize)
            hold on
            plotCI(THETA(:,ndxmid), THETA(:,ndxtail90), dates, -1, 'k-.', 'linewidth', 2)
            plot(dates, THETAsmoother(:,ndxmid), 'r-.', 'linewidth', 2);
            plot(dates, THETAsmoother(:,ndxtail90), 'r-', 'linewidth', 1);
            ylim([-1 1])
            set(gca, 'ytick', -1 : .1 : 1)
            set(gca, 'xtick', xtickdates)
            title({'Filtered and Smoothed: \theta_{t|t} and \theta_{t|T}'})
            nberlines(xtickdates)
                        
            % panels complete, now store the figure
            orient landscape
            print(figname, '-fillpage', '-dpdf')
            
        end
        
    end
end

%% finish
finishwrap
dockAllFigures
finishscript
