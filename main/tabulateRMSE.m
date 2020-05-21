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
    };


vintageChoices = {'GDPD', 'CPI'};

dofAPFfattail = 0;
lambdaN       = 1;


%% setting

doWrap        = true;

if doWrap
    initwrap
end

%% quantile settings
ndxmean     = 1;
ndxmedian   = 2;

ndxmid      = ndxmedian;

ndxtailIQR    = 2 + [5 6];
% ndxtails    = 2 + [3 8];
% fractiles = [0.005, 0.025, 0.05, .1586, .25, .75, .8413, 0.95, 0.975, 0.995] * 100;

% for tables

ndxtail90   = 2 + [3 8];
ndxtails    = ndxtail90;
ndxtails9068  = 2 + [3 4 7 8];

Nsurveys = 5;

SurveyLabels = cell(Nsurveys,1);
for n = 1 : Nsurveys
    SurveyLabels{n} = sprintf('h = %d', n);
end

sam(1).start = [];
sam(1).end   = [];

s = 2;
sam(s).start = [];
sam(s).end   = datenum(1984,12,1);

s = 3;
sam(s).start = datenum(1985,1,1);
sam(s).end   = datenum(2000,1,1);

s = 4;
sam(s).start = datenum(2000,1,1);
sam(s).end   = [];

panelNums = {'a', 'b', 'c', 'd', 'e'};

%% loop over vintageChoices

for v = 1 : 2

    %% parameters
    switch v
        case 1
            datadir   = 'datGDPD';
        case 2
            datadir   = 'datCPI';
        otherwise
            datadir = pwd;
    end



    %% loop over models, vintageChoices and doInflationNoise

    for doInflationNoise = [true false]

        %% set up table
        if doInflationNoise
            filename = 'relrmseREvsSPF';
        else
            filename = 'relrmseREvsSPF-nonoise';
        end
        if v == 2 % CPI
            filename = strcat(filename, '-CPI');
        end

        filename = strcat(filename, '.tex');
        fid = fopen(fullfile(wrap.dir, filename), 'wt');

        Ntabcol = 1 + 1 + 2 * length(MODELS);

        fprintf(fid, '\\begin{center}\n');
        fprintf(fid, '\\begin{tabular}{c%s}\n', repmat('.4', 1, Ntabcol - 1));
        fprintf(fid, '\\toprule\n');
        % table header
        fprintf(fid, ' & & \\multicolumn{%d}{c}{RE model forecasts} & \\multicolumn{%d}{c}{SI model forecasts} \\\\  \n', length(MODELS), length(MODELS));
        fprintf(fid, ' & \\ccol{SPF} & \\multicolumn{%d}{c}{(rel. RMSE)} & \\multicolumn{%d}{c}{(rel. RMSE)} \\\\ \\cmidrule(r){3-%d}\\cmidrule(l){%d-%d} \n', length(MODELS),length(MODELS), 2+length(MODELS), 2+length(MODELS)+1, 2+length(MODELS)+length(MODELS));

        fprintf(fid, '\\ccol{horizon} ');
        fprintf(fid, '& \\ccol{(RMSE)}');
        for i = 1 : length(MODELS)
            fprintf(fid, '& \\ccol{$\\mathcal{M}_{%d}$}', i-1);
        end
        for i = 1 : length(MODELS)
            fprintf(fid, '& \\ccol{$\\mathcal{M}_{%d}$}', i-1);
        end
        fprintf(fid, '\\\\\n');

        for s = 1 : length(sam)

            %% prepare data matrices for table output
            tabledata.relRMSEuc  = NaN(Nsurveys,length(MODELS));
            tabledata.relRMSEsi  = NaN(Nsurveys,length(MODELS));

            for m = 1 : length(MODELS)

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


                timestamp   = [];


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

                if isempty(sam(s).start)
                    samNdx = true(length(dates), 1);
                else
                    samNdx = (dates >= sam(s).start);
                end
                if ~isempty(sam(s).end)
                    samNdx = samNdx & (dates <= sam(s).end);
                end

                theseDates = dates(samNdx);
                sampleTxt = sprintf('%s to %s', datestr(theseDates(1), 'yyyy:qq'), datestr(theseDates(end), 'yyyy:qq'));

                %% read results
                T   = size(y,1);
                Ny  = size(y,2);
                Nstates = 4;
                Nsv     = 2;


                type(fullfile(datadir, strcat('settings.', filext)));
                %     if ~isempty(wrap)
                %         copyfile(fullfile(datadir, strcat('settings.', filext)), fullfile(wrap.dir, strcat('settings.', filext)))
                %         latexwrapper(wrap, 'add', 'listing', strcat('settings.', filext))
                %     end

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



                % inflation forecasts
                PIhatRE = loaddat(fullfile(datadir, sprintf('PIHAT.%s', filext)));
                PIhatSI = loaddat(fullfile(datadir, sprintf('PIFCST.%s', filext)));



                %% Forecast comparison
                % compute forecast errors
                fcsterr = NaN(T,Nsurveys);
                for t = 1 : T
                    for n = 1 : Nsurveys
                        if t + n <= T
                            fcsterr(t,n) = PIhatRE(t,n) - y(t+n,1);
                        end
                    end
                end


                SPFerr  = NaN(T,Nsurveys);
                for t = 1 : T
                    for n = 1 : Nsurveys
                        if t + n <= T
                            SPFerr(t,n) = y(t,1+n) - y(t+n,1);
                        end
                    end
                end

                SIerr  = NaN(T,Nsurveys);
                for t = 1 : T
                    for n = 1 : Nsurveys
                        if t + n <= T
                            SIerr(t,n) = PIhatSI(t,n) - y(t+n,1);
                        end
                    end
                end


                % cut sample
                fcsterr(~samNdx,:) = NaN;
                SPFerr(~samNdx,:)  = NaN;
                SIerr(~samNdx,:)   = NaN;

                ucRMSE  = sqrt(nanmean(fcsterr.^2));
                spfRMSE = sqrt(nanmean(SPFerr.^2));
                siRMSE  = sqrt(nanmean(SIerr.^2));

                ucRMSEcum  = sqrt(nancummean(fcsterr.^2));
                spfRMSEcum = sqrt(nancummean(SPFerr.^2));
                siRMSEcum  = sqrt(nancummean(SIerr.^2));

                %% run DM regressions



                for n = 1 : Nsurveys
                    if m == 1
                        tabledata.spfRMSE(n)  = spfRMSE(n);
                    end

                    tabledata.relRMSEuc(n, m)  = ucRMSE(n) / spfRMSE(n);
                    tabledata.relRMSEsi(n, m)  = siRMSE(n) / spfRMSE(n);

                end

            end % model m


            %% write results into table

            fprintf(fid, '\\midrule\n');

            fprintf(fid, '\\multicolumn{%d}{c}{\\textsc{Panel (%s)}: %s}', Ntabcol, panelNums{s}, sampleTxt);

            fprintf(fid, '\\\\\n');
            fprintf(fid, '\\midrule\n');



            for n = 1 : Nsurveys
                fprintf(fid, '%d ', n);
                fprintf(fid, ' & \\ccol{ %5.2f} ', tabledata.spfRMSE(n));
                fprintf(fid, ' & \\ccol{ %5.2f} ', tabledata.relRMSEuc(n,:));
                fprintf(fid, ' & \\ccol{ %5.2f} ', tabledata.relRMSEsi(n,:));
                %         for m = 1 : length(MODELS)
                %             fprintf(fid, ' & %5.2f%s ', tabledata.relRMSE(n,m), Zstar(tabledata.dmTstat(n,m)));
                %         end
                fprintf(fid, '\\\\\n');
            end


        end % sam

        %% close table
        fprintf(fid, '\\bottomrule\n');
        fprintf(fid, '\\end{tabular}\n');
        fprintf(fid, '\\end{center}\n');
        fprintf(fid, '\n');
        fprintf(fid, '\n');

        fprintf(fid, 'Note: RMSE of SPF forecasts and relative RMSE of RE and SI predictions generated by a given model compared to the SPF forecast. ');
        fprintf(fid, '(Numbers below one indicate a lower RMSE of the model forecasts). ');
        fprintf(fid, 'In each panel, model forecasts used are based on filtered estimates using data since %s. Forecast errors are then collected over the (sub)periods indicated in each panel.\n', datestr(dates(1), 'yyyy:qq'));



        fclose(fid);
        type(fullfile(wrap.dir, filename));
        latexwrapper(wrap, 'add', 'tab', filename)

    end % inflationNoise
end % vintageChoice

%% finish
finishwrap
finishscript
dockAllFigures
