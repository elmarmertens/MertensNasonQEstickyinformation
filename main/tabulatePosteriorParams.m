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

matdir   = pwd;

%% setting
vintageChoices = {'GDPD', 'CPI'};
v = 2;

for doInflationNoise = [true false]
    
    %% allocate memory for results
    dummy  = struct('siga', [], 'siglambda', [], 'logSVvol', [], 'lambda', [], 'a', [], 'noisevar', [], 'ndxmid', []);
    params = repmat(dummy, length(MODELS), 1);
    
    
    %% collect MDD
    switch v 
        case 1
            MDDdat = load('cambridgeMDD.mat');
        case 2
            MDDdat = load('cambridgeMDD-CPI.mat');
        otherwise
            error houston
    end
    
    MDDmeans  = NaN(1, length(MODELS));
    MDDstderr = NaN(1, length(MODELS));
    
    % select appropriate MDD
    p = 3;
    %#ok<*UNRCH>
    if doInflationNoise
        i = 1;
    else
        i = 2;
    end
    
    for m = 1 : 4
        modellabel = MODELS{m};
        thisModel = ismember(MDDdat.MODELS, modellabel);
        MDDmeans(m) = MDDdat.MDD(thisModel,p,i).mean(end);
        MDDstderr(m) = MDDdat.MDD(thisModel,p,i).stderr(end);
        
    end
    
    %% loop over models, vintageChoices and doInflationNoise
    
    for m = 1 : 4 % : length(MODELS)
        
        modellabel = MODELS{m};
        vintagechoice = vintageChoices{v};
        
        
        %% process settings
        datasetlabel = sprintf('cambridge2018%s', vintagechoice);
        
        
        datalabel = sprintf('%s.%s', datasetlabel, modellabel);
        if ~doInflationNoise
            datalabel = strcat(datalabel, '.nonoise');
        end
        
        matfilename = strcat('params-', strrep(datalabel, '.', '-'));
        
        
        %         tvpA      = strcmpi(modellabel, 'thetaTVP');
        %         tvpLambda = strcmpi(modellabel, 'lambdaTVP');
        
        jack = load(fullfile(matdir, matfilename));
        params(m) = jack;
        clear jack
        
    end
    
    %% create table
    if doInflationNoise
        tabname = 'posteriorMoments';
    else
        tabname = 'posteriorMomentsNoNoise';
    end
    if v == 2
        tabname = strcat(tabname, '-CPI');
    end    
    tabname = strcat(tabname, '.tex');
    fid = fopen(tabname, 'wt');
    fprintf(fid, '\\begin{tabular}{c%s}\n', repmat('c', 1, length(MODELS)));
    fprintf(fid, '\\toprule\n');
    fprintf(fid, ' & \\multicolumn{%d}{c}{Models} \\\\ \\cmidrule{2-%d} \n', length(MODELS), 1+length(MODELS));
    fprintf(fid, 'Parameter ');
    for i = 1 : length(MODELS)
        fprintf(fid, '& $\\mathcal{M}_{%d}$', i-1);
    end
    fprintf(fid, '\\\\\n');
    
    fprintf(fid, '\\midrule\n');
    fprintf(fid, '\\multicolumn{5}{c}{Variances of shocks to SV processes}\\\\\n');
    fprintf(fid, '\\midrule\n');

    % sigSVtrend
    ndxSV = 1;
    fprintf(fid, '$\\sigma_{\\eta}^2$ (Trend SV)');
    for i = 1 : length(MODELS)
        fprintf(fid, ' & %6.3f ', params(i).logSVvol.mid(ndxSV));
    end
    fprintf(fid, '\\\\\n');
    for i = 1 : length(MODELS)
        fprintf(fid, ' & $\\bigl[ %6.3f,%6.3f \\bigr]$ ', params(i).logSVvol.tails(1,ndxSV), params(i).logSVvol.tails(2,ndxSV));
    end
    fprintf(fid, '\\vspace{\\skipper}\\\\\n');
    
    % sigSVgap
    ndxSV = 2;
    fprintf(fid, '$\\sigma_{\\upsilon}^2$ (Gap SV)');
    for i = 1 : length(MODELS)
        fprintf(fid, ' & %6.3f ', params(i).logSVvol.mid(ndxSV));
    end
    fprintf(fid, '\\\\\n');
    for i = 1 : length(MODELS)
        fprintf(fid, ' & $\\bigl[ %6.3f,%6.3f \\bigr]$ ', params(i).logSVvol.tails(1,ndxSV), params(i).logSVvol.tails(2,ndxSV));
    end
    fprintf(fid, '\\vspace{\\skipper}\\\\\n');
    
    fprintf(fid, '\\midrule\n');
    fprintf(fid, '\\multicolumn{5}{c}{Persistence of inflation gap}\\\\\n');
    fprintf(fid, '\\midrule\n');

    % THETA
    fprintf(fid, '$\\theta$');
    for i = 1 : length(MODELS)
        this = params(i).a;
        if isempty(this)
            fprintf(fid, ' & -- ');
        else
            fprintf(fid, ' & %6.3f ', this.mid);
        end
    end
    fprintf(fid, '\\\\\n');
    
    for i = 1 : length(MODELS)
        this = params(i).a;
        if isempty(this)
            fprintf(fid, ' & ');
        else
            fprintf(fid, ' & $\\bigl[ %6.3f,%6.3f \\bigr]$ ', this.tails(1), this.tails(2));
        end
    end
    fprintf(fid, '\\vspace{\\skipper}\\\\\n');
    
    % SIGTHETA
    fprintf(fid, '$\\sigma_{\\phi}^2$');
    for i = 1 : length(MODELS)
        this = params(i).siga;
        if isempty(this)
            fprintf(fid, ' & -- ');
        else
            fprintf(fid, ' & %6.3f ', this.mid);
        end
    end
    fprintf(fid, '\\\\\n');
    for i = 1 : length(MODELS)
        this = params(i).siga;
        if isempty(this)
            fprintf(fid, ' & ');
        else
            fprintf(fid, ' & $\\bigl[ %6.3f,%6.3f \\bigr]$ ', this.tails(1), this.tails(2));
        end
    end
    fprintf(fid, '\\vspace{\\skipper}\\\\\n');
    
    fprintf(fid, '\\midrule\n');
    fprintf(fid, '\\multicolumn{5}{c}{Forecast stickiness}\\\\\n');
    fprintf(fid, '\\midrule\n');

    fprintf(fid, '$\\lambda$');
    for i = 1 : length(MODELS)
        this = params(i).lambda;
        if isempty(this)
            fprintf(fid, ' & -- ');
        else
            fprintf(fid, ' & %6.3f ', this.mid);
        end
    end
    fprintf(fid, '\\\\\n');
    for i = 1 : length(MODELS)
        this = params(i).lambda;
        if isempty(this)
            fprintf(fid, ' & ');
        else
            fprintf(fid, ' & $\\bigl[ %6.3f,%6.3f \\bigr]$ ', this.tails(1), this.tails(2));
        end
    end
    fprintf(fid, '\\vspace{\\skipper}\\\\\n');
    
    
    % SIGLAMBDA
    fprintf(fid, '$\\sigma_{\\kappa}^2$');
    for i = 1 : length(MODELS)
        this = params(i).siglambda;
        if isempty(this)
            fprintf(fid, ' & -- ');
        else
            fprintf(fid, ' & %6.3f ', this.mid);
        end
    end
    fprintf(fid, '\\\\\n');
    for i = 1 : length(MODELS)
        this = params(i).siglambda;
        if isempty(this)
            fprintf(fid, ' & ');
        else
            fprintf(fid, ' & $\\bigl[ %6.3f,%6.3f \\bigr]$ ', this.tails(1), this.tails(2));
        end
    end
    fprintf(fid, '\\vspace{\\skipper}\\\\\n');
    
    fprintf(fid, '\\midrule\n');
    fprintf(fid, '\\multicolumn{5}{c}{Measurement error variances}\\\\\n');
    fprintf(fid, '\\midrule\n');

    % NOISEVAR
    for j = 1 : 6
        
        if j == 1
            fprintf(fid, '$\\sigma^2_{\\zeta,\\pi}$');
        else
            fprintf(fid, '$\\sigma^2_{\\zeta,%d}$', j-1);
        end
        for i = 1 : length(MODELS)
            if ~doInflationNoise && j == 1
                fprintf(fid, ' & --');
            else
                fprintf(fid, ' & %6.3f ', params(i).noisevar.mid(j));
            end
        end
        fprintf(fid, '\\\\\n');
        if doInflationNoise || j > 1
            for i = 1 : length(MODELS)
                fprintf(fid, ' & $\\bigl[ %6.3f,%6.3f \\bigr]$ ', params(i).noisevar.tails(1,j), params(i).noisevar.tails(2,j));
            end
        end
        fprintf(fid, '\\vspace{\\skipper}\\\\\n');
        
    end
    
    fprintf(fid, '\\midrule\n');
    % MDD
    fprintf(fid, '$\\ln \\text{MDD} \\Bigl ( \\mathcal{M}_i \\Bigl | \\mathcal{Y}_{1:T} \\Bigr )$ ');
    for i = 1 : length(MODELS)
        fprintf(fid, ' & $%6.3f \\qquad $ ', MDDmeans(i));
    end
    fprintf(fid, '\\\\\n');
    for i = 1 : length(MODELS)
        fprintf(fid, ' & $\\bigl ( %6.3f \\bigr )$ ', MDDstderr(i));
    end
    fprintf(fid, '\\\\\n');
    
    
    % wrap up table
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    if params(1).ndxmid == 1
        midTxt = 'mean';
    else
        midTxt = 'median';
    end
    
    % fprintf(fid, '\\\\[\\baselineskip]\nNote: Posterior %ss for parameters and 90\\%% confidence bands in square brackets. logMDD reported as average over 250 repetitions of the particle filter with 100''000 particles, standard deviations across the 250 repetitions in brackets.\n', midTxt);
    
    fclose(fid);
    
    type(tabname)
    
end
%% finish
% !pdflatex list-parameterTables.tex
finishscript
