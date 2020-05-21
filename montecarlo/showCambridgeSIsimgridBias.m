%% plot results

%% load toolboxes
path(pathdef)

addpath ../matlabbox/emtools/
addpath ../matlabbox/emtexbox/
addpath ../matlabbox/emgibbsbox/
addpath ../matlabbox/emeconometrics/
addpath ../matlabbox/emstatespace/

%% settings
initscript

modellabel = 'SIthetaTVPlambdaTVP';

doInflationNoise = true;

T          = 200;
Ngrid      = 320;
Nparticles = 1e4;

doSmoother         = false;
Nsmoother          = 1e2;
smootherNparticles = 1e3;

datasetlabel = sprintf('simdataT%d', T);
datalabel    = sprintf('%s.%s', datasetlabel, modellabel);
if ~doInflationNoise
    datalabel = strcat(datalabel, '.nonoise');
end
if doSmoother
    datalabel    = sprintf('%s.Nsmoother%d.smootherNparticles%d', datalabel, Nsmoother, smootherNparticles); %#ok<*UNRCH>
end
datalabel    = sprintf('%s.Nparticles%d.Ngrid%d', datalabel, Nparticles, Ngrid);

titlename = strrep(datalabel, '.', '_');
wrap = [];
initwrap

timestamp   = [];

datadir = pwd;

%% prep objects
filext = sprintf('%s.dat', datalabel);

dates = 1:T;

Nx      = 2;
Nsv     = 2;

percy = [5 25 75 95];

% define whether to use mean or median
mittel = @(x) nanmedian(x,1)';


%% select biastype
biastype    = 'Bias';
medbiastype = 'Medbias';

%% print settings
type(fullfile(datadir, strcat('settings.', filext)));
if ~isempty(wrap)
    copyfile(fullfile(datadir, strcat('settings.', filext)), fullfile(wrap.dir, strcat('settings.', filext)))
    latexwrapper(wrap, 'add', 'listing', strcat('settings.', filext))
end

%% read filter results
filterX = NaN(Ngrid,T,Nx);
for n = 1 : Nx
    filterX(:,:,n) = importdata(fullfile(datadir, sprintf('filter%sX%d.%s', biastype, n, filext)));
end


filterSV = NaN(Ngrid,T,Nsv);
for n = 1 : Nsv
    filterSV(:,:,n) = importdata(fullfile(datadir, sprintf('filter%sSV%d.%s', biastype, n, filext)));
end

filterHinno = loaddat(fullfile(datadir, sprintf('filter%sHinno.%s', biastype, filext)));

filterTheta  = importdata(fullfile(datadir, sprintf('filter%sTheta.%s', biastype, filext)));
filterLambda = importdata(fullfile(datadir, sprintf('filter%sLambda.%s', biastype, filext)));

filterSigTheta  = loaddat(fullfile(datadir, sprintf('filter%sSigTheta.%s', biastype, filext)));
filterSigLambda = loaddat(fullfile(datadir, sprintf('filter%sSigLambda.%s', biastype, filext)));

filterNoisevol = importdata(fullfile(datadir, sprintf('filter%sNoisevol.%s', biastype, filext)));
Ny = size(filterNoisevol,2);

%% read smoother results
if doSmoother
    smootherX = NaN(Ngrid,T,Nx);
    for n = 1 : Nx
        smootherX(:,:,n) = importdata(fullfile(datadir, sprintf('smoother%sX%d.%s', biastype, n, filext)));
    end


    smootherSV = NaN(Ngrid,T,Nsv);
    for n = 1 : Nsv
        smootherSV(:,:,n) = importdata(fullfile(datadir, sprintf('smoother%sSV%d.%s', biastype, n, filext)));
    end

    smootherTheta  = importdata(fullfile(datadir, sprintf('smoother%sTheta.%s', biastype, filext)));
    smootherLambda = importdata(fullfile(datadir, sprintf('smoother%sLambda.%s', biastype, filext)));


end

%% check results
if any(isnan(filterX(:)))
    warning('There are %d NaN in filterX', sum(isnan(filterX(:))))
    find(any(isnan(filterX(:,:,1)),2))
    find(any(isnan(filterX(:,:,2)),2))
end
if any(isnan(filterSV(:)))
    warning('There are %d NaN in filterSV', sum(isnan(filterSV(:))))
    find(any(isnan(filterSV(:,:,1)),2))
    find(any(isnan(filterSV(:,:,2)),2))
end

if any(isnan(filterTheta(:)))
    warning('There are %d NaN in filterTheta', sum(isnan(filterTheta(:))))
    find(any(isnan(filterTheta),2))
end

if any(isnan(filterLambda(:)))
    warning('There are %d NaN in filterLambda', sum(isnan(filterLambda(:))))
    find(any(isnan(filterLambda),2))
end

% coild also check NaN for other parameters ...

%% report filter results

% X
for n = 1 : Nx
    plotthis(filterX(:,:,n), sprintf('filterX%d', n), dates, [], [], wrap, mittel, percy)
end



% SV
for n = 1 : Nsv
    plotthis(filterSV(:,:,n), sprintf('filterSV%d', n), dates, [], [], wrap, mittel, percy)
end

% Hinno
if ~isempty(filterHinno)
    for n = 1 : Nsv
        plotthisparameter(filterHinno(:,n), sprintf('filterHinno%d', n), wrap, mittel, biastype)
    end
end

%% noisevol
for n = 1 : Ny
    if n == 1 && ~doInflationNoise
        continue
    end
    plotthisparameter(filterNoisevol(:,n), sprintf('filterNoisevol%d', n), wrap, mittel, biastype)
end

%% Theta / Lambda
plotthis(filterTheta, [], dates, [], [], wrap, mittel, percy)
plotthis(filterLambda, [], dates, [], [], wrap, mittel, percy)

%% Sigtheta/lambda
if ~isempty(filterSigTheta)
    plotthisparameter(filterSigTheta, [], wrap, mittel, biastype)
end

if ~isempty(filterSigLambda)
    plotthisparameter(filterSigLambda, [], wrap, mittel, biastype)
end


%% plot smoother results
if doSmoother
    % X
    for n = 1 : Nx
        plotthis(smootherX(:,:,n), sprintf('smootherX%d', n), dates, filterX(:,:,n), [], wrap, mittel, percy)
    end



    % SV
    for n = 1 : Nsv
        plotthis(smootherSV(:,:,n), sprintf('smootherSV%d', n), dates, filterSV(:,:,n), [], wrap, mittel, percy)
    end


    % Theta / Lambda
    plotthis(smootherTheta, [], dates, filterTheta, [], wrap, mittel, percy)
    plotthis(smootherLambda, [], dates, filterLambda, [], wrap, mittel, percy)

end


%% finish
finishwrap
dockAllFigures
finishscript

%% function plotthisparameter
function plotthisparameter(this, thatlabel, wrap, mittel, biastype)

if isempty(thatlabel)
    thatlabel = inputname(1);
end

mid   = mittel(this);
if strcmpi(biastype, 'absbias')
    [pdf, x] = ksdensity(this, 'support', 'positive');
else
    [pdf, x] = ksdensity(this);
end
thisfig = figure;
plot(x,pdf, 'k-')
plotvertline(mid, [], 'k--');
plotvertline(0, [], 'r--')
wrapthisfigure(thisfig,thatlabel, wrap)
end


%% function plotthis
function plotthis(this, thislabel, dates, that, thisthatlabel, wrap, mittel, percy)

if isempty(thislabel)
    thislabel = inputname(1);
end
if isempty(thisthatlabel)
    thisthatlabel = {'smoother', 'filter'};
end


mid   = mittel(this);
tails = prctile(this, percy, 1)';
thisfig = figure;
h1 = plotCIlines(mid, tails, dates, [], 'k');

if ~isempty(that)
    mid   = mittel(that);
    tails = prctile(that, percy, 1)';
    h2 = plotCIlines(mid, tails, dates, [], 'r');
end
xlim(dates([1 end]))
plotOrigin('r--',[],[],1)
if ~isempty(that)
    legend([h1, h2], thisthatlabel{:})
end
wrapthisfigure(thisfig,thislabel, wrap)
end
