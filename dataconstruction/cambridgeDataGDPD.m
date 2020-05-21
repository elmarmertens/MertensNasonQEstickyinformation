%% construct data set for project "cambridge"
% data sources obtained from FRB-PHIL
% 1) Mean responses of SPF for GNP/GDP deflator
%    https://www.phil.frb.org/research-and-data/real-time-center/survey-of-professional-forecasters/data-files/files/Mean_PGDP_Level.xls
% 2) The corresponding real-time data for (second revision)
%    https://www.philadelphiafed.org/-/media/research-and-data/real-time-center/real-time-data/data-files/files/xlsx/p_first_second_third.xlsx?la=en
%
% This code lags the forecasts by one period (i.e. the nowcast is
% actually F(t-1) \pi(t) (the adjustement occurs vis differences in the variables spfdates and dates, see below)
%
%
% BEFORE RUNNNING THIS SCRIPT:
% - download data (xlsx format is fine); note: coversion via csv can result in loss of some decimals ...
% - adjust date vectors in the code, see [ UPDATE ] tags in the code

%% clear workspace
clear variables
clear global
close all
fclose all;
clc

%% prepare variables 

% quarters are time-stamped on the first day of the quarter (FRED convention)
% note: the following could also be automated, but might be good to set a
% few things manually (also forces us to check things when updating the data)


% [ UPDATE ]
dates    = genrQdates(1968,2018);
spfdates = dates(4:end);     % these must be the dates as stated in the SPF file (1968:Q4 through 2018:Q4)
dates    = dates(3:end-1);   % these are the model dates (after lagging the forecast horizon, 1968:Q3 through 2018:Q3)


T        = length(dates);
Ny       = 1 + 5 + 1; % actual inflation plus five quarterly forecasts and one annual forecast
horizons = 1 : 5;     % SPF forecast horizons (adjusted)
Ydata    = NaN(T,Ny); % this will be the output of this program

vintagechoice = 'GDPD';

%% load SPF responses for the level of PGDP and convert into inflation rates
% the results of this cell are to be stored in Ydata(:,2:6)

SPFimport = importdata('Mean_PGDP_Level.xlsx');
% col 1: year
% col 2: quarter
% col 3: "forecasts" for previous quarter
% col 4: nowcast
% col 5: forecast for next quarter
% col 6: 2-quarter forecast
% col 7: 3-quarter forecast Y
% col 8: 4-quarter forecast
% col 9: Annual forecast for current year
% col 10: Annual forecast for next year


% convert -999 into NaN
SPFimport.data(SPFimport.data == -999) = NaN;

% check dates
Y = SPFimport.data(:,1);
if ~isequal(Y, year(spfdates))
    error('inconsistent date vector in the SPF file (years)');
end
Q = SPFimport.data(:,2);
if ~isequal(Q, quarter(spfdates))
    error('inconsistent date vector in the SPF file (quarters)');
end

% convert SPF level forecasts into inflation forecasts
% note: the base level is given by the first data column of the SPF file,
% which reports the SPF "forecasts" for the previous quarter


Ydata(:,2:6) = 100 * bsxfun(@minus, log(SPFimport.data(:,4:8)), ...
    log(SPFimport.data(:,3)));
% annualize inflation forecasts
Ydata(:,2:6) = bsxfun(@times, Ydata(:,2:6), 4 ./ horizons);

% import B forecast, notice: this is the average of the quarterly price
% levels expected for next year (will be normalized by current price level,
% notice that this average forecast does not get properly annualized here)
Ydata(:,7)   = 100 * (log(SPFimport.data(:,10)) - log(SPFimport.data(:,3)));

% Bhorizons    = 5 - quarter(dates); % need to check

%% load second revision data
majorreleases = importdata('p_first_second_third.xlsx');

% [ UPDATE ]
infdates = genrQdates(1965,2018);
infdates = infdates(3:end-1); % 1965:Q3 through 2018:Q3;

% check dates
checkdates = datenum(majorreleases.textdata.DATA(6:end,1), 'YYYY:QQ');

if ~isequal(checkdates, infdates)
    error('mismatch in release dates')
end


secondrelease = majorreleases.data.DATA(:,2);
secondrelease(secondrelease == -999) = NaN;

% patch in third release value for 2003Q3 (Govt shutdown)
ndxShutdown2003 = infdates == datenum(2003, 7, 1);
if isnan(secondrelease(ndxShutdown2003,1))
    secondrelease(ndxShutdown2003) = majorreleases.data.DATA(ndxShutdown2003,3);
else
    error('problem identifying missing obs related to 2003Q3 shutdown')
end

basendx        = ismember(dates,infdates);
ndx            = ismember(infdates, dates);

inflation2nd = NaN(T,1);
inflation2nd(basendx) = secondrelease(ndx);

% convert from simple rate of change into log-difference
inflation2nd = log(1 + inflation2nd / 100) * 100;

Ydata(:,1) = inflation2nd;
Ydata      = Ydata(:,1:end-1); % disregard the B forecast
Ny         = size(Ydata, 2);   % number of observables

%% store data

Ylabel    = {'Inflation', 'Surv1Q', 'Surv2Q', 'Surv3Q', 'Surv4Q', 'Surv5Q'};
Ynames    = {'Inflation', 'Survey 1Q', 'Survey 2Q', 'Survey 3Q', 'Survey 4Q', 'Survey 5Q'};
if length(Ylabel) ~= Ny
    error('dimension mismatch')
end


%% Transform data into forecasts of quarterly inflation
yNaNndx  = isnan(Ydata);
Nsurveys = length(horizons);


% note: observations are cumulative forecasts
w = diag(1./horizons);
for n = 2 : Nsurveys
    w(n,1:n-1) = w(n,n);
end
ww = blockdiag(1,w);

% rotate the data
for t = 1 : T
    
    thisY = Ydata(t,2:end)';
    thatY = NaN(Nsurveys,1);
    nanny = ~yNaNndx(t,2:end);
    thatY(nanny) = w(nanny,nanny) \ thisY(nanny);
    
    % alternative
    thatY2 = NaN(Nsurveys,1);
    thatY2(1) = thisY(1);
    for n = 2 : Nsurveys
        thatY2(n) = n * thisY(n) - (n-1) * thisY(n-1);
    end
    checkdiff(thatY,thatY2);
    
    Ydata(t,2:end) = thatY;
end

% remove missing data
Ydata(yNaNndx) = 0;

% store data
mat2fortran(sprintf('cambridge2018%s.yData.txt', vintagechoice), Ydata)
logical2fortran(sprintf('cambridge2018%s.yNaN.txt', vintagechoice), yNaNndx)
mat2fortran(sprintf('cambridge2018%s.dates.txt', vintagechoice), dates)

%% finish
dir *.txt

dockAllFigures
