%% construct data set for project "cambridge"
% data sources obtained from FRB-PHIL
% 1) Mean responses of SPF for CPI inflation
%    https://www.philadelphiafed.org/-/media/research-and-data/real-time-center/survey-of-professional-forecasters/data-files/files/mean_cpi_level.xlsx?la=en
% 2) FRED data for CPI seasonally adjusted
% 
%
% This code lags the forecasts by one period (i.e. the nowcast is
% actually F(t-1) \pi(t) (the adjustement occurs vis differences in the variables spfdates and dates, see below)
%
%
% BEFORE RUNNNING THIS SCRIPT:
% - download data (xlsx format is fine); note: coversion via csv can result in loss of some decimals ...
% - adjust date vectors in the code, see [ UPDATE ] tags in the code

% Note: all simple rates of change are converted into log-differences

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

vintagechoice = 'CPI';

%% load SPF responses for the level of PGDP and convert into inflation rates
% the results of this cell are to be stored in Ydata(:,2:6)

SPFimport = importdata('Mean_CPI_Level.xlsx'); 
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


Ydata(:,2:6) = SPFimport.data(:,4:8);

% import B forecast, notice: this is the average of the quarterly price
% levels expected for next year (will be normalized by current price level,
% notice that this average forecast does not get properly annualized here)
Ydata(:,7)   = SPFimport.data(:,10); 

% Bhorizons    = 5 - quarter(dates); % need to check

%% collect actual inflation data from FRED file and store in Ydata(:,1)
freddata = fredreadcsv('CPIAUCSL', 'cpi', 'm', 'q', 'avg');

infdates = freddata.DATES;
inflation = [NaN; (freddata.VALUES(2:end) ./ freddata.VALUES(1:end-1)).^4 - 1 ] * 100;
% check = [NaN; diff(log(freddata.VALUES)) * 400];

% fill inflation series into Ydata
% (fill in NaN)
ndx1             = ismember(infdates, dates);
ndx2             = ismember(dates, infdates);
Ydata(ndx2,1)    = inflation(ndx1);

%% convert simple rates of change into log-differences
simpleYdata = Ydata;
Ydata = log(1 + simpleYdata / 100) * 100;

%% compare actual inflation against SPF forecasts
   
for n = 1 : 4
    figure
    hold on
    plot(dates, Ydata(:,1), 'b-', 'linewidth', 1)
    plot(dates, Ydata(:,2), 'r-', 'linewidth', 1)
    plot(dates, Ydata(:,2+n), 'k-', 'linewidth', 3)
    legend('Inflation', 'SPF nowcast', sprintf('SPF %d quarter forecast', n))
    grid on
    xlim(dates([1 end]))
    datetick('x', 10, 'keeplimits')
    % print('-djpeg', sprintf('forecast%d', n))
end

%% cut sample
trimSample = true;
if trimSample
    firstobs = find(~isnan(Ydata(:,2)), 1);
    dates = dates(firstobs:end);
    Ydata = Ydata(firstobs:end,:);
    T     = length(dates);
end

% disregard the B forecast
Ydata    = Ydata(:,1:end-1); 
Nsurveys = length(horizons);

%% finish: store data

Ylabel    = {'Inflation', 'Surv1Q', 'Surv2Q', 'Surv3Q', 'Surv4Q', 'Surv5Q'};
Ynames    = {'Inflation', 'Survey 1Q', 'Survey 2Q', 'Survey 3Q', 'Survey 4Q', 'Survey 5Q'};


% treat missing data 
yNaNndx         = isnan(Ydata);
Ydata(yNaNndx) = 0;

% store data
mat2fortran(sprintf('cambridge2018%s.yData.txt', vintagechoice), Ydata)
logical2fortran(sprintf('cambridge2018%s.yNaN.txt', vintagechoice), yNaNndx)
mat2fortran(sprintf('cambridge2018%s.dates.txt', vintagechoice), dates)


%% finish
dir *.txt

dockAllFigures
