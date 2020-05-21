%% plot draws stored in grid-debug

%% load toolboxes
path(pathdef)

addpath ../matlabbox/emtools/
addpath ../matlabbox/emtexbox/
addpath ../matlabbox/emgibbsbox/
addpath ../matlabbox/emeconometrics/
addpath ../matlabbox/emstatespace/

%% settings
initscript

% thisname = 'THETA';
% xsupport = [-1 1];

thisname = 'LAMBDA';
xsupport = [0 1];

Ngrid     = 16;
Nparticles= 1e2;

doSmoother         = false;
Nsmoother          = 1e3;
smootherNparticles = 1e3;

T         = 200;
datalabel = 'simdata';


datadir = pwd;

%% collect draws
these = NaN(T+1,Ngrid);
for n = 1 : Ngrid
    if doSmoother
        filename   = fullfile(datadir, sprintf('true%s.grid%d.%sT%d.SIthetaTVPlambdaTVP.Nsmoother%d.smootherNparticles%d.Nparticles%d.Ngrid%d.dat', ...
            thisname, n, datalabel, T, Nsmoother, smootherNparticles, Nparticles, Ngrid)); %#ok<*UNRCH>
    else
        filename   = fullfile(datadir, sprintf('true%s.grid%d.%sT%d.SIthetaTVPlambdaTVP.Nparticles%d.Ngrid%d.dat', ...
            thisname, n, datalabel, T, Nparticles, Ngrid));
    end
    these(:,n) = importdata(filename);
end

%% plot
figure
for n = 1 : 4
    subplot(2,2,n)
    plot(0:T, these(:,n))
    switch upper(thisname)
        case 'THETA'
            ylim([-1 1])
        case 'LAMBDA'
            ylim([0 1])
    end
end

figure
hold on
plot(0:T,these)
plot(0:T,mean(these,2), 'k-', 'linewidth', 2)
title(thisname)


%% prepare 3D surf plot
xgrid = xsupport(1) : .01: xsupport(2);
pdf   = NaN(T + 1,length(xgrid));
for t = 1 : T + 1
    pdf(t,:) = ksdensity(these(t,:), xgrid, 'support', xsupport, 'bandwidth', .2, 'kernel', 'triangle');
end
figure
surf(xgrid, 0:T, pdf)
ylabel('time')
xlabel(thisname)
shading interp
