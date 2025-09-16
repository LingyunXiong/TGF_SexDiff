%Clear
clear all
close all

set(0,                          ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

%% Load sample data 
table =  readtable('../Data/IVM_M3G_Baseline_Region4_Data.csv');
time = table.Time; % unit: second
dt_dex = table.Dextran_Mean; 
dt_dex_MA = movmean(dt_dex,15);

%% Plots
figure(1)
plot(time,dt_dex,'k','LineWidth',2)
xlabel("Time (second)",FontSize=20)
ylabel("Fluorescence Intensity (a.u.)",FontSize=20)
axis([0 400 2.5 9])

%% Fitting Sinewaves 
X = time;
Y = dt_dex;

B0 = mean(Y);               % Vertical shift
B1 = (max(Y) - min(Y))/2;   % Amplitude
B2 = 2*pi*0.032;            % Phase (Number of peaks)
B3 = 0;                     % Phase shift (eyeball the Curve)
myFit = NonLinearModel.fit(X,Y, 'y ~ b0 + b1*sin(b2*x1 + b3)', [B0, B1, B2, B3])

figure(2);
plot(time,dt_dex,'-k','LineWidth',2); hold on
plot(time(210:end),myFit.Fitted(210:end),'-r','LineWidth',2)
xlabel('Time (second)','fontsize',20)
ylabel("Fluorescence Intensity (a.u.)",FontSize=20)
axis([0 400 2.5 9])
xticks(0:100:400)
lgd = legend('Raw','Ref: 0.032 Hz');
lgd.FontSize = 20;
hold off
