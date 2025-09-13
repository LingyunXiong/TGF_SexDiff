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
axis([0 500 2.5 9])

%% FFT
Y = fft(detrend(dt_dex));

Fr = 0.74;              % frame rate
Fs = 1/Fr;              % sampling frequency
L = length(time);       % number of samples 
f = transpose(Fs/L*(0:L/2-1));

P1 = Y(1:L/2);
power = abs(P1).^2/L;    % power of the DFT

idx1 = find(f<0.01);
idx2 = find(f<0.02);
idx3 = find(f<0.08);

figure(2)
plot(f((idx1(end)+1):end),power((idx1(end)+1):end),"-o")
title("Single-Sided Spectrum of Raw Signal")
xlabel("Frequency (Hz)")
ylabel("Power")
xlim([0 0.5])

%% Inverse FFT
Y = fft(detrend(dt_dex));
Y2 = [zeros(1,idx2(end)),transpose(Y((idx2(end)+1):idx3(end))),zeros(1,(L-idx3(end)))]; % filter frequency: 0.02-0.08Hz
X2 = real(ifft(Y2));

figure(3)
plot(time,X2,'r','LineWidth',2); hold on
xlabel("Time (second)",FontSize=20)
ylabel("Detrended Signal (a.u.)",FontSize=20)
xlim([0 500])

figure(4)
plot(time,dt_dex,'k','LineWidth',2); hold on
plot(time,X2'+dt_dex_MA,'r','LineWidth',2)
xlabel("Time (second)",FontSize=20)
ylabel("Fluorescence Intensity (a.u.)",FontSize=20)
axis([0 500 2.5 9])
lgd = legend('Raw','FFT Smoothing');
lgd.FontSize = 20;
