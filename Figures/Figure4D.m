%% code to visualize simulation results
close all
clear all

set(0,                          ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);
  
%% Parameter Sensitivity Test
load('../Simulations/TGF_Compartmental_Model_ParamSweep2_Refine_N3_Output.mat')

% Set the range of parameter values (Refine)
param1=flip(0.4:0.005:1); % parameter range for k
param2=flip(0.005:0.00025:0.04); % parameter range for lambda

% simulation time step and duration
tstep = 1;  % time step interval
t_end = 6000; % total time to run simulation
tspan = 0:tstep:t_end;
T_cutoff = 2000;
vec_AveGFR = zeros(length(param1),length(param2));
vec_Amplitude = NaN(length(param1),length(param2));
vec_AvePeriod = zeros(length(param1),length(param2));
vec_AveFreq = NaN(length(param1),length(param2));

% calculate amplitude and period
for i = 1:length(param1)
    for j = 1:length(param2)
        G = mtx_G(:,i,j); % GFR
        vec_AveGFR(i,j) = mean(G(end-T_cutoff:end,end));

        if length(G)<6001
	        vec_Amplitude(i,j) = -1;
	        vec_AvePeriod(i,j) = -1;
            continue
        end

        [pks,plocs] = findpeaks(G(end-T_cutoff:end,end));
        [vks,vlocs] = findpeaks(-G(end-T_cutoff:end,end));
        minimum = min(G(end-T_cutoff:end));

        if ~isempty(pks) && ~isempty(vks) && minimum>0
	        Amp = pks(end)+vks(end);
            if Amp > 0.5 && length(plocs)>1
                vec_Amplitude(i,j) = Amp;
                vec_AvePeriod(i,j) = mean(diff(plocs))/10; % unit: sec
                vec_AveFreq(i,j) = 10/mean(diff(plocs)); % unit: Hz
            end
        end
    end
end

%% Cluster Plots of Simulation Results %%
t=tspan/10;

G1=squeeze(mtx_G(:,17,121));   %k=0.92,lambda=0.01
G2=squeeze(mtx_G(:,17,109));    %k=0.92,lambda=0.013
G3=squeeze(mtx_G(:,17,65));   %k=0.92,lambda=0.024
G4=squeeze(mtx_G(:,29,49));   %k=0.86,lambda=0.026
G5=squeeze(mtx_G(:,101,129));   %k=0.5,lambda=0.008
G6=squeeze(mtx_G(:,113,121));   %k=0.44,lambda=0.01

figure(8)
subplot(6,1,1)
plot(t,0.3*G1(:),'.k','LineWidth',3); axis([0 t(5001) 0 10]); xticks([])
subplot(6,1,2)
plot(t,0.3*G2(:),'.','LineWidth',3); axis([0 t(5001) 0 10]); xticks([])
subplot(6,1,3)
plot(t,0.3*G3(:),'.','LineWidth',3); axis([0 t(5001) 0 10]); xticks([])
subplot(6,1,4)
plot(t,0.3*G4(:),'.k','LineWidth',3); axis([0 t(5001) 0 10]); xticks([])
subplot(6,1,5)
plot(t,0.3*G5(:),'.r','LineWidth',3); axis([0 t(5001) 0 10]); xticks([])
subplot(6,1,6)
plot(t,0.3*G6(:),'.r','LineWidth',3); axis([0 t(5001) 0 10])
xlabel('Time (second)','fontsize',20)
% ylabel('SNGFR(t) (nl/min)','fontsize',20)
