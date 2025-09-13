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

%% Heatmap and Contour Plots
figure(1);
ax = axes();
hold(ax);
p = pcolor(param2,param1,vec_Amplitude);
p.FaceAlpha = 0.5;
p.EdgeAlpha=0;
c = contour(param2,param1,vec_Amplitude,[5 10 15 20 25],'ShowText','on');
xlabel('Parameter Value (\lambda)'); ylabel('Parameter Value (\kappa)');
colormap gray; colorbar

vec_RelativeAmp = 100 * vec_Amplitude ./ vec_AveGFR;
figure(2);
ax = axes();
hold(ax);
p = pcolor(param2,param1,vec_RelativeAmp);
p.FaceAlpha = 0.75;
p.EdgeAlpha=0;
set(gca, 'clim', [0 200]);
c = contour(param2,param1,vec_RelativeAmp,[5 25 50 100 150],'ShowText','on');
xlabel('Parameter Value (\lambda)'); ylabel('Parameter Value (\kappa)');
colormap jet; colorbar
