%% code to simulate advection process via Method of lines 
close all
clear all

set(0,                          ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);
  
%% Parameter Sensitivity Test
% Set the range of parameter values (Fine)
param1=flip(0.5:0.01:0.98); % parameter range for k
param2=flip(0.005:0.0005:0.02); % parameter range for lambda

% Set the range of parameter values (Coarse)
% param1=flip(0.5:0.02:0.98); % parameter range for k
% param2=flip(0.007:0.001:0.02); % parameter range for lambda
 
% simulation time step and duration
tstep = 1;  % time step interval
t_end = 6000; % total time to run simulation
tspan = 0:tstep:t_end;
T_cutoff = 2000;
vec_Amplitude = zeros(length(param1),length(param2));
vec_AvePeriod = zeros(length(param1),length(param2));
vec_AveFreq = NaN(length(param1),length(param2));

% load('Simulations/TGF_Model_v1_ParamSweep2_Coarse.mat')
% load('Simulations/TGF_Model_v1_ParamSweep2_Coarse_MaleSigmaN.mat')
% load('Simulations/TGF_Model_v1_ParamSweep2_Fine.mat')
load('Simulations/TGF_Model_v1_ParamSweep2_Fine_MaleSigmaN.mat')

% calculate amplitude and period
for i = 1:length(param1)
    for j = 1:length(param2)
        G = mtx_GFR(:,i,j); % GFR
        
        if length(G)<6001
	        vec_Amplitude(i,j) = -1;
	        vec_AvePeriod(i,j) = -1;
            continue
        end

        [pks,plocs] = findpeaks(G(end-T_cutoff:end,end));
        [vks,vlocs] = findpeaks(-G(end-T_cutoff:end,end));
        
        if ~isempty(pks) && ~isempty(vks)
	        Amp = pks(end)+vks(end);
            if Amp > 0.5 && length(plocs)>1
                vec_Amplitude(i,j) = Amp;
                vec_AvePeriod(i,j) = mean(diff(plocs))/10; % unit: sec
                vec_AveFreq(i,j) = 1000*10/mean(diff(plocs)); % unit: mHz
            end
        end
    end
end

% figure(1);
% h=heatmap(param1,param2,vec_Amplitude','GridVisible','off');
% h.Colormap=jet;
% title('Parameter Sweep for Amplitude (k and \lambda)');
% xlabel('Parameter Value (k)'); ylabel('Parameter Value (\lambda)');
% XLabels=1:length(param1);
% CustomXLabels = string(param1);
% CustomXLabels(mod(XLabels,4) ~= 1) = " ";
% h.XDisplayLabels = CustomXLabels;
% YLabels=1:length(param2);
% CustomYLabels = string(param2);
% CustomYLabels(mod(YLabels,5) ~= 1) = " ";
% h.YDisplayLabels = CustomYLabels;
% clim(h,[0 28]);
% exportgraphics(gcf,'ParamSweep_Codim2_Amplitude-K&Lambda.png','Resolution',200)

% figure(2);
% h=heatmap(param1,param2,vec_AvePeriod','GridVisible','off');
% h.Colormap=jet;
% title('Parameter Sweep for Period (k and \lambda)');
% xlabel('Parameter Value (k)'); ylabel('Parameter Value (\lambda)');
% XLabels=1:length(param1);
% CustomXLabels = string(param1);
% CustomXLabels(mod(XLabels,4) ~= 1) = " ";
% h.XDisplayLabels = CustomXLabels;
% YLabels=1:length(param2);
% CustomYLabels = string(param2);
% CustomYLabels(mod(YLabels,5) ~= 1) = " ";
% h.YDisplayLabels = CustomYLabels;
% clim(h,[0 58]);
% exportgraphics(gcf,'ParamSweep_Codim2_Period-K&Lambda.png','Resolution',200)

% figure(3);
% h=heatmap(param1,param2,vec_AveFreq','GridVisible','off');
% h.Colormap=jet;
% title('Parameter Sweep for Frequency (k and \lambda)');
% xlabel('Parameter Value (k)'); ylabel('Parameter Value (\lambda)');
% XLabels=1:length(param1);
% CustomXLabels = string(param1);
% CustomXLabels(mod(XLabels,4) ~= 1) = " ";
% h.XDisplayLabels = CustomXLabels;
% YLabels=1:length(param2);
% CustomYLabels = string(param2);
% CustomYLabels(mod(YLabels,5) ~= 1) = " ";
% h.YDisplayLabels = CustomYLabels;
% clim(h,[15 35]);
% exportgraphics(gcf,'ParamSweep_Codim2_Frequency-K&Lambda.png','Resolution',200)

figure(4);
h=heatmap(param2,param1,vec_AveFreq,'GridVisible','off','CellLabelColor','none');
h.Colormap=jet;
title('Parameter Sweep for Frequency (k and \lambda)');
xlabel('Parameter Value (\lambda)'); ylabel('Parameter Value (k)');
XLabels=1:length(param2);
CustomXLabels = string(param2);
CustomXLabels(mod(XLabels,5) ~= 1) = " ";
h.XDisplayLabels = CustomXLabels;
% h.XDisplayData = flipud(h.XDisplayData);
YLabels=1:length(param1);
CustomYLabels = string(param1);
CustomYLabels(mod(YLabels,8) ~= 1) = " ";
h.YDisplayLabels = CustomYLabels;
% h.YDisplayData = flipud(h.YDisplayData);
clim(h,[15 35]);
% exportgraphics(gcf,'ParamSweep_Codim2Fine_Frequency-Lambda&K_MaleSigmaN_v2.png','Resolution',200)

 
