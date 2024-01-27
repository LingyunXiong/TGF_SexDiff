%% Load BP-RBF data (Male WT CTL)
% 19D16 (500Hz)
load('Cupples_SimonFraser/WT/19D16_W_M_CTL.mat');
BP_m1 = p1;
RBF_m1 = q1;
tspan = 0:.002:1600; % [second] 
tspan(end)=[];
% RBF_m1_MA = movmean(RBF_m1,500); % 1Hz
RBF_m1_MA = movmean(RBF_m1,5000); % 0.1Hz

% Plot the time-series
m1=figure;
plot(tspan,RBF_m1_MA,'LineWidth',1); hold on
title('19D16 WT Male CTL (0.1Hz)')
xlabel('Time (s)'); ylabel('RBF (mL/min)') 
xlim([0 1600])
hold off

% 19E15 (500Hz)
load('Cupples_SimonFraser/WT/19E15_W_M_CTL.mat');
BP_m2 = p1;
RBF_m2 = q1;
tspan = 0:.002:1600; % [second]
tspan(end)=[];
% RBF_m2_MA = movmean(RBF_m2,500); % 1Hz
RBF_m2_MA = movmean(RBF_m2,5000); % 0.1Hz

% Plot the time-series
m2=figure;
plot(tspan,RBF_m2_MA,'LineWidth',1); hold on
title('19E15 WT Male CTL (0.1Hz)')
xlabel('Time (s)'); ylabel('RBF (mL/min)') 
xlim([0 1600])
hold off

% 19K14 (1000Hz)
load('Cupples_SimonFraser/WT/19K14_W_M_CTL.mat');
BP_m3 = p1;
RBF_m3 = q1;
tspan = 0:.001:1600; % [second]
tspan(end)=[];
% RBF_m3_MA = movmean(RBF_m3,1000); % 1Hz
RBF_m3_MA = movmean(RBF_m3,10000); % 0.1Hz

% Plot the time-series
m3=figure;
plot(tspan,RBF_m3_MA,'LineWidth',1); hold on
title('19K14 WT Male CTL (0.1Hz)')
xlabel('Time (s)'); ylabel('RBF (mL/min)') 
xlim([0 1600])
hold off

% 19K26 (1000Hz)
load('Cupples_SimonFraser/WT/19K26_W_M_CTL.mat');
BP_m4 = p1;
RBF_m4 = q1;
tspan = 0:.001:1600; % [second]
tspan(end)=[];
% BP_m4_MA = movmean(BP_m4,1000); % 1Hz
BP_m4_MA = movmean(BP_m4,10000); % 0.1Hz
% RBF_m4_MA = movmean(RBF_m4,1000); % 1Hz
RBF_m4_MA = movmean(RBF_m4,10000); % 0.1Hz

% Plot the time-series
m4=figure;
plot(tspan,RBF_m4_MA,'LineWidth',1); hold on
title('19K26 WT Male CTL (0.1Hz)')
xlabel('Time (s)'); ylabel('RBF (mL/min)') 
xlim([0 1600])
hold off

% Subset time-series:
idx_vec=1000*(1000:1500);
m4s=figure;
plot(tspan(idx_vec),RBF_m4_MA(idx_vec),'LineWidth',1); hold on
title('Male WT (19K26)')
xlabel('Time (s)'); ylabel('RBF (mL/min)') 
xlim([1000 1500])
xticks([1000 1100 1200 1300 1400 1500])
xticklabels({'0','100','200','300','400','500'})
hold off

m4=figure;
plot(tspan,BP_m4_MA,'LineWidth',1); hold on
title('19K26 WT Male CTL (0.1Hz)')
xlabel('Time (s)'); ylabel('BP (mmHg)') 
xlim([0 1600])
hold off

figure(1)
plot(BP_m4_MA,RBF_m4_MA,'LineWidth',1);
title('19K26 WT Male CTL (0.1Hz) BP-RBF Crossplot')
xlabel('BP (mmHg)','fontsize',20)
ylabel('RBF (mL/min)','fontsize',20) 

figure(2)
for j = 1:length(tspan)
    if mod(j,500)==0
        plot(BP_m4_MA(j),RBF_m4_MA(j),'*','MarkerSize',3); hold on
        xlabel('BP (mmHg)','fontsize',20)
        ylabel('RBF (mL/min)','fontsize',20) 
        title('19K26 WT Male CTL (0.1Hz) BP-RBF Crossplot', sprintfc('time step = %i',j))
        axis([50 130 6 12])
        pause(.05)
    end
end
hold off

% 19L05 (1000Hz)
load('Cupples_SimonFraser/WT/19L05_W_M_CTL.mat');
BP_m5 = p1;
RBF_m5 = q1;
tspan = 0:.001:1600; % [second]
tspan(end)=[];
% RBF_m5_MA = movmean(RBF_m5,1000); % 1Hz
RBF_m5_MA = movmean(RBF_m5,10000); % 0.1Hz

% Plot the time-series
m5=figure;
plot(tspan,RBF_m5_MA,'LineWidth',1); hold on
title('19L05 WT Male CTL (0.1Hz)')
xlabel('Time (s)'); ylabel('RBF (mL/min)') 
xlim([0 1600])
hold off

%% Load BP-RBF data (Female WT CTL)
% 19F13 (500Hz)
load('Cupples_SimonFraser/WT/19F13_W_F_CTL.mat');
BP_f1 = p1;
RBF_f1 = q1;
tspan = 0:.002:1600; % [second]
tspan(end)=[];
BP_f1_MA = movmean(BP_f1,5000); % 0.1Hz
% RBF_f1_MA = movmean(RBF_f1,500); % 1Hz
RBF_f1_MA = movmean(RBF_f1,5000); % 0.1Hz

% Plot the time-series
f1=figure;
plot(tspan,RBF_f1_MA,'LineWidth',1); hold on
title('19F13 WT Female CTL (0.1Hz)')
xlabel('Time (s)'); ylabel('RBF (mL/min)') 
xlim([0 1600])
hold off

% Subset time-series:
idx_vec=500*(200:700);
f1s=figure;
plot(tspan(idx_vec),RBF_f1_MA(idx_vec),'LineWidth',1); hold on
title('Female WT (19F13)')
xlabel('Time (s)'); ylabel('RBF (mL/min)') 
xlim([200 700])
xticks([200 300 400 500 600 700])
xticklabels({'0','100','200','300','400','500'})
hold off

figure(3)
plot(BP_f1_MA,RBF_f1_MA,'LineWidth',1);
title('19F13 WT Female CTL (0.1Hz) BP-RBF Crossplot')
xlabel('BP (mmHg)','fontsize',20)
ylabel('RBF (mL/min)','fontsize',20) 

h = figure;
for j = 1:length(tspan)
    if mod(j,500)==0
        plot(BP_f1_MA(j),RBF_f1_MA(j),'*','MarkerSize',3); hold on
        xlabel('BP (mmHg)','fontsize',20)
        ylabel('RBF (mL/min)','fontsize',20) 
        title('19F13 WT Female CTL (0.1Hz) BP-RBF Crossplot', sprintfc('time step = %i',j))
        axis([85 120 2.6 4])
        pause(.05)
    end
end
hold off

% 19D10 (500Hz)
load('Cupples_SimonFraser/WT/19D10_W_F_CTL.mat');
BP_f2 = p1;
RBF_f2 = q1;
tspan = 0:.002:1600; % [second]
tspan(end)=[];
RBF_f2_MA = movmean(RBF_f2,500); % 1Hz
RBF_f2_MA = movmean(RBF_f2,5000); % 0.1Hz

% Plot the time-series
f2=figure;
plot(tspan,RBF_f2_MA,'LineWidth',1); hold on
title('19D10 WT Female CTL (0.1Hz)')
xlabel('Time (s)'); ylabel('RBF (mL/min)') 
xlim([0 1600])
hold off

% 19E16 (500Hz)
load('Cupples_SimonFraser/WT/19E16_W_F_CTL.mat');
BP_f3 = p1;
RBF_f3 = q1;
tspan = 0:.002:1600; % [second]
tspan(end)=[];
RBF_f3_MA = movmean(RBF_f3,500); % 1Hz
RBF_f3_MA = movmean(RBF_f3,5000); % 0.1Hz

% Plot the time-series
f3=figure;
plot(tspan,RBF_f3_MA,'LineWidth',1); hold on
title('19E16 WT Female CTL (0.1Hz)')
xlabel('Time (s)'); ylabel('RBF (mL/min)') 
xlim([0 1600])
hold off

% 19K20 (1000Hz)
load('Cupples_SimonFraser/WT/19K20_W_F_CTL.mat');
BP_f4 = p1;
RBF_f4 = q1;
tspan = 0:.001:1600; % [second]
tspan(end)=[];
BP_f4_MA = movmean(BP_f4,10000); % 0.1Hz
% RBF_f4_MA = movmean(RBF_f4,1000); % 1Hz
RBF_f4_MA = movmean(RBF_f4,10000); % 0.1Hz

% Plot the time-series
f4=figure;
plot(tspan,RBF_f4_MA,'LineWidth',1); hold on
title('19K20 WT Female CTL (0.1Hz)')
xlabel('Time (s)'); ylabel('RBF (mL/min)') 
xlim([0 1600])
hold off

f4=figure;
plot(tspan,BP_f4_MA,'LineWidth',1); hold on
title('19K20 WT Female CTL (0.1Hz)')
xlabel('Time (s)'); ylabel('BP (mmHg)') 
xlim([0 1600])
hold off

% 19L11 (1000Hz)
load('Cupples_SimonFraser/WT/19L11_W_F_CTL.mat');
BP_f5 = p1;
RBF_f5 = q1;
tspan = 0:.001:1600; % [second]
tspan(end)=[];
% RBF_f5_MA = movmean(RBF_f5,1000); % 1Hz
RBF_f5_MA = movmean(RBF_f5,10000); % 0.1Hz

% Plot the time-series
f5=figure;
plot(tspan,RBF_f5_MA,'LineWidth',1); hold on
title('19L11 WT Female CTL (0.1Hz)')
xlabel('Time (s)'); ylabel('RBF (mL/min)') 
xlim([0 1600])
hold off
