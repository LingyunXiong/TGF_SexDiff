%% Load BP-RBF data (Male WT CTL)
%% 19K26 (1000Hz)
load('Cupples_SimonFraser/WT/19K26_W_M_CTL.mat');
BP_m4 = p1;
RBF_m4 = q1;
tspan = 0:.001:1600; % [second]
tspan(end)=[];
BP_m4_MA = movmean(BP_m4,10000); % 0.1Hz
RBF_m4_MA = movmean(RBF_m4,10000); % 0.1Hz

figure(1)
plot(BP_m4_MA,RBF_m4_MA,'LineWidth',1);
title('19K26 WT Male CTL (0.1Hz) BP-RBF Crossplot')
xlabel('BP (mmHg)','fontsize',20)
ylabel('RBF (mL/min)','fontsize',20) 

v = VideoWriter('CrossPlot_BP_RBF_19K26_WT_M_CTL_0.1Hz','MPEG-4');
open(v)
figure(2)
for j = 1:length(tspan)
    if mod(j,1000)==0
        plot(BP_m4_MA(j),RBF_m4_MA(j),'*','MarkerSize',3); hold on
        xlabel('BP (mmHg)','fontsize',20)
        ylabel('RBF (mL/min)','fontsize',20) 
        title('19K26 WT Male CTL (0.1Hz) BP-RBF Crossplot', sprintfc('time (s) = %i',j/1000))
        axis([50 130 6.5 11])
        frame = getframe(gcf);
        writeVideo(v,frame)
        pause(.05)
    end
end
close(v)
hold off

%% Load BP-RBF data (Female WT CTL)
%% 19F13 (500Hz)
load('Cupples_SimonFraser/WT/19F13_W_F_CTL.mat');
BP_f1 = p1;
RBF_f1 = q1;
tspan = 0:.002:1600; % [second]
tspan(end)=[];
BP_f1_MA = movmean(BP_f1,5000); % 0.1Hz
RBF_f1_MA = movmean(RBF_f1,5000); % 0.1Hz

figure(3)
plot(BP_f1_MA,RBF_f1_MA,'LineWidth',1);
title('19F13 WT Female CTL (0.1Hz) BP-RBF Crossplot')
xlabel('BP (mmHg)','fontsize',20)
ylabel('RBF (mL/min)','fontsize',20) 

v = VideoWriter('CrossPlot_BP_RBF_19F13_WT_F_CTL_0.1Hz','MPEG-4');
open(v)
figure(4)
for j = 1:length(tspan)
    if mod(j,500)==0
        plot(BP_f1_MA(j),RBF_f1_MA(j),'*','MarkerSize',3); hold on
        xlabel('BP (mmHg)','fontsize',20)
        ylabel('RBF (mL/min)','fontsize',20) 
        title('19F13 WT Female CTL (0.1Hz) BP-RBF Crossplot', sprintfc('time (s) = %i',j/500))
        axis([85 120 2.6 4])
        frame = getframe(gcf);
        writeVideo(v,frame)
        pause(.05)
    end
end
close(v)
hold off
