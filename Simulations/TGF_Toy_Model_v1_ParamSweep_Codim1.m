%% code to simulate advection process via Method of lines 
function TGF_Toy_Model_v1_ParamSweep_Codim1
global k lambda sigma alpha beta cop n
 
set(0,                          ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);
  
% parameters
k = 0.72;  % time constant (Female)
lambda = 0.0072; % reabsorption rate in PT (Female): lower fractional resorption
sigma = 0.1; % reabsorption rate in TAL (Female)
alpha = 0.5; % sensitivity of the negative feedback
beta = 0.5; % shift of the sigmoid function
cop = 0.3; % target [Na+] at MD
n = 3; % slope of the sigmoid function (Female)

% initial condition
N = 68; % number of cells in the 1-D domain (PT and TAL)  
g0 = 30;
u0 = g0*ones(N,1); % uniform fluid volume in PT
ci = 1;
c0 = ci*ones(N,1); % uniform sodium concentration in TAL
s0 = [u0; c0; g0]; 

% modify parameters
k = 0.98; % time constant (Male)
sigma = 0.12; % reabsorption rate in TAL (Male)
n = 4; % slope of the sigmoid function (Male)

%% Parameter Sensitivity Test
% Set the range of parameter values
% param=flip(0.5:0.01:0.97); % parameter range for k
param=flip(0.007:0.0005:0.02); % parameter range for lambda
 
% simulation time step and duration
tstep = 1;  % time step interval
t_end = 6000; % total time to run simulation
tspan = 0:tstep:t_end;
% T_cutoff = 2000;
mtx_GFR = zeros(length(tspan),length(param));
% vec_Amplitude = zeros(1,length(param));
% vec_AvePeriod = zeros(1,length(param));

% integrate the system of ode's:
for i = 1:length(param)
    % k = param(i);
    lambda = param(i);
    [T,S] = ode23s(@deRHS, tspan, s0, odeset('maxstep',1));  
    G = k*(S(:,2*N+1)+0); % GFR

    if length(G)<6001
        continue
    end
    mtx_GFR(:,i)=G;
end

save('TGF_Model_v1_ParamSweep.mat','mtx_GFR')

%% Visualizations
% load('TGF_Model_v1_ParamSweep.mat')
tspan = tspan/10; 
% param=flip(0.7:0.01:0.97); % parameter range for k
% param=flip(0.007:0.0005:0.02); % parameter range for lambda

figure;
h=heatmap(tspan,param,mtx_GFR','GridVisible','off');
h.Colormap=jet;
title('Parameter Sweep (\lambda)');
xlabel('Time (s)'); ylabel('Parameter Value');
XLabels=1:length(tspan);
CustomXLabels = string(tspan);
CustomXLabels(mod(XLabels,600) ~= 1) = " ";
h.XDisplayLabels = CustomXLabels;
YLabels=1:length(param);
CustomYLabels = string(param);
CustomYLabels(mod(YLabels,5) ~= 1) = " ";
h.YDisplayLabels = CustomYLabels;
% clim(h,[0 46]);
exportgraphics(gcf,'ParamSweep_Codim1_GFR-Lambda_wMaleK_n4_Sigma.png','Resolution',200)

%% plot GFR
% figure(1)
% plot(tspan,G(:),'*'); hold on
% ylabel('GFR(t) (nl/min)','fontsize',20)
% xlabel('t (s)','fontsize',20)
% axis([0 tspan(end) 0 40])
% title('SNGFR')
% hold off

%% find out about time delay
% [G_pks,G_locs] = findpeaks(G(end-2500:end));
% [G_vks,G_vlocs] = findpeaks(-G(end-2500:end));
% [Fex_pks,Fex_locs] = findpeaks(V(end-2500:end,end));
% [Cmd_pks,Cmd_locs] = findpeaks(C(end-2500:end,end));
% 
% tau1 = Fex_locs - G_locs
% tau2 = Cmd_locs - G_locs
% % tau3 = G_vlocs(2:end) - Cmd_locs
% 
% FRatio = mean(V(end-2500:end,end))./mean(G(end-2500:end))
% CRatio = mean(C(end-2500:end,end))/ci

%% calculate cumulative salt reabsorption
% T_cutoff = 2000;
% tspan=T/10;
% Q = cumtrapz(tspan(end-T_cutoff:end,end),G(end-T_cutoff:end,end));
% Sin=ci*Q(end);
% 
% A = ci*k*V(:,end);
% Asum = cumtrapz(tspan(end-T_cutoff:end,end),A(end-T_cutoff:end,end));
% Sout_PT = Asum(end);
% 
% B = k*V(:,end).*C(:,end);
% Bsum = cumtrapz(tspan(end-T_cutoff:end,end),B(end-T_cutoff:end,end));
% Sout_TAL = Bsum(end);
% 
% Frac_PT = 1-Sout_PT/Sin
% Frac_TAL = 1-Sout_TAL/Sin

%% parameter sensitivity test (periodicity)
% [pks,plocs] = findpeaks(G(end-T_cutoff:end,end));
% [vks,vlocs] = findpeaks(-G(end-T_cutoff:end,end));
% AveAmp = mean(pks+vks);
% AvePeriod = mean(diff(plocs))/10; % unit: sec

end
 
%% the right hand side for ode simulation:
function s_prime = deRHS(t,s)
    global k lambda sigma alpha beta cop n

    nc = (length(s)-1)/2;
    v = s(1:nc);
    c = s(nc+1:2*nc);
    g = s(end);
 
    % fluid volume in the PT compartment
    Fv = k*([g;v(1:end-1)]-v) - lambda*v; % with linear reabsorption

    % sodium concentration in the TAL compartment
    Fc = k*v(end)*([1;c(1:end-1)]-c) - sigma*c; % with linear reabsorption

    % negative feedback from MD to AA
    Fg = alpha*(cop^n/(cop^n+c(end)^n)-beta); % downward sigmoid function
    
    s_prime = [Fv; Fc; Fg];
end

