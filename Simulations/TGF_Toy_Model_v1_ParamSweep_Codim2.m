%% code to simulate advection process via Method of lines 
function TGF_Toy_Model_v1_ParamSweep_Codim2
global k lambda sigma alpha beta cop n
 
set(0,                          ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);
  
% parameters
k = 0.72;  % time constant (female)
lambda = 0.0072; % reabsorption rate in PT (female): lower fractional resorption
sigma = 0.1; % reabsorption rate in TAL (neutral)
alpha = 0.5; % sensitivity of the negative feedback
beta = 0.5; % shift of the sigmoid function
cop = 0.3; % target [Na+] at MD
n = 3; % slope of the sigmoid function

% initial condition
N = 68; % number of cells in the 1-D domain (PT and TAL)  
g0 = 30;
u0 = g0*ones(N,1); % uniform fluid volume in PT
ci = 1;
c0 = ci*ones(N,1); % uniform sodium concentration in TAL
s0 = [u0; c0; g0]; 

%% Parameter Sensitivity Test
% Set the range of parameter values
param1=flip(0.5:0.01:0.97); % parameter range for k
param2=flip(0.007:0.0005:0.02); % parameter range for lambda
 
% simulation time step and duration
tstep = 1;  % time step interval
t_end = 6000; % total time to run simulation
tspan = 0:tstep:t_end;
T_cutoff = 2000;
% mtx_GFR = zeros(length(tspan),length(param1),length(param2));
vec_Amplitude = zeros(length(param1),length(param2));
vec_AvePeriod = zeros(length(param1),length(param2));

% integrate the system of ode's:
for i = 1:length(param1)
    k = param1(i);
    for j = 1:length(param2)
        lambda = param2(j);
        [T,S] = ode23s(@deRHS, tspan, s0, odeset('maxstep',1));  
        G = k*(S(:,2*N+1)+0); % GFR
        
        if length(G)<6001
            continue
        end

        % mtx_GFR(:,i,j)=G;
        [pks,plocs] = findpeaks(G(end-T_cutoff:end,end));
        [vks,vlocs] = findpeaks(-G(end-T_cutoff:end,end));
        vec_Amplitude(i,j) = pks(end)+vks(end);
        vec_AvePeriod(i,j) = mean(diff(plocs))/10; % unit: sec
    end
end

save('TGF_Model_v1_ParamSweep2.mat','mtx_GFR','vec_Amplitude','vec_AvePeriod')

%% Visualizations
% load('TGF_Model_v1_ParamSweep.mat')
% tspan = tspan/10; 
% param=flip(0.7:0.01:0.97); % parameter range for k
% param=flip(0.007:0.0005:0.02); % parameter range for lambda

figure;
h=heatmap(param1,param2,vec_Amplitude,'GridVisible','off');
h.Colormap=jet;
title('Parameter Sweep (k and \lambda)');
xlabel('Parameter Value (k)'); ylabel('Parameter Value (\lambda)');
XLabels=1:length(param1);
CustomXLabels = string(param1);
CustomXLabels(mod(XLabels,5) ~= 1) = " ";
h.XDisplayLabels = CustomXLabels;
YLabels=1:length(param2);
CustomYLabels = string(param2);
CustomYLabels(mod(YLabels,5) ~= 1) = " ";
h.YDisplayLabels = CustomYLabels;
% clim(h,[0 46]);
exportgraphics(gcf,'ParamSweep_Amplitude-K&Lambda.png','Resolution',200)

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

