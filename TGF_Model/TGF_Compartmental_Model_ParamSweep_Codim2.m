%% Simulation the effect of TGF on SNGFR in the kidney  
function TGF_Compartmental_Model_ParamSweep_Codim2
global k lambda sigma alpha beta cop n v_unit N M
 
% baseline parameters
k = 0.92;  % PT fluid flow rate constant 
lambda = 0.024; % PT sodium reabsorption rate
sigma = 0.1; % TAL sodium reabsorption rate
alpha = 0.4; % Scaling factor for the sigmoid function
beta = 0.5; % Vertical shift of the sigmoid function (dimensionless)
cop = 0.2; % Operational [Na+] at MD (relative to plasma [Na+])
n = 3; % Maximal slope of the sigmoid function
N = 50; % Number of PT compartments
M = 60; % Number of TAL compartments
v_unit = 1; % Unit volume (5e-3 nl)

% initial condition
v0 = 10; % initial value of fluid volume (relative to unit volume)
u0 = v0*ones(N,1); 
ci = 1; % initial value of [Na+] (relative to plasma [Na+])
c0 = ci*ones(M,1); 
s0 = [u0; c0; v0]; 

%% Parameter Sensitivity Test
% set the range of parameter values (Refine)
param1=flip(0.4:0.005:1); % parameter range for k
param2=flip(0.005:0.00025:0.04); % parameter range for lambda

% simulation time step and duration
tstep = 1;  % unit time: 0.1 sec
t_end = 6000; % simulate for a total of 10 minutes
tspan = 0:tstep:t_end;

mtx_G = zeros(length(tspan),length(param1),length(param2));
vec_Frac_PT = zeros(length(param1),length(param2));
vec_Frac_TAL = zeros(length(param1),length(param2));

% integrate the system of ode's:
for i = 1:length(param1)
    k = param1(i);
    for j = 1:length(param2)
        lambda = param2(j);
        [T,S] = ode15s(@deRHS, tspan, s0, odeset('maxstep',1));  
	    V = S(:,1:N); % fluid volume in all PT segments (relative to unit volume)
	    C = S(:,N+1:N+M); % sodium concentration in TAL (relative to plasma [Na+])
	    G = k*S(:,N+M+1); % GFR (relative to unit volume)
        
        if length(G)<6001
            continue
        end

        mtx_G(:,i,j)=G;
	    T_cutoff = 1000;
	    t=T/10;
	    Q = cumtrapz(t(end-T_cutoff:end,end),G(end-T_cutoff:end,end));
	    Sin=ci*Q(end);

	    A = ci*k*V(:,end);
	    Asum = cumtrapz(t(end-T_cutoff:end,end),A(end-T_cutoff:end,end));
	    Sout_PT = Asum(end);

	    B = k*V(:,end).*C(:,end);
	    Bsum = cumtrapz(t(end-T_cutoff:end,end),B(end-T_cutoff:end,end));
	    Sout_TAL = Bsum(end);

	    vec_Frac_PT(i,j) = 1-Sout_PT/Sin;
	    vec_Frac_TAL(i,j) = Sout_PT/Sin - Sout_TAL/Sin;
    end
end

save('TGF_Compartmental_Model_ParamSweep2_Refine_N3_Output.mat','mtx_G','vec_Frac_PT','vec_Frac_TAL')

end
 
%% the right hand side for ode simulation:
function s_prime = deRHS(t,s)
    global k lambda sigma alpha beta cop n v_unit N M 

    v = s(1:N);
    c = s(N+1:N+M);
    g = s(end);

    % changes in fluid volume in the PT compartment
    Fv = k*([g;v(1:end-1)]-v) - lambda*v; 

    % changes in sodium concentration in the TAL compartment
    Fc = k*v(end)*([1;c(1:end-1)]-c)/v_unit - sigma*c; 

    % negative feedback from MD to AA
    Fg = alpha*(cop^n/(cop^n+c(end)^n)-beta); % sigmoid function
    
    s_prime = [Fv; Fc; Fg];
end

