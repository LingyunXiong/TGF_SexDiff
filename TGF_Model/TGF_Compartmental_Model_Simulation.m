%% Simulation the effect of TGF on SNGFR in the kidney 
function TGF_Compartmental_Model_Simulation
global k lambda sigma alpha beta cop n v_unit N M 
 
set(0,                          ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);
  
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
v0 = 10; % Initial value of fluid volume (relative to unit volume)
u0 = v0*ones(N,1); 
ci = 1; % Initial value of [Na+] (relative to plasma [Na+])
c0 = ci*ones(M,1); 
s0 = [u0; c0; v0]; 
 
% simulation time step and duration
tstep = 1;  % Unit time: 0.1 second
t_end = 6000; % Simulate for a total of 600 seconds (10 minutes)
tspan = 0:tstep:t_end;

% integrate the system of ode's:
[T,S] = ode15s(@deRHS, tspan, s0, odeset('maxstep',1));  
V = S(:,1:N); % Fluid volume in all PT segments (relative to unit volume)
C = S(:,N+1:N+M); % Sodium concentration in TAL (relative to plasma [Na+])
G = k*S(:,N+M+1); % SNGFR (relative to unit volume)
GFR = 0.3*k*S(:,N+M+1); % SNGFR (nl/min)

time = T/10; % Unit: 1 second

%% plots
figure(1)
plot(time,GFR(:),'.k',LineWidth=5); hold on
ylabel('GFR (nl/min)','fontsize',20)
xlabel('Time (second)','fontsize',20)
axis([0 600 0 12])
hold off

% calculate fractional salt reabsorption
T_cutoff = 1000;
Q = cumtrapz(time(end-T_cutoff:end,end),G(end-T_cutoff:end,end));
Sin=ci*Q(end);

A = ci*k*V(:,end);
Asum = cumtrapz(time(end-T_cutoff:end,end),A(end-T_cutoff:end,end));
Sout_PT = Asum(end);

B = k*V(:,end).*C(:,end);
Bsum = cumtrapz(time(end-T_cutoff:end,end),B(end-T_cutoff:end,end));
Sout_TAL = Bsum(end);

Frac_PT = 1-Sout_PT/Sin
Frac_Total = 1-Sout_TAL/Sin;
Frac_TAL = Frac_Total - Frac_PT;

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

