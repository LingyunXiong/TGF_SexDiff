%% Simulation the effect of TGF on SNGFR in the kidney 
function Figure3B
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
c=0:0.01:1;
conc=150*c; % unit: mM 

n = 2; % slope of the sigmoid function
Fg1 = alpha*(cop^n./(cop^n+c.^n)-beta);

n = 3; % slope of the sigmoid function
Fg2 = alpha*(cop^n./(cop^n+c.^n)-beta);

n = 4; % slope of the sigmoid function
Fg3 = alpha*(cop^n./(cop^n+c.^n)-beta);

figure(1)
plot(conc,Fg1(:),"LineWidth",3); hold on
plot(conc,Fg2(:),"LineWidth",3)
plot(conc,Fg3(:),"LineWidth",3)
ylabel('\DeltaGFR / \Deltat','fontsize',16)
xlabel('[Na+]_{MD} (mM)','fontsize',16)
axis([0 60 -0.25 0.25])
yline(0)
legend('n=2','n=3','n=4')
hold off

figure(2)
plot(time,GFR(:),'.k',LineWidth=5); hold on
ylabel('GFR (nl/min)','fontsize',20)
xlabel('Time (second)','fontsize',20)
axis([0 600 0 10])
hold off

end
 
%% the right hand side for ode simulation:
function s_prime = deRHS(t,s)
    global k lambda sigma alpha beta cop n v_unit N M 

    % Hill Coef:
    if t>3000 && t<6000
        n=2;
    else
        n=4;
    end

    v = s(1:N);
    c = s(N+1:N+M);
    g = s(end);

 
    % Changes in fluid volume in the PT compartment
    Fv = k*([g;v(1:end-1)]-v) - lambda*v; 

    % Changes in sodium concentration in the TAL compartment
    Fc = k*v(end)*([1;c(1:end-1)]-c)/v_unit - sigma*c; 

    % Negative feedback from MD to AA
    Fg = alpha*(cop^n/(cop^n+c(end)^n)-beta); % sigmoid function
    
    s_prime = [Fv; Fc; Fg];
end

