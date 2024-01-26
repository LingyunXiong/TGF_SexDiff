%% code to simulate advection process via Method of lines 
% function TwoCompartment_Sim_Sinusoidal
global k lambda sigma
 
set(0,                          ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);
  
N = 100; % number of cells   

k = 0.8;  % time constant
% p = 0.004;  % permeability coefficient
% sc = p*ones(N,1); % reabsorption rate in PT
lambda = 0.008; % reabsorption rate in PT
sigma = 0.005; % reabsorption rate in TAL
 
% initial condition
% u0 = zeros(N,1); u0(1) = 1; % a single bolus
% u0 = zeros(N,1); % uniform zero
u0 = ones(N,1); % uniform fluid volume in PT
c0 = ones(N,1); % uniform sodium concentration in TAL
s0 = [u0; c0]; 
 
tstep = 1;  % time step interval
t_end = 500; % total time to run simulation
 
% specify the output points
tspan = 0:tstep:t_end;

% integrate the system of ode's:
[T,S] = ode23s(@deRHS,tspan, s0, odeset('maxstep',1));  
V = S(:,1:N); % fluid volume in PT
C = S(:,N+1:2*N); % sodium concentration in TAL

%% plot the solution at each time step in PT
figure(1)
for j = 1:length(T)
    plot(V(j,:),'*')
    ylabel('V(n,t)','fontsize',20)
    xlabel('n','fontsize',20)
    axis([0 N 0 1.5])
    title('Simulation of fluid volume along PT over time', sprintfc('time step = %i',j))
    pause(.1)
end

%% plot the solution at each time step in TAL
figure(2)
for j = 1:length(T)
    plot(C(j,:),'*')
    ylabel('C(n,t)','fontsize',20)
    xlabel('n','fontsize',20)
    axis([0 N 0 1])
    title('Simulation of Na+ conc. along TAL over time', sprintfc('time step = %i',j))
    pause(.1)
end

%% export as gif files
% gifFile = 'PT_Flow_Simulation_Sinusoidal_GFRx2_n100_lambda0.8%_sigma0.5%.gif';
% obj = figure;
% exportgraphics(obj, gifFile);
% 
% for j = 1:length(T)
%     if mod(j,2)==0
%         plot(V(j,:),'*')
%         ylabel('V(n,t)','fontsize',20)
%         xlabel('n','fontsize',20)
%         axis([0 N 0 2.5])
%         title('Simulation of fluid volume along PT over time', sprintfc('time step = %i',j))
%         exportgraphics(obj, gifFile, Append=true);
%     end
% end

% gifFile = 'TAL_Conc_Simulation_Sinusoidal_GFRx2_n100_lambda0.8%_sigma0.5%.gif';
% obj = figure;
% exportgraphics(obj, gifFile);
% 
% for j = 1:length(T)
%     if mod(j,2)==0
%         plot(C(j,:),'*')
%         ylabel('C(n,t)','fontsize',20)
%         xlabel('n','fontsize',20)
%         axis([0 N 0 1])
%         title('Simulation of Na+ conc. along TAL over time', sprintfc('time step = %i',j))
%         exportgraphics(obj, gifFile, Append=true);
%     end
% end

%% plot the timeseries across PT segment
% figure(3)
% plot(V(:,N/2),'*')
% ylabel('V(n,t)','fontsize',20)
% xlabel('time step','fontsize',20)
% axis([0 t_end 0 1])
% title('Fluid volume at mid-PT segment')
% 
% figure(4)
% plot(V(:,N),'*')
% ylabel('V(n,t)','fontsize',20)
% xlabel('time step','fontsize',20)
% axis([0 t_end 0 1])
% title('Fluid volume at end-PT segment')

%% plot the timeseries of Na+ conc. at the end of TAL
figure(5)
plot(C(:,N),'*')
ylabel('C(n,t)','fontsize',20)
xlabel('time step','fontsize',20)
axis([0 t_end 0 1])
title('Sodium conc. at TAL terminal')

% return
 
%% the right hand side for ode simulation:
function s_prime=deRHS(t,s)
    global k lambda sigma

    n = length(s)/2;
    v = s(1:n);
    c = s(n+1:2*n);
 
    % fluid volume in PT compartment
    % Fv = k*([(1+0.5*sin(0.1*t));v(1:end-1)]-v); % without reabsorption
    % Fv = k*([(1+0.5*sin(0.1*t));v(1:end-1)]-v) - sc; % with constant reabsorption
    Fv = k*([(1+0.5*sin(0.1*t));v(1:end-1)]-v) - lambda*v; % with linear reabsorption

    % sodium concentration in TAL compartment
    Fc = k*v(n)*([1;c(1:end-1)]-c) - sigma*c; % with linear reabsorption
  
    s_prime = [Fv; Fc];
end

