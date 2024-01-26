%% code to simulate advection process via Method of lines 
function PT_Flow_Sim_Sinusoidal
global alpha lambda sc
 
set(0,                          ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);
  
N = 100; %number of cells   

alpha = 0.8;  % time constant
p = 0.004;  % permeability coefficient
sc = p*ones(N,1); % constant reabsorption rate
lambda = 0.005; % proportional reabsorption rate
 
% initial condition
% u0 = zeros(N,1); u0(1) = 1; % a single bolus
% u0 = zeros(N,1); % uniform zero
u0 = ones(N,1); % uniform one
 
tstep = 1;  % time step interval
t_end = 500; % total time to run simulation
 
% specify the output points
tspan = 0:tstep:t_end;

% integrate the system of ode's:
[T,S] = ode23s(@deRHS,tspan, u0, odeset('maxstep',1));  

%% plot the solution at each time step
figure(1)
for j = 1:length(T)
    plot(S(j,:),'*')
    ylabel('V(n,t)','fontsize',20)
    xlabel('n','fontsize',20)
    axis([0 N 0 2])
    title('Simulation of fluid volume along PT over time', sprintfc('time step = %i',j))
    pause(.1)
end

%% export as gif files
% gifFile = 'PT_Flow_Simulation_Sinusoidal_GFR_n100_lambda_0.75%.gif';
% obj = figure;
% exportgraphics(obj, gifFile);
% 
% for j = 1:length(T)
%     plot(S(j,:),'*')
%     ylabel('V(n,t)','fontsize',20)
%     xlabel('n','fontsize',20)
%     axis([0 N 0 1.5])
%     title('Simulation of fluid volume along PT over time', sprintfc('time step = %i',j))
%     exportgraphics(obj, gifFile, Append=true);
% end

%% plot the timeseries across PT segment
% figure(2)
% plot(S(:,N/2),'*')
% ylabel('V(n,t)','fontsize',20)
% xlabel('time step','fontsize',20)
% axis([0 t_end 0 1])
% title('Fluid volume at mid-PT segment')
% 
figure(3)
plot(S(:,N),'*')
ylabel('V(n,t)','fontsize',20)
xlabel('time step','fontsize',20)
axis([0 t_end 0 1])
title('Fluid volume at end-PT segment')

%% plot the flow rate across PT segment
figure(4)
plot(gradient(S(:,N)),'*')
ylabel('v(n,t)','fontsize',20)
xlabel('time step','fontsize',20)
% axis([0 t_end 0 1])
title('Fluid velocity at end-PT segment')

return
 
%% the right hand side for ode simulation:
function s_prime=deRHS(t,s)
    global alpha lambda sc
 
    % Fu = alpha*([(1+0.5*sin(0.1*t));s(1:end-1)]-s); % without leak
    % Fu = alpha*([(1+0.5*sin(0.1*t));s(1:end-1)]-s) - lambda*s; % with linear leaky term
    Fu = alpha*([(1+0.5*sin(0.1*t));s(1:end-1)]-s) - sc; % with constant leaky term
  
    s_prime = Fu;
return

