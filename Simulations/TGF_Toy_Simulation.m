%% code to simulate advection process via Method of lines 
function TGF_Toy_Simulation
global k lambda sigma alpha beta cop n
 
set(0,                          ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);
  
% parameters
k = 0.9;  % time constant
% lambda = 0.01; % reabsorption rate in PT
sigma = 0.1; % reabsorption rate in TAL
alpha = 0.5; % sensitivity of the negative feedback
beta = 0.5; % shift of the sigmoid function
cop = 0.3; % target [Na+] at MD
n = 4; % slope of the sigmoid function

% modify parameters
lambda = 0.01*0.58;
 
% initial condition
N = 100; % number of cells in the 1-D domain (PT and TAL)  
g0 = 30;
% u0 = zeros(N,1); u0(1) = 1; % a single bolus
u0 = g0*ones(N,1); % uniform fluid volume in PT
c0 = ones(N,1); % uniform sodium concentration in TAL
s0 = [u0; c0; g0*1.01]; 
 
% simulation time step and duration
tstep = 1;  % time step interval
t_end = 20000; % total time to run simulation
tspan = 0:tstep:t_end;

% integrate the system of ode's:
[T,S] = ode23s(@deRHS, tspan, s0, odeset('maxstep',1));  
V = S(:,1:N); % fluid volume in PT
C = S(:,N+1:2*N); % sodium concentration in TAL
G = S(:,2*N+1); % GFR

save('TGF_Toy_Sim_lambda_0.58_g0_101%.mat','T','S','G')

%% plot the solution at each time step in PT
% figure(1)
% for j = 1:length(T)
%     plot(V(j,:),'*')
%     ylabel('V(n,t)','fontsize',20)
%     xlabel('n','fontsize',20)
%     axis([0 N 0 60])
%     title('Simulation of fluid volume along PT over time', sprintfc('time step = %i',j))
%     pause(.02)
% end

%% plot the solution at each time step in TAL
% figure(2)
% for j = 1:length(T)
%     plot(C(j,:),'*')
%     ylabel('C(n,t)','fontsize',20)
%     xlabel('n','fontsize',20)
%     axis([0 N 0 1])
%     title('Simulation of Na+ conc. along TAL over time', sprintfc('time step = %i',j))
%     pause(.02)
% end

%% export as gif files
% gifFile = 'TGF_Toy_PT_Flow_k0.9_lambda0.01_sigma0.1_alpha0.5_beta0.5_cop0.3_n4.gif';
% obj = figure;
% exportgraphics(obj, gifFile);
% 
% for j = 1:length(T)
%     if mod(j,10)==0
%         plot(V(j,:),'*')
%         ylabel('V(n,t)','fontsize',20)
%         xlabel('n','fontsize',20)
%         axis([0 N 0 60])
%         title('Simulation of fluid volume along PT over time', sprintfc('time step = %i',j))
%         exportgraphics(obj, gifFile, Append=true);
%     end
% end

% gifFile = 'TGF_Toy_TAL_Conc_k0.9_lambda0.01_sigma0.1_alpha0.5_beta0.5_cop0.3_n4.gif';
% obj = figure;
% exportgraphics(obj, gifFile);
% 
% for j = 1:length(T)
%     if mod(j,10)==0
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
% figure(5)
% plot(C(:,N),'*')
% ylabel('C(n,t)','fontsize',20)
% xlabel('time step','fontsize',20)
% axis([0 t_end 0 1])
% title('Sodium conc. at TAL terminal')

%% plot GFR
figure(6)
plot(G(:),'*')
ylabel('GFR (nl/min)','fontsize',20)
xlabel('time step','fontsize',20)
axis([0 t_end 0 70])
title('GFR [\lambda(fold) = 0.58]')

% (furosemide affects C_{op} & n)

%% plot GFR with delay embedding
% G_1=G(1:end-125);
% G_2=G(126:end);
% figure(7)
% plot(G_2,G_1)
% ylabel('GFR(t+\tau) (nl/min)','fontsize',20);
% xlabel('GFR(t) (nl/min)','fontsize',20)
% title('GFR delay embedding [\lambda(fold) = 0.58]')
end
 
%% the right hand side for ode simulation:
function s_prime = deRHS(t,s)
    global k lambda sigma alpha beta cop n

    % furosemide treatment:
    % if t>2000 && t<4000
    %     n=2;
    %     cop = 0.6;
    % else
    %     n=4;
    %     cop = 0.3;
    % end

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

