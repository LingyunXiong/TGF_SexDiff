%% code to simulate advection process via Method of lines 
function TGF_Toy_Model_v1
global k lambda sigma alpha beta cop n
 
set(0,                          ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);
  
% parameters
k = 0.72;  % time constant (female)
% lambda = 0.015; % reabsorption rate in PT (female)
lambda = 0.006; % reabsorption rate in PT (female): lower fractional resorption
sigma = 0.1; % reabsorption rate in TAL (neutral)
alpha = 0.5; % sensitivity of the negative feedback
beta = 0.5; % shift of the sigmoid function
cop = 0.3; % target [Na+] at MD
n = 3; % slope of the sigmoid function

% modify parameters
% k = 0.98; % time constant (male)
% lambda = 0.02; % reabsorption rate in PT (male)
% sigma = 0.12; % reabsorption rate in TAL (male)
% n = 4; % slope of the sigmoid function (male)
% alpha = 1;
% cop = 0.15;

% initial condition
N = 68; % number of cells in the 1-D domain (PT and TAL)  
% N = 100; % Original N %
g0 = 30;
u0 = g0*ones(N,1); % uniform fluid volume in PT
ci = 1;
% high-salt condition
% ci = 1*1.5;
% low-salt condition
% ci = 1*0.5;
c0 = ci*ones(N,1); % uniform sodium concentration in TAL
s0 = [u0; c0; g0]; 
 
% simulation time step and duration
tstep = 1;  % time step interval
t_end = 2000; % total time to run simulation
tspan = 0:tstep:t_end;

% integrate the system of ode's:
[T,S] = ode23s(@deRHS, tspan, s0, odeset('maxstep',1));  
V = S(:,1:N); % fluid volume in PT
C = S(:,N+1:2*N); % sodium concentration in TAL
G = k*(S(:,2*N+1)+0); % GFR

% save('TGF_Toy_v1_Sim.mat','T','S','G')

tspan=T/10;

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
%     axis([0 N 0 1.5])
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
end

%% plot the timeseries across PT segment
% figure(3)
% plot(V(:,N/2),'*')
% ylabel('V(n,t)','fontsize',20)
% xlabel('time step','fontsize',20)
% axis([0 t_end 0 60])
% title('Fluid volume at mid-PT segment')

% figure(4)
% plot(V(:,end),'*')
% ylabel('V(n,t)','fontsize',20)
% xlabel('time step','fontsize',20)
% axis([0 t_end 0 30])
% title('Fluid volume at end-PT segment')

%% plot GFR
figure(5)
plot(tspan,G(:),'*'); hold on
ylabel('GFR(t) (nl/min)','fontsize',20)
xlabel('t (s)','fontsize',20)
axis([0 tspan(end) 0 40])
title('SNGFR')
hold off

% (furosemide affects C_{op} & n)
% [\lambda(fold) = 0.58]
% (n: 4 \rightarrow 2)
% [k=0.98 \lambda=0.02 p=30s]
% [k=0.98 \lambda=0.02 \sigma=0.12 p=30s]

%% plot the timeseries of Na+ conc. at the end of TAL
% figure(6)
% plot(tspan,C(:,end),'*')
% ylabel('C(n,t) (a.u.)','fontsize',20)
% xlabel('t (s)','fontsize',20)
% axis([0 tspan(end) 0 1])
% title('[Na+] at MD')

%% plot GFR with delay embedding
% G_1=G(1:end-125);
% G_2=G(126:end);
% figure(7)
% plot(G_2,G_1)
% ylabel('GFR(t+\tau) (nl/min)','fontsize',20);
% xlabel('GFR(t) (nl/min)','fontsize',20)
% title('GFR delay embedding')

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

%%% calculate cumulative salt reabsorption
T_cutoff = 1000;
tspan=T/10;
Q = cumtrapz(tspan(end-T_cutoff:end,end),G(end-T_cutoff:end,end));
Sin=ci*Q(end);

A = ci*k*V(:,end);
Asum = cumtrapz(tspan(end-T_cutoff:end,end),A(end-T_cutoff:end,end));
Sout_PT = Asum(end);

B = k*V(:,end).*C(:,end);
Bsum = cumtrapz(tspan(end-T_cutoff:end,end),B(end-T_cutoff:end,end));
Sout_TAL = Bsum(end);

Frac_PT = 1-Sout_PT/Sin
Frac_TAL = 1-Sout_TAL/Sin
% end
 
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

