function feature_detection()

% Constants
C_m  = 1.0; % membrane capacitance, in uF/cm^2

% System definition
N = 2;
sim_time = 300;
step = 0.05;
time = 0:step:sim_time;

%topo = [0, 1; 1, 0];

V = zeros(N, length(time));
m = zeros(N, length(time));
h = zeros(N, length(time));
n = zeros(N, length(time));
V(:, 1) = -70;
m(:, 1) = 0.053;
h(:, 1) = 0.596;
n(:, 1) = 0.317;
dmdt = zeros(N, 1);
dhdt = zeros(N, 1);
dndt = zeros(N, 1);
I_self = zeros(N, 1);

r = zeros(N, 1);
drdt = zeros(N, 1);
I_syn = zeros(N, 1);

%I_ext = (10+10*rand(N,1))*ones(1, length(time));
%I_ext = [20; 0] * ones(1, length(time));
%I_ext(1, 1:1000) = 0;
I_ext = zeros(N, length(time));
I_ext(1, 1001:1200) = 50;

% Parameters
g_gaba = 3;  % mS/cm^2
g_glu = 1;

for t = 1:length(time)-1
%    system('pause');
	[I_self, dmdt, dhdt, dndt] = HH_burst(V(:,t), m(:,t), h(:,t), n(:,t));
	%[I_syn(1), drdt(2)] = inhibitory(V(2,t), V(1,t), g_gaba, r(2));
	[I_syn(2), drdt(1)] = inhibitory(V(1,t), V(2,t), g_gaba, r(1));
	%[I_syn(2), drdt(1)] = excitatory(V(1,t), V(2,t), g_glu, r(1));
	%[I_syn(3), drdt(2)] = excitatory(V(1,t), V(3,t), g_glu, r(2));
	I_total = I_ext(:,t) + I_self; 
	%I_total = I_ext(:,t) + I_self + I_syn; 
	dVdt = I_total / C_m;
	V(:,t+1) = V(:,t) + step*dVdt;
	m(:,t+1) = m(:,t) + step*dmdt;
	h(:,t+1) = h(:,t) + step*dhdt;
	n(:,t+1) = n(:,t) + step*dndt;
    
    % synaptic parameters update
    r = r + step*drdt;
end

figure; plot(time, V);



function [I_syn, drdt] = inhibitory(V_pre, V_post, g_gaba, r)
% Constants
E_Cl = -80;  % mV
alpha_r = 5;  % mM^-1ms^-1
beta_r = 0.18;  % ms^-1
T_max = 1.5;  % mM
K_p = 5;  % mV
V_p = 7;  %mV

I_syn = g_gaba * r * (E_Cl - V_post);
T = T_max / (1 + exp(-(V_pre - V_p)/K_p));
drdt = alpha_r * T * (1-r) - beta_r * r;

function [I_syn, drdt] = excitatory(V_pre, V_post, g_glu, r)
% Constants
E = -38;  % mV
alpha_r = 2.4;  % mM^-1ms^-1
beta_r = 0.56;  % ms^-1
T_max = 1.0;  % mM
K_p = 5;  % mV
V_p = 7;  %mV

I_syn = g_glu * r * (E - V_post);
T = T_max / (1 + exp(-(V_pre - V_p)/K_p));
drdt = alpha_r * T * (1-r) - beta_r * r;

function [I, dmdt, dhdt, dndt] = HH_burst(V, m, h, n)

% Constants
C_m  = 1.0; % membrane capacitance, in uF/cm^2
g_Na = 35.0; % maximum conducances, in mS/cm^2
g_K  = 9.0;
g_L  = 0.1;
E_Na = 55.0; % Nernst reversal potentials, in mV
E_K  = -90.0;
E_L  = -65;
phi = 5;

%Functionalize Single HH Neuron
alpha_m = 0.1*(V+35) ./ (1-exp(-(V+35)/10));
alpha_h = 0.07 * exp(-(V+58)/20);
alpha_n = 0.01*(V+34) ./ (1-exp(-(V+34)/10));
beta_m = 4 * exp(-(V+60)/18);
beta_h = 1 ./ (1+exp(-(V+28)/10));
beta_n = 0.125 * exp(-(V+44)/80);

m = alpha_m ./ (alpha_m + beta_m);

% Membrane currents (in uA/cm^2)
I_Na = g_Na * (m.^3) .* h .* (E_Na - V);
I_K = g_K  * (n.^4) .* (E_K - V);
I_L = g_L * (E_L - V);
I = I_Na + I_K + I_L;

% Diffeq's for gating variables
%dmdt = alpha_m.*(1-m) - beta_m.*m;
dmdt = 0;
dhdt = phi * (alpha_h.*(1-h) - beta_h.*h);
dndt = phi * (alpha_n.*(1-n) - beta_n.*n);

