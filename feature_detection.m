function feature_detection()
% Constants
C_m  = 1.0; % membrane capacitance, in uF/cm^2

% System definition
N = 5;
sim_time = 100;
step = 0.05;
time = 0:step:sim_time;

V = zeros(N, length(time));
m = zeros(N, length(time));
h = zeros(N, length(time));
n = zeros(N, length(time));
slow = zeros(N, length(time));
V(:, 1) = -70;
m(:, 1) = 0.053;
h(:, 1) = 0.596;
n(:, 1) = 0.318;
slow(:, 1) = 0;
dmdt = zeros(N, 1);
dhdt = zeros(N, 1);
dndt = zeros(N, 1);
dslowdt = zeros(N, 1);
I_self = zeros(N, 1);

r = zeros(N);
%drdt = zeros(N);
topo_e = zeros(N);
topo_i = zeros(N);
topo_e(1, 2:3) = 1;
topo_e(5, 3) = 1;
topo_e(3, 4) = 1;
topo_i(2, 4:5) = 1;

%I_ext = (10+10*rand(N,1))*ones(1, length(time));
%I_ext = [20; 0] * ones(1, length(time));
%I_ext(1, 1:1000) = 0;
I_ext = zeros(N, length(time));
I_ext(1, 1001:1200) = 40;
% I_ext(3, 1001:1300) = [linspace(1, 10, 150), linspace(10, 1, 150)];
% I_ext(1, 1501:1700) = 40;
% I_ext(5, 1001:1200) = -10;

% Parameters
%g_gaba = 20;  % mS/cm^2
%g_glu = 8;
%g = [10, 8, 10, 20, 1, 5];
g = zeros(N);
g(1, 2:3) = [3, 1];
g(2, 4:5) = [1, 5];
g(3, 4) = 3;
g(5, 3) = 8;

for t = 1:length(time)-1
%    system('pause');
	%[I_self, dmdt, dhdt, dndt, dslowdt] = HH_burst(V(:,t), m(:,t), h(:,t), n(:,t), slow(:,t));
% 	[I_self(1), dmdt(1), dhdt(1), dndt(1)] = HH(V(1,t), m(1,t), h(1,t), n(1,t));
	[I_self(1), dmdt(1), dhdt(1), dndt(1), dslowdt(1)] = AN1(V(1,t), m(1,t), h(1,t), n(1,t), slow(1,t));
% 	[I_self(2), dmdt(2), dhdt(2), dndt(2)] = HH(V(2,t), m(2,t), h(2,t), n(2,t));
	[I_self(2), dmdt(2), dhdt(2), dndt(2), dslowdt(2)] = LN2(V(2,t), m(2,t), h(2,t), n(2,t), slow(2,t));
	[I_self(3), dmdt(3), dhdt(3), dndt(3), dslowdt(3)] = LN3(V(3,t), m(3,t), h(3,t), n(3,t), slow(3,t));
	[I_self(4), dmdt(4), dhdt(4), dndt(4), dslowdt(4)] = HH_burst(V(4,t), m(4,t), h(4,t), n(4,t), slow(4,t));
	[I_self(5), dmdt(5), dhdt(5), dndt(5), dslowdt(5)] = LN5(V(5,t), m(5,t), h(5,t), n(5,t), slow(5, t));
	%[I_self, dmdt, dhdt, dndt] = HH(V(:,t), m(:,t), h(:,t), n(:,t));
	%[I_self, dmdt, dhdt, dndt] = Oscillate(V(:,t), m(:,t), h(:,t), n(:,t));
	%[I_syn(1), drdt(2)] = inhibitory(V(2,t), V(1,t), g_gaba, r(2));
	[I_syn_e, drdt_e] = excitatory(topo_e, V(:,t), g, r);
	[I_syn_i, drdt_i] = inhibitory(topo_i, V(:,t), g, r);
    I_syn = sum(I_syn_e + I_syn_i, 1)';
    drdt = drdt_e + drdt_i;
% 	[I_syn(1), drdt(1)] = excitatory(V(1,t), V(2,t), g(1), r(1));
% 	[I_syn(2), drdt(2)] = excitatory(V(1,t), V(3,t), g(2), r(2));
% 	[I_syn(3), drdt(3)] = inhibitory(V(2,t), V(5,t), g(3), r(3));
% 	[I_syn(4), drdt(4)] = inhibitory(V(2,t), V(4,t), g(4), r(4));
% 	[I_syn(5), drdt(5)] = excitatory(V(5,t), V(3,t), g(5), r(5));
% 	[I_syn(6), drdt(6)] = excitatory(V(3,t), V(4,t), g(6), r(6));
	%I_total = I_ext(:,t) + I_self; 
	I_total = I_ext(:,t) + I_self + I_syn; 
	dVdt = I_total / C_m;
	V(:,t+1) = V(:,t) + step*dVdt;
	m(:,t+1) = m(:,t) + step*dmdt;
	h(:,t+1) = h(:,t) + step*dhdt;
	n(:,t+1) = n(:,t) + step*dndt;
	slow(:,t+1) = slow(:,t) + step*dslowdt;
    
    % synaptic parameters update
    r = r + step*drdt;
end

figure; plot(time, V);
%figure; plot(time, m(5,:))
%figure; plot(time, slow);



function [I_syn, drdt] = inhibitory(topo, V, g_gaba, r)
% Constants
E_Cl = -80;  % mV
alpha_r = 5;  % mM^-1ms^-1
beta_r = 0.18;  % ms^-1
T_max = 1.5;  % mM
K_p = 5;  % mV
V_p = 7;  %mV

V_post = ones(length(V),1) * V';
V_pre = V * ones(1,length(V));

I_syn = g_gaba .* r .* (E_Cl - V_post) .* topo;
T = T_max ./ (1 + exp(-(V_pre - V_p)/K_p)) .* topo;
drdt = (alpha_r * T * (1-r) - beta_r * r) .* topo;

function [I_syn, drdt] = excitatory(topo, V, g_glu, r)
% Constants
E = -38;  % mV
alpha_r = 2.4;  % mM^-1ms^-1
beta_r = 0.56;  % ms^-1
T_max = 1.0;  % mM
K_p = 5;  % mV
V_p = 7;  %mV

V_post = ones(length(V),1) * V';
V_pre = V * ones(1,length(V));

I_syn = g_glu .* r .* (E - V_post) .* topo;
T = T_max ./ (1 + exp(-(V_pre - V_p)/K_p));
drdt = (alpha_r * T * (1-r) - beta_r * r) .* topo;


function [I, dmdt, dhdt, dndt, dslowdt] = HH_burst(V, m, h, n, slow)

% Constants
C_m  = 1.0; % membrane capacitance, in uF/cm^2
g_Na = 20.0; % maximum conducances, in mS/cm^2
g_K  = 9.0;
g_L  = 8;
g_slow = 5;
E_Na = 60.0; % Nernst reversal potentials, in mV
E_K  = -90.0;
E_L  = -80;
E_slow  = -90.0;

%Functionalize Single HH Neuron
%alpha_m = 0.1*(V+45) ./ (1-exp(-(V+45)/10));
alpha_h = 0.07 * exp(-(V+70)/20);
%alpha_n = 0.01*(V+60) ./ (1-exp(-(V+60)/10));
%beta_m = 4 * exp(-(V+70)/18);
beta_h = 1 ./ (1+exp(-(V+40)/10));
%beta_n = 0.125 * exp(-(V+70)/80);

%m = alpha_m ./ (alpha_m + beta_m);
m_inf = 1 ./ (1 + exp(-(V+20)/15));
m = m_inf;
n_inf = 1 ./ (1 + exp(-(V+25)/5));
tau_n = 0.152;
slow_inf = 1 ./ (1 + exp(-(V+50)/5));
tau_slow = 20;

% Membrane currents (in uA/cm^2)
%I_Na = g_Na * (m.^3) .* h .* (E_Na - V);
%I_K = g_K  * (n.^4) .* (E_K - V);
I_Na = g_Na * (m.^1) .* (E_Na - V);
I_K = g_K  * (n.^1) .* (E_K - V);
I_slow = g_slow  * (slow.^1) .* (E_slow - V);
I_L = g_L * (E_L - V);
I = I_Na + I_K + I_slow + I_L;

% Diffeq's for gating variables
%dmdt = alpha_m.*(1-m) - beta_m.*m;
dmdt = 0;
dhdt= alpha_h.*(1-h) - beta_h.*h;
%dndt = alpha_n.*(1-n) - beta_n.*n;
dndt = (n_inf - n) ./ tau_n;
dslowdt = (slow_inf - slow) ./ tau_slow;


function [I, dmdt, dhdt, dndt] = Oscillate(V, m, h, n)

% Constants
C_m  = 1.0; % membrane capacitance, in uF/cm^2
g_Na = 30.0; % maximum conducances, in mS/cm^2
g_K  = 24.0;
g_L  = 0.3;
E_Na = 45.0; % Nernst reversal potentials, in mV
E_K  = -82.0;
E_L  = -59.387;

%Functionalize Single HH Neuron
alpha_m = 0.1*(V+45) ./ (1-exp(-(V+45)/10));
alpha_h = 0.07 * exp(-(V+70)/20);
alpha_n = 0.01*(V+60) ./ (1-exp(-(V+60)/10));
beta_m = 4 * exp(-(V+70)/18);
beta_h = 1 ./ (1+exp(-(V+40)/10));
beta_n = 0.125 * exp(-(V+70)/80);


% Membrane currents (in uA/cm^2)
I_Na = g_Na * (m.^3) .* h .* (E_Na - V);
I_K = g_K  * (n.^4) .* (E_K - V);
I_L = g_L * (E_L - V);
I = I_Na + I_K + I_L;

% Diffeq's for gating variables
dmdt = alpha_m.*(1-m) - beta_m.*m;
dhdt= alpha_h.*(1-h) - beta_h.*h;
dndt = alpha_n.*(1-n) - beta_n.*n;

