function [I_syn, dVdt_nmda, dVdt_synapse, dCadt, dgdt_ampa, drdt] = synapse(topo, V, g, g_ampa, r, V_nmda, Ca)
% Constants
C_m  = 1.0; % membrane capacitance, in uF/cm^2
tau_ampa = 2.5;  % in ms
tau_nmda = 30;
E_L  = -59.387;
g_L = 0.3;

% Constants that need to look up (or adjust)
tau_Ca = 20;
tau_stdp = 20;
E_ampa = 0;
E_nmda = 0;
E_Ca = 110;
Mg = 1.2 * ones(size(Ca)) .* topo;  % mM

V_post = ones(length(V),1) * V';
%V_pre = V * ones(1,length(V));
%V_vgcc = ones(length(V),1) * (V + V_nmda)';  % V_vgcc = V + V_nmda
V_vgcc = V_post + V_nmda;  % V_vgcc = V + V_nmda
g_fb = g(1);
g_nmda = g(2);
g_vgcc = g(3);
r_ampa = r(:,:,1);  % receptor channel kinetics, r, for each synapse
r_nmda = r(:,:,2);

rho_nmda = 1 ./ (1 + Mg/3.57.*exp(-0.062*V_post));
rho_vgcc = 1 ./ (1 + exp(-V_vgcc/7.2));

% % Does V_soma mean V_pre? check
% V_soma = V_pre;
% I_soma = g_fb * (V_soma - V_post) .* topo;
I_soma = zeros(size(topo));
I_ampa = g_ampa .* r_ampa .* (E_ampa - V_post) .* topo;
I_nmda = g_nmda * rho_nmda * r_nmda .* (E_nmda - V_post) .* topo;
I_vgcc = g_vgcc * rho_vgcc .* (E_Ca - V_vgcc) .* topo;
I_L = g_L * (E_L - V_post) .* topo;
% Total synaptic current for each synapse
I_synapse = I_soma + I_ampa + I_nmda + I_vgcc + I_L;
% Sum synaptic current for each neuron
I_syn = sum(I_synapse, 1)'; 

drdt_ampa = -r_ampa / tau_ampa;
drdt_nmda = -r_nmda / tau_nmda;
drdt = cat(3, drdt_ampa, drdt_nmda);

% I_nmda_total = sum(I_nmda, 1)';
% I_vgcc_total = sum(I_vgcc, 1)';
dVdt_nmda = 30 * I_nmda / C_m - V_nmda / tau_nmda;
dCadt = ((I_nmda + I_vgcc) - Ca / tau_Ca) .* topo;
dVdt_synapse = I_synapse / C_m;

% Calcium based STDP
D_rate = -g_ampa;
P_rate = 2 * exp(-g_ampa.^2);
theta_D = 0.215 * g_ampa + 0.895;
theta_P = -1.363 * g_ampa + 2.626;

omega = D_rate .* (Ca*1000 > theta_D) .* (Ca*1000 <= theta_P) + P_rate .* (Ca*1000 > theta_P);
eta = 1 ./ (1 + exp(-5 * (Ca*1000 - 0.5)));
dgdt_ampa = (eta * omega / tau_stdp) .* topo;
