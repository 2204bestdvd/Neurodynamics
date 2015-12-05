function [I, dmdt, dhdt, dndt, dslowdt] = LN2(V, m, h, n, slow)
% Constants
g_Na = 20.0; % maximum conducances, in mS/cm^2
g_K  = 9.0;
g_L  = 8;
g_slow = 10;
E_Na = 60.0; % Nernst reversal potentials, in mV
E_K  = -80.0;
E_L  = -80;
E_slow  = -80.0;

%Functionalize Single HH Neuron
% alpha_m = 0.1*(V+45) ./ (1-exp(-(V+45)/10));
alpha_h = 0.07 * exp(-(V+70)/20);
%alpha_n = 0.01*(V+60) ./ (1-exp(-(V+60)/10));
% beta_m = 4 * exp(-(V+70)/18);
beta_h = 1 ./ (1+exp(-(V+40)/10));
%beta_n = 0.125 * exp(-(V+70)/80);

% m = alpha_m ./ (alpha_m + beta_m);
m_inf = 1 ./ (1 + exp(-(V+20)/15));
m = m_inf;
n_inf = 1 ./ (1 + exp(-(V+25)/5));
tau_n = 0.152;
slow_inf = 1 ./ (1 + exp(-(V+10)/5));
tau_slow = 40;

% Membrane currents (in uA/cm^2)
%I_Na = g_Na * (m.^3) .* h .* (E_Na - V);
%I_K = g_K  * (n.^4) .* (E_K - V);
I_Na = g_Na * (m.^1) .* (E_Na - V);
I_K = g_K  * (n.^1) .* (E_K - V);
I_slow = g_slow  * (slow.^1) .* (E_slow - V);
I_L = g_L * (E_L - V);
I = I_Na + I_K + I_slow + I_L;

% Diffeq's for gating variables
% dmdt = alpha_m.*(1-m) - beta_m.*m;
dmdt = 0;
dhdt= alpha_h.*(1-h) - beta_h.*h;
%dndt = alpha_n.*(1-n) - beta_n.*n;
dndt = (n_inf - n) ./ tau_n;
dslowdt = (slow_inf - slow) ./ tau_slow;
