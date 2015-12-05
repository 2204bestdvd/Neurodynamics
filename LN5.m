function [I, dmdt, dhdt, dndt, dslowdt] = LN5(V, m, h, n, slow)

% Constants
C_m  = 1.0; % membrane capacitance, in uF/cm^2
g_Na = 50.0; % maximum conducances, in mS/cm^2
g_K  = 30.0;
g_L  = 0.3;
g_slow = 10;
E_Na = 45.0; % Nernst reversal potentials, in mV
E_K  = -82.0;
E_L  = -59.387;
E_slow  = -80.0;


%Functionalize Single HH Neuron
alpha_m = 0.1*(V+45) ./ (1-exp(-(V+45)/10));
alpha_h = 0.07 * exp(-(V+70)/20);
alpha_n = 0.01*(V+60) ./ (1-exp(-(V+60)/10));
beta_m = 4 * exp(-(V+70)/18);
beta_h = 1 ./ (1+exp(-(V+40)/10));
beta_n = 0.125 * exp(-(V+70)/80);

slow_inf = 1 ./ (1 + exp(-(V+20)/5));
tau_slow = 20;

% Membrane currents (in uA/cm^2)
I_Na = g_Na * (m.^3) .* h .* (E_Na - V);
I_K = g_K  * (n.^4) .* (E_K - V);
I_L = g_L * (E_L - V);
I_slow = g_slow  * (slow.^1) .* (E_slow - V);
I = I_Na + I_K + I_slow + I_L;

% Diffeq's for gating variables
dmdt = alpha_m.*(1-m) - beta_m.*m;
dhdt= alpha_h.*(1-h) - beta_h.*h;
dndt = alpha_n.*(1-n) - beta_n.*n;
dslowdt = (slow_inf - slow) ./ tau_slow;
