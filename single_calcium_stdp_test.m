clear; clc;
% Constants
C_m  = 1.0; % membrane capacitance, in uF/cm^2

% System definition
N = 2;
sim_time = 600;
step = 0.05;
time = 0:step:sim_time;
topo = eye(N);  % N x N matrix indicating network topology, ex. topo(i,j)=1 for i -> j
topo = topo(circshift((1:N)', -1), :);  % test: circular 1 -> 2-> ... -> N -> 1
% topo = zeros(N);
topo = [0, 1; 0, 0];

V = zeros(N, length(time));
V(:, 1) = -70;
% V_nmda = zeros(N, length(time));
V_synapse = zeros(N);
V_nmda = zeros(N);
Ca = zeros(N, N, length(time));
Ca(:,:,1) = 1e-6 * ones(N) .* topo;  % Assume Ca concentration independent across synapses
m = zeros(N, length(time));
h = zeros(N, length(time));
n = zeros(N, length(time));
m(:, 1) = 0.053;
h(:, 1) = 0.596;
n(:, 1) = 0.317;
I_ext = zeros(N, length(time));
%I_ext = (10+10*rand(N,1))*ones(1, length(time));
%I_ext = [20; 0; 0; 0; 0] * ones(1, length(time));
I_ext(1, 4000:4200) = 20;
% I_ext(1, 8000:9000) = 20;
I_ext(2, 8000:8200) = 20;

% Parameters
% Assume peak synaptic conductance, g, fixed (over time and across synapses) at this point
g = [1e-4, 1e-6, 3e-5];  % g_fb, g_nmda, g_vgcc
g_ampa = zeros(N, N, length(time));
g_ampa(:,:,1) = 1e-2 * ones(N) .* topo;
r_ampa = zeros(N) .* topo;
r_nmda = zeros(N) .* topo;
r = cat(3, r_ampa, r_nmda);
thresh = -20;  % threshold beyond which neuron is considered to have fired

debug = zeros(1, length(time));
debug2 = zeros(1, length(time));
for t = 1:length(time)-1
%    system('pause');
	[I_self, dmdt, dhdt, dndt] = HH(V(:,t), m(:,t), h(:,t), n(:,t));
	[I_syn, dVdt_nmda, dVdt_synapse, dCadt, dgdt_ampa, drdt] = synapse(topo, V(:,t), g, g_ampa(:,:,t), r, V_nmda, Ca(:,:,t));
    
    %I_total = I_ext(:,t) + I_self; 
	I_total = I_ext(:,t) + I_self + I_syn; 
	dVdt = I_total / C_m;
	V(:,t+1) = V(:,t) + step*dVdt;
	m(:,t+1) = m(:,t) + step*dmdt;
	h(:,t+1) = h(:,t) + step*dhdt;
	n(:,t+1) = n(:,t) + step*dndt;
    
    % synaptic parameters update
	%V_nmda(:,t+1) = V_nmda(:,t) + step*dVdt_nmda;
	V_nmda = V_nmda + step*dVdt_nmda;
	V_synapse = V_synapse + step*dVdt_synapse;
	Ca(:,:,t+1) = Ca(:,:,t) + step*dCadt;
	g_ampa(:,:,t+1) = g_ampa(:,:,t) + step*dgdt_ampa;

	% If presynaptic spike occurs, set r_ampa and r_nmda to unity, otherwise update
	drdt_ampa = drdt(:,:,1);
	drdt_nmda = drdt(:,:,2);
	check_spike = ((V(:,t+1) > thresh) .* (V(:,t) < thresh)) * ones(1,N);
	r_ampa = (check_spike + (r_ampa + step*drdt_ampa) .* (1 - check_spike)) .* topo;
	r_nmda = (check_spike + (r_nmda + step*drdt_nmda) .* (1 - check_spike)) .* topo;
    r = cat(3, r_ampa, r_nmda);
debug(t) = r_ampa(1,2);
debug2(t) = r_nmda(1,2);
end

figure; plot(time, V);
figure;
subplot(4,1,1);  plot(time, squeeze(sum(sum(Ca))));
subplot(4,1,2);  plot(time, squeeze(sum(sum(g_ampa))));
subplot(4,1,3);  plot(time, debug);
subplot(4,1,4);  plot(time, debug2);

