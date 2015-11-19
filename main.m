% Constants
C_m  = 1.0; % membrane capacitance, in uF/cm^2

% System definition
N = 5;
sim_time = 300;
step = 0.05;
time = 0:step:sim_time;
V = zeros(N, length(time));
% V_nmda = zeros(N, length(time));
V_synapse = zeros(N);
V_nmda = zeros(N);
Ca = zeros(N);
%Ca = 1e-9 * ones(N) .* topo;  % Assume Ca concentration independent across synapses
m = zeros(N, length(time));
h = zeros(N, length(time));
n = zeros(N, length(time));
I_ext = (10+10*rand(N,1))*ones(1, length(time));
topo = eye(N);  % N x N matrix indicating network topology, ex. topo(i,j)=1 for i -> j
topo = topo(circshift((1:N)', -1), :);  % test: circular 1 -> 2-> ... -> N -> 1

% Parameters
% Assume peak synaptic conductance, g, fixed (over time and across synapses) at this point
g = [0.01, 0.005, 0.005];  % g_fb, g_nmda, g_vgcc
g_ampa = 0.1 * ones(N) .* topo;
r_ampa = ones(N) .* topo;
r_nmda = ones(N) .* topo;
r = cat(3, r_ampa, r_nmda);
thresh = -20;  % threshold beyond which neuron is considered to have fired

for t = 1:length(time)-1
%    system('pause');
	[I_self, dmdt, dhdt, dndt] = HH(V(:,t), m(:,t), h(:,t), n(:,t));
	[I_syn, dVdt_nmda, dVdt_synapse, dCadt, dgdt_ampa, drdt] = synapse(topo, V(:,t), g, g_ampa, r, V_nmda, Ca);
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
	Ca = Ca + step*dCadt;
	g_ampa = g_ampa + step*dgdt_ampa;

	% If presynaptic spike occurs, set r_ampa and r_nmda to unity, otherwise update
	drdt_ampa = drdt(:,:,1);
	drdt_nmda = drdt(:,:,2);
	check_spike = ((V(:,t+1) > thresh) .* (V(:,t) < thresh)) * ones(1,N);
	r_ampa = (check_spike + (r_ampa + step*drdt_ampa) .* (1 - check_spike)) .* topo;
	r_nmda = (check_spike + (r_nmda + step*drdt_nmda) .* (1 - check_spike)) .* topo;
end

figure; plot(time, V);

