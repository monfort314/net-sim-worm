function out = simulateNetwork(Iext, T, Tw, s, h, x0, resetIdx, perturbationIdx, synapseAdjustIdx)
% out = simulateNetwork(Iext, T, Tw, s, h)
% runs simulation with external input applied, 
% outputs voltages as function
% of time and plots voltages and modes of motorneurons
% 
% INPUT:
% Iext:         [279 x N] input current over the time to all the neurons,
% T:            Total simulation time, 
% s:            sampling period,
% h:            simulation timestep, 
% Tw:           and time to wait before recording
%
% OUTPUTS:
% 'out':        struct, fields are motor neuron voltages 'V' 
%               and timestep vector 't'

arguments
    Iext = 100 * ones(279, 1);
    T = 20;
    Tw = 15;
    s = 0.001;
    h = 1e-4;
    x0 = [];
    resetIdx = -1;
    perturbationIdx = -1;
    synapseAdjustIdx = -1;
end
%% Add folder for subfunctions
addpath 'Varshney';
addpath 'CEresponsefunctions';
addpath 'subfunctions_simulation';
%make s integer multiple of h, T integer mult of s
T = s * floor(T / s);
h = s / floor(s / h);

% number of timesteps kept
nts = int32((T - Tw) / s);
%initialize time, voltage vector
t = 0;
ts = linspace(Tw, T, nts);

%% Define Input/Stimulus
%Amplitude of input:
%Amp = 2e4;
%Vector for Neurons receiving input:
if size(Iext, 2) == 1
    Iext = repmat(Iext, 1, nts);
end



%indices of neurons receiving random input and amplitude:
nosN = 1:279;
nosA = 0.00;%1;

InpF=@(vf,tf) Iext(:, tf);

%% Import Model Parameters
%Model parameters Structure
modcon = open('modcon.mat');

modcon.tau = 0.01;
modcon.gelec = 1;
modcon.gchem = 1;
modcon.memG = 0.1;
modcon.memV = -35;
modcon.EchemEx = 0;
modcon.EchemInh = -48;
modcon.beta = 0.125;
modcon.ar = 1;
modcon.ad = 5.0;
sig = 0.5 * ones(279,1);



%Chemical and Gap Adjacency Matrices from Varshney

Ag = datareader('gap','weighted');
Ag = full(Ag);
Ac = datareader('chem','weighted');
Ac = full(Ac);
Dg = diag(sum(Ag));
L = Dg - Ag;
%treat GABAergic synapses in chemical network as inhibitory
GABAergic = GABA;
Ec = ones(size(Ac,1),1) * modcon.EchemEx;
Ec(find(GABAergic)) = modcon.EchemInh;

N = size(Ag, 1);


%Calculate Equilibrium/vmean
eqS = diag((modcon.ar) ./ ((modcon.ar) * (sig) + modcon.ad)) * (sig);
C = modcon.memG * eye(size(Ac)) + modcon.gelec * L + modcon.gchem * diag(Ac' * eqS);
b = modcon.memG * modcon.memV + modcon.gchem * Ac' * (eqS .* Ec) + Iext(:, 1);
eqV = C \ b;
xin = [eqV; eqS];
vmean = eqV + 1 ./ modcon.beta .* log(1 ./ sig - 1);
modcon.vmean = vmean;

calculated = 0;

function y = resetValues
    if calculated
        y = xin .* (1 + 0.01 * randn(size(xin)));
    else
        y0 = load('xCycle_PLM2e4.mat');
        y = y0.xCycle;
        vmean = eqV + 1 ./ modcon.beta .* log(1 ./ sig - 1);
        modcon.vmean = vmean;
    end
end

function params_new = adjustSynapses(I, params)
    eqS = diag((params.ar) ./ ((params.ar) * (sig) + params.ad)) * (sig);
    C = params.memG * eye(size(Ac)) + params.gelec * L + params.gchem * diag(Ac' * eqS);
    b = params.memG * params.memV + params.gchem * Ac' * (eqS .* Ec) + I;
    eqV = C \ b;
    vmean = eqV + 1 ./ modcon.beta .* log(1 ./ sig - 1);
    params_new = params;
    params_new.vmean = vmean;
end

function y = perturbation(y)
    y = y .* (1 + 0.02 * rand(size(y)));
end

if isempty(x0)
   x0 = resetValues; 
end

%if isempty(x0)
%    x0 = load('xCycle_PLM2e4.mat');
%    x0 = x0.xCycle;
%end
x = x0;

%% Define EM1 Timestep Function
F = @(t, v, Ink, modcon) Model2(v, Ac, Ec, L, Ink, modcon, Ag);
    function [tnn, vnn, dv, out] = timestepEM(tn, vn, sigma, ind, Ink, i, modcon)
        dW = zeros(size(vn));
        dW(ind) = sqrt(h) * randn(size(ind));
        [dv, out] = F(tn, vn, Ink(vn, i), modcon);
        vnn = vn + h * dv + sigma * dW;
        tnn = tn + h;
end

%% Perform Simulation keeping every s/h timesteps
X = zeros(N*2, nts);
dX = X;


i = 0;
tic
fprintf('\n%03.2f', 0);
while i < nts
    [t, x, dx, ~] = timestepEM(t, x, nosA, nosN, InpF, i + 1, modcon);
    if (mod(t, s) < h)
        fprintf('\b\b\b\b%03.2f', (t / T));
        if (t >= Tw)
            i = i + 1;
            X(:, i) = x;
            dX(:, i) = dx;
            if resetIdx == i
                x = resetValues;
                disp('resetting values')
            end
            if perturbationIdx == i
                x = perturbation(x);
                disp('perturbing values')
            end
            if synapseAdjustIdx == i
                modcon = adjustSynapses(Iext(:, i), modcon);
            end
            
        end
    end
end
fprintf('\n');
disp('Calculation Complete');
toc;

%Keep only Voltages
V = X(1 : N, :);
% Subtract out mean
V = V - repmat(xin(1 : N), 1, size(V, 2));
XX = X - repmat(xin, 1, size(X, 2));

%% Define Output

out = struct;
out.V = V;
out.t = ts;
out.X = X;
out.XX = XX;
out.Xs = XX;
end