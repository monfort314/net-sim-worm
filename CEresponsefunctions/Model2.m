function [dv, out] = Model2(v,Achem,Echem,Lgap,In,constants,Gel)
%Non-linear model of C. elegans Neuron from Feree 1999
%as on p. 8
%v is the voltage vector
%A is the adjacency matrix
%AEchem is adjancency times reversal potentials matrix
%L is the Laplacian matrix
%In is current injected into neuron
%constants is structure of constants with inputs shown below
%%Electrical Properties:
%%with standard values shown
%Cell Time Constant
%tau=0.01s as in Varshney
tau=constants.tau;
%linear membrane term
memG=constants.memG;
memV=constants.memV;
%Electrical/Chemical Synapse Ratio (assumed equal for diff. connections)
%Ratio is 1 in Varshney
gelec=constants.gelec;
gchem=constants.gchem;
%Steepness of chemical conductance Beta (assumed equal)
%beta=~0.125 for mV;
beta=constants.beta;
%Center of chemical conductance vmean
%assumed half open at eq?
vmean=constants.vmean;

%For Synaptic Dynamics:
%corresponding to rise and decay times of synaptic activity
ar = constants.ar;
ad = constants.ad;

%
%
N=length(v);N=N/2;
s=v((N+1):(2*N));
v=v(1:N);

%Synaptic Dynamics
sig=1./(1+exp(-beta*(v-vmean)));
ds = ar*(sig.*(1-s))-ad*s;

%First term:
Iohm=memG*(v-memV);


%Second term: Electrical Synapses
Ielec=gelec*Lgap*v;


%Third term: Chemical Synapses
Ichem=gchem*(v.*(Achem'*s)-Achem'*(s.*Echem));


%Function for dV/dt
dv=-1/tau*(Iohm+Ichem+Ielec-In);
dv=[dv;ds];

%% Ola added this
debug = 0;
out = {};
if debug
    % input current to PLM in every step
    % find connectivity for this neuron
    n = [277, 279];
    Agn = Gel(:, n); % receiving from gap junctions
    AcC = Achem(:,n); % receiving from synapses
    
    input_el = zeros(N, 2);
    input_syn = zeros(N, 2);
    
    for i = 1 : N
        input_syn(i, :) = gchem * AcC(i, :)' .* s(i) .* (v(n) - Echem(i));
    
        input_el(i, :) = gelec .* Agn(i, :)' .* (v(n) - v(i));
    end
    out = {};
    out.Ielec = input_el;
    out.Ichem = input_syn;
end
end












