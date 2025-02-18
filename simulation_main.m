

close all; 
clear all;

addpath 'Varshney';
addpath 'CEresponsefunctions';

normalised = 0;
[~, labels, class] = datareader('chem', 'weighted');
savedir = pwd;
%%
% the code to investigate the impact of removing the input current from the TRNs
% when the network is already driven by stimulation of the PLM and the
% neuron of interest

T = 10;
s = 0.001;
h = 1e-4;
Tw = 0;

% number of timesteps kept
nts = int32((T - Tw) / s);
% initialize time, voltage vector
t = 0;
ts = linspace(Tw, T, nts);
crop = nts/2;

driver = 'PLM';
idriver = find(strncmpi(driver, labels, length(driver)) == 1);

amp = 2e4; 

%%
% find SVM given joint PLM and PVM inputs

MDB = open('MDB.mat');
DB = MDB.DB; DD = MDB.DD; VB = MDB.VB; VD = MDB.VD;
nDB = length(DB); nDD = length(DD);
nVB = length(VB); nVD = length(VD);
Nplot = [DB, DD, VB, VD];

INPUTS = {'AVM', 'ALM', 'PVM'};
INPUTS = {'PVM'};

AMP_ADD = {-0.01, -0.005, -0.001, 0, 0.001, 0.005, 0.01};
AMP_ADD = {0, 0.001, 0.005, 0.01, 0.05};

RES0 = {};
plot_figs = 1;

for iinput = 1 : length(INPUTS)
    input = INPUTS{iinput};
    i_add = find(strncmpi(input, labels, length(input)) == 1);
    RES0.(input) = {};
    for iamp = 1 : length(AMP_ADD)
        factor = AMP_ADD{iamp};
        amp_add = amp * factor;

        Iext = zeros(279, nts);
        Iext(idriver, :) = amp;
        
        Iext(i_add, :) = Iext(i_add, :) + amp_add;
        matfile = fullfile(savedir, ['baseline_' input '_' strrep(num2str(amp_add), '-', '_') '.mat']);
        if ~exist(matfile, 'file')
            out0 = simulateNetwork(Iext, T, Tw, s, h);
            
            V = out0.V(Nplot, crop:end);
            mu = mean(V, 2);
            stdev = std(V, [], 2) + 1e-3;
            if normalised
                stdev_array = repmat(stdev,  1, size(V, 2));
                V = (V - mu) ./ stdev_array;
            end
            
            [U, S, P] = svd(V);
        
            proj = U(:, 1:3);
            
            save(matfile, "proj", "mu", "stdev")
            if plot_figs
                fig = figure;
                subplot(3, 2, [1 2])
                plot(P(:, 1), P(:, 2))
                
                subplot(3, 2, [3 4])
                Vproj0 = proj' * V;
                plot(Vproj0(1,:), Vproj0(2,:), 'LineWidth', 2);
                
                subplot(3, 1, 3)
                plot(ts(crop:end), P(:, 1)); hold on;
                plot(ts(crop:end), P(:, 2))
                title('modes')
                saveas(fig, fullfile('figures', ['basic_svd_input_into_' input '_amp_' num2str(factor) '.png']))
                close(fig)
            end
        else
            load(matfile)
        end
        RES0.(input).(['amp' strrep(num2str(amp_add), '-', '_')]) = {};
        RES0.(input).(['amp' strrep(num2str(amp_add), '-', '_')]).('proj') = proj;
        RES0.(input).(['amp' strrep(num2str(amp_add), '-', '_')]).('mu') = mu;
        RES0.(input).(['amp' strrep(num2str(amp_add), '-', '_')]).('stdev') = stdev;
        
    end
end
%%
% switch of the additional (PVM) current for the time tstart - tstop

T = 60;
s = 0.001;
h = 1e-4;
Tw = 0;

% number of timesteps kept
nts = int32((T - Tw) / s);
% initialize time, voltage vector
t = 0;
ts = linspace(Tw, T, nts);

tstart = 15;
tstop = 35;
RES = {};

for iinput = 1 : length(INPUTS)
    input = INPUTS{iinput};
    i_add = find(strncmpi(input, labels, length(input)) == 1);
    RES.(input) = {};
    for iamp = 1 : length(AMP_ADD)
        factor = AMP_ADD{iamp};
        amp_add = amp * factor;
        matfile = fullfile(savedir, ['current_' input '_' strrep(num2str(amp_add), '-', '_') '.mat']);
        if ~exist(matfile, 'file')
            Iext = zeros(279, nts);
            Iext(idriver, :) = amp;
            
            Iext(i_add, :) = Iext(i_add, :) + amp_add;
            Iext(i_add, ts > tstart & ts < tstop) = 0;
        
            out = simulateNetwork(Iext, T, Tw, s, h);
            
            save(matfile, "out")
        else
            disp(['loading', matfile])
            load(matfile)
        end
        RES.(input).(['amp' strrep(num2str(amp_add), '-', '_')]) = out;
    end
end



%%

a = 15;
fig = figure('Position', [60, 60, length(INPUTS) * 100 + 60, length(AMP_ADD) * 100]);
[AX, pos] = tight_subplot(length(AMP_ADD), length(INPUTS), [.01 .01], [.1 .06], [.23 .01]);
ax0 = AX(end);
hold(ax0, 'on')
for iinput = 1 : length(INPUTS)
    input = INPUTS{iinput};
    for iamp = 1 : length(AMP_ADD)
        factor = AMP_ADD{iamp};
        amp_add = amp * factor;
        out = RES.(input).(['amp' strrep(num2str(amp_add), '-', '_')]);
        
        proj = RES0.(input).(['amp' strrep(num2str(amp_add), '-', '_')]).('proj');
        disp([input, num2str(factor)])
        
        ax = AX((iamp - 1) * length(INPUTS) + iinput);
        hold(ax, 'on')
        
        V = out.V;
        V = V(Nplot, :); % MNs only
        if normalised
            V = (V - mu) ./ stdev_array;
        end
        
        Vproj = proj' * V;
        %idxs = (ts < tstart & ts > 5);
        %plot(ax, Vproj(1, idxs), Vproj(2, idxs), 'LineWidth', 2);
        
        idxs = ts > (tstop + 5);
        plot(ax, Vproj(1, idxs), Vproj(2, idxs), 'LineWidth', 2);

        idxs = ts > (tstart + 5) & ts < (tstop);
        plot(ax, Vproj(1, idxs), Vproj(2, idxs), 'LineWidth', 2);
        
        
        %X = out.X;
        %X = X(Nplot, :); % MNs only
        fig2 = figure('Position', [100, 100, 500, 200]);
        ax2 = subplot(1, 1, 1);
        plot(ax2, ts, V)
        box(ax2, 'off')
        title(ax2, [input ' ' num2str(factor)])
        xlabel(ax2, 'Time (sec)')
        %saveas(fig2, fullfile('figures', ['driver_' driver '_' num2str(factor) input '_currentOFF_traces.png']))
    
        box(ax, 'off')
        axis(ax, 'square')
        if iinput == 1
            if iamp == 3
                ylabel(ax, {'\it{a_2}'; ['\it{' num2str(factor * 100) '%}']})
            else
                ylabel(ax, ['\it{' num2str(factor * 100) '% }'])
            end
            yticklabels(ax, 'auto');
        else
            yticklabels(ax, []);
        end
   
        if iamp == 1
            title(ax, input)
        end
        if iamp == length(AMP_ADD) 
            if iinput == 2
                xlabel(ax, '\it{a_1}')
            end
            xticklabels(ax, 'auto')
        else
            xticklabels(ax, [])
        end
        fontsize(ax, 12, 'points')
        xlim(ax, [-a, a]); ylim(ax, [-a, a])
        
    end

end
%ax0.Visible = 'off';
lg = legend(ax0, {'baseline', 'no current'});
lg.Orientation = 'horizontal';
lg.Position = [0.3665    0.2104    0.5278    0.0257];
lg.FontSize = 8; %lg.Box = 'off';
uistack(lg, 'top')
linkaxes(AX, 'xy')
%sgtitle('Removing input into')
saveas(fig, fullfile('figures', ['driver_' driver '_input' '.png']))
exportgraphics(fig, fullfile('figures', ['driver_' driver '_input' '.pdf']), 'ContentType', 'vector')

%%  
% adjust synaptic equlibrium while current is off

T = 60;
s = 0.001;
h = 1e-4;
Tw = 0;

% number of timesteps kept
nts = int32((T - Tw) / s);
% initialize time, voltage vector
t = 0;
ts = linspace(Tw, T, nts);

tstart = 15;
tstop = 45;
RESET = {};
INPUTS = {'AVM', 'ALM', 'PVM'};
for iinput = 1 : length(INPUTS)
    input = INPUTS{iinput};
    disp(input)
    i_add = find(strncmpi(input, labels, length(input)) == 1);
    RESET.(input) = {};
    for iamp = 1 : length(AMP_ADD)
        factor = AMP_ADD{iamp};
        amp_add = amp * factor;
        matfile = fullfile(savedir, ['current_reset_' input '_' strrep(num2str(amp_add), '-', '_') '.mat']);
        if ~exist(matfile, 'file')
            Iext = zeros(279, nts);
            Iext(idriver, :) = amp;
            
            Iext(i_add, :) = Iext(i_add, :) + amp_add;
            Iext(i_add, ts > tstart & ts < tstop) = 0;
            resetIdx = -1;
            perturbationIdx = -1;
            synapseAdjustIdx = find(ts >= (tstop + tstart) / 2, 1);
            out = simulateNetwork(Iext, T, Tw, s, h, [], resetIdx, perturbationIdx, synapseAdjustIdx);
            
            save(matfile, "out")
        else
            load(matfile)
        end
        RESET.(input).(['amp' strrep(num2str(amp_add), '-', '_')]) = out;
    end
end

%%
tmid = (tstart + tstop) / 2;
sections = [5, tstart; tstart + 5, tmid; tmid + 5, tstop; tstop + 5, T];
labs = {'I_{MS}', 'no I_{MS}', 'no I_{MS}, \phi(V_{eq})', 'I_{MS}, \phi(V_{eq})'};

sections = [5, tstart; tstart + 5, tmid];
labs = {'I_{MS}', 'no I_{MS}'};

a = 15;
fig = figure('Position', [60, 60, length(INPUTS) * 100 + 60, length(AMP_ADD) * 100]);
[AX, pos] = tight_subplot(length(AMP_ADD), length(INPUTS), [.01 .01], [.1 .06], [.23 .01]);
ax0 = AX(end);
hold(ax0, 'on')
for iinput = 1 : length(INPUTS)
    input = INPUTS{iinput};
    for iamp = 1 : length(AMP_ADD)
        factor = AMP_ADD{iamp};
        amp_add = amp * factor;
        out = RESET.(input).(['amp' strrep(num2str(amp_add), '-', '_')]);
        
        proj = RES0.(input).(['amp' strrep(num2str(amp_add), '-', '_')]).('proj');
        disp([input, num2str(factor)])
        
        ax = AX((iamp - 1) * length(INPUTS) + iinput);
        hold(ax, 'on')
        
        V = out.V;
        V = V(Nplot, :); % MNs only
        if normalised
            V = (V - mu) ./ stdev_array;
        end
        
        Vproj = proj' * V;
        %idxs = (ts < tstart & ts > 5);
        %plot(ax, Vproj(1, idxs), Vproj(2, idxs), 'LineWidth', 2);
        for isection = 1 : size(sections, 1)
            t1 = sections(isection, 1);
            t2 = sections(isection, 2);
        
            idxs = ts > t1 & ts < t2;
            plot(ax, Vproj(1, idxs), Vproj(2, idxs), 'LineWidth', 2);
        end
        fig2 = figure('Position', [100, 100, 500, 200]);
        ax2 = subplot(1, 1, 1);
        plot(ax2, ts, V)
        box(ax2, 'off')
        title(ax2, [input ' ' num2str(factor * 100) '%'])
        xlabel(ax2, 'Time (sec)')
        saveas(fig2, fullfile('figures', ['driver_' driver '_' num2str(factor) input '_currentOFF_reset_synapse2.png']))
    
        box(ax, 'off')
        axis(ax, 'square')
        if iinput == 1
            if iamp == 4
                ylabel(ax, {'\it{a_2}'; ['\it{' num2str(factor * 100) '%}']})
            else
                ylabel(ax, ['\it{' num2str(factor * 100) '% }'])
            end
            yticklabels(ax, 'auto');
        else
            yticklabels(ax, []);
        end
   
        if iamp == 1
            title(ax, input)
        end
        if iamp == length(AMP_ADD) 
            if iinput == 2
                xlabel(ax, '\it{a_1}')
            end
            xticklabels(ax, 'auto')
        else
            xticklabels(ax, [])
        end
        fontsize(ax, 12, 'points')
        xlim(ax, [-a, a]); ylim(ax, [-a, a])
        
    end

end
%ax0.Visible = 'off';
lg = legend(ax0, labs);
%lg.Position(1) = 0.65 * ax0.Position(1);
%lg.Position(2) = 2.1 * ax0.Position(2);
lg.NumColumns = 2;
lg.Position = [0.3290    0.4385    0.5611    0.0579];

lg.FontSize = 8; %lg.Box = 'off';
uistack(lg, 'top')
linkaxes(AX, 'xy')
%sgtitle('Removing input into')
saveas(fig, fullfile('figures', ['driver_' driver '_reset_synapse2' '.png']))
exportgraphics(fig, fullfile('figures', ['driver_' driver 'reset_synapse2' '.pdf']), 'ContentType', 'vector')


%% make a plot for a paper

input = 'PVM';
RESET_PVM = {};
AMP_ADD = {0, 0.001, 0.005, 0.01, 0.05};

for iamp = 1 : length(AMP_ADD)
    factor = AMP_ADD{iamp};
    amp_add = amp * factor;
    matfile = fullfile(savedir, ['current_reset_' input '_' strrep(num2str(amp_add), '-', '_') '.mat']);
    disp(['loading' matfile])
    load(matfile)
    RESET_PVM.(['amp' strrep(num2str(amp_add), '-', '_')]) = out;
end

%%
AMP_ADD = {0.005, 0.01, 0.05};
n = numel(AMP_ADD); 

fig0 = figure('Position', [200, 200, 450 + 100, n * 150 - 40]);
[AX, POS] = tight_subplot(n, 3, [0.01 0.08], [0.15, 0.02], [0.15, 0.02]);
AXS = [];
for i = 1 : 3 * n
    if mod(i, 3) == 2
        AX(i).Position(3) = 2 * AX(i).Position(3) + 0.1;
    elseif mod(i, 3) == 0
        AX(i).Visible = 'off';
    end
    if mod(i, 3) ~= 0
        AXS = [AXS, AX(i)];
    end
end
T = 60;
tstart = 15;
tstop = 45;
tmid = (tstart + tstop) / 2;
sections = [5, tstart; tstart + 5, tmid; tmid + 5, tstop];
labs = {'$I_{PVM}$', 'no $I_{PVM}$', 'no $I_{PVM}$, $V_{eq}$'};

a = 13;

for iamp = 1 : length(AMP_ADD)
    factor = AMP_ADD{iamp};
    amp_add = amp * factor;
    out = RESET_PVM.(['amp' strrep(num2str(amp_add), '-', '_')]);
        
    proj = RES0.(input).(['amp' strrep(num2str(amp_add), '-', '_')]).('proj');
    disp([input, num2str(factor)])
        
    ax = AXS((iamp - 1) * 2 + 1);
    hold(ax, 'on')
    ax2 = AXS((iamp - 1) * 2 + 2);
    hold(ax2, 'on')

    V = out.V;
    V = V(Nplot, :); % MNs only
        
    Vproj = proj' * V;
    for isection = 1 : size(sections, 1)
        t1 = sections(isection, 1);
        t2 = sections(isection, 2);
        
        idxs = ts > t1 & ts < t2;
        p = plot(ax, Vproj(1, idxs), Vproj(2, idxs), 'LineWidth', 2.5);
        set(p, 'Color', [get(p, 'Color') 0.6]); 
        p = plot(ax2, [t1 t2], [2 2], 'LineWidth', 3);
        set(p, 'Color', [get(p, 'Color') 0.6]); 
        if iamp == 1
            text(ax2, (t1 + t2) /2, 2.5, labs{isection}, "FontSize", 7,...
                'HorizontalAlignment','center', 'VerticalAlignment','baseline', 'Interpreter','latex')
        end
    end
    box(ax, 'off')
    axis(ax, 'square')
    ylabel(ax, [num2str(factor * 100) '%', num2str(amp_add * 10) '(mV)'])
    if iamp == 2
        yticks(ax, 'auto')
        yticklabels(ax, 'auto')
        ylabel(ax2, '$\overline{\Delta V}$', 'Interpreter', 'latex');
        yticks(ax2, 'auto')
        yticklabels(ax2, 'auto')
    else
        set(ax2, 'YColor', 'none');
        set(ax, 'YColor', 'none');
    end

    if iamp == 2
        ylabel(ax, {'\it{a_2}'; ['\it{' num2str(factor * 100) '% }', num2str(amp_add * 1e-2) '(V)']},'Color', 'k')
    else
        ylabel(ax, ['\it{' num2str(factor * 100) '% }', num2str(amp_add * 1e-2) '(V)'], 'Color', 'k')
    end
    yticklabels(ax, 'auto');
    xticklabels(ax, [])
    fontsize(ax, 12, 'points')
    xlim(ax, [-a, a]); ylim(ax, [-a, a])
    if iamp ~= length(AMP_ADD)  
        set(ax2, 'XColor', 'none');
        set(ax, 'XColor', 'none');
    end
    
    
    plot(ax2, ts, mean(V, 1), 'Color', [0 0 0 0.6], LineWidth=2)
    box(ax2, 'off')
    
    
    xticklabels(ax2, [])
    fontsize(ax2, 12, 'points')
    ylim(ax2, [-3, 3])
    xlim(ax2, [3, tstop+3])

end
xticks(ax, 'auto')
xticklabels(ax, 'auto')
xlabel(ax, '\it{a_1}')
xticks(ax2, 'auto')
xticklabels(ax2, 'auto')
xlabel(ax2, 'Time (sec)')
%lg = legend(AXS(end-1), labs);
%lg.Orientation = 'horizontal';
%lg.Box = 'off'
%lg.Position = [0.2648    0.7876    0.6000    0.0518];
saveas(fig0, fullfile('figures', ['driver_' driver '_' input '_final.png']))
exportgraphics(fig0, fullfile('figures', ['driver_' driver '_' input '_final.pdf']), 'ContentType', 'vector')
