function run_sd2_adcp_dataset_grouped_split()
% Grouped generation with group-wise Train/Val/Test split and parfor progress (CLI)

%% ===== Path fix: client + workers (put at very top) =====
currentDir = fileparts(mfilename('fullpath'));

oldDir = 'D:\MatlabProject\MyADCPproject_2order2.0';
if exist(oldDir,'dir')
    rmpath(genpath(oldDir));
end

addpath(genpath(currentDir));
rehash toolboxcache;

p = gcp('nocreate');
if isempty(p)
    parpool('local'); % or parpool('local', 16)
end

pctRunOnAll(sprintf( ...
    "if exist('%s','dir'), rmpath(genpath('%s')); end; addpath(genpath('%s')); rehash toolboxcache;", ...
    oldDir, oldDir, currentDir));

clc;

%% ===== Scale params =====
M = 5;   % number of sea-state groups
R = 2;   % realizations per group
N = M * R;

%% ===== Output dirs =====
outRoot  = fullfile(currentDir, 'out_seastates_grouped');
outTrain = fullfile(outRoot,'train');
outVal   = fullfile(outRoot,'val');
outTest  = fullfile(outRoot,'test');
if ~exist(outTrain,'dir'), mkdir(outTrain); end
if ~exist(outVal,'dir'), mkdir(outVal); end
if ~exist(outTest,'dir'), mkdir(outTest); end

%% ===== baseCfg & opts =====
baseCfg = struct();
baseCfg.h = 80;
baseCfg.t = 0:0.5:1200;
baseCfg.beam_angle = 20;
baseCfg.n_bins = 30;
baseCfg.bin_size = 1.0;
baseCfg.adcp_xy = [20,20];
baseCfg.adcp_z0 = -45;
baseCfg.freq = 0.05:0.01:0.5;
baseCfg.dir  = 0:5:355;
baseCfg.Hs = 2.0; baseCfg.Tp = 8.0; baseCfg.gamma = 3.3;
baseCfg.theta0 = 60; baseCfg.spread_s = 25;
baseCfg.use_sd2 = false;
baseCfg.amp_threshold = 0;
baseCfg.pair_amp_threshold = 0;
baseCfg.sd2 = struct();
baseCfg.Hs_window_sec = 600;
baseCfg.Hs_update_sec = 10;
baseCfg.Hs_fmin = 0.04;
baseCfg.Hs_fmax = 0.50;

opts = struct();
opts.steepHL = [0.03 0.10];
opts.TpRange = [6 12];
opts.HsClip  = [0.5 4.0];
opts.gammaRange  = [1.5 5.0];
opts.spreadRange = [10 35];
opts.maxTries = 5000;

%% ===== Label downsample grid =====
freq2 = 0.05:0.02:0.49;   % 23
dir2  = 0:10:350;         % 36

%% ===== 1) Sample M groups =====
groups = cell(M,1);
for g = 1:M
    cfg = sea_state_sampler_for_training(baseCfg, opts);
    cfg.verbose = false;       % IMPORTANT: reduce parfor spam
    cfg.use_sd2 = false;       % IMPORTANT: avoid build_pairs memory explosion
    groups{g} = cfg;
end

%% ===== 2) Group-wise split 80/10/10 (min 1 group for val/test) =====
rng(12345);
perm = randperm(M);

nVal   = max(1, round(0.1*M));
nTest  = max(1, round(0.1*M));
nTrain = M - nVal - nTest;
assert(nTrain >= 1, 'M is too small (%d). Need at least 3 groups.', M);

trainG = perm(1:nTrain);
valG   = perm(nTrain+1:nTrain+nVal);
testG  = perm(nTrain+nVal+1:end);

splitOf = strings(M,1);
splitOf(trainG) = "train";
splitOf(valG)   = "val";
splitOf(testG)  = "test";

%% ===== 3) CLI progress for parfor (print every 1%) =====
dq = parallel.pool.DataQueue;

progressCount = 0;
startT = tic;

nextPct = 1;  % next percent to print (1..100)

afterEach(dq, @onTick);

    function onTick(~)
        progressCount = progressCount + 1;

        pct = floor(100 * progressCount / N);
        if pct >= nextPct
            elapsed = toc(startT);
            if progressCount > 0
                eta = elapsed * (N / progressCount - 1);
            else
                eta = NaN;
            end
            fprintf('[%3d%%] %d/%d | elapsed %.1fs | ETA %.1fs\n', ...
                pct, progressCount, N, elapsed, eta);
            nextPct = pct + 1;
        end
    end

fprintf('Start generating %d samples (M=%d, R=%d)...\n', N, M, R);

%% ===== 4) Generate N samples =====
parfor k = 1:N
    g = ceil(k / R);
    r = k - (g-1)*R;

    cfg = groups{g};
    cfg.phase_seed = 1000000*g + r;
    cfg.platform_mode = 'passive_drift';
    cfg.platform = struct();
    cfg.platform.current_speed_mean_mps = 0.08 + 0.08*rand();
    cfg.platform.current_speed_std_mps = 0.02 + 0.03*rand();
    cfg.platform.current_dir_deg = [];
    cfg.platform.drift_corr_time_sec = 120 + 240*rand();
    cfg.platform.drift_rms_mps = 0.015 + 0.03*rand();
    cfg.platform.include_stokes = true;
    cfg.platform.stokes_scale = 0.5 + rand();

    out = synthesize_adcp_sd2_dataset(cfg);

    [S2, f2, d2] = downsample_directional_spectrum_energy_conserving( ...
        out.label.S, out.label.freq, out.label.dir, freq2, dir2);

    dataset = struct();
    dataset.adcp_radial = out.adcp.radial4;
    dataset.Hs_ts       = out.adcp.Hs_ts;
    dataset.eta_ts      = out.adcp.eta_ts;
    dataset.time        = out.adcp.t;

    dataset.label_S     = S2;
    dataset.label_freq  = f2;
    dataset.label_dir   = d2;

    dataset.meta        = out.meta;
    dataset.cfg_used    = cfg;

    dataset.meta.group_id = g;
    dataset.meta.realization_id = r;
    dataset.meta.split = splitOf(g);

    if splitOf(g) == "train"
        outDir = outTrain;
    elseif splitOf(g) == "val"
        outDir = outVal;
    else
        outDir = outTest;
    end

    fname = sprintf('g%05d_r%02d_Hs%.2f_Tp%.2f_th%.0f_g%.2f.mat', ...
        g, r, cfg.Hs, cfg.Tp, cfg.theta0, cfg.gamma);

    Ssave = struct();
    Ssave.dataset = dataset;
    save(fullfile(outDir, fname), '-fromstruct', Ssave, '-v7');

    send(dq, 1); % progress tick
end

elapsedAll = toc(startT);
fprintf('All done in %.1fs.\n', elapsedAll);

fprintf('Done: %s\nTrain=%d, Val=%d, Test=%d samples (Total=%d)\n', ...
    outRoot, numel(trainG)*R, numel(valG)*R, numel(testG)*R, N);
end
