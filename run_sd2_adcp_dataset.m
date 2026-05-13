function run_sd2_adcp_dataset()
clc;
% 添加当前脚本所在目录到 MATLAB 路径
currentDir = fileparts(mfilename('fullpath'));
addpath(currentDir);
N = 2;
outDir = fullfile(fileparts(mfilename('fullpath')), 'out_seastates');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% baseCfg（固定）
baseCfg = struct();
baseCfg.h = 30;
baseCfg.t = 0:0.5:1200;
baseCfg.beam_angle = 20;
baseCfg.n_bins = 30;
baseCfg.bin_size = 1.0;
baseCfg.adcp_xy = [20,20];
baseCfg.adcp_z0 = -28.5;
baseCfg.freq = 0.05:0.01:0.5;  % 46
baseCfg.dir  = 0:5:355;        % 72

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

% 海况采样 opts
opts = struct();
opts.steepHL = [0.03 0.10];
opts.TpRange = [6 12];
opts.HsClip  = [0.5 4.0];
opts.gammaRange  = [1.5 5.0];
opts.spreadRange = [10 35];
opts.maxTries = 5000;

% ===== 目标粗网格：23×36 =====
freq2 = 0.05:0.02:0.49;   % 23 points
dir2  = 0:10:350;         % 36 points

% ===== 能量检查设置 =====
enable_energy_check = true;
energy_check_first_n = 3;     % 只检查前 n 个样本，避免并行输出太多
energy_relerr_tol = 1e-2;     % 相对误差容差：1%

% 并行池
p = gcp('nocreate');
if isempty(p)
    parpool('local');   % 进程池支持并行保存
end
% 在所有工作进程上添加当前目录
currentDir = fileparts(mfilename('fullpath'));
pctRunOnAll(sprintf('addpath(''%s'')', currentDir));

parfor i = 1:N
    cfg = sea_state_sampler_for_training(baseCfg, opts);
    cfg.platform_mode = 'auv_const_vel';

    % 合理速度范围：0.5~2.0 m/s（你可以按需要调大/调小）
    cfg.auv_speed = 0.5 + (2.0-0.5)*rand();

    % 随机水平航向（如果你想固定某个方向就给具体角度）
    cfg.auv_dir_deg = 360*rand();
    out = synthesize_adcp_sd2_dataset(cfg);

    % ===== 方案 B：能量守恒重采样 =====
    [S2, f2, d2] = downsample_directional_spectrum_energy_conserving( ...
        out.label.S, out.label.freq, out.label.dir, freq2, dir2);

    % ===== 重采样前后能量一致性验证（m0）=====
    % m0 = sum(S) * df * dtheta  （注意 dtheta 用“度”，与你的 dir 网格一致）
    df  = mean(diff(out.label.freq));
    dd  = mean(diff(out.label.dir));
    m0_fine = sum(out.label.S, 'all') * df * dd;

    df2 = mean(diff(f2));
    dd2 = mean(diff(d2));
    m0_coarse = sum(S2, 'all') * df2 * dd2;

    rel_err = (m0_coarse - m0_fine) / m0_fine;

    if enable_energy_check && (i <= energy_check_first_n)
        fprintf('[%d] m0 fine=%.6g, coarse=%.6g, rel_err=%.3g\n', i, m0_fine, m0_coarse, rel_err);
        assert(abs(rel_err) < energy_relerr_tol, ...
            '[%d] Resample energy mismatch too large: rel_err=%.3g', i, rel_err);
    end

    dataset = struct();
    dataset.adcp_radial = out.adcp.radial4;
    dataset.Hs_ts       = out.adcp.Hs_ts;
    dataset.eta_ts      = out.adcp.eta_ts;
    dataset.time        = out.adcp.t;

    dataset.label_S     = S2;  % [23×36]
    dataset.label_freq  = f2;
    dataset.label_dir   = d2;

    dataset.meta        = out.meta;
    dataset.cfg_used    = cfg;

    % 把能量检查结果也存下来（方便后处理统计）
    dataset.meta.resample = struct();
    dataset.meta.resample.m0_fine   = m0_fine;
    dataset.meta.resample.m0_coarse = m0_coarse;
    dataset.meta.resample.rel_err   = rel_err;
    dataset.meta.resample.freq2     = f2;
    dataset.meta.resample.dir2      = d2;
    
    fname = sprintf('sd2_%03d_Hs%.2f_Tp%.2f_s%.3f_th%.0f_g%.2f.mat', ...
        i, cfg.Hs, cfg.Tp, cfg.meta.steep_H_over_L, cfg.theta0, cfg.gamma);
    outFile = fullfile(outDir, fname);

    % parfor 中保存：-fromstruct
    S = struct();
    S.dataset = dataset;
    save(outFile, '-fromstruct', S, '-v7');
end

fprintf('Done. Output dir: %s\n', outDir);
end