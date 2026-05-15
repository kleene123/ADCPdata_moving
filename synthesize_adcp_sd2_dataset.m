function out = synthesize_adcp_sd2_dataset(cfg)
%SYNTHESIZE_ADCP_SD2_DATASET
% 合成：方向谱标签 + （一阶+二阶）速度场在ADCP采样点 + 四波束投影 + 第5束Hs时间序列
%
% 主要流程：
%   1) 生成方向谱标签 S(f,theta)
%   2) 从方向谱离散抽样波成分 comp（含随机相位）
%   3) 求解色散关系得到 k，并分解得到 kx/ky
%   4) 构建4束几何与采样bin坐标
%   5) 合成表面升沉 eta(t)，并“缩放 comp.a”使该 realization 的 Hs 匹配 cfg.Hs
%      然后用 pwelch 实时估计 Hs_ts（短窗口热启动版本）
%   6) 合成(线性 + 可选二阶)速度场，并投影到4束径向速度
%
% 可选 cfg 字段（新增，方便分组+realization 生成）：
%   cfg.phase_seed : (int)  每个 realization 的随机种子（控制相位抽样等随机过程）
%   cfg.verbose    : (bool) 是否打印进度/校准信息（parfor 生成时建议 false）
%   cfg.use_sd2    : (bool) 是否启用二阶 SD2（二阶会导致 build_pairs 占用极大内存）
%
% 输出 out 结构：
%   out.label.freq / dir / S
%   out.adcp.radial4  [4×nBins×Nt] single
%   out.adcp.Hs_ts    [1×Nt]       single
%   out.adcp.eta_ts   [1×Nt]       single
%   out.adcp.t        [1×Nt]
%   out.meta.comp / cfg

%% 0) Optional: seed & verbosity
if isfield(cfg,'phase_seed') && ~isempty(cfg.phase_seed)
    rng(cfg.phase_seed, 'twister');
end
if ~isfield(cfg,'verbose') || isempty(cfg.verbose)
    cfg.verbose = true;
end

g = 9.81;

if ~isfield(cfg,'adcp_z0') || isempty(cfg.adcp_z0)
    cfg.adcp_z0 = -45;
end

%% 1) 方向谱标签
[fgrid, thgrid, Sfd] = make_directional_spectrum(cfg);
out.label.freq = fgrid;
out.label.dir  = thgrid;
out.label.S    = Sfd;

% 从已生成方向谱提取主频/主方向（用于被动漂流中的 Stokes 漂移）
cfg.wave_drift = infer_wave_drift_from_spectrum(cfg, fgrid, thgrid, Sfd, g);

%% 2) 采样线性成分（离散波成分）
comp = sample_components_from_spectrum(fgrid, thgrid, Sfd, cfg);

%% 3) 色散求 k，并分解到 kx/ky
Nc = numel(comp.f);
comp.omega = 2*pi*comp.f;
comp.k = zeros(1,Nc);
for i = 1:Nc
    comp.k(i) = solve_dispersion(comp.omega(i), cfg.h, g);
end
comp.kx = comp.k .* cosd(comp.theta);
comp.ky = comp.k .* sind(comp.theta);

%% 4) 4斜束几何（机体系）
[cfg, nbin_adjust] = enforce_bin_depth_constraints(cfg);
[beam_vec_body, bin_offset] = make_4beam_geometry(cfg);

%% 4b) 平台运动
[platform, cfg] = build_platform_motion(cfg);

%% 5) 第5束：合成 eta(t)，缩放 comp.a 使 Hs 匹配 cfg.Hs，然后计算 Hs_ts
Nt = numel(cfg.t);

eta_ts = zeros(1,Nt);
for it = 1:Nt
    tt = cfg.t(it);
    xt = platform.x(it);
    yt = platform.y(it);
    eta_ts(it) = surface_elevation_at_point(comp, xt, yt, tt);
end

eta0 = eta_ts - mean(eta_ts);
Hs_eta = 4*sqrt(var(eta0,1));
scale  = cfg.Hs / Hs_eta;

if cfg.verbose
    fprintf('[Hs calibration] Hs_eta=%.4f m, target=%.4f m, scale=%.6f\n', ...
        Hs_eta, cfg.Hs, scale);
end

comp.a = comp.a * scale;

eta_ts = zeros(1,Nt);
for it = 1:Nt
    tt = cfg.t(it);
    xt = platform.x(it);
    yt = platform.y(it);
    eta_ts(it) = surface_elevation_at_point(comp, xt, yt, tt);
end

fs = 1/median(diff(cfg.t));
Hs_ts = hs_welch_realtime(eta_ts, fs, cfg.Hs_window_sec, cfg.Hs_update_sec, cfg.Hs_fmin, cfg.Hs_fmax);

%% 6) 合成速度并投影（仅4束）
nB = size(bin_offset,1);
nZ = size(bin_offset,2);

radial4 = zeros(nB, nZ, Nt, 'single');

% ===== 关键修复：只有在 use_sd2=true 时才 build_pairs =====
use_sd2 = isfield(cfg,'use_sd2') && cfg.use_sd2;
pair = [];
if use_sd2
    pair = build_pairs(comp, cfg);  % 二阶才需要；否则会无意义占内存
end

for it = 1:Nt
    if cfg.verbose && mod(it,50) == 1
        fprintf('t step %d/%d (t=%.1fs)\n', it, Nt, cfg.t(it));
    end
    tt = cfg.t(it);
    R = rpy_to_rotmat(platform.roll_deg(it), platform.pitch_deg(it), platform.yaw_deg(it));
    beam_vec = (R * beam_vec_body')';
    pos = [platform.x(it), platform.y(it), platform.z(it)];
    vplat = [platform.vx(it), platform.vy(it), platform.vz(it)];

    % 一阶速度
    U1 = zeros(nB, nZ); V1 = zeros(nB, nZ); W1 = zeros(nB, nZ);
    for b = 1:nB
        for iz = 1:nZ
            offset = R * squeeze(bin_offset(b,iz,:));
            x = pos(1) + offset(1);
            y = pos(2) + offset(2);
            z = pos(3) + offset(3);
            [u1,v1,w1] = linear_velocity_at_point(comp, x, y, z, tt, cfg.h, g);
            U1(b,iz) = u1; V1(b,iz) = v1; W1(b,iz) = w1;
        end
    end

    % 二阶速度（近似核）
    U2 = zeros(nB, nZ); V2 = zeros(nB, nZ); W2 = zeros(nB, nZ);
    if use_sd2
        for b = 1:nB
            for iz = 1:nZ
                offset = R * squeeze(bin_offset(b,iz,:));
                x = pos(1) + offset(1);
                y = pos(2) + offset(2);
                z = pos(3) + offset(3);
                [u2,v2,w2] = sd2_velocity_at_point(comp, pair, x, y, z, tt, cfg.h, g, cfg.sd2);
                U2(b,iz) = u2; V2(b,iz) = v2; W2(b,iz) = w2;
            end
        end
    end

    % 投影为4束径向速度
    U = U1 + U2;
    V = V1 + V2;
    W = W1 + W2;

    for b = 1:nB
        bv = beam_vec(b,:);
        vproj = dot(bv, vplat);
        radial4(b,:,it) = single( ...
            bv(1).*U(b,:) + bv(2).*V(b,:) + bv(3).*W(b,:) - vproj );
    end
end

%% 7) 输出打包
out.adcp.radial4 = radial4;
out.adcp.Hs_ts   = single(Hs_ts);
out.adcp.eta_ts  = single(eta_ts);
out.adcp.t       = cfg.t;

out.meta.comp = comp;
out.meta.cfg  = cfg;
% --- platform meta for downstream use ---
out.meta.platform = struct();
out.meta.platform.mode = cfg.platform_mode;

out.meta.platform.t = platform.t;
out.meta.platform.x = platform.x;
out.meta.platform.y = platform.y;
out.meta.platform.z = platform.z;
out.meta.platform.vx = platform.vx;
out.meta.platform.vy = platform.vy;
out.meta.platform.vz = platform.vz;
out.meta.platform.speed_mps = hypot(platform.vx, platform.vy);
out.meta.platform.roll_deg  = platform.roll_deg;
out.meta.platform.pitch_deg = platform.pitch_deg;
out.meta.platform.yaw_deg   = platform.yaw_deg;
out.meta.platform.params    = platform.params;
out.meta.platform.bin_adjust = nbin_adjust;
end

function [platform, cfg] = build_platform_motion(cfg)
if ~isfield(cfg,'platform_mode') || isempty(cfg.platform_mode)
    cfg.platform_mode = 'passive_drift';
end
if ~isfield(cfg,'platform') || isempty(cfg.platform)
    cfg.platform = struct();
end

t = cfg.t(:)';
Nt = numel(t);
x0 = cfg.adcp_xy(1);
y0 = cfg.adcp_xy(2);
z0 = cfg.adcp_z0;

params = cfg.platform;

switch lower(cfg.platform_mode)
    case 'passive_drift'
        params.current_speed_mean_mps = get_or_default(params, 'current_speed_mean_mps', 0.12);
        params.current_speed_std_mps = get_or_default(params, 'current_speed_std_mps', 0.04);
        params.current_dir_deg = get_or_default(params, 'current_dir_deg', []);
        params.drift_corr_time_sec = get_or_default(params, 'drift_corr_time_sec', 180);
        params.drift_rms_mps = get_or_default(params, 'drift_rms_mps', 0.03);
        params.include_stokes = get_or_default(params, 'include_stokes', true);
        params.stokes_scale = get_or_default(params, 'stokes_scale', 1.0);

        if isempty(params.current_dir_deg)
            params.current_dir_deg = 360*rand();
        end
        current_speed = max(0, params.current_speed_mean_mps + params.current_speed_std_mps*randn());
        current_vx = current_speed * cosd(params.current_dir_deg);
        current_vy = current_speed * sind(params.current_dir_deg);

        if params.include_stokes && isfield(cfg,'wave_drift') && isstruct(cfg.wave_drift)
            stokes_speed = params.stokes_scale * cfg.wave_drift.stokes_speed_mps;
            stokes_dir_deg = cfg.wave_drift.dom_dir_deg;
        else
            stokes_speed = 0;
            stokes_dir_deg = NaN;
        end
        stokes_vx = stokes_speed * cosd(stokes_dir_deg);
        stokes_vy = stokes_speed * sind(stokes_dir_deg);

        dt = median(diff(t));
        tau = max(params.drift_corr_time_sec, dt);
        alpha = exp(-dt / tau);
        sigma = params.drift_rms_mps * sqrt(1 - alpha^2);

        ou_vx = zeros(1,Nt);
        ou_vy = zeros(1,Nt);
        for it = 2:Nt
            ou_vx(it) = alpha * ou_vx(it-1) + sigma * randn();
            ou_vy(it) = alpha * ou_vy(it-1) + sigma * randn();
        end

        vx = current_vx + stokes_vx + ou_vx;
        vy = current_vy + stokes_vy + ou_vy;
        vz_vec = zeros(1,Nt);

        x = zeros(1,Nt); y = zeros(1,Nt);
        x(1) = x0; y(1) = y0;
        for it = 2:Nt
            dti = t(it) - t(it-1);
            x(it) = x(it-1) + 0.5*(vx(it) + vx(it-1))*dti;
            y(it) = y(it-1) + 0.5*(vy(it) + vy(it-1))*dti;
        end
        z = z0 * ones(1,Nt);

        roll = zeros(1,Nt);
        pitch = zeros(1,Nt);
        yaw = zeros(1,Nt);

        params.current_speed_sampled_mps = current_speed;
        params.current_vx_mps = current_vx;
        params.current_vy_mps = current_vy;
        params.stokes_speed_mps = stokes_speed;
        params.stokes_dir_deg = stokes_dir_deg;
        params.stokes_vx_mps = stokes_vx;
        params.stokes_vy_mps = stokes_vy;
        params.wave_peak_freq_hz = cfg.wave_drift.peak_freq_hz;
        params.wave_peak_period_sec = cfg.wave_drift.peak_period_sec;
        params.wave_dom_dir_deg = cfg.wave_drift.dom_dir_deg;

        vx0 = mean(vx);
        vy0 = mean(vy);

    case 'fixed'
        vx0 = 0; vy0 = 0; vz_vec = zeros(1,Nt);
        x = x0 * ones(1,Nt);
        y = y0 * ones(1,Nt);
        z = z0 * ones(1,Nt);
        roll = zeros(1,Nt);
        pitch = zeros(1,Nt);
        yaw = zeros(1,Nt);
        params.speed_mps = 0;
        params.dir_deg = NaN;

    case 'auv_const_vel'
        if ~isfield(cfg,'auv_speed') || isempty(cfg.auv_speed)
            cfg.auv_speed = 0.5 + (2.0-0.5)*rand();
        end
        if ~isfield(cfg,'auv_dir_deg')
            cfg.auv_dir_deg = [];
        end
        [vx0, vy0, dir_deg] = sample_auv_constant_velocity(cfg.auv_speed, cfg.auv_dir_deg);
        cfg.auv_dir_deg = dir_deg;
        params.speed_mps = cfg.auv_speed;
        params.dir_deg = dir_deg;

        x = x0 + vx0 * t;
        y = y0 + vy0 * t;
        z = z0 * ones(1,Nt);
        vz_vec = zeros(1,Nt);
        roll = zeros(1,Nt);
        pitch = zeros(1,Nt);
        yaw = zeros(1,Nt);

    case 'moving_full'
        if ~isfield(params,'speed_mps') || isempty(params.speed_mps)
            if isfield(cfg,'auv_speed') && ~isempty(cfg.auv_speed)
                params.speed_mps = cfg.auv_speed;
            else
                params.speed_mps = 0.5 + (2.0-0.5)*rand();
            end
        end
        if ~isfield(params,'dir_deg') || isempty(params.dir_deg)
            if isfield(cfg,'auv_dir_deg') && ~isempty(cfg.auv_dir_deg)
                params.dir_deg = cfg.auv_dir_deg;
            else
                params.dir_deg = 360*rand();
            end
        end

        [vx0, vy0, dir_deg] = sample_auv_constant_velocity(params.speed_mps, params.dir_deg);
        params.dir_deg = dir_deg;

        params.yaw0_deg = get_or_default(params, 'yaw0_deg', dir_deg);
        params.yaw_rate_deg_s = get_or_default(params, 'yaw_rate_deg_s', 0);
        params.yaw_amp_deg = get_or_default(params, 'yaw_amp_deg', 5);
        params.yaw_period_sec = get_or_default(params, 'yaw_period_sec', 60);
        params.yaw_phase_deg = get_or_default(params, 'yaw_phase_deg', 360*rand());

        params.roll_amp_deg = get_or_default(params, 'roll_amp_deg', 3);
        params.roll_period_sec = get_or_default(params, 'roll_period_sec', 12);
        params.roll_phase_deg = get_or_default(params, 'roll_phase_deg', 360*rand());

        params.pitch_amp_deg = get_or_default(params, 'pitch_amp_deg', 2);
        params.pitch_period_sec = get_or_default(params, 'pitch_period_sec', 10);
        params.pitch_phase_deg = get_or_default(params, 'pitch_phase_deg', 360*rand());

        params.heave_amp_m = get_or_default(params, 'heave_amp_m', 0.3);
        params.heave_period_sec = get_or_default(params, 'heave_period_sec', 8);
        params.heave_phase_deg = get_or_default(params, 'heave_phase_deg', 360*rand());

        x = x0 + vx0 * t;
        y = y0 + vy0 * t;

        heave = sinusoid_series(t, params.heave_amp_m, params.heave_period_sec, params.heave_phase_deg);
        z = z0 + heave;
        vz_vec = sinusoid_derivative(t, params.heave_amp_m, params.heave_period_sec, params.heave_phase_deg);

        roll  = sinusoid_series(t, params.roll_amp_deg, params.roll_period_sec, params.roll_phase_deg);
        pitch = sinusoid_series(t, params.pitch_amp_deg, params.pitch_period_sec, params.pitch_phase_deg);
        yaw = params.yaw0_deg + params.yaw_rate_deg_s * t + ...
            sinusoid_series(t, params.yaw_amp_deg, params.yaw_period_sec, params.yaw_phase_deg);

    otherwise
        error('Unknown platform_mode: %s', cfg.platform_mode);
end

platform = struct();
platform.t = t;
platform.x = x;
platform.y = y;
platform.z = z;
if exist('vx','var') && exist('vy','var')
    platform.vx = vx;
    platform.vy = vy;
else
    platform.vx = vx0 * ones(1,Nt);
    platform.vy = vy0 * ones(1,Nt);
end
platform.vz = vz_vec;
platform.roll_deg = roll;
platform.pitch_deg = pitch;
platform.yaw_deg = yaw;
platform.params = params;

cfg.platform = params;
cfg.auv_vx = vx0;
cfg.auv_vy = vy0;
cfg.auv_vz = mean(vz_vec);
cfg.auv_speed = hypot(vx0, vy0);
if isfield(params,'dir_deg')
    cfg.auv_dir_deg = params.dir_deg;
elseif hypot(vx0,vy0) > 0
    cfg.auv_dir_deg = mod(atan2d(vy0, vx0), 360);
else
    cfg.auv_dir_deg = NaN;
end
end

function val = get_or_default(s, field, defaultVal)
if isfield(s, field) && ~isempty(s.(field))
    val = s.(field);
else
    val = defaultVal;
end
end

function y = sinusoid_series(t, amp, period, phase_deg)
if amp == 0 || period <= 0
    y = zeros(size(t));
else
    phase = deg2rad(phase_deg);
    y = amp * sin(2*pi*t/period + phase);
end
end

function y = sinusoid_derivative(t, amp, period, phase_deg)
if amp == 0 || period <= 0
    y = zeros(size(t));
else
    phase = deg2rad(phase_deg);
    y = amp * (2*pi/period) * cos(2*pi*t/period + phase);
end
end

function wave_drift = infer_wave_drift_from_spectrum(cfg, fgrid, thgrid, Sfd, g)
df = abs(fgrid(2)-fgrid(1));
dd_rad = abs(thgrid(2)-thgrid(1))*pi/180;

Sf = sum(Sfd, 2);                       % 频率方向积分前的频谱切片
[~, ifp] = max(Sf);
peak_freq_hz = fgrid(ifp);
peak_period_sec = 1 / peak_freq_hz;

[~, ith] = max(Sfd(ifp,:));
dom_dir_deg = thgrid(ith);

zref = cfg.adcp_z0;
omega_p = 2*pi*peak_freq_hz;
kp = solve_dispersion(omega_p, cfg.h, g);

a_peak = sqrt(max(2*Sfd(ifp,ith)*df*dd_rad, 0));
stokes_speed = (a_peak^2) * omega_p * kp * exp(2*kp*zref);

wave_drift = struct();
wave_drift.peak_freq_hz = peak_freq_hz;
wave_drift.peak_period_sec = peak_period_sec;
wave_drift.dom_dir_deg = dom_dir_deg;
wave_drift.k_peak = kp;
wave_drift.a_peak_m = a_peak;
wave_drift.stokes_speed_mps = stokes_speed;
end

function [cfg, info] = enforce_bin_depth_constraints(cfg)
surface_margin_m = 0.5;
z0 = cfg.adcp_z0;

if z0 >= 0
    error('cfg.adcp_z0 must be negative (underwater), got %.3f', z0);
end
if cfg.h <= abs(z0)
    error('cfg.h (%.3f m) must exceed |cfg.adcp_z0| (%.3f m).', cfg.h, abs(z0));
end

max_upward_range = (-z0 - surface_margin_m) / cosd(cfg.beam_angle);
max_bins = floor(max_upward_range / cfg.bin_size);
if max_bins < 1
    error('No valid bins remain below surface with current adcp_z0/bin_size/beam_angle.');
end

info = struct();
info.surface_margin_m = surface_margin_m;
info.original_n_bins = cfg.n_bins;
info.adjusted = false;

if cfg.n_bins > max_bins
    cfg.n_bins = max_bins;
    info.adjusted = true;
end
info.applied_n_bins = cfg.n_bins;
info.max_n_bins_allowed = max_bins;
end
