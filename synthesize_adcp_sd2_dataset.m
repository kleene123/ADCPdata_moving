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

%% 1) 方向谱标签
[fgrid, thgrid, Sfd] = make_directional_spectrum(cfg);
out.label.freq = fgrid;
out.label.dir  = thgrid;
out.label.S    = Sfd;

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

%% 4) 4斜束几何
[beam_vec, bin_pos] = make_4beam_geometry(cfg);


% --- AUV motion (MVP: horizontal constant velocity, no attitude) ---
if ~isfield(cfg,'platform_mode') || isempty(cfg.platform_mode)
    cfg.platform_mode = 'fixed';
end

vx = 0; vy = 0; vz = 0;  % vz 固定为 0（你要求先不加垂向）
if strcmpi(cfg.platform_mode, 'auv_const_vel')
    if ~isfield(cfg,'auv_speed') || isempty(cfg.auv_speed)
        % 合理 AUV 速度范围（可改）：0.5~2.0 m/s
        cfg.auv_speed = 0.5 + (2.0-0.5)*rand();
    end
    if ~isfield(cfg,'auv_dir_deg')
        cfg.auv_dir_deg = []; % 允许随机
    end

    [vx, vy, dir_deg] = sample_auv_constant_velocity(cfg.auv_speed, cfg.auv_dir_deg);
    cfg.auv_dir_deg = dir_deg; % 记录下来，便于复现/存档
end

% 统一回填：把“最终真实使用的平台运动”写入 cfg，确保会进入 .mat
cfg.auv_vx = vx;
cfg.auv_vy = vy;
cfg.auv_vz = vz;
if ~isfield(cfg,'auv_speed') || isempty(cfg.auv_speed)
    cfg.auv_speed = hypot(vx, vy); % fixed 模式下为 0
end

%% 5) 第5束：合成 eta(t)，缩放 comp.a 使 Hs 匹配 cfg.Hs，然后计算 Hs_ts
Nt = numel(cfg.t);
x0 = cfg.adcp_xy(1);
y0 = cfg.adcp_xy(2);

eta_ts = zeros(1,Nt);
for it = 1:Nt
    tt = cfg.t(it);
    xt = x0 + vx*tt;
    yt = y0 + vy*tt;
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
    xt = x0 + vx*tt;
    yt = y0 + vy*tt;
    eta_ts(it) = surface_elevation_at_point(comp, xt, yt, tt);
end

fs = 1/median(diff(cfg.t));
Hs_ts = hs_welch_realtime(eta_ts, fs, cfg.Hs_window_sec, cfg.Hs_update_sec, cfg.Hs_fmin, cfg.Hs_fmax);

%% 6) 合成速度并投影（仅4束）
nB = size(bin_pos,1);
nZ = size(bin_pos,2);

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

    % 一阶速度
    U1 = zeros(nB, nZ); V1 = zeros(nB, nZ); W1 = zeros(nB, nZ);
    for b = 1:nB
        for iz = 1:nZ
            x = bin_pos(b,iz,1) + vx*tt;
            y = bin_pos(b,iz,2) + vy*tt;
            z = bin_pos(b,iz,3);   % vz=0，所以不动
            [u1,v1,w1] = linear_velocity_at_point(comp, x, y, z, tt, cfg.h, g);
            U1(b,iz) = u1; V1(b,iz) = v1; W1(b,iz) = w1;
        end
    end

    % 二阶速度（近似核）
    U2 = zeros(nB, nZ); V2 = zeros(nB, nZ); W2 = zeros(nB, nZ);
    if use_sd2
        for b = 1:nB
            for iz = 1:nZ
                x = bin_pos(b,iz,1) + vx*tt;
                y = bin_pos(b,iz,2) + vy*tt;
                z = bin_pos(b,iz,3);
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
        vproj = bv(1)*vx + bv(2)*vy + bv(3)*vz; % 这里 vz=0
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

% 水平常量速度（MVP）
out.meta.platform.vx = cfg.auv_vx;
out.meta.platform.vy = cfg.auv_vy;
out.meta.platform.vz = cfg.auv_vz;          % 现在固定 0
out.meta.platform.speed_mps = hypot(cfg.auv_vx, cfg.auv_vy);

% 方向（如果有）
if isfield(cfg,'auv_dir_deg') && ~isempty(cfg.auv_dir_deg)
    out.meta.platform.dir_deg = cfg.auv_dir_deg;
else
    out.meta.platform.dir_deg = NaN;
end
end