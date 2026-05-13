function cfg = sea_state_sampler_for_training(baseCfg, opts)
% 面向训练的随机海况采样器（H=Hs，波陡用 H/L）
% 约束：
%   steep = Hs/Lp in [0.03, 0.10]
%   Lp 由有限水深色散关系 ω^2=gk*tanh(kh) 求得
% 同时避免谱峰落在 cfg.freq 频带边缘，减少谱截断

if nargin < 2, opts = struct(); end
if ~isfield(opts,'maxTries'), opts.maxTries = 5000; end

if ~isfield(opts,'steepHL'), opts.steepHL = [0.03 0.08]; end
if ~isfield(opts,'TpRange'),  opts.TpRange = [6 12]; end

% 为了“现实 + 训练稳定”，给个 Hs 额外裁剪（可按需改/去掉）
if ~isfield(opts,'HsClip'),  opts.HsClip = [0.5 3.5]; end

% 方向谱形参数随机范围（可按需改）
if ~isfield(opts,'gammaRange'),  opts.gammaRange  = [1.5 5.0]; end
if ~isfield(opts,'spreadRange'), opts.spreadRange = [10 35]; end

% 谱峰离频带边界留白（避免峰值落边上被截断）
if ~isfield(opts,'fpMarginLow'),  opts.fpMarginLow = 0.02; end
if ~isfield(opts,'fpMarginHigh'), opts.fpMarginHigh = 0.05; end

g = 9.81;
h = baseCfg.h;

f_band_min = min(baseCfg.freq);
f_band_max = max(baseCfg.freq);

smin = opts.steepHL(1);
smax = opts.steepHL(2);

for ntry = 1:opts.maxTries
    % 1) 随机 Tp（峰值周期）
    Tp = opts.TpRange(1) + (opts.TpRange(2)-opts.TpRange(1))*rand;
    fp = 1/Tp;

    % 频带约束：谱峰不要太靠近边界
    if fp < (f_band_min + opts.fpMarginLow) || fp > (f_band_max - opts.fpMarginHigh)
        continue;
    end

    % 2) 色散 -> kp, Lp
    omega_p = 2*pi*fp;
    kp = solve_dispersion(omega_p, h, g);
    Lp = 2*pi/kp;

    % 3) 随机波陡 steep=Hs/Lp，并得到 Hs
    steep = smin + (smax-smin)*rand;
    Hs = steep * Lp;

    % 4) 额外裁剪 Hs（训练更稳定、更“现实”）
    if Hs < opts.HsClip(1) || Hs > opts.HsClip(2)
        continue;
    end

    % 5) 输出 cfg（覆盖 baseCfg 中的海况参数）
    cfg = baseCfg;
    cfg.Tp = Tp;
    cfg.Hs = Hs;

    % 方向/扩散/JONSWAP峰度随机（利于训练泛化）
    cfg.theta0 = 360*rand;
    cfg.spread_s = opts.spreadRange(1) + (opts.spreadRange(2)-opts.spreadRange(1))*rand;
    cfg.gamma = opts.gammaRange(1) + (opts.gammaRange(2)-opts.gammaRange(1))*rand;

    % 记录诊断信息（可选）
    cfg.meta.fp = fp;
    cfg.meta.kp = kp;
    cfg.meta.Lp = Lp;
    cfg.meta.steep_H_over_L = Hs/Lp;
    cfg.meta.kph = kp*h;

    return;
end

error('Failed to sample a valid sea state within %d tries. Relax TpRange/HsClip/steepHL.', opts.maxTries);
end