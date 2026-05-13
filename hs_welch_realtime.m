function Hs_ts = hs_welch_realtime(eta_ts, fs, Tw_sec, update_sec, fmin, fmax)
%HS_WELCH_REALTIME  用 Welch 方法（pwelch）从表面升沉 eta(t) 实时估计 Hs(t)
%
% 【功能】
%   输入表面升沉时间序列 eta_ts，按“实时更新”的方式输出 Hs_ts：
%     - 目标窗长 Tw_sec（例如 600 s = 10 分钟）
%     - 在窗口未满 Tw_sec 之前，采用“短窗口热启动”：用当前可用数据长度作为窗口
%     - 每隔 update_sec 秒更新一次 Hs，并在两次更新之间用零阶保持（ZOH）
%     - 频带 [fmin, fmax] 内积分得到 m0，再由 Hs = 4*sqrt(m0)
%
% 【输入参数】
%   eta_ts     : [1×Nt] 或 [Nt×1]，表面升沉 (m)
%   fs         : 采样率 (Hz)，例如 2
%   Tw_sec     : 目标窗口长度 (秒)，例如 600（10分钟）
%   update_sec : 更新间隔 (秒)，例如 10（每10秒更新一次）
%   fmin,fmax  : 积分频带 (Hz)，例如 0.04 ~ 0.50
%
% 【输出】
%   Hs_ts      : [1×Nt]，实时 Hs 时间序列（在更新间隔内保持常值）
%
% 【说明/注意】
%   1) “短窗口热启动”含义：
%      当数据还不足 Tw_sec 时，用更短的窗口 W_eff 来估计 Hs；
%      数据越多，W_eff 越大，直到达到目标窗长 W_target。
%   2) 为避免极短窗口导致 Hs 乱跳，可以设定最短可信窗口 minWin_sec：
%      当可用数据不足 minWin_sec 时，输出仍为 NaN（安全策略）。
%      如果你坚持“从第一个点开始永远不 NaN”，可把 minWin_sec 设为 0（不推荐）。
%   3) PSD 使用 MATLAB 自带 pwelch（Signal Processing Toolbox），避免手写归一化出错。

eta_ts = double(eta_ts(:));   % 统一成列向量 double
Nt = numel(eta_ts);

% ---- 目标窗长（点数）
W_target = max(16, round(Tw_sec*fs));       % 目标窗口点数
% ---- 更新步长（点数）
U        = max(1,  round(update_sec*fs));   % 每次更新跨越的点数

% ---- 热启动：最短可信窗口（秒/点）
minWin_sec = 60;                             % 可调：30~120s较常用
W_min      = max(16, round(minWin_sec*fs));  % 最短可信窗口点数

Hs_ts = nan(Nt,1); % 先全部置 NaN

for n = 1:U:Nt
    % ============================================================
    % 1) 决定本次估计使用的“有效窗口长度” W_eff（短窗热启动）
    %    - 若 n < W_target：W_eff = n（窗口随着时间增长）
    %    - 若 n >= W_target：W_eff = W_target（固定10分钟窗）
    % ============================================================
    W_eff = min(W_target, n);

    % 可用数据太短：不输出（保持 NaN）
    if W_eff < W_min
        continue;
    end

    % 取最近 W_eff 点作为当前窗口段
    i1 = n - W_eff + 1;
    seg = eta_ts(i1:n);

    % ============================================================
    % 2) 预处理：去均值（稳健）
    %    如果你有明显的低频漂移，可以改成 detrend(seg,1) 去线性趋势
    % ============================================================
    seg = seg - mean(seg);

    % ============================================================
    % 3) Welch PSD 参数（在 W_eff 这个大窗口内部再分段）
    %    nperseg 不宜超过 W_eff；W_eff 很短时也要能工作
    % ============================================================
    nperseg = min(W_eff, 256);     % 每段长度（点）
    if nperseg < 16
        continue;                 % 太短就不算（避免数值不稳定）
    end
    if mod(nperseg,2)==1
        nperseg = nperseg - 1;    % 保证偶数（方便）
    end
    noverlap = round(0.5*nperseg);% 50% 重叠
    win = hann(nperseg);          % Hann窗

    % ============================================================
    % 4) 用 pwelch 估计 PSD：单位 m^2/Hz
    % ============================================================
    [Pxx,f] = pwelch(seg, win, noverlap, [], fs);

    % ============================================================
    % 5) 频带积分得到 m0，再换算 Hs = 4*sqrt(m0)
    % ============================================================
    band = (f >= fmin) & (f <= fmax);
    if ~any(band)
        Hs_now = NaN;
    else
        m0 = trapz(f(band), Pxx(band));      % m0 = ∫ S(f) df
        Hs_now = 4*sqrt(max(m0,0));
    end

    % ============================================================
    % 6) 零阶保持（ZOH）：把本次 Hs 填到下一次更新之前
    % ============================================================
    n2 = min(Nt, n + U - 1);
    Hs_ts(n:n2) = Hs_now;
end

Hs_ts = Hs_ts.';  % 按你数据集习惯，输出为 1×Nt 行向量
end