function [u2,v2,w2] = sd2_velocity_at_point(comp, pair, x, y, z, t, h, g, opts)
% 二阶（Sharma–Dean近似束缚波）在点(x,y,z,t)的速度
% 说明：Q± 使用近似核，包含深度形状与方向耦合。若需精确，请替换为文献公式或WAFO核。

u2=0; v2=0; w2=0;

% 可调权重（近似核缩放），根据海况调优
if ~isfield(opts,'alpha_plus'),  opts.alpha_plus  = 0.25; end   % 和频幅度缩放
if ~isfield(opts,'alpha_minus'), opts.alpha_minus = 0.35; end   % 差频幅度缩放
if ~isfield(opts,'min_km'),      opts.min_km      = 1e-6; end   % 差频稳定阈

for n=1:numel(pair.i)
    ai = pair.ai(n); aj = pair.aj(n);

    % -------- 和频 (+) --------
    kp = pair.kp(n); thp = pair.thp(n); wsp = pair.omega_p(n);
    if kp > 0
        phase_p = kp*(cosd(thp)*x + sind(thp)*y) - wsp*t + pair.phi_sum(n);
        [Tu_p,Tv_p,Tw_p] = sd_kernel_sum(kp, thp, z, h, wsp, ai, aj, comp, pair, n, opts);
        amp_p = opts.alpha_plus * ai * aj;
        u2 = u2 + amp_p * Tu_p * cos(phase_p);
        v2 = v2 + amp_p * Tv_p * cos(phase_p);
        w2 = w2 + amp_p * Tw_p * sin(phase_p);
    end

    % -------- 差频 (−) --------
    km = pair.km(n); thm = pair.thm(n); wdm = pair.omega_m(n);
    if km > opts.min_km && wdm > 0
        phase_m = km*(cosd(thm)*x + sind(thm)*y) - wdm*t + pair.phi_dif(n);
        [Tu_m,Tv_m,Tw_m] = sd_kernel_diff(km, thm, z, h, wdm, ai, aj, comp, pair, n, opts);
        amp_m = opts.alpha_minus * ai * aj;
        u2 = u2 + amp_m * Tu_m * cos(phase_m);
        v2 = v2 + amp_m * Tv_m * cos(phase_m);
        w2 = w2 + amp_m * Tw_m * sin(phase_m);
    else
        % 近零差频：可在此加入 set-down 的近定常项（略）
    end
end
end