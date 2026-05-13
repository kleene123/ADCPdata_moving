function [fgrid, thgrid, Sfd] = make_directional_spectrum(cfg)
% 二维方向谱 S(f,θ) = JONSWAP(f) × cos^(2s)(θ-θ0)（±90度有效）

g = 9.81;
fgrid = cfg.freq(:)';% 行vectors46个频率点
thgrid = cfg.dir(:)'; % 行vectors72个方向点
Nf = numel(fgrid); Nd = numel(thgrid);% 频率和方向网格大小

S_f = jonswap_spectrum(fgrid, cfg.Hs, cfg.Tp, cfg.gamma, g); % [1×Nf]
Gth = spreading_cos2s(thgrid, cfg.theta0, cfg.spread_s);     % [1×Nd], normalized
Sfd = S_f(:) * Gth(:)'; % [Nf×Nd]
end