function pair = build_pairs(comp, cfg)
% 构建所有(i,j)配对；可加筛选
Nc = numel(comp.f);
[idx_i, idx_j] = ndgrid(1:Nc, 1:Nc);
ii = idx_i(:); jj = idx_j(:);

% 能量阈值筛选（可选）
if isfield(cfg,'pair_amp_threshold') && cfg.pair_amp_threshold>0
    keep = comp.a(ii).*comp.a(jj) >= cfg.pair_amp_threshold;
    ii = ii(keep); jj = jj(keep);
end

pair.i = ii; pair.j = jj;

pair.omega_p = comp.omega(ii) + comp.omega(jj);
pair.omega_m = abs(comp.omega(ii) - comp.omega(jj));

kx_i = comp.kx(ii); ky_i = comp.ky(ii);
kx_j = comp.kx(jj); ky_j = comp.ky(jj);
kpx = kx_i + kx_j; kpy = ky_i + ky_j;
kmx = kx_i - kx_j; kmy = ky_i - ky_j;

pair.kp = hypot(kpx,kpy);
pair.km = hypot(kmx,kmy);
pair.thp = atan2d(kpy,kpx);
pair.thm = atan2d(kmy,kmx);

pair.phi_sum = comp.phi(ii) + comp.phi(jj);
pair.phi_dif = comp.phi(ii) - comp.phi(jj);

pair.ai = comp.a(ii); pair.aj = comp.a(jj);
pair.ki = comp.k(ii); pair.kj = comp.k(jj);
pair.thetai = comp.theta(ii); pair.thetaj = comp.theta(jj);
end