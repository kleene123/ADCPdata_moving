function [Tu,Tv,Tw] = sd_kernel_diff(km, thm, z, h, wdm, ai, aj, comp, pair, n, opts)
% 二阶差频近似核：长波调制，深度形状 × 方向耦合
Ch = cosh(km*(z+h))/max(sinh(km*h), 1e-12); % 防止数值发散
Sh = sinh(km*(z+h))/max(sinh(km*h), 1e-12);

theta_i = pair.thetai(n); theta_j = pair.thetaj(n);
w_dir = abs(cosd(theta_i - theta_j));

Tu = wdm * Ch * cosd(thm) * w_dir;
Tv = wdm * Ch * sind(thm) * w_dir;
Tw = wdm * Sh * w_dir;
end