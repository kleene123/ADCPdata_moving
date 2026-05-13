function [Tu,Tv,Tw] = sd_kernel_sum(kp, thp, z, h, wsp, ai, aj, comp, pair, n, opts)
% 二阶和频近似核：深度形状 × 方向耦合 × 频率尺度
% 物理动机：
% - 水平速度 ~ ω C(k,z)，竖向 ~ ω S(k,z)
% - 和频方向取 k⃗⁺ 的方向；耦合权重与 cosΔθ 相关
% - 以 ω⁺ 为频率尺度，幅度与 ai·aj 成正比
Ch = cosh(kp*(z+h))/sinh(kp*h);
Sh = sinh(kp*(z+h))/sinh(kp*h);

% 耦合权重：与成分夹角相关（近似）
theta_i = pair.thetai(n); theta_j = pair.thetaj(n);
w_dir = abs(cosd(theta_i - theta_j)); % 0~1

% 水平分解到 x/y（沿 thp）
Tu = wsp * Ch * cosd(thp) * w_dir;
Tv = wsp * Ch * sind(thp) * w_dir;

% 竖向
Tw = wsp * Sh * w_dir;
end