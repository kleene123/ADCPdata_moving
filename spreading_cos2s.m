function G = spreading_cos2s(dir_deg, theta0, s)%方向角向量 主波向 方向扩散参数s越大，波浪越集中在主方向
% 余弦幂方向扩散，±90°有效；归一化到 ∫ G dθ = 1（θ为度）
d = wrapTo180(dir_deg - theta0);%计算每个方向角与主波向的偏移角
G = zeros(size(dir_deg));%72个0的数组
mask = abs(d) <= 90;%标记 ±90° 范围内的有效方向 物理意义: 波浪能量主要集中在主波向 ±90° 范围内
G(mask) = cosd(d(mask)).^(2*s);
dth = abs(dir_deg(2)-dir_deg(1));
G = G / (sum(G)*dth*pi/180); % 度→弧度的积分修正
end
%G 的值越大，表示波浪能量在该方向上的分布越集中。
%G 的分布形状由参数 s 控制：
%当 s 较小时，G 的分布较宽，波浪能量分布在更大的方向范围内。
%当 s 较大时，G 的分布较窄，波浪能量更集中在主波向附近。