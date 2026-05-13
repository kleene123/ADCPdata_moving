function S = jonswap_spectrum(f, Hs, Tp, gamma, g)%S 是一个数组，表示每个频率对应的能量密度。
% 简式JONSWAP，归一化到给定Hs
fp = 1/Tp; %峰值频率 (Hz) = 1/峰值周期
alpha = 0.0081;% JONSWAP谱的标准化常数
S = zeros(size(f));%46个0的数组
for i=1:numel(f)
    fi = f(i);
    sigma = (fi<=fp)*0.07 + (fi>fp)*0.09;%频率宽度参数
    r = exp(-(fi - fp)^2/(2*sigma^2*fp^2));%谱峰增强因子，表示谱峰的尖锐程度
    S(i) = alpha*g^2/(2*pi)^4/fi^5 * exp(-1.25*(fp/fi)^4) * gamma^r;
end
m0 = trapz(f, S);%总能量
S = S * (Hs^2/(16*m0)); % Hs^2 = 16 m0，缩放归一化，能量水平可能与期望的有效波高 Hs 不匹配
end