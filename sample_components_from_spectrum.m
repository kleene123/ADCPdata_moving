function comp = sample_components_from_spectrum(fgrid, thgrid, Sfd, cfg)
% 等网格抽样；幅值 a_i = sqrt(2 S Δf Δθ)；相位随机
df  = abs(fgrid(2)-fgrid(1));%频率间隔
dth = abs(thgrid(2)-thgrid(1))*pi/180; % 方向间隔（转为弧度）
[Nf,Nd] = size(Sfd); % Nf=46频率点, Nd=72方向点  
Nc = Nf*Nd; % 总成分数 = 46×72 = 3312
%comp：结构体，包含抽样得到的波浪分量信息：(初始化的)
comp.f = zeros(1,Nc);%频率分量(一个数组，单位为Hz)
comp.theta = zeros(1,Nc);%方向分量（单位为度）
comp.a = zeros(1,Nc);%波浪分量的幅值。物理意义每个离散波浪分量的波幅，它反映了该分量在波浪场中的能量大小
comp.phi = 2*pi*rand(1,Nc);%波浪分量的随机相位，范围为 [0, 2π)

k = 0;
for i=1:Nf
    for j=1:Nd
        k = k+1;
        comp.f(k) = fgrid(i);
        comp.theta(k) = thgrid(j);
        comp.a(k) = sqrt(2*Sfd(i,j)*df*dth);
    end
end

% 裁剪小幅值
if isfield(cfg,'amp_threshold') && cfg.amp_threshold>0
    keep = comp.a >= cfg.amp_threshold;
    fields = fieldnames(comp);
    for m=1:numel(fields), comp.(fields{m}) = comp.(fields{m})(keep); end
end
end