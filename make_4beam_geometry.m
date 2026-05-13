function [beam_vec, bin_pos] = make_4beam_geometry(cfg)
% 四条斜束的波束方向与采样点（不包含第5束）
nB = 4;
beam_vec = zeros(nB,3);
ba = cfg.beam_angle;%波束角度20

for b=1:4
    phi = (b-1)*90; % 0,90,180,270，波束在水平面上的方位角
    beam_vec(b,:) = [sind(ba)*cosd(phi), sind(ba)*sind(phi), cosd(ba)];
end

bin_pos = zeros(nB, cfg.n_bins, 3);
for b=1:nB
    for iz=1:cfg.n_bins
        r = iz*cfg.bin_size;
        bin_pos(b,iz,1) = cfg.adcp_xy(1) + r*beam_vec(b,1);%x
        bin_pos(b,iz,2) = cfg.adcp_xy(2) + r*beam_vec(b,2);%y
        bin_pos(b,iz,3) = cfg.adcp_z0   + r*beam_vec(b,3);%z
    end
end
end