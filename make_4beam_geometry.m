function [beam_vec_body, bin_offset] = make_4beam_geometry(cfg)
% 四条斜束的波束方向与采样点（不包含第5束）
% 输出为机体系（body frame）下的波束方向与bin偏移量
nB = 4;
beam_vec_body = zeros(nB,3);
ba = cfg.beam_angle;%波束角度20

for b=1:4
    phi = (b-1)*90; % 0,90,180,270，波束在水平面上的方位角
    beam_vec_body(b,:) = [sind(ba)*cosd(phi), sind(ba)*sind(phi), cosd(ba)];
end

bin_offset = zeros(nB, cfg.n_bins, 3);
for b=1:nB
    for iz=1:cfg.n_bins
        r = iz*cfg.bin_size;
        bin_offset(b,iz,1) = r*beam_vec_body(b,1);%x
        bin_offset(b,iz,2) = r*beam_vec_body(b,2);%y
        bin_offset(b,iz,3) = r*beam_vec_body(b,3);%z
    end
end
end
