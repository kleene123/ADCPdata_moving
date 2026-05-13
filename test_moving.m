S = load('D:\MatlabProject\MyADCPproject_moving1.0\out_seastates\sd2_001_Hs3.41_Tp7.93_s0.036_th26_g4.96.mat');
rad = S.dataset.adcp_radial;     % [4 x nBins x Nt]
plat = S.dataset.meta.platform;

% beam geometry（用 cfg_used 生成，确保一致）
cfg = S.dataset.cfg_used;
[beam_vec, ~] = make_4beam_geometry(cfg);

v = [plat.vx, plat.vy, plat.vz];

pred = zeros(4,1);
meas = zeros(4,1);
for b=1:4
    pred(b) = -dot(beam_vec(b,:), v);               % 理论平台偏置
    meas(b) = mean(rad(b,:,:), 'all');              % 测到的总体均值（粗略）
end

table((1:4)', pred, meas, meas-pred, ...
    'VariableNames', {'beam','pred_bias','meas_mean','差'})