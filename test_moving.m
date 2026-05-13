S = load('D:\MatlabProject\MyADCPproject_moving1.0\out_seastates\sd2_001_Hs3.41_Tp7.93_s0.036_th26_g4.96.mat');
rad = S.dataset.adcp_radial;     % [4 x nBins x Nt]
plat = S.dataset.meta.platform;

% beam geometry（用 cfg_used 生成，确保一致）
cfg = S.dataset.cfg_used;
[beam_vec_body, ~] = make_4beam_geometry(cfg);

if isfield(plat, 't')
    t = plat.t;
else
    t = 1:size(rad,3);
end
Nt = numel(t);

roll = get_field_or_default(plat, 'roll_deg', zeros(1,Nt));
pitch = get_field_or_default(plat, 'pitch_deg', zeros(1,Nt));
yaw = get_field_or_default(plat, 'yaw_deg', zeros(1,Nt));
vx = expand_to_length(get_field_or_default(plat, 'vx', 0), Nt);
vy = expand_to_length(get_field_or_default(plat, 'vy', 0), Nt);
vz = expand_to_length(get_field_or_default(plat, 'vz', 0), Nt);

pred = zeros(4,Nt);
meas = zeros(4,1);
for it = 1:Nt
    R = rpy_to_rotmat(roll(it), pitch(it), yaw(it));
    beam_vec = (R * beam_vec_body')';
    v = [vx(it), vy(it), vz(it)];
    for b = 1:4
        pred(b,it) = -dot(beam_vec(b,:), v);               % 理论平台偏置
    end
end

pred_mean = mean(pred, 2);
for b=1:4
    meas(b) = mean(rad(b,:,:), 'all');              % 测到的总体均值（粗略）
end

function out = expand_to_length(val, Nt)
if isscalar(val)
    out = repmat(val, 1, Nt);
else
    out = val;
end
end

table((1:4)', pred_mean, meas, meas-pred_mean, ...
    'VariableNames', {'beam','pred_bias','meas_mean','差'})

function v = get_field_or_default(s, field, defaultVal)
if isfield(s, field) && ~isempty(s.(field))
    v = s.(field);
else
    v = defaultVal;
end
end
