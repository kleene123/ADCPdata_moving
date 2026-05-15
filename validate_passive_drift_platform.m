function R = validate_passive_drift_platform(inputObj, opts)
%VALIDATE_PASSIVE_DRIFT_PLATFORM Validate passive-drift platform assumptions.
% Usage:
%   R = validate_passive_drift_platform('out_file.mat');
%   R = validate_passive_drift_platform(outStruct);

if nargin < 2, opts = struct(); end
if ~isfield(opts,'z_tol_m'), opts.z_tol_m = 1e-10; end
if ~isfield(opts,'rpy_tol_deg'), opts.rpy_tol_deg = 1e-10; end
if ~isfield(opts,'drift_speed_max_mps'), opts.drift_speed_max_mps = 0.5; end
if ~isfield(opts,'lowfreq_ratio_max'), opts.lowfreq_ratio_max = 0.35; end
if ~isfield(opts,'stokes_dir_tol_deg'), opts.stokes_dir_tol_deg = 5; end

if ischar(inputObj) || isstring(inputObj)
    S = load(inputObj);
    if isfield(S,'dataset')
        meta = S.dataset.meta;
        cfg_used = S.dataset.cfg_used;
    elseif isfield(S,'out')
        meta = S.out.meta;
        cfg_used = S.out.meta.cfg;
    else
        error('Input .mat must contain dataset or out struct.');
    end
elseif isstruct(inputObj) && isfield(inputObj,'meta')
    meta = inputObj.meta;
    cfg_used = inputObj.meta.cfg;
else
    error('Unsupported input type for inputObj.');
end

plat = meta.platform;
t = plat.t(:);
dt = median(diff(t));

roll = plat.roll_deg(:);
pitch = plat.pitch_deg(:);
yaw = plat.yaw_deg(:);
z = plat.z(:);
vx = plat.vx(:);
vy = plat.vy(:);
speed = hypot(vx, vy);

z_ref = cfg_used.adcp_z0;
z_const_err = max(abs(z - z_ref));
rpy_max_abs = max(abs([roll; pitch; yaw]));

dvx = [0; diff(vx)] / dt;
dvy = [0; diff(vy)] / dt;
accel_to_vel_ratio = rms([dvx; dvy]) / max(rms([vx; vy]), eps);

stokes_ok = true;
stokes_dir_err_deg = NaN;
if isfield(plat,'params') && isfield(plat.params,'include_stokes') && plat.params.include_stokes ...
        && isfield(plat.params,'stokes_vx_mps') && isfield(plat.params,'stokes_vy_mps') ...
        && isfield(plat.params,'wave_dom_dir_deg')
    stokes_dir = mod(atan2d(plat.params.stokes_vy_mps, plat.params.stokes_vx_mps), 360);
    dom_dir = mod(plat.params.wave_dom_dir_deg, 360);
    stokes_dir_err_deg = abs(mod((stokes_dir - dom_dir) + 180, 360) - 180);
    stokes_ok = stokes_dir_err_deg <= opts.stokes_dir_tol_deg;
end

R = struct();
R.mode = meta.platform.mode;
R.z_ref_m = z_ref;
R.z_const_max_err_m = z_const_err;
R.rpy_max_abs_deg = rpy_max_abs;
R.speed_mean_mps = mean(speed);
R.speed_rms_mps = rms(speed);
R.speed_max_mps = max(speed);
R.drift_smoothness_metric = accel_to_vel_ratio;
R.stokes_dir_error_deg = stokes_dir_err_deg;
R.flags = struct( ...
    'mode_passive_drift', strcmpi(meta.platform.mode, 'passive_drift'), ...
    'z_constant', z_const_err <= opts.z_tol_m, ...
    'rpy_zero', rpy_max_abs <= opts.rpy_tol_deg, ...
    'drift_speed_reasonable', max(speed) <= opts.drift_speed_max_mps, ...
    'drift_low_frequency', accel_to_vel_ratio <= opts.lowfreq_ratio_max, ...
    'stokes_alignment_ok', stokes_ok);

disp('=== Passive drift validation ===');
disp(R);
end
