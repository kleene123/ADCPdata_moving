function [vx, vy, dir_deg] = sample_auv_constant_velocity(speed_mps, dir_deg)
% speed_mps: 标量速度 (m/s)
% dir_deg:   航向角 (deg)，空则随机 [0,360)

    if nargin < 2 || isempty(dir_deg)
        dir_deg = 360*rand();  % 随机方向
    end
    vx = speed_mps * cosd(dir_deg);
    vy = speed_mps * sind(dir_deg);
end