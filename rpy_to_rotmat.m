function R = rpy_to_rotmat(roll_deg, pitch_deg, yaw_deg)
%RPY_TO_ROTMAT Body (roll/pitch/yaw) to ENU rotation matrix
% roll  : rotation about body x (deg)
% pitch : rotation about body y (deg)
% yaw   : rotation about body z (deg)

cr = cosd(roll_deg);  sr = sind(roll_deg);
cp = cosd(pitch_deg); sp = sind(pitch_deg);
cy = cosd(yaw_deg);   sy = sind(yaw_deg);

Rz = [cy, -sy, 0;
      sy,  cy, 0;
       0,   0, 1];

Ry = [ cp, 0, sp;
        0, 1,  0;
      -sp, 0, cp];

Rx = [1,  0,   0;
      0, cr, -sr;
      0, sr,  cr];

R = Rz * Ry * Rx;
end
