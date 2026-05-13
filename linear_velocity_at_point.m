function [u1,v1,w1] = linear_velocity_at_point(comp, x, y, z, t, h, g)
% 一阶（Airy）速度，有限水深
Nc = numel(comp.f);
u1=0; v1=0; w1=0;
for i=1:Nc
    omega=comp.omega(i); ki=comp.k(i);
    cx = cosd(comp.theta(i)); sx = sind(comp.theta(i));
    phase = comp.kx(i)*x + comp.ky(i)*y - omega*t + comp.phi(i);
    Ch = cosh(ki*(z+h))/sinh(ki*h);
    Sh = sinh(ki*(z+h))/sinh(ki*h);
    aomega = comp.a(i)*omega;
    u1 = u1 + aomega*Ch*cos(phase)*cx;
    v1 = v1 + aomega*Ch*cos(phase)*sx;
    w1 = w1 + aomega*Sh*sin(phase);
end
end