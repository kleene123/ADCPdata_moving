function eta = surface_elevation_at_point(comp, x, y, t)
% 线性表面升沉 eta(t) at z=0
Nc = numel(comp.f);
eta = 0;
for i=1:Nc
    phase = comp.kx(i)*x + comp.ky(i)*y - comp.omega(i)*t + comp.phi(i);
    eta = eta + comp.a(i)*cos(phase);
end
end