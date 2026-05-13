function k = solve_dispersion(omega, h, g)
% 有限水深色散，牛顿迭代
k = max(omega^2/g, 1e-6);
for it=1:40
    th = tanh(k*h);
    f  = g*k*th - omega^2;
    fp = g*th + g*k*h*(sech(k*h))^2;
    k_new = k - f/fp;
    if abs(k_new-k) < 1e-12, k = k_new; break; end
    k = k_new;
end
end

function s = sech(x), s = 1./cosh(x); end