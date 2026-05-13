function [S2, f2, dir2] = downsample_directional_spectrum_energy_conserving(S, f, dir, f2, dir2)
% 能量守恒思路的方向谱重采样（分箱积分再除以粗 bin 面积）
% 方向采用周期边界 [0,360)。

assert(isvector(f) && isvector(dir) && isvector(f2) && isvector(dir2), 'f/dir/f2/dir2 must be vectors.');
f = f(:)'; dir = dir(:)'; f2 = f2(:)'; dir2 = dir2(:)';

[Nf, Nd] = size(S);
assert(Nf == numel(f) && Nd == numel(dir), 'Size of S must match f and dir.');

df  = mean(diff(f));
dd  = mean(diff(dir));
df2 = mean(diff(f2));
dd2 = mean(diff(dir2));

% Frequency edges
f_edges  = [f - df/2,  f(end) + df/2];
f2_edges = [f2 - df2/2, f2(end) + df2/2];

% Direction edges (periodic)
dir  = mod(dir, 360);
dir2 = mod(dir2, 360);
dir_edges  = mod([dir - dd/2,  dir(end) + dd/2], 360);
dir2_edges = mod([dir2 - dd2/2, dir2(end) + dd2/2], 360);

overlap_1d = @(a1,a2,b1,b2) max(0, min(a2,b2) - max(a1,b1));

% Wf: [Nf2 x Nf]
Wf = zeros(numel(f2), numel(f));
for i = 1:numel(f2)
    a1 = f2_edges(i);
    a2 = f2_edges(i+1);
    for j = 1:numel(f)
        b1 = f_edges(j);
        b2 = f_edges(j+1);
        Wf(i,j) = overlap_1d(a1,a2,b1,b2);  % Hz
    end
end

% Wd: [Nd2 x Nd]
Wd = zeros(numel(dir2), numel(dir));
for i = 1:numel(dir2)
    a1 = dir2_edges(i);
    a2 = dir2_edges(i+1);
    for j = 1:numel(dir)
        b1 = dir_edges(j);
        b2 = dir_edges(j+1);
        Wd(i,j) = circular_overlap_deg(a1,a2,b1,b2); % deg
    end
end

% integrate then normalize by coarse bin area
S2 = (Wf * S * Wd.') / (df2 * dd2);
end

function ov = circular_overlap_deg(a1,a2,b1,b2)
ov = overlap_wrapped(a1,a2,b1,b2);

    function o = overlap_wrapped(x1,x2,y1,y2)
        xs = split_interval(x1,x2);
        ys = split_interval(y1,y2);
        o = 0;
        for ii = 1:size(xs,1)
            for jj = 1:size(ys,1)
                o = o + max(0, min(xs(ii,2), ys(jj,2)) - max(xs(ii,1), ys(jj,1)));
            end
        end
    end

    function segs = split_interval(p1,p2)
        if p2 >= p1
            segs = [p1 p2];
        else
            segs = [p1 360; 0 p2];
        end
    end
end