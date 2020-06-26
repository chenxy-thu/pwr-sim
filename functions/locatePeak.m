function out = locatePeak(vec, kMin)
h = size(vec,1);
w = size(vec,2);

xx = [(1:1:h/2)-1, (-h/2:1:-1)];
yy = [-(1:1:h/2)+1, (h/2:-1:1)];
[xx,yy] = meshgrid(xx,yy);

rad = sqrt(xx.^2 + yy.^2);
mask = rad>kMin;

vec_mask = vec .* mask;
[~,ii] = max(abs(vec_mask(:)));

iPos = fix((ii-1)/w)+1; % col
jPos = mod(ii-1,w)+1; % row



vMax = abs(vec(iPos,jPos));
vPhase = phase(vec(iPos,jPos));

% dx
if iPos < w/2
    iPos = iPos - 1;
else
    iPos = iPos-h-1;
end

if jPos < h/2
    jPos = -jPos+1;
else
    jPos = h - jPos+1;
end

% dx, dy, v, phase
out = [iPos,jPos, vMax, vPhase];
end