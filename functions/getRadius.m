function radius = getRadius(h, w, kx, ky, dk)
xx = [ (1:1:w/2)-1, -(w/2:-1:1)];
yy = [-(1:1:h/2)+1, (h/2:-1:1)];
[xx,yy] = meshgrid(xx,yy);
radius = sqrt((xx-kx).^2+(yy-ky).^2)*dk;
end