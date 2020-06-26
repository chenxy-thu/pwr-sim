function [out,cntrl] = fitPeak(band0, band1, bn0, bn1, otfPr, kx,ky, weightLimit, search)
resPhase =0;
resMag =0;
cntrl = zeros(30,10);
for iter = 1: 1: 10
    b0 = band0;
    b1 = band1;
    [b0_out,b1_out] = commonRegion( b0, b1, bn0, bn1, otfPr, kx, ky, 0.15, weightLimit, true);
%     commonRegion( band0, band1, bn0, bn1, otf, kx, ky, dist, weightLimit, divideByOtf ) 
    b0 = ifft2(b0_out);
    b1 = ifft2(b1_out);
    
    corr = zeros(10,10);
    scal = 1 ./ norm(b0,2);
%     maxV = 0;
%     minV = -inf;
%     newKx = 0;
%     newKy = 0;
    
    
    tkx = kx;
    tky = ky;
    ts = search;
    for p = 1: 1: 100
        xi = mod(p-1,10)+1;
        yi = fix((p-1)/10)+1;
        xpos = tkx + ((xi-4.5)/4.5)*ts;
        ypos = tky + ((yi-4.5)/4.5)*ts;
        b1s = fourierShift(b1, xpos, -ypos);
        b1s = b1s.*conj(b0);
        corr(xi,yi) = sum(b1s(:)) * scal;
    end
    [vv_max,ii_max] = max(abs(corr(:)));
    [vv_min,~] = min(abs(corr(:)));
    
    xi = fix((ii_max-1)/10)+1;
    yi = mod(ii_max-1,10)+1;
    newKx = tkx + ((xi-4.5)/4.5)*ts;
    newKy = tky + ((yi-4.5)/4.5)*ts;
    resPhase = phase(corr(xi,yi));
    resMag   = abs(corr(xi,yi));
    cntrl(((iter-1)*10+1):(iter*10),:) = (abs(corr)-vv_min)/(vv_max-vv_min);
    kx = newKx;
    ky = newKy;
    search = search/5; 
end
out = [ kx, ky, resPhase, resMag];
end