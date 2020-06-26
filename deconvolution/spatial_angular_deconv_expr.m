function sr = spatial_angular_deconv_expr(pm, pm_guess, psf, theta, n_iter, lk, saveDir)
if ~exist([saveDir,'stack\'],'dir')
    mkdir([saveDir,'stack\']);
end
if ~exist([saveDir,'avg\'],'dir')
    mkdir([saveDir,'avg\']);
end

%% parameters
global zero_tol
zero_tol = 1e-10;
%% initialization
% image shape
img_shape = size(pm);
psf_shape = size(psf);
xk_shape = img_shape+[psf_shape(1),psf_shape(2),1]-1;
slice{1} = (psf_shape(1)+1)/2:((psf_shape(1)+1)/2+img_shape(1)-1);
slice{2} = (psf_shape(2)+1)/2:((psf_shape(2)+1)/2+img_shape(2)-1);
% polarization modulation
f = cos(2*theta)+1;
ft = fft(f);
ft = reshape(ft, 1, 1, length(ft));
% xk0
xk = zeros(xk_shape);
if size(pm_guess,3)==3
    for ii = 1: 1: xk_shape(3)
        xk(slice{1},slice{2},ii) = pm_guess(:,:,ii);
    end
else
    for ii = 1: 1: xk_shape(3)
        xk(slice{1},slice{2},ii) = pm_guess;
    end
end


% tk, yk
tkp1 = 1;
ykp1 = xk;
% df
mu = forward( xk, psf, ft, slice);
func1 = sqrt(lsqError( mu, pm) / length(pm(:)));
df = 0;
%% iteration
for kk = 1 : n_iter
    xkm1 = xk;
    yk = ykp1;
    tk = tkp1;
    %% find ltest and slove xk
    mu = forward( yk, psf, ft, slice);
    grad = gradient(mu, pm, psf, ft, size(yk));
    lsqY = lsqError(mu, pm);
    for jj = 0:1000
        ltest = lk*1.1^jj;
        xtest = step(yk, grad, ltest);
        new_mu = forward(xtest, psf, ft, slice);
        newLsq = lsqError(new_mu, pm);
        quadratic = lsqY + quadraticApprox(xtest, yk, grad, ltest);
        if newLsq < quadratic
            xk = xtest;
            lk = ltest;
            break;            
        end
        disp( ['ltest=', num2str(ltest), '  iter_', num2str(kk, '%.3d'), '  df=', num2str(df)])
    end
    if newLsq >= quadratic
        disp( 'cannot find ltest')
        break;            
    end
    %% iteration and analysis
    tkp1 = 1+sqrt(1+4*tk*tk)/2;
    ykp1 = xk + (tk-1)/tkp1*(xk-xkm1);
    % target function
    func0 = func1;
    func_tmp = mu-pm.*log(mu);
    func1 = mean(func_tmp(:));
    df = func0 - func1;
    disp(['iter_', num2str(kk, '%.3d'), '  df=', num2str(df)])
    if mod(kk,2) == 1
%         display3d(xk, 21, ft);
        % save image
        imgsave = xk(slice{1},slice{2},:);        
        bfsave(uint16(imgsave),[saveDir,'stack\', 'stack_iter_', num2str(kk, '%.4d'), '.tif']);
        imgsave_avg = sum(imgsave,3)/3;
        imwrite(uint16(imgsave_avg), [saveDir,'avg\', 'avg_iter_', num2str(kk, '%.4d'), '.tif']);     
    end   
end
%
sr = xk(slice{1},slice{2},:);
%%
% figure(3)
% for kk = 1 : 6
%     subplot(2, 6, kk)
%     imshow(pm(:,:,kk*3),[])
% end
% for kk = 1 : 6
%     subplot(2, 6, kk+6)
%     imshow(g_tmp(:,:,kk*3),[])
% end
void = 0;

function xk = step( yk, grad, lk)
xk = yk-grad/lk;
xk = max(xk, 0);
% sigmoid

function mu = forward(yk, psf, ft, slice)
global zero_tol
for kk = 1 : size(yk,3)
    mu_tmp(:,:,kk) = conv2(yk(:,:,kk), psf(:,:,kk), 'same');
end
mu_f = fft(mu_tmp, [],3);
mu = ifft(mu_f.*repmat(ft, size(mu_f,1), size(mu_f,2)), [], 3);
mu = max(mu(slice{1},slice{2},:), zero_tol);


function grad = gradient(mu, img, psf, ft, g_size)
tmp = 2*(mu - img);
grad = backward(tmp, psf, ft, g_size);

function grad = backward(h, psf, ft, g_size)
grad_g = zeros(g_size);
for kk = 1 : size(h,3)
    grad_g(:,:,kk) = conv2(h(:,:,kk), psf(:,:,kk), 'full');
end
grad_f = fft(grad_g, [],3);
grad = ifft(grad_f.*repmat(ft, size(grad_f,1), size(grad_f,2)), [], 3);

function lsq = lsqError( mu, img)
tmp = (mu-img).^2;
lsq = sum(tmp(:));

function quadratic = quadraticApprox( xk, yk, grad, lk)
delta = xk - yk;
tmp1 = delta.*grad+lk/2*delta.*delta;
quadratic = sum(tmp1(:));
