[ny,nx] = size(img);
[px,py] = meshgrid(1:nx,1:ny);

% starting values:
% 1 - constant offset
% 2 - mode_x
% 3 - mode_y
% 4 - amplitude
% 5,6,7 - cholesky parameters for cov matrix
params0 = [mean(img(:)),nx/2,ny/2,1,5,5,0];
LB = [-inf,1,1,0,-inf,-inf-inf];
UB = [inf,nx,ny,inf,inf,inf,inf];

options = optimset('lsqnonlin');
options.Display = 'iter';
params = lsqnonlin(@(P) objfun(P,px,py,img),params0,LB,UB,options);

cov = [params(5),0;params(6),params(7)];
cov = cov*cov';


