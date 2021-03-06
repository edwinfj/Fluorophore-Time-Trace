%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To fit a 2-D gaussian 
% m = Image 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cx,cy,sx,sy,PeakOD] = Gaussian2D(m,tol);

[sizey sizex] = size(m);
[x,y] = MeshGrid(1:sizex,1:sizey);
fit = abs(PeakOD)*(exp(-0.5*(x-cx).^2./(sx^2)-0.5*(y-cy).^2./(sy^2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%