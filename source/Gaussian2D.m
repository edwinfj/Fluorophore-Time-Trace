%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% a function to fit a thermal cloud 2-D
function [cx,cy,sx,sy,PeakOD] = Gaussian2D(m,tol);
%% m = image
%% tol = fitting tolerance

options = optimset('Display','off','TolFun',tol,'LargeScale','off');

[sizey sizex] = size(m);
[cx,cy,sx,sy] = centerofmass(m);
pOD = max(max(m));

mx = m(round(cy),:);
x1D = 1:sizex;
ip1D = [cx,sx,pOD];
fp1D = fminunc(@fitgaussian1D,ip1D,options,mx,x1D);

cx = fp1D(1);
sx = fp1D(2);
PeakOD = fp1D(3);

my = m(:,round(cx))';
y1D = 1:sizey;
ip1D = [cy,sy,pOD];
fp1D = fminunc(@fitgaussian1D,ip1D,options,my,y1D);

cy = fp1D(1);
sy = fp1D(2);
PeakOD = fp1D(3);
[X,Y] = meshgrid(1:sizex,1:sizey);

initpar = [cx,cy,sx,sy,PeakOD];
fp = fminunc(@fitgaussian2D,initpar,options,m,X,Y);
cx = fp(1);
cy = fp(2);
sx = fp(3);
sy = fp(4);
PeakOD = fp(5);