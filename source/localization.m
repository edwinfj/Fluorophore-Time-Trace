%%
%theoretical localization accuracy -Thompson et al
%Nyquist resolution
wvlnth=519;%set emission wavelength in nm -parameter
psf_scale=1.22;%set scale of psf -parameter
NA=1.49;%set NA -parameter
bkgn=5.32;%set background (std photons/pixel) -parameter
q=160;%set the pixel size in nm -parameter
N=1069.3+bkgn;%set the total photons detected
dc=0.0011;%duty cycle

fwhm=psf_scale*0.55*wvlnth/NA;%FWHM of psf
psf_w0 = fwhm/1.1774; % 1/e2 radius of PSF in nm
psf_std=psf_w0/2;% std of psf
%theoretical localization precision in nm
lp2=((psf_std^2)+(q^2)/12)*1./N+8*pi*(psf_std^4)*(bkgn^2)/(q^2)*1./(N.*N);
lp=sqrt(lp2);

%Nyquist resolution
nyqr=2*sqrt(fwhm*fwhm*dc);
%%
%draw localization distribution of one molecule
%10nm scale
[y,x]=size(coord_distr);
y=((1:y)-ceil(y/2))*10;x=((1:x)-ceil(x/2))*10; % -input ?nm/pixel
figure
surf(x,y,coord_distr),axis tight,xlabel('x (nm)'),ylabel('y (nm)');
zlabel(['points (total localizations:' num2str(sum(sum(coord_distr))) ')']);
%title(['Localization distribution of One Molecule (total localizations:' num2str(sum(sum(npoints2))) ')']);
colorbar;
%%
npoints2_sumx=sum(npoints2,1); npoints2_sumy=sum(npoints2,2);
cftool(x,npoints2_sumx);%Gaussian fit for x
cftool(y,npoints2_sumy);%Gaussian fit for y

