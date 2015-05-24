% ==============================
% The objective function will be vaguely like this:
function resids = objfun(params,px,py,img);
covchol = [params(5),0;params(6),params(7)];
temp = [px(:) - params(2),py(:)-params(3)]*covchol;
pred = params(1) + params(4)*exp(-sum(temp.*temp,2)/2);
resids = pred - img(:);