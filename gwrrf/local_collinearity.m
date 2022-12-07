function [ci,vdf] = local_collinearity(x,w)

w = sqrt(w);
sw = sum(w);
w = bsxfun(@rdivide,w,sw);
xw = bsxfun(@times,x,w);
sxw = sqrt(sum(xw.^2));
sxw = bsxfun(@rdivide,xw,sxw);
[u,s,v] = svd(sxw,'econ');
ci = s(1,1)/s(end,end);
phi = v*diag(1./diag(s));
phi = (phi.^2)';
pi_ij = bsxfun(@rdivide,phi,sum(phi));
vdf = pi_ij(end,:);
