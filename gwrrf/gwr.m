function [pred,betas,se_pred,se_betas,cn,vdf,local_r2,fm] = gwr(x,y,coors,bandwidth,varargin)
%Inference function for GWR model at sampled/unsampled points
%using the bandwidth
%
% USAGE:
%   [pred,betas,ci,se_pred,se_betas] = gwr(x,y,coors,bandwidth,'Name1',Value1,'Name2',Value2,...);
%
% INPUT:
%   x                        -   n*m covariate matrix containing m covariates for n observations 
%   y                        -   n*1 responce vector containing n observations
%   coors                -   n*2 coordinate matrix containing horizontal axis and
%                                 vertical axis for n observations.
%   bandwidth      -   sclar value indicating bandwidth parameter
%   varargin -   Name/Value paris, 
%                       adapt : logical value indicating the adaptive bandwidth
%                       gweight: char value indicating the bandwith 
%                                        kernel type, the supported kernel 
%                                        function is 'guass' and 'bisquare'
%                       disttype: char value indicating distance type such
%                                       as 'euclidean' and 'great circle' type
%                       parfor: logical value indicating whether to use parpool
%                       predx: n_pred*m covariate matrix containing 
%                                   m covariates for n_pred predicted points
%                       predcoor: n_pred*m coordinate matrix containing horizontal axis and
%                                 vertical axis for n_pred predicted points
%
%
% OUTPUT:
%   pred       - n_pred*1 vector containing predicted value
%   betas     - n_pred*m matrix containing regression coefficients for m
%                    covaraites at n_pred points
%   ci            - n_pred*1 vector containing condition number for n_pred
%                    local gwr regression
%   se_pred - n_pred*1 vector containing standard errors for pred
%   se_betas - n_pred*m matrix containing standard errors for betas

warning('off','MATLAB:nearlySingularMatrix')

% error checking
if size(y,1) ~= size(x,1) || size(coors,1) ~= size(x,1)
    error('x , y and coors must have the same number of rows')
end

if ~isscalar(bandwidth)
    error('bandwidth must be a scalar')
end

% check input using PARSEARGS
params.predx                = x;
params.predcoor           = coors;
params.gweight             = 'bisquare';
params.adapt                = true;
params.seiga_v1             = true;
params.disttype             = 'euclidean';
params.parfor             = false;
params = parseargs(params,varargin{:});

[n,m] = size(x);

if size(params.predcoor,1) ~= size(params.predx,1)
    error('predicted x and coors must have the same number of rows')
end

if size(params.predcoor,2) ~= size(coors,2)
    error('predicted coors and coors must have the same number of cols')
end

if size(params.predx,2) ~= size(x,2)
    error('predicted predx and x must have the same number of cols')
end

n_pred = size(params.predx,1);

%%%%%%%%%%%%%% fitting %%%%%%%%%%%%%
pool = gcp('nocreate');
if params.parfor && ~isempty(pool)
parfor i=1:n
    if strcmp(params.disttype,'greatcircle')
        dist = gc_distance(coors,coors(i,:));
    else
        dist = sqrt(sum(bsxfun(@minus,coors,coors(i,:)).^2,2));
    end
    
    if params.adapt
        sort_dist = sort(dist);
        bw = sort_dist(bandwidth);
    else
        bw = bandwidth;
    end
    w = feval(params.gweight,dist,bw);
    xw = bsxfun(@times,x,w);
    lhati = x(i,:)/(x'*xw)*xw';
    v2(i) = sum(lhati.^2);
    v1(i) = lhati(1,i);
    lhati = -lhati;
    lhati(1,i) = lhati(1,i) + 1;
    sigma2(i) = lhati*y;
end
v2 = sum(v2);
v1 = sum(v1);
else
v1 = 0;
v2 = 0;
sigma2 = zeros(n,1);
for i=1:n
    if strcmp(params.disttype,'greatcircle')
        dist = gc_distance(coors,coors(i,:));
    else
        dist = sqrt(sum(bsxfun(@minus,coors,coors(i,:)).^2,2));
    end
    
    if params.adapt
        sort_dist = sort(dist);
        bw = sort_dist(bandwidth);
    else
        bw = bandwidth;
    end
    w = feval(params.gweight,dist,bw);
    xw = bsxfun(@times,x,w);
    lhati = x(i,:)/(x'*xw)*xw';
    v2 = v2 + sum(lhati.^2);
    v1 = v1 + lhati(1,i);
    lhati = -lhati;
    lhati(1,i) = lhati(1,i) + 1;
    sigma2(i) = lhati*y;
end
end
ssr = sum(sigma2.^2);
if params.seiga_v1
    sei = sqrt(ssr/(n-v1));
else
    sei = sqrt(ssr/(n-2*v1+v2));
end

fm.bw = bandwidth;
fm.adapt = params.adapt;
fm.trS = v1;
fm.trSS = v2;
fm.sei = sei;
fm.seiga_v1 = params.seiga_v1;
fm.disttype = params.disttype;
fm.gweight = params.gweight;


%%%%%%%% predict %%%%%%%%%%
if params.parfor && ~isempty(pool)
parfor i=1:n_pred
    if strcmp(params.disttype,'greatcircle')
        dist = gc_distance(coors,params.predcoor(i,:));
    else
        dist = sqrt(sum(bsxfun(@minus,coors,params.predcoor(i,:)).^2,2));
    end
    if ~(params.adapt)
        bw = bandwidth;
    else
        sort_dist = sort(dist);
        bw = sort_dist(bandwidth);
    end
    w = feval(params.gweight,dist,bw);
	
	[cni,vdfi] = local_collinearity(x,w);
    cn(i,:) = cni;
	vdf(i,:) = vdfi;
	
	wy = sum(y.*w)./sum(w);
	tss = sum(w.*((y-wy).^2));
	rss = sum(w.*(sigma2).^2);
	local_r2(i,:) = (tss-rss)/tss;
    
    xw = bsxfun(@times,x,w);
    invnbb = eye(m)/(x'*xw);
    Sp = params.predx(i,:)*invnbb*xw';

    betas(i,:) = (invnbb*x'*(y.*w))';
    pred(i,:) = Sp*y;
    se_pred(i,:) = sqrt(Sp*Sp')*sei;
    Cp = invnbb*xw';
    se_betas(i,:) = sqrt(diag(Cp*Cp')')*sei;

end
else
betas = nan(n_pred,m);
pred = nan(n_pred,1);
cn = nan(n_pred,1);
local_r2 = nan(n_pred,1);
vdf = nan(n_pred,m);
se_pred = nan(n_pred,1);
se_betas = nan(n_pred,m);
for i=1:n_pred

    if strcmp(params.disttype,'greatcircle')
        dist = gc_distance(coors,params.predcoor(i,:));
    else
        dist = sqrt(sum(bsxfun(@minus,coors,params.predcoor(i,:)).^2,2));
    end
    if ~(params.adapt)
        bw = bandwidth;
    else
        sort_dist = sort(dist);
        bw = sort_dist(bandwidth);
    end
    w = feval(params.gweight,dist,bw);
	
    [cni,vdfi] = local_collinearity(x,w);
    cn(i,:) = cni;
	vdf(i,:) = vdfi;
	
	wy = sum(y.*w)./sum(w);
	tss = sum(w.*((y-wy).^2));
	rss = sum(w.*(sigma2).^2);
	local_r2(i,:) = (tss-rss)/tss;

    xw = bsxfun(@times,x,w);
    invnbb = eye(m)/(x'*xw);
    betas(i,:) = (invnbb*x'*(y.*w))';
    Sp = params.predx(i,:)*invnbb*xw';

    pred(i,:) = Sp*y;
    se_pred(i,:) = sqrt(Sp*Sp')*sei;
    Cp = invnbb*xw';
    se_betas(i,:) = sqrt(diag(Cp*Cp')')*sei;

end
end
end
