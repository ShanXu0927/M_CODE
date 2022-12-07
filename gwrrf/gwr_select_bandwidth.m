function [bandwidth,score] = gwr_select_bandwidth(x,y,coors,varargin)
%Search the optimal bandwith for the gwr model
%using Golden Section optimization method 
%
% USAGE:
%   [bw,score] = gwr_select_bandwidth(x,y,coors,'Name1',Value1,'Name2',Value2,...);
%
% INPUT:
%   x              -   n*m covariate matrix containing m covariates for n observations 
%   y              -   n*1 responce vector containing n observations
%   coors      -   n*2 coordinate matrix containing horizontal axis and
%                     vertical axis for n observations.
%   varargin -   Name/Value paris, 
%                       adapt : logical value indicating the adaptive bandwidth
%                       gweight: char value indicating the bandwith 
%                                        kernel type, the supported kernel 
%                                        function is 'guass' and 'bisquare'
%                       method: char value indicating the bandwith search
%                                       method, such as 'cv' (cross validation) 
%                                       and 'aic' (Akaike information criterion)
%                       disttype: char value indicating distance type such
%                                       as 'euclidean' and 'great circle' type
%                       parfor: logical value indicating whether to use parpool
%                       min_beta: sclar value containg minmum bandwidth
%                       max_beta: sclar value containg maximum bandwidth
%
% OUTPUT:
%   bandwidth - sclar value containg the optimal bandwidth result
%   score           - sclar value containg the minimal cv or aic result, if
%                          this value is inf, there is serious model design
%                          for establing gwr model


% error checking
if size(y,1) ~= size(x,1) || size(coors,1) ~= size(x,1)
    error('x , y and coors must have the same number of rows')
end

% check input using PARSEARGS
params.adapt           = true;
params.gweight        = 'bisquare';
params.method        = 'cv';
params.cnflag           = true;
params.disttype        = 'euclidean';
params.parfor        = false;
params.rmseflag          = true;
params.min_beta     = [];
params.max_beta    = [];
params = parseargs(params,varargin{:});

[n,m] = size(x);

if ~islogical(params.adapt)
    warning('adapt parameter should be logical number!')
end

pool = gcp('nocreate');
if ~params.adapt

    min_coors_x = min(coors(:,1));
    max_coors_x = max(coors(:,1));
    min_coors_y = min(coors(:,2));
    max_coors_y = max(coors(:,2));
    
    if strcmp(params.disttype,'greatcircle')
        difmin = gc_distance([min_coors_x,min_coors_y],...
            [max_coors_x,max_coors_y]);
    else
        difmin = sqrt((min_coors_x-max_coors_x).^2 ...
            + (max_coors_y-min_coors_y).^2);
    end
    
    if isempty(params.min_beta)
        beta1 = difmin/5000;
    else
        beta1 = params.min_beta;
    end
    if isempty(params.max_beta)
        beta2 = difmin;
    else
        beta2 = params.max_beta;
    end
    if strcmp(params.method, 'cv')
        if params.parfor && ~isempty(pool)
            funs = @(bw) par_cv_f(bw,x,y,coors,...
                params.gweight,params.disttype,params.rmseflag,params.cnflag);
            [bandwidth,score] = golden_section(beta1, beta2, 0.38197, funs, 1.0e-6, 200, false);
        else
            funs = @(bw) cv_f(bw,x,y,coors,...
                params.gweight,params.disttype,params.rmseflag,params.cnflag);
            [bandwidth,score] = golden_section(beta1, beta2, 0.38197, funs, 1.0e-6, 200, false);
        end
    else
        if params.parfor && ~isempty(pool)
            funs = @(bw) par_aic_f(bw,x,y,coors,...
                params.gweight,params.disttype,params.cnflag);
            [bandwidth,score] = golden_section(beta1, beta2, 0.38197, funs, 1.0e-6, 200, false);
        else
            funs = @(bw) aic_f(bw,x,y,coors,...
                params.gweight,params.disttype,params.cnflag);
            [bandwidth,score] = golden_section(beta1, beta2, 0.38197, funs, 1.0e-6, 200, false);
        end
    end
else
    if isempty(params.min_beta)
        beta1 = 20+2*m;
    else
        beta1 = params.min_beta;
    end
    if isempty(params.max_beta)
        beta2 = n;
    else
        beta2 = params.max_beta;
    end
    if strcmp(params.method, 'cv')
        if params.parfor && ~isempty(pool)
            funs = @(bw) par_cv_adapt_f(bw,x,y,coors,...
                params.gweight,params.disttype,params.rmseflag,params.cnflag);
            [bandwidth,score] = golden_section(beta1, beta2, 0.38197, funs, 1.0e-6, 200, true);
        else
            funs = @(bw) cv_adapt_f(bw,x,y,coors,...
                params.gweight,params.disttype,params.rmseflag,params.cnflag);
%             [bandwidth,score] = fminbnd(funs,beta1, beta2);
            [bandwidth,score] = golden_section(beta1, beta2, 0.38197, funs, 1.0e-6, 200, true);
        end
    else
        if params.parfor && ~isempty(pool)
            funs = @(bw) par_aic_adapt_f(bw,x,y,coors,...
                params.gweight,params.disttype,params.cnflag);
            [bandwidth,score] = golden_section(beta1, beta2, 0.38197, funs, 1.0e-6, 200, true);
        else
            funs = @(bw) aic_adapt_f(bw,x,y,coors,...
                params.gweight,params.disttype,params.cnflag);
            [bandwidth,score] = golden_section(beta1, beta2, 0.38197, funs, 1.0e-6, 200, true);
        end
    end
end
end


function score = cv_f(bandwidth,x,y,coors,gweight,disttype,...
            rmseflag,cnflag)
warning('off','MATLAB:nearlySingularMatrix')

n = length(y);
cv = zeros(n,1);
for i=1:n
    if strcmp(disttype,'greatcircle')
        dist = gc_distance(coors,coors(i,:));
    else
        dist = sqrt(sum(bsxfun(@minus,coors,coors(i,:)).^2,2));
    end
    w = feval(gweight,dist,bandwidth);
    w(i) = 0;
    if cnflag
        cn = condnum(x,w);
    else
        cn = 0;
    end
            
    if cn<30
        xw = bsxfun(@times,x,w);
        lhati = x(i,:)/(x'*xw)*(xw)';
        cv(i) = y(i,:) - lhati*y;
    else
        score = Inf;
        return;
    end
end

n = length(cv);
ssr = nansum(cv.^2);
sei = sqrt(ssr/(n));
if rmseflag
    score = sei;
else
    score = ssr;
end
end


function score = par_cv_f(bandwidth,x,y,coors,gweight,disttype,...
            rmseflag,cnflag)
warning('off','MATLAB:nearlySingularMatrix')

n = length(y);
parfor i=1:n
    if strcmp(disttype,'greatcircle')
        dist = gc_distance(coors,coors(i,:));
    else
        dist = sqrt(sum(bsxfun(@minus,coors,coors(i,:)).^2,2));
    end
    w = feval(gweight,dist,bandwidth);
    w(i) = 0;
    if cnflag
        cn = condnum(x,w);
    else
        cn = 0;
    end
    
    if cn <30
        xw = bsxfun(@times,x,w);
        lhati = x(i,:)/(x'*xw)*xw';
        cv(i) = y(i,:) - lhati*y;
    else
        cv(i) = Inf;
    end
end

if any(isinf(cv))
    score = Inf;
    return;
end

n = length(cv);
ssr = nansum(cv.^2);
sei = sqrt(ssr/(n));
if rmseflag
    score = sei;
else
    score = ssr;
end
end

function score = cv_adapt_f(q,x,y,coors,gweight,disttype,...
            rmseflag,cnflag)
warning('off','MATLAB:nearlySingularMatrix')

q = round(q);
[n,m] = size(x);
cv = zeros(n,1);
for i=1:n
    if strcmp(disttype,'greatcircle')
        dist = gc_distance(coors,coors(i,:));
    else
        dist = sqrt(sum(bsxfun(@minus,coors,coors(i,:)).^2,2));
    end
    sort_dist = sort(dist);
    bandwidth = sort_dist(q);
    w = feval(gweight,dist,bandwidth);
    w(i) = 0;
    if cnflag
        cn = condnum(x,w);
    else
        cn = 0;
    end
    
    if cn<30
        xw = bsxfun(@times,x,w);
        lhati = x(i,:)/(x'*xw)*xw';
        cv(i) = y(i,:) - lhati*y;
    else
%             sw = sqrt(w);
%             Xw = bsxfun(@times,x,sw);
%             yw = y.*sw;
%             Xsd = [std(x(:,1:end-1)) 1];
%             Xws = bsxfun(@rdivide,Xw,Xsd);
%             ysd = std(y);
%             yws = yw / ysd;
%             beta = (((Xws'*Xws+la*eye(m))\(Xws'*yws))*ysd)./Xsd';
%             cv(i) = y(i,:) - x(i,:)*beta;
            
        score = Inf;
        return;
    end
end

n = length(cv);
ssr = nansum(cv.^2);
sei = sqrt(ssr/(n));
if rmseflag
    score = sei;
else
    score = ssr;
end
end

function score = par_cv_adapt_f(q,x,y,coors,gweight,disttype,...
            rmseflag,cnflag)
warning('off','MATLAB:nearlySingularMatrix')

q = round(q);
n = length(y);
%     cv = zeros(n,1);
parfor i=1:n
    if strcmp(disttype,'greatcircle')
        dist = gc_distance(coors,coors(i,:));
    else
        dist = sqrt(sum(bsxfun(@minus,coors,coors(i,:)).^2,2));
    end
    sort_dist = sort(dist);
    bandwidth = sort_dist(q);
    w = feval(gweight,dist,bandwidth);
    w(i) = 0;
    if cnflag
        cn = condnum(x,w);
    else
        cn = 0;
    end
    if cn<30
        xw = bsxfun(@times,x,w);
        lhati = x(i,:)/(x'*xw)*xw';
        cv(i) = y(i,:) - lhati*y;
    else
        cv(i) = Inf;
    end
end

if any(isinf(cv))
    score = Inf;
    return;
end

n = length(cv);
ssr = nansum(cv.^2);
sei = sqrt(ssr/(n));
if rmseflag
    score = sei;
else
    score = ssr;
end
end

function score = aic_f(bandwidth,x,y,coors,gweight,disttype,cnflag)
warning('off','MATLAB:nearlySingularMatrix')

[n,m]=size(x);
v1 = 0;
v2 = 0;
sigma2 = zeros(n,1);
for i=1:n
    if strcmp(disttype,'greatcircle')
        dist = gc_distance(coors,coors(i,:));
    else
        dist = sqrt(sum(bsxfun(@minus,coors,coors(i,:)).^2,2));
    end
    w = feval(gweight,dist,bandwidth);
    if cnflag
        cn = condnum(x,w);
    else
        cn = 0;
    end
    if cn<30
        xw = bsxfun(@times,x,w);
        lhati = x(i,:)/(x'*xw)*xw';
        v2 = v2 + sum(lhati.^2);
        v1 = v1 + lhati(1,i);
        lhati = -lhati;
        lhati(1,i) = lhati(1,i) + 1;
        sigma2(i) = lhati*y;
    else
        score = Inf;
        return;
    end
    
end

ssr = sum(sigma2.^2);
sei = sqrt(ssr/n);
% nobs2 = n/2;
% llf = -log(ssr)*nobs2;
% llf = llf - (1+log(pi/nobs2))*nobs2;
% score = -2*llf+2*n*(v1+1)/(n-v1-2);
score = 2*n*log(sei)  + n*log(2*pi) + n*(v1+n)/(n-v1-2);
end

function score = par_aic_f(bandwidth,x,y,coors,gweight,disttype,cnflag)
warning('off','MATLAB:nearlySingularMatrix')

[n,m]=size(x);
% v1 = 0;
% v2 = 0;
% sigma2 = zeros(n,1);
parfor i=1:n
    if strcmp(disttype,'greatcircle')
        dist = gc_distance(coors,coors(i,:));
    else
        dist = sqrt(sum(bsxfun(@minus,coors,coors(i,:)).^2,2));
    end
    w = feval(gweight,dist,bandwidth);
    if cnflag
        cn = condnum(x,w);
    else
        cn = 0;
    end
    if cn<30
        xw = bsxfun(@times,x,w);
        lhati = x(i,:)/(x'*xw)*xw';
        v2(i) = sum(lhati.^2);
        v1(i) = lhati(1,i);
        lhati = -lhati;
        lhati(1,i) = lhati(1,i) + 1;
        sigma2(i) = lhati*y;
    else
        sigma2(i) = Inf;
    end
    
end

if any(isinf(sigma2))
    score = Inf;
    return;
end
v1 = sum(v1);
v2 = sum(v2);

ssr = sum(sigma2.^2);
sei = sqrt(ssr/n);
% nobs2 = n/2;
% llf = -log(ssr)*nobs2;
% llf = llf - (1+log(pi/nobs2))*nobs2;
% score = -2*llf+2*n*(v1+1)/(n-v1-2);
score = 2*n*log(sei)  + n*log(2*pi) + n*(v1+n)/(n-v1-2);
end

function score = aic_adapt_f(q,x,y,coors,gweight,disttype,cnflag)
warning('off','MATLAB:nearlySingularMatrix')

q = round(q);
[n,m]=size(x);
v2 = 0;
v1 = 0;
sigma2 = zeros(n,1);
for i=1:n
    if strcmp(disttype,'greatcircle')
        dist = gc_distance(coors,coors(i,:));
    else
        dist = sqrt(sum(bsxfun(@minus,coors,coors(i,:)).^2,2));
    end
    sort_dist = sort(dist);
    bw = sort_dist(q);
    w = feval(gweight,dist,bw);
    if cnflag
        cn = condnum(x,w);
    else
        cn = 0;
    end
    if cn<30
        xw = bsxfun(@times,x,w);
        lhati = x(i,:)/(x'*xw)*xw';
        v2 = v2 + sum(lhati.^2);
        v1 = v1 + lhati(1,i);
        lhati = -lhati;
        lhati(1,i) = lhati(1,i) + 1;
        sigma2(i) = lhati*y;
    else
        score = Inf;
        return;
    end
    
end

ssr = sum(sigma2.^2);
sei = sqrt(ssr/n);
% nobs2 = n/2;
% llf = -log(ssr)*nobs2;
% llf = llf - (1+log(pi/nobs2))*nobs2;
% score = -2*llf+2*n*(v1+1)/(n-v1-2);
score = 2*n*log(sei)  + n*log(2*pi) + n*(v1+n)/(n-v1-2);
end


function score = par_aic_adapt_f(q,x,y,coors,gweight,disttype,cnflag)
warning('off','MATLAB:nearlySingularMatrix')

q = round(q);
[n,m]=size(x);
% v2 = 0;
% v1 = 0;
% sigma2 = zeros(n,1);
parfor i=1:n
    if strcmp(disttype,'greatcircle')
        dist = gc_distance(coors,coors(i,:));
    else
        dist = sqrt(sum(bsxfun(@minus,coors,coors(i,:)).^2,2));
    end
    sort_dist = sort(dist);
    bw = sort_dist(q);
    w = feval(gweight,dist,bw);
    if cnflag
        cn = condnum(x,w);
    else
        cn = 0;
    end
    if cn<30
        xw = bsxfun(@times,x,w);
        lhati = x(i,:)/(x'*xw)*xw';
        v2(i) = sum(lhati.^2);
        v1(i) = lhati(1,i);
        lhati = -lhati;
        lhati(1,i) = lhati(1,i) + 1;
        sigma2(i) = lhati*y;
    else
        sigma2(i) = Inf;
    end
    
end

if any(isinf(sigma2))
    score = Inf;
    return;
end
v1 = sum(v1);
v2 = sum(v2);

ssr = sum(sigma2.^2);
sei = sqrt(ssr/n);
% nobs2 = n/2;
% llf = -log(ssr)*nobs2;
% llf = llf - (1+log(pi/nobs2))*nobs2;
% score = -2*llf+2*n*(v1+1)/(n-v1-2);
score = 2*n*log(sei)  + n*log(2*pi) + n*(v1+n)/(n-v1-2);
end


