function [dd,xc,yc] = datadensitymat(x,y,varargin)
%Calculate the number of the data pairs under the specified regular grid 
%
% USAGE:
%   [dd,xc,yc] = datadensitymat(x,y,'Name1',Value1,'Name2',Value2,...);
%
% INPUT:
%   x              -   n*1 verctor containing one type observation at horizontal axis
%   y              -   n*1 vector containing another type observation at vertical axis
%   varargin -   Name/Value paris, used to specify the regular grid
%                       xmin: the minimum limit of the regualr grid at horizontal axis
%                       xmax: the maximum limit of the regualr grid at horizontal axis
%                       ymin: the minimum limit of the regualr grid at vertical axis
%                       ymax: the maximum limit of the regualr grid at vertical axis
%                       M: the grid number at horizontal axis for the regular grid
%                       N: the grid number at vertical axis for the regular grid
%
% OUTPUT:
%   dd - the output data pairs number at each regular grid
%   xc  - the regular grid x limit at horizontal axis
%   yc  - the regular grid y limit at vertical axis

if nargin<2
    error('you must input two variable')
end

% error checking
if size(y,1) ~= size(x,1)
    error('x , y must have the same number of rows')
end

nandex = isnan(x) | isnan(y);
x=x(~nandex);
y = y(~nandex);

% check input using PARSEARGS
params.M          = 500;
params.N           = 500;
params.xmin     = min(x);
params.xmax     = max(x);
params.ymin     = min(y);
params.ymax     = max(y);
params = parseargs(params,varargin{:});

if params.M<=0
    params.M=500;
end

if  params.N<=0
    params.N=500;
end

xedge = linspace(params.xmin,params.xmax,params.N+1);
xc = (xedge + [xedge(2:end) NaN])/2;
xc(end)=[];

yedge = linspace(params.ymin,params.ymax,params.M+1);
yc = (yedge + [yedge(2:end) NaN])/2;
yc(end)=[];

[~,~,xbins] = histcounts(x,xedge);
[~,~,ybins] = histcounts(y,yedge);
dex = xbins==0 | ybins==0;
xbins = xbins(~dex);
ybins = ybins(~dex);
subs = sub2ind([params.M,params.N],ybins,xbins);
tbl = tabulate(subs);
dd = zeros(params.M,params.N);
dd(tbl(:,1)) = tbl(:,2);

end