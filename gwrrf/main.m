clear
clc

startmatlabpool(8);

T = readtable('PM1320_modelling.xlsx');
txt = T.Properties.VariableNames ;
num = table2array(T);

y = num(:,1);
coors = num(:,2:3);
x = num(:,5:end);
x = [x ones(length(y),1)];
names = [txt(5:end) 'Intercept'];
[n,q]=size(x);

inmodel = false(0,q);
for i=0:q-1
    tmp1 = 1:q-1;
    tmp2 = nchoosek(tmp1,i);
    len_tmp2 = size(tmp2,1);
    flag = false(len_tmp2,q-1);
    for j=1:len_tmp2
        flag(j,tmp2(j,:)) = 1;
    end
    inmodel = [inmodel;[true(len_tmp2,1) flag]];
end

in_model([8,9,10] )= true;

% gwr fitting
[bw,score] = gwr_select_bandwidth(x(:,in_model),y,coors(:,1:2));
[gwr_fit,gwr_betas,ci,se_pred,se_betas] = gwr(x(:,in_model),y,coors(:,1:2),bw,...
    'predx',x(:,in_model),'predcoor',coors(:,1:2));
save Fit_result bw score gwr_fit x y coors in_model
rsquare(y,gwr_fit)

% ten fold cross validation
rng(100);
m = length(y);
id = randperm(m);
len = floor(m/10);
gwr_validate = zeros(m,1);
for i=1:10
    start_id = len*i-len+1;
    if i<10
        end_id = len*i;
    else
        end_id = m;
    end
    xt = x(id(start_id:end_id),:);
    yt = y(id(start_id:end_id),:);
    coorst = coors(id(start_id:end_id),:);
    idm = id;
    idm(start_id:end_id) = [];
    xm = x(idm,:);
    ym = y(idm,:);
    coorsm = coors(idm,:);
    
    [bw,score] = gwr_select_bandwidth(xm(:,in_model),ym,coorsm(:,1:2));
    if isinf(score)
        gwr_validate(id(start_id:end_id),:) = Inf;
        break;
    end
    [predm,betas,ci,se_pred,se_betas] = gwr(xm(:,in_model),ym,coorsm(:,1:2),bw);
    [predt,betas,ci,se_pred,se_betas] = gwr(xm(:,in_model),ym,coorsm(:,1:2),bw,'predx',xt(:,in_model),'predcoor',coorst(:,1:2));
    gwr_validate(id(start_id:end_id),:) = predt;
    
    str = sprintf('Validate_result_%d',i);
    save(str,'ym','xm','coorsm','predm','yt','xt','coorst','predt','in_model');
end
rsquare(y,gwr_validate)
