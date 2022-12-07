function [rf_fit, rf_pred] = rf_model(x,y,in_model,predx,names,flag,workpath)

rf_var_names = names(in_model);
len_rf_var_names = length(rf_var_names);
in_table = array2table([y,x(:,in_model)],...
    'VariableNames', ['res',rf_var_names]);
out_table = array2table(predx(:,in_model),...
    'VariableNames', [rf_var_names]);

% flag = num2str(randi(10000000000));
str_intable = [workpath filesep 'in_table_' flag '.csv'];
if exist(str_intable,'file')
    flag = num2str(randi(10000000000));
    str_intable = [workpath filesep 'in_table_' flag '.csv'];
end
writetable(in_table, str_intable);

str_outtable = [workpath filesep 'out_table_' flag '.csv'];
if exist(str_outtable,'file')
    flag = num2str(randi(10000000000));
    str_outtable = [workpath filesep 'out_table_' flag '.csv'];
end
writetable(out_table, str_outtable);
% toc

%% 写RF代码
% tic
str_R_code = [workpath filesep 'RFcode_' flag '.R'];
if exist(str_R_code,'file')
    flag = num2str(randi(10000000000));
    str_R_code = [workpath filesep 'RFcode_' flag '.R'];
end
fid = fopen(str_R_code,'w');
fprintf(fid,'library(randomForest)\r\n');
fprintf(fid,'\r\n');
fprintf(fid,'indata<-read.table("%s",header=TRUE,sep=",")\r\n',strrep(str_intable,'\','\\'));
fprintf(fid,'rf<-randomForest(res~');
for i=1:len_rf_var_names
    fprintf(fid,'+%s',rf_var_names{i});
end
fprintf(fid,', data=indata, ntree=100, important=TRUE, proximity=TRUE)\r\n');
fprintf(fid,'RF_fit<-rf$predicted\r\n');

fprintf(fid,'\r\n');
fprintf(fid,'outdata<-read.table("%s",header=TRUE,sep=",")\r\n',strrep(str_outtable,'\','\\'));
fprintf(fid,'RF_pred<-predict(rf,outdata)\r\n');

str_outtable_rf_fit = [workpath filesep 'out_table_rf_fit' flag '.csv'];
if exist(str_outtable_rf_fit,'file')
    flag = num2str(randi(10000000000));
    str_outtable_rf_fit = [workpath filesep 'out_table_rf_fit' flag '.csv'];
end
str_outtable_rf_pred = [workpath filesep 'out_table_rf_pred' flag '.csv'];
if exist(str_outtable_rf_pred,'file')
    flag = num2str(randi(10000000000));
    str_outtable_rf_pred = [workpath filesep 'out_table_rf_fit' flag '.csv'];
end

fprintf(fid,'write.table(RF_fit,"%s",sep=",")\r\n',strrep(str_outtable_rf_fit,'\','\\'));
fprintf(fid,'write.table(RF_pred,"%s",sep=",")\r\n',strrep(str_outtable_rf_pred,'\','\\'));
fclose(fid);

% str_cmd = ['"C:\Program Files\R\R-3.5.1\bin\R.exe" CMD BATCH --vanilla --slave ' str_R_code];
str_cmd = ['Rscript ' str_R_code];
[status,cmdout] = system(str_cmd);
% toc

%% 计算�?��结果
% tic
T = readtable(str_outtable_rf_fit);
rf_fit = table2array(T(:,2));

T = readtable(str_outtable_rf_pred);
rf_pred = table2array(T(:,2));

% toc

%% 删除文件
% tic
delete(str_intable);
delete(str_outtable);
delete(str_R_code);
delete(str_outtable_rf_fit);
delete(str_outtable_rf_pred);
% toc
end
