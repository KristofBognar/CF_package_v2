function [output1,output2] = shrink_table_for_merge(table1,table2)
% this function is a general one, which can take two tables, and compare
% the variables in them. the output tables from this function will only
% content common variables

allvars1 = table1.Properties.VariableNames; % get all variables in table 1
allvars2 = table2.Properties.VariableNames; % get all variables in table 2

for i = 1:numel(allvars1)
    TF1 = strcmp(allvars1(i),allvars2);
    if i == 1
        TF1_all = TF1;
    else
        TF1_all = TF1_all | TF1;
    end
end

for i = 1:numel(allvars2)
    TF2 = strcmp(allvars2(i),allvars1);
    if i == 1
        TF2_all = TF2;
    else
        TF2_all = TF2_all | TF2;
    end
end

output1 = table1(:,TF2_all);
output2 = table2(:,TF1_all);