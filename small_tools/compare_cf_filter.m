function [data, data_cfonly, data_clearonly] = compare_cf_filter()
% this function loads two EWS info included final table data, then find the
% inner join part, left only , and right only part of the two tables.
% outputs: 
% data --> inner join part
% data_cfonly --> left only
% data_clearonly --> right only
% the output tables, can be used in X_days_low_f_test function, to perform
% uncertainty estimation

load('E:\H\work\Eureka\GBS\CI\weather_impact_v2_cf\weather_impact.mat'); % load the left data table
data_cf = final_table_concat;

load('E:\H\work\Eureka\GBS\CI\weather_impact_clear_test\weather_impact.mat');% load the right data table
data_clear = final_table_concat;

first_match = true;
first_match_cfonly = true;
data_clearonly = data_clear;
for i = height(data_cf):-1:1
    
    TF = data_clear.mean_datetime_ampm == data_cf.mean_datetime_ampm(i);
    
    if sum(TF) > 0
        data_clearonly(TF,:) = [];
        if first_match == true
            data = data_cf(i,:);
            first_match = false;
        else
            data = [data;data_cf(i,:)];
        end
    end
    
    if sum(TF) == 0
        if first_match_cfonly == true
            data_cfonly = data_cf(i,:);
            first_match_cfonly = false;
        else
            data_cfonly = [data_cfonly;data_cf(i,:)];
        end
    end
end
        
        
    
