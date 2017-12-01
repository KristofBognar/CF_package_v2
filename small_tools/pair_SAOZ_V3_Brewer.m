function output = pair_SAOZ_V3_Brewer()
% note that when using our GBS CF package to process SAOZ data, a
% SAOZ-Brewer paired table will be given, in which Brewer data have been
% resampled to semi-daily f. 
% This function won't resample Brewer data, but just grab Brewer data in those SAOZ-Brewer
% paired data table, and then pair those data with SAOZ-V3 data.

load('E:\H\work\Eureka\GBS\CI\archive\gbs_saoz_brewer_merra2_ews.mat'); % load archive datasets, which has EWS

data1 = SAOZ_V3_reformat;
data2 = SAOZ_Brewer_EWS;

% add date for SAOZ_V3 table
data1.date = datetime(data1.UTC_str.Year,data1.UTC_str.Month,data1.UTC_str.Day);
data1(data1.date.Year <2010,:) = [];% now we only want SAOZ data start from 2010
data1(isnan(data1.mean_vcd),:) = [];% filter all NaN in SAOZ-V3 data
% let's rename SAOZ_V3 table, to make sure we won't have repeated column names
varnames1 = data1.Properties.VariableNames;
for i = 1:numel(varnames1)
    varnames1_new(i) = cellstr(['SAOZ_V3_' cell2mat(varnames1(i))]);
end
data1.Properties.VariableNames = varnames1_new;

output = table;
j = 1;
for i = 1:height(data1)
    TF = (data1.SAOZ_V3_date(i) == data2.date) & (data1.SAOZ_V3_ampm(i) == data2.ampm);
    if sum(TF) == 1
        output(j,:) = [data1(i,:),data2(TF,:)];
        j=j+1;
    elseif sum(TF) == 0
        disp('Warning: one data point within SAOZ_V3 is not paired with SAOZ_Brewer_EWS!');
        disp(data1.SAOZ_V3_UTC_str(i));
    else
        disp('Warning: one data point within SAOZ_V3 can be paired with more than one SAOZ_Brewer_EWS data point! Pls check if SAOZ_Brewer_EWS has only unique value!');
    end
end


