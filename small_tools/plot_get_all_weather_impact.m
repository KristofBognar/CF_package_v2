function mean_delta_o3_table = plot_get_all_weather_impact(instrument,input_table,save_fig,labels)
% this function read in one data table, which should have GBS VCD, MERRA2
% VCD, and EWS weather information.
% eg.: 
% mean_delta_o3_table = plot_get_all_weather_impact('GBS',GBS_VCD_MERRA2_EWS,1,'_daily');
% mean_delta_o3_table = plot_get_all_weather_impact('GBS',GBS_VCD_MERRA2_EWS,1,'_ampm');

DU = 2.6870e+16;
%%
data = input_table;
varnames = data.Properties.VariableNames;
if sum(strcmp(varnames,'MERRA2_Ozone')) % check if we are working with GBS/MERRA-2 paired table
    data(isnan(data.MERRA2_Ozone),:) = [];
    VCD1_column_nm = 'mean_vcd';
    VCD2_column_nm = 'MERRA2_Ozone';
    benchmark_nm = 'MERRA-2';
elseif sum(strcmp(varnames,'mean_ColumnO3')) % check if we are working with GBS/Brewer paired table
    data(isnan(data.mean_ColumnO3),:) = [];
    VCD1_column_nm = 'mean_vcd';
    VCD2_column_nm = 'mean_ColumnO3';
    benchmark_nm = 'Brewer';
end
    
eval(['VCD1 = data.' VCD1_column_nm ';']);
eval(['VCD2 = data.' VCD2_column_nm ';']);

if strcmp(instrument, 'GBS') | strcmp(instrument, 'SAOZ')
    VCD1 = VCD1./DU; % if VCD1 is GBS data, need divide by DU! If this is SAOZ (processed by using GBS algorithm), we also need divide DU!
elseif strcmp(instrument, 'SAOZ-V3')
    VCD1 = VCD1; % if VCD1 is GBS data, no need divide by DU!
end
data.delta_o3 = VCD1 - VCD2;
data.p_delta_o3 = (VCD1 - VCD2)./VCD2.*100;

data(isnan(data.delta_o3),:) = [];
% weathers = unique(delta_o3_table.weather);
% N_weathers = size(weathers);
% N_weathers = N_weathers(1);

if strcmp(labels,'_daily')
    weathers = unique(data.weather_median);
    mean_delta_o3_table = grpstats(data,'weather_median',{'mean','std','sem'},'DataVars',{VCD2_column_nm,'delta_o3'});
    mean_delta_o3_table = sortrows(mean_delta_o3_table,'weather_median');
elseif strcmp(labels,'_ampm')
    weathers = unique(data.weather_median_ampm);
    mean_delta_o3_table = grpstats(data,'weather_median_ampm',{'mean','std','sem'},'DataVars',{VCD2_column_nm,'delta_o3'});
    mean_delta_o3_table = sortrows(mean_delta_o3_table,'weather_median_ampm');
elseif strcmp(labels,'_ref')
    weathers = unique(final_table.ref_weather);
    mean_delta_o3_table = grpstats(data,'ref_weather',{'mean','std','sem'},'DataVars',{VCD2_column_nm,'delta_o3'});
    mean_delta_o3_table = sortrows(mean_delta_o3_table,'ref_weather');
end
N_weathers = size(weathers);

%% 
% group data by weather 
%mean_delta_o3_table = grpstats(delta_o3_table,'weather');
%mean_delta_o3_table = grpstats(data,'weather_median_ampm',{'mean','std'},'DataVars',{VCD2_column_nm,'delta_o3'});
new_table = data;
%% figure 1, 
figure;
ax = gca;
index = 1:1:N_weathers(1);
x = index;
%y = mean_delta_o3_table.mean_mean_delta_o3;
y = mean_delta_o3_table.mean_delta_o3;
%y_err = mean_delta_o3_table.mean_std_delta_o3;
%y_err = mean_delta_o3_table.std_delta_o3;
y_err = mean_delta_o3_table.sem_delta_o3;
errorbar(x,y,y_err,'.');
set(gca,'XTick',index);
set(gca,'XTickLabel',str2mat(mean_delta_o3_table.Properties.RowNames));
xmax = N_weathers(1) + 1;
xlim([0 xmax]);
xlabel('EWS reported weather');
ylabel(['Delta (' instrument '-' benchmark_nm ') Ozone VCD [DU]']);
rotateXLabels( gca(), 45);
print_setting('narrow2',save_fig,['Delta_o3_vcd' labels]);
%% figure 1.1, 
figure;
ax = gca;
index = 1:1:N_weathers(1);
%plot(index,delta_o3_table.freq,'.');
%bar(index,mean_delta_o3_table.mean_freq.*mean_delta_o3_table.GroupCount);
bar(index,mean_delta_o3_table.GroupCount);
set(gca,'XTick',index);

%set(gca,'XTickLabel',str2mat(mean_delta_o3_table.weather));
set(gca,'XTickLabel',str2mat(mean_delta_o3_table.Properties.RowNames));
rotateXLabels( gca(), 45);
xlabel('EWS reported weather');
ylabel('Number of measurements [half day]');
print_setting('narrow2',save_fig,['Delta_o3_freq' labels]);
%% figure 1.2, 
figure;hold all;
ax = gca;
index = 1:1:N_weathers(1);
x = index;

eval(['mean_VCD2 = mean_delta_o3_table.mean_' VCD2_column_nm ';']);
y = mean_delta_o3_table.mean_delta_o3./mean_VCD2.*100;
%y_err = mean_delta_o3_table.std_delta_o3./mean_VCD2.*100;
y_err = mean_delta_o3_table.sem_delta_o3./mean_VCD2.*100;
bar(x,y);
errorbar(x,y,y_err,'.');
set(gca,'XTick',index);
%set(gca,'XTickLabel',str2mat(mean_delta_o3_table.weather));
set(gca,'XTickLabel',str2mat(mean_delta_o3_table.Properties.RowNames));
xmax = N_weathers(1) + 1;
xlim([0 xmax]);
xlabel('EWS reported weather');
ylabel(['Delta (' instrument '-' benchmark_nm ') Ozone VCD [%]']);
rotateXLabels( gca(), 45);
print_setting('narrow2',save_fig,['Delta_percentage_o3_vcd' labels]);

%% figure 2
figure;
boxplot(data.delta_o3,data.weather_median_ampm);
ylabel(['Delta (' instrument '-' benchmark_nm ') Ozone VCD [DU]']);
rotateXLabels(gca,45);
print_setting('narrow2',save_fig,['Delta_o3_vcd_boxplot' labels]);
%% figure 2.1
figure;
boxplot(data.p_delta_o3,data.weather_median_ampm);
ylabel(['Delta (' instrument '-' benchmark_nm ') Ozone VCD [%]']);
rotateXLabels(gca,45);
print_setting('narrow2',save_fig,['Delta_percentage_o3_vcd_boxplot' labels]);


