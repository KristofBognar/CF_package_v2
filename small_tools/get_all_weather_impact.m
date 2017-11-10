function get_all_weather_impact(mode)
%% by Xiaoyi 2017-11-08
% this function can concatinate all weather impact results from "weather_impact.m"
% 1) if want manually input file pathes, fill the mode 1 part
% then just run "get_all_weather_impact(1)"
% 2) if provide inputs in the "CF_input_file_smalltools.txt"
% then run "get_all_weather_impact(2)"

% mode = 1, use local inputs from mode 1
% mode = 2, use inputs from "CF_input_file_smalltools.txt"

if nargin == 0
    mode = 1;
end

if mode == 1
    % use local inputs here
    save_fig =1; % save fig or not, 1 = yes, 0 = no
    instrument = 'SAOZ'; % instrument name str
    table_nm = 'gbs_brewer'; % the table we want work on  
    if strcmp(instrument , 'GBS')
        concat_data_plot_path = 'E:\H\work\Eureka\GBS\CI\weather_impact\';% file path used to save final concatinated data and figrues
        % load data from the following folders
        % 2016-2017 EWS record has some issue, all weather column shows 'NA', need
        % figure out what happend
        CF_temp_files{1,:} = ['E:\H\work\Eureka\GBS\CI\2010\CF_450_550_minCI_v2_VCDcodev2_rerun\'];
        CF_temp_files{2,:} = ['E:\H\work\Eureka\GBS\CI\2011\CF_450_550_minCI_v2_VCDcodev2\'];
        CF_temp_files{3,:} = ['E:\H\work\Eureka\GBS\CI\2013\CF_450_550_minCI_v2\'];
        CF_temp_files{4,:} = ['E:\H\work\Eureka\GBS\CI\2014\CF_450_550_minCI\'];
        CF_temp_files{5,:} = ['E:\H\work\Eureka\GBS\CI\2015\CF_450_550_minCI\'];
        % CF_temp_files{6,:} = ['E:\H\work\Eureka\GBS\CI\2016\CF_450_550_minCI_v2_VCDcodev2\'];
    elseif strcmp(instrument , 'SAOZ')
        concat_data_plot_path = 'E:\H\work\Eureka\SAOZ\weather_impact\'; % file path used to save final concatinated data and figrues
        CF_temp_files{1,:} = ['E:\H\work\Eureka\SAOZ\2011\CF_450_550_minCI_fixref\'];
        CF_temp_files{2,:} = ['E:\H\work\Eureka\SAOZ\2013\CF_450_550_minCI_v2_fixref\'];
        CF_temp_files{3,:} = ['E:\H\work\Eureka\SAOZ\2014\CF_450_550_minCI_v2_fixref\'];
        CF_temp_files{4,:} = ['E:\H\work\Eureka\SAOZ\2015\CF_450_550_minCI_v2_fixref\'];
    end
else
    % use settings from input file
    % Change the current folder to the folder of this m-file.
    if(~isdeployed)
      cd(fileparts(which(mfilename)));
      cd ..
      addpath(pwd);
    end
    input_table = read_input_file_smalltools();
    
    % interpret inputs from input_table
    save_fig = input_table.save_fig; 
    instrument = input_table.instrument;
    table_nm = input_table.table_nm;
    concat_data_plot_path = input_table.concat_data_plot_path;
    general_file_path = input_table.general_file_path;
    temp_file_folder = input_table.temp_file_folder;
    years = input_table.years;
    
    N = size(years);
    for i = 1:1:N(2)
        CF_temp_files{i,:} = [general_file_path num2str(years(i)) '/' temp_file_folder];
    end
end

mkdir(concat_data_plot_path);

%%
DU = 2.6870e+16;
N_files = size(CF_temp_files);
for i = 1:1:N_files(1)
    cd(CF_temp_files{i});
    weather_impact_plot_path = [pwd '/weather_impact'];
    try
        %[mean_weather,delta_o3_table,final_table] = weather_impact(2,CF_temp_files{i},weather_impact_plot_path,save_fig);
        disp(['Working on year: ' num2str(years(i))]);
        [mean_weather,mean_weather_ampm,delta_o3_table,delta_o3_table_ampm,final_table] = weather_impact(2,table_nm,CF_temp_files{i},weather_impact_plot_path,save_fig);
    catch 
        disp(['no results from file: ' num2str(i)]);
    end
    if i == 1
        mean_weather_concat = mean_weather;
        delta_o3_table_concat = delta_o3_table;
        mean_weather_concat_ampm = mean_weather_ampm;
        delta_o3_table_concat_ampm = delta_o3_table_ampm;
        final_table_concat = final_table;
    else
        mean_weather_concat = [mean_weather_concat;mean_weather];
        delta_o3_table_concat = [delta_o3_table_concat;delta_o3_table];
        mean_weather_concat_ampm = [mean_weather_concat_ampm;mean_weather_ampm];
        delta_o3_table_concat_ampm = [delta_o3_table_concat_ampm;delta_o3_table_ampm];
        final_table_concat = [final_table_concat;final_table];
    end
end

cd(concat_data_plot_path);


mean_delta_o3_table = plot_get_all_weather_impact(instrument,final_table_concat,save_fig,'_daily');
mean_delta_o3_table_ampm = plot_get_all_weather_impact(instrument,final_table_concat,save_fig,'_ampm');
try
    mean_delta_o3_table_ref = plot_get_all_weather_impact(final_table_concat,save_fig,'_ref');
catch
end
close all;
%% save data
try
    save('weather_impact.mat','mean_delta_o3_table','mean_delta_o3_table_ampm','mean_delta_o3_table_ref','mean_weather_concat','delta_o3_table_concat','mean_weather_concat_ampm','delta_o3_table_concat_ampm','final_table_concat');
catch
    save('weather_impact.mat','mean_delta_o3_table','mean_delta_o3_table_ampm','mean_weather_concat','delta_o3_table_concat','mean_weather_concat_ampm','delta_o3_table_concat_ampm','final_table_concat');
end


function mean_delta_o3_table = plot_get_all_weather_impact(instrument,final_table_concat,save_fig,labels)
DU = 2.6870e+16;
%%

final_table = final_table_concat;

% weathers = unique(delta_o3_table.weather);
% N_weathers = size(weathers);
% N_weathers = N_weathers(1);

if strcmp(labels,'_daily')
    weathers = unique(final_table.weather_median);
    mean_delta_o3_table = grpstats(final_table_concat,'weather_median',{'mean','std'},'DataVars',{'mean_ColumnO3','delta_o3'});
    mean_delta_o3_table = sortrows(mean_delta_o3_table,'weather_median');
elseif strcmp(labels,'_ampm')
    weathers = unique(final_table.weather_median_ampm);
    mean_delta_o3_table = grpstats(final_table_concat,'weather_median_ampm',{'mean','std'},'DataVars',{'mean_ColumnO3','delta_o3'});
    mean_delta_o3_table = sortrows(mean_delta_o3_table,'weather_median_ampm');
elseif strcmp(labels,'_ref')
    weathers = unique(final_table.ref_weather);
    mean_delta_o3_table = grpstats(final_table_concat,'ref_weather',{'mean','std'},'DataVars',{'mean_ColumnO3','delta_o3'});
    mean_delta_o3_table = sortrows(mean_delta_o3_table,'ref_weather');
end
N_weathers = size(weathers);

%% 
% group data by weather 
%mean_delta_o3_table = grpstats(delta_o3_table,'weather');
%mean_delta_o3_table = grpstats(final_table_concat,'weather_median_ampm',{'mean','std'},'DataVars',{'mean_ColumnO3','delta_o3'});
new_table = final_table_concat;
%% figure 1, 
figure;
ax = gca;
index = 1:1:N_weathers(1);
x = index;
%y = mean_delta_o3_table.mean_mean_delta_o3;
y = mean_delta_o3_table.mean_delta_o3;
%y_err = mean_delta_o3_table.mean_std_delta_o3;
y_err = mean_delta_o3_table.std_delta_o3;
errorbar(x,y,y_err,'.');
set(gca,'XTick',index);
set(gca,'XTickLabel',str2mat(mean_delta_o3_table.Properties.RowNames));
xmax = N_weathers(1) + 1;
xlim([0 xmax]);
xlabel('EWS reported weather');
ylabel(['Delta (' instrument '-Brewer) Ozone VCD [DU]']);
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
xlabel('EWS reported weather');
ylabel('Number of measurements [half day]');
print_setting('narrow2',save_fig,['Delta_o3_freq' labels]);
%% figure 1.2, 
figure;hold all;
ax = gca;
index = 1:1:N_weathers(1);
x = index;
%y = mean_delta_o3_table.mean_mean_delta_o3./mean_delta_o3_table.mean_mean_Brewer_o3.*100;
y = mean_delta_o3_table.mean_delta_o3./mean_delta_o3_table.mean_mean_ColumnO3.*100;
y_err = mean_delta_o3_table.std_delta_o3./mean_delta_o3_table.mean_mean_ColumnO3.*100;
bar(x,y);
errorbar(x,y,y_err,'.');
set(gca,'XTick',index);
%set(gca,'XTickLabel',str2mat(mean_delta_o3_table.weather));
set(gca,'XTickLabel',str2mat(mean_delta_o3_table.Properties.RowNames));
xmax = N_weathers(1) + 1;
xlim([0 xmax]);
xlabel('EWS reported weather');
ylabel(['Delta (' instrument '-Brewer) Ozone VCD [%]']);
rotateXLabels( gca(), 45);
print_setting('narrow2',save_fig,['Delta_percentage_o3_vcd' labels]);
%% figure 2,
figure;
hold all;
gscatter(new_table.fd,new_table.mean_vcd./DU,final_table.weather_idx);
gscatter(new_table.fd,new_table.mean_ColumnO3,final_table.weather_idx);
auto_legend = [str2mat(mean_delta_o3_table.Properties.RowNames)];
legend(auto_legend);
plot(new_table.fd,new_table.mean_vcd./DU,'s');
plot(new_table.fd,new_table.mean_ColumnO3,'x');
textbp(['s: ' instrument '; x: Brewer']);
xlabel('fractional day of the year');
ylabel('Ozone VCD [DU]');
print_setting('narrow2',save_fig,['O3_timeserise_with_weather_info' labels]);
%% figure 3,
figure;
hold all;
%plot(new_table.fd,new_table.mean_vcd,'.');
%plot(new_table.fd,new_table.mean_ColumnO3,'.');
gscatter(new_table.fd,new_table.delta_o3,final_table.weather_idx);
auto_legend = [str2mat(mean_delta_o3_table.Properties.RowNames)];
legend(auto_legend);
xlabel('fractional day of the year');
ylabel(['Delta (' instrument '-Brewer) Ozone VCD [DU]']);
print_setting('narrow2',save_fig,['Delta_o3_timeserise_with_weather_info' labels]);
%% figure 3.1,
figure;
hold all;
%plot(new_table.fd,new_table.mean_vcd,'.');
%plot(new_table.fd,new_table.mean_ColumnO3,'.');
y = new_table.delta_o3./new_table.mean_ColumnO3.*100;
gscatter(new_table.fd,y,final_table.weather_idx);
auto_legend = [str2mat(mean_delta_o3_table.Properties.RowNames)];
legend(auto_legend);
xlabel('fractional day of the year');
ylabel(['Delta (' instrument '-Brewer) Ozone VCD [%]']);
print_setting('narrow2',save_fig,['Delta_percentage_o3_timeserise_with_weather_info' labels]);




