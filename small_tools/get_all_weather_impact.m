function get_all_weather_impact()
save_fig =1;
DU = 2.6870e+16;
instrument = 'SAOZ';

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

if strcmp(instrument , 'GBS')
    concat_data_plot_path = 'E:\H\work\Eureka\GBS\CI\weather_impact\';
    mkdir(concat_data_plot_path);
    % load data from the following folders
    % 2016-2017 EWS record has some issue, all weather column shows 'NA', need
    % figure out what happend
    CF_temp_files{1,:} = ['E:\H\work\Eureka\GBS\CI\2010\CF_450_550_minCI_v2_VCDcodev2_rerun'];
    CF_temp_files{2,:} = ['E:\H\work\Eureka\GBS\CI\2011\CF_450_550_minCI_v2_VCDcodev2'];
    CF_temp_files{3,:} = ['E:\H\work\Eureka\GBS\CI\2013\CF_450_550_minCI_v2'];
    CF_temp_files{4,:} = ['E:\H\work\Eureka\GBS\CI\2014\CF_450_550_minCI'];
    CF_temp_files{5,:} = ['E:\H\work\Eureka\GBS\CI\2015\CF_450_550_minCI'];
    % CF_temp_files{6,:} = ['E:\H\work\Eureka\GBS\CI\2016\CF_450_550_minCI_v2_VCDcodev2'];
elseif strcmp(instrument , 'SAOZ')
    concat_data_plot_path = 'E:\H\work\Eureka\SAOZ\weather_impact\';
    mkdir(concat_data_plot_path);
    CF_temp_files{1,:} = ['E:\H\work\Eureka\SAOZ\2011\CF_450_550_minCI_fixref'];
    CF_temp_files{2,:} = ['E:\H\work\Eureka\SAOZ\2013\CF_450_550_minCI_v2_fixref'];
    CF_temp_files{3,:} = ['E:\H\work\Eureka\SAOZ\2014\CF_450_550_minCI_v2_fixref'];
    CF_temp_files{4,:} = ['E:\H\work\Eureka\SAOZ\2015\CF_450_550_minCI_v2_fixref'];
    
    %CF_temp_files{6,:} = ['E:\H\work\Eureka\SAOZ\2016\CF_450_550_minCI_v2_VCDcodev2_rerun_fixref']
end
%%
N_files = size(CF_temp_files);
for i = 1:1:N_files(1)
    cd(CF_temp_files{i});
    %cd(cell2mat(CF_temp_files{i}))
    cd ..
    weather_impact_plot_path = [pwd '\weather_impact'];
    try
        [mean_weather,delta_o3_table,final_table] = weather_impact(2,CF_temp_files{i},weather_impact_plot_path,save_fig);
    catch 
        disp(['no results from file: ' num2str(i)]);
    end
    if i == 1
        mean_weather_concat = mean_weather;
        delta_o3_table_concat = delta_o3_table;
        final_table_concat = final_table;
    else
        mean_weather_concat = [mean_weather_concat;mean_weather];
        delta_o3_table_concat = [delta_o3_table_concat;delta_o3_table];
        final_table_concat = [final_table_concat;final_table];
    end
end

cd(concat_data_plot_path);


%%
mean_weather = mean_weather_concat;
delta_o3_table = delta_o3_table_concat;
final_table = final_table_concat;

weathers = unique(delta_o3_table.weather);
N_weathers = size(weathers);
N_weathers = N_weathers(1);

%% 
% group data by weather 
mean_delta_o3_table = grpstats(delta_o3_table,'weather');
new_table = final_table_concat;
%% figure 1, 
figure;
ax = gca;
index = 1:1:N_weathers(1);
x = index;
y = mean_delta_o3_table.mean_mean_delta_o3;
y_err = mean_delta_o3_table.mean_std_delta_o3
errorbar(x,y,y_err,'.');
set(gca,'XTick',index);
set(gca,'XTickLabel',str2mat(weathers));
xmax = N_weathers(1) + 1;
xlim([0 xmax]);
xlabel('EWS reported weather');
ylabel('Delta (GBS-Brewer) Ozone VCD [DU]');
print_setting('narrow2',save_fig,'Delta_o3_vcd');
%% figure 1.1, 
figure;
ax = gca;
index = 1:1:N_weathers(1);
%plot(index,delta_o3_table.freq,'.');
bar(index,mean_delta_o3_table.mean_freq);
set(gca,'XTick',index);
set(gca,'XTickLabel',str2mat(weathers));
xlabel('EWS reported weather');
ylabel('Frequence [days]');
print_setting('narrow2',save_fig,'Delta_o3_freq');
%% figure 1.2, 
figure;
ax = gca;
index = 1:1:N_weathers(1);
x = index;
y = mean_delta_o3_table.mean_mean_delta_o3./mean_delta_o3_table.mean_mean_Brewer_o3.*100;
bar(x,y);
set(gca,'XTick',index);
set(gca,'XTickLabel',str2mat(weathers));
xmax = N_weathers(1) + 1;
xlim([0 xmax]);
xlabel('EWS reported weather');
ylabel('Delta (GBS-Brewer) Ozone VCD [%]');
print_setting('narrow2',save_fig,'Delta_percentage_o3_vcd');
%% figure 2,
figure;
hold all;
gscatter(new_table.fd,new_table.mean_vcd./DU,final_table.weather_idx);
gscatter(new_table.fd,new_table.mean_ColumnO3,final_table.weather_idx);
auto_legend = [str2mat(weathers)];
legend(auto_legend);
plot(new_table.fd,new_table.mean_vcd./DU,'s');
plot(new_table.fd,new_table.mean_ColumnO3,'x');
textbp('s: GBS; x: Brewer');
xlabel('fractional day of the year');
ylabel('Ozone VCD [DU]');
print_setting('narrow2',save_fig,'O3_timeserise_with_weather_info');
%% figure 3,
figure;
hold all;
%plot(new_table.fd,new_table.mean_vcd,'.');
%plot(new_table.fd,new_table.mean_ColumnO3,'.');
gscatter(new_table.fd,new_table.delta_o3,final_table.weather_idx);
auto_legend = [str2mat(weathers)];
legend(auto_legend);
xlabel('fractional day of the year');
ylabel('Delta (GBS-Brewer) Ozone VCD [DU]');
print_setting('narrow2',save_fig,'Delta_o3_timeserise_with_weather_info');
%% figure 3.1,
figure;
hold all;
%plot(new_table.fd,new_table.mean_vcd,'.');
%plot(new_table.fd,new_table.mean_ColumnO3,'.');
y = new_table.delta_o3./new_table.mean_ColumnO3.*100;
gscatter(new_table.fd,y,final_table.weather_idx);
auto_legend = [str2mat(weathers)];
legend(auto_legend);
xlabel('fractional day of the year');
ylabel('Delta (GBS-Brewer) Ozone VCD [%]');
print_setting('narrow2',save_fig,'Delta_percentage_o3_timeserise_with_weather_info');



%% save data
save('weather_impact.mat','mean_weather','delta_o3_table','final_table','mean_delta_o3_table');
close all;
