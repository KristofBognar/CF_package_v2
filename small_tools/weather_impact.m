function [mean_weather,delta_o3_table,final_table] =weather_impact(mode,CF_temp_file,weather_impact_plot_path,save_fig)
addpath('E:\F\Work\MatlabCode');
if nargin == 0
    mode = 1;
end
if mode == 1 % use code locally
    %load('E:\H\work\Eureka\GBS\CI\2010\CF_450_550_minCI_v2_VCDcodev2_rerun\temp.mat');
    load('E:\H\work\Eureka\GBS\CI\2011\CF_450_550_minCI_v2_VCDcodev2\temp.mat');

    save_fig = 1;
    
    %plot_path = 'E:\H\work\Eureka\GBS\CI\2010\weather_inpact\'
    weather_impact_plot_path = 'E:\H\work\Eureka\GBS\CI\2011\weather_inpact\'
elseif mode == 2 % call by other function
    try
        load([CF_temp_file '\temp.mat']);
        C = who;
        k = strfind(C,'gbs_brewer');       
        if isempty(cell2mat(k))
            try
                load([CF_temp_file '\vs_Brewer\temp.mat']);
            catch
                disp('Still no gbs_brewer table found');
            end
        end
     catch
        disp('No temp.mat found, will try to load temp_fixref.mat');
    end
    
    try
        load([CF_temp_file '\vs_Brewer\temp_fixref.mat']);
     catch
        disp('No temp_fixref.mat found, please check file path.');
    end
end

C = who;
k = strfind(C,'gbs_brewer');  
if isempty(cell2mat(k))
    disp('Error: Did not find temp file content gbs_brewer, please check file path.');
end
k = strfind(C,'EWS');  
if isempty(cell2mat(k))
    disp('Error: Did not find temp file content EWS, please check file path.');
end

mkdir(weather_impact_plot_path);
cd(weather_impact_plot_path);
%% calculate gbs - brewer VCD difference
DU = 2.6870e+16;
gbs_brewer.delta_o3 = gbs_brewer.mean_vcd./DU - gbs_brewer.mean_ColumnO3;

%% make EWS daily report
EWS.datetime = datetime(EWS.DateTime);
EWS.DoY = day(EWS.datetime,'dayofyear');

days = unique(EWS.DoY);
N_days = size(days);
for i = 1:1:N_days(1)
    TF = EWS.DoY == days(i);
    EWS_1day = EWS(TF,:); % slice 1 day out of the EWS data
    
    weathers_in_1day = unique(EWS_1day.Weather);% find all weather type reported in that day
    N_weathers = size(weathers_in_1day);
    for j = 1:1:N_weathers(1)
        freq(j) = sum(strcmp(weathers_in_1day(j),EWS_1day.Weather));
    end
    [v,ind] = max(freq);
       
    weather_median(i) = weathers_in_1day(ind);
    median_obs_hrs(i) = v;
    total_obs_hrs(i) = sum(freq);
    mean_datetime(i) = mean(EWS_1day.datetime);
    DoY(i) = days(i);
    
%     TF_localnoon = EWS_1day.datetime.Hour == 12;
%     if sum(TF_localnoon) == 1
%         weather_localnoon(i) = EWS_1day.Weather(TF_localnoon);
%     elseif sum(EWS_1day.datetime.Hour == 13) == 1
%         TF_localnoon = EWS_1day.datetime.Hour == 13;
%         weather_localnoon(i) = EWS_1day.Weather(TF_localnoon);
%     else
%         weather_localnoon(i) = NaN;
%     end
    
    freq=[];

end
mean_weather = table;
mean_weather.DoY = DoY';
mean_weather.mean_datetime = mean_datetime';
mean_weather.weather_median = weather_median';
mean_weather.median_obs_hrs = median_obs_hrs';
mean_weather.total_obs_hrs = total_obs_hrs';


%% merge gbs-brewer obs table with EWS daily report
new_table = table;
j = 1;
N_gbs_brewer = size(gbs_brewer);
for i=1:1:N_gbs_brewer
    TF = gbs_brewer(i,:).day == mean_weather.DoY;
    if sum(TF) > 0
        new_table(j,:) = [gbs_brewer(i,:),mean_weather(TF,:),];
        j = j + 1;
    end
end


%% prepair statistical table 
weathers = unique(new_table.weather_median);
N_weathers = size(weathers);

for i = 1:1:N_weathers(1)
    TF = strcmp(weathers(i), new_table.weather_median);
    new_subtable = new_table(TF,:);
    new_subtable.weather_idx = repmat(i,[sum(TF),1]);
    if i == 1
        final_table = [new_subtable];
    else
        final_table = [final_table;new_subtable];
    end
    
    weather_freq(i) = sum(TF);
    mean_delta_o3(i) = mean(new_table.delta_o3(TF,:));
    std_delta_o3(i) = std(new_table.delta_o3(TF,:));
    mean_Brewer_o3(i) = mean(gbs_brewer.mean_ColumnO3(TF,:));
end
final_table = sortrows(final_table,'UTC');

delta_o3_table = table();
delta_o3_table.weather = weathers;
delta_o3_table.freq = weather_freq';
delta_o3_table.mean_delta_o3 = mean_delta_o3';
delta_o3_table.std_delta_o3 = std_delta_o3';
delta_o3_table.mean_Brewer_o3 = mean_Brewer_o3';

%% figure 1, 
figure;
ax = gca;
index = 1:1:N_weathers(1);
errorbar(index,delta_o3_table.mean_delta_o3,delta_o3_table.std_delta_o3,'.');
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
bar(index,delta_o3_table.freq);
set(gca,'XTick',index);
set(gca,'XTickLabel',str2mat(weathers));
xlabel('EWS reported weather');
ylabel('Frequence [days]');
print_setting('narrow2',save_fig,'Delta_o3_freq');
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
save('weather_impact.mat','mean_weather','delta_o3_table','final_table');
close all;
    
    



