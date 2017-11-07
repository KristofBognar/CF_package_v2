function [mean_weather,mean_weather_ampm,delta_o3_table,delta_o3_table2,final_table2] =weather_impact(mode,CF_temp_file,weather_impact_plot_path,save_fig)

% Change the current folder to the folder of this m-file.
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end

addpath('E:\F\Work\MatlabCode');
if nargin == 0
    mode = 1;
end
if mode == 1 % use code locally
    %load('E:\H\work\Eureka\GBS\CI\2010\CF_450_550_minCI_v2_VCDcodev2_rerun\temp.mat');
    %load('E:\H\work\Eureka\GBS\CI\2011\CF_450_550_minCI_v2_VCDcodev2\temp.mat');
    %load('E:\H\work\Eureka\GBS\CI\2014\CF_450_550_minCI_v2_VCDcodeonGit_test\temp.mat');
    %load('E:\H\work\Eureka\GBS\CI\2011\CF_VCD_CF_onGit_test2\temp.mat');
    temp_file_path = 'E:\H\work\Eureka\GBS\CI\2010\CF_newLangely\';
    load([temp_file_path 'temp.mat']);
    cd(temp_file_path);
    save_fig = 1;
    
    %plot_path = 'E:\H\work\Eureka\GBS\CI\2010\weather_impact\'
    weather_impact_plot_path = [temp_file_path 'weather_impact\'];
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

%% make EWS ampm report
%EWS.datetime = datetime(EWS.DateTime);
%EWS.DoY = day(EWS.datetime,'dayofyear');

%days = unique(EWS.DoY);
%N_days = size(days);
k = 1;
for i = 1:1:N_days(1)
    TF = EWS.DoY == days(i);
    EWS_1day = EWS(TF,:); % slice 1 day out of the EWS data
    
    TF_am = EWS_1day.Time.Hour <= 12;
    TF_pm = EWS_1day.Time.Hour >= 12;
    
    if sum(TF_am) > 0
        EWS_1day_am = EWS_1day(TF_am,:); % slice am values from 1 day out of the EWS data    
        weathers_in_1day_am = unique(EWS_1day_am.Weather);% find all weather type reported in that day
        N_weathers = size(weathers_in_1day_am);
        for j = 1:1:N_weathers(1)
            freq(j) = sum(strcmp(weathers_in_1day_am(j),EWS_1day_am.Weather));
        end
        [v,ind] = max(freq);       
        weather_median_ampm(k) = weathers_in_1day_am(ind);
        median_obs_hrs_ampm(k) = v;
        total_obs_hrs_ampm(k) = sum(freq);
        mean_datetime_ampm(k) = mean(EWS_1day_am.datetime);
        DoY(k) = days(i);
        ampm(k) = 0;
        freq=[];
        k= k + 1;
    else
        disp('No morning EWS records found!');
    end
    
    if sum(TF_pm) > 0
        EWS_1day_pm = EWS_1day(TF_pm,:); % slice pm values from 1 day out of the EWS data
        weathers_in_1day_pm = unique(EWS_1day_pm.Weather);% find all weather type reported in that day
        N_weathers = size(weathers_in_1day_pm);
        for j = 1:1:N_weathers(1)
            freq(j) = sum(strcmp(weathers_in_1day_pm(j),EWS_1day_pm.Weather));
        end
        [v,ind] = max(freq);       
        weather_median_ampm(k) = weathers_in_1day_pm(ind);
        median_obs_hrs_ampm(k) = v;
        total_obs_hrs_ampm(k) = sum(freq);
        mean_datetime_ampm(k) = mean(EWS_1day_pm.datetime);
        DoY(k) = days(i);
        ampm(k) = 1;
        freq=[];
        k = k + 1;
    else
        disp('No afternoon EWS records found!');
    end

end
mean_weather_ampm = table;
mean_weather_ampm.DoY_ampm = DoY';
mean_weather_ampm.EWS_ampm = ampm';
mean_weather_ampm.mean_datetime_ampm = mean_datetime_ampm';
mean_weather_ampm.weather_median_ampm = weather_median_ampm';
mean_weather_ampm.median_obs_hrs_ampm = median_obs_hrs_ampm';
mean_weather_ampm.total_obs_hrs_ampm = total_obs_hrs_ampm';

%% merge gbs-brewer obs table with EWS daily report
new_table = table;
j = 1;
N_gbs_brewer = size(gbs_brewer);
TF_ref_time_exist = strcmp('ref_sza_utc1',gbs_brewer.Properties.VariableNames);
if sum(TF_ref_time_exist) == 0
    disp('No reference spec time was found! Will not search weather record for ref spec.');
end
for i=1:1:N_gbs_brewer
    TF = gbs_brewer(i,:).day == mean_weather.DoY;
    TF_ampm = (gbs_brewer(i,:).day == mean_weather_ampm.DoY_ampm) & (gbs_brewer(i,:).ampm == mean_weather_ampm.EWS_ampm);
    
    
    if sum(TF_ref_time_exist) > 0
        ref_time1 = datetime(datevec(gbs_brewer.ref_sza_utc1));
        ref_time2 = datetime(datevec(gbs_brewer.ref_sza_utc2));
        TF_ref1_isnat = isnat(ref_time1);
        ref_time1(TF_ref1_isnat,:) = ref_time2(TF_ref1_isnat,:);

        UTC_offset = 5;
        EWS_ref_weather = table;
        TF_ref = (gbs_brewer(i,:).day == EWS.DoY) & (ref_time1(i,:).Hour -UTC_offset == EWS.datetime.Hour);
        if sum(TF_ref) == 1
            EWS_ref_weather.ref_weather = EWS.Weather(TF_ref,:);
            EWS_ref_weather.ref_datetime_in_EWS =  EWS.datetime(TF_ref,:);
        elseif sum(TF_ref) == 0
            disp('Warning: no weather record was found match with time of ref sepc');
            disp(['Year:' num2str(gbs_brewer(i,:).year) '; Day:' num2str(gbs_brewer(i,:).day)]);
            disp('Warning: will try to find record in +1 hr');
            TF_ref = (gbs_brewer(i,:).day == EWS.DoY) & (ref_time1(i,:).Hour -UTC_offset +1 == EWS.datetime.Hour);
            if sum(TF_ref) == 1
                EWS_ref_weather.ref_weather = EWS.Weather(TF_ref,:);
                EWS_ref_weather.ref_datetime_in_EWS =  EWS.datetime(TF_ref,:);
            elseif sum(TF_ref) == 0
                disp('Warning: no weather record was found match with time of ref sepc + 1hr');
                disp('Warning: will try to find record in -1 hr');
                TF_ref = (gbs_brewer(i,:).day == EWS.DoY) & (ref_time1(i,:).Hour -UTC_offset -1 == EWS.datetime.Hour);
                if sum(TF_ref) == 1
                    EWS_ref_weather.ref_weather = EWS.Weather(TF_ref,:); 
                    EWS_ref_weather.ref_datetime_in_EWS =  EWS.datetime(TF_ref,:);
                elseif sum(TF_ref) == 0
                    disp('Warning: no weather record was found match with time of ref sepc +/- 1hr');
                    disp('Warning: will just place NaN in there!');
                    EWS_ref_weather.ref_weather = 'NaN'; 
                    EWS_ref_weather.ref_datetime_in_EWS =  'NaT';
                    
                end
            end
        end
    else
        %disp('No reference spec time was found! Will not search weather record for ref spec.');
    end

        
    if sum(TF_ref_time_exist) > 0
        if (sum(TF) > 0) & (sum(TF_ampm) > 0) & (sum(TF_ref) > 0) 
            new_table(j,:) = [gbs_brewer(i,:),mean_weather(TF,:),mean_weather_ampm(TF_ampm,:), EWS_ref_weather];
            j = j + 1;
        end
    else
        if (sum(TF) > 0) & (sum(TF_ampm) > 0)
            new_table(j,:) = [gbs_brewer(i,:),mean_weather(TF,:),mean_weather_ampm(TF_ampm,:)];
            j = j + 1;
        end
    end
end


%% prepair statistical table (daily median)
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


%% prepair statistical table (ampm median)
weathers = unique(final_table.weather_median_ampm);
N_weathers = size(weathers);

for i = 1:1:N_weathers(1)
    TF = strcmp(weathers(i), final_table.weather_median_ampm);
    new_subtable = final_table(TF,:);
    new_subtable.weather_idx = repmat(i,[sum(TF),1]);
    if i == 1
        final_table2 = [new_subtable];
    else
        final_table2 = [final_table;new_subtable];
    end
    
    weather_freq_ampm(i) = sum(TF);
    mean_delta_o3_ampm(i) = mean(final_table.delta_o3(TF,:));
    std_delta_o3_ampm(i) = std(final_table.delta_o3(TF,:));
    mean_Brewer_o3_ampm(i) = mean(gbs_brewer.mean_ColumnO3(TF,:));
end
final_table2 = sortrows(final_table2,'UTC');

delta_o3_table2 = table();
delta_o3_table2.weather = weathers;
delta_o3_table2.freq = weather_freq_ampm';
delta_o3_table2.mean_delta_o3 = mean_delta_o3_ampm';
delta_o3_table2.std_delta_o3 = std_delta_o3_ampm';
delta_o3_table2.mean_Brewer_o3 = mean_Brewer_o3_ampm';

%% save data
save('weather_impact.mat','mean_weather','mean_weather_ampm','delta_o3_table','delta_o3_table2','final_table2');


weather_inpact_plots(delta_o3_table,final_table2,save_fig,'_daily');
weather_inpact_plots(delta_o3_table2,final_table2,save_fig,'_ampm') ;
%% 
function weather_inpact_plots(delta_o3_table,final_table,save_fig,labels)
DU = 2.6870e+16;
% figure 1,
figure;
ax = gca;
if strcmp(labels,'_daily')
    weathers = unique(final_table.weather_median);
elseif strcmp(labels,'_ampm')
    weathers = unique(final_table.weather_median_ampm);
end
N_weathers = size(weathers);

index = 1:1:N_weathers(1);
errorbar(index,delta_o3_table.mean_delta_o3,delta_o3_table.std_delta_o3,'.');
set(gca,'XTick',index);
set(gca,'XTickLabel',str2mat(delta_o3_table.weather));
xmax = N_weathers(1) + 1;
xlim([0 xmax]);
xlabel('EWS reported weather');
ylabel('Delta (GBS-Brewer) Ozone VCD [DU]');
print_setting('narrow2',save_fig,['Delta_o3_vcd' labels]);
%% figure 1.1, 
figure;
ax = gca;
index = 1:1:N_weathers(1);
%plot(index,delta_o3_table.freq,'.');
bar(index,delta_o3_table.freq);
set(gca,'XTick',index);
set(gca,'XTickLabel',str2mat(delta_o3_table.weather));
xlabel('EWS reported weather');
ylabel('Frequence [days]');
print_setting('narrow2',save_fig,['Delta_o3_freq' labels]);
%% figure 2,
figure;
hold all;
gscatter(final_table.fd,final_table.mean_vcd./DU,final_table.weather_idx);
gscatter(final_table.fd,final_table.mean_ColumnO3,final_table.weather_idx);
auto_legend = [str2mat(delta_o3_table.weather)];
legend(auto_legend);
plot(final_table.fd,final_table.mean_vcd./DU,'s');
plot(final_table.fd,final_table.mean_ColumnO3,'x');
textbp('s: GBS; x: Brewer');
xlabel('fractional day of the year');
ylabel('Ozone VCD [DU]');
print_setting('narrow2',save_fig,['O3_timeserise_with_weather_info' labels]);
%% figure 3,
figure;
hold all;
%plot(final_table.fd,final_table.mean_vcd,'.');
%plot(final_table.fd,final_table.mean_ColumnO3,'.');
gscatter(final_table.fd,final_table.delta_o3,final_table.weather_idx);
auto_legend = [str2mat(delta_o3_table.weather)];
legend(auto_legend);
xlabel('fractional day of the year');
ylabel('Delta (GBS-Brewer) Ozone VCD [DU]');
print_setting('narrow2',save_fig,['Delta_o3_timeserise_with_weather_info' labels]);
%% figure 3.1,
figure;
hold all;
%plot(final_table.fd,final_table.mean_vcd,'.');
%plot(final_table.fd,final_table.mean_ColumnO3,'.');
y = final_table.delta_o3./final_table.mean_ColumnO3.*100;
gscatter(final_table.fd,y,final_table.weather_idx);
auto_legend = [str2mat(delta_o3_table.weather)];
legend(auto_legend);
xlabel('fractional day of the year');
ylabel('Delta (GBS-Brewer) Ozone VCD [%]');
print_setting('narrow2',save_fig,['Delta_percentage_o3_timeserise_with_weather_info' labels]);


close all;

    
    



