function new_table = weather_impact_MERRA2()
% this function read in GBS-MERRA2 data, and EWS data, then merge them

% load GBS VCD paired with MERRA data table
load('E:\H\work\Eureka\GBS\CI\MERRA2\GBS_VCD_2010_2017_MERRA2_2010_2015.mat');
% load merged multi-year EWS data
load('E:\H\work\Eureka\Eureka_weather_station\EWS_1999_2017.mat');

data = GBS_VCD_MERRA2; % just rename it
data.date = datetime(data.UTC_str.Year,data.UTC_str.Month,data.UTC_str.Day);

% give the name of instrument
input_table = table;
input_table.instrument = 'GBS';
save_fig = 0;
labels = 'test';

% filter EWS
EWS((EWS.Year < 2010),:) = [];
%% make EWS daily report
EWS.datetime = datetime(EWS.DateTime);
EWS.DoY = day(EWS.datetime,'dayofyear');
EWS.date = datetime(EWS.datetime.Year,EWS.datetime.Month,EWS.datetime.Day);

days = unique(EWS.date);
N_days = size(days);
for i = 1:1:N_days(1)
    TF = EWS.date == days(i);
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
    date(i) = days(i);
    
    freq=[];

end
mean_weather = table;
mean_weather.ews_date = date';
mean_weather.mean_datetime = mean_datetime';
mean_weather.weather_median = weather_median';
mean_weather.median_obs_hrs = median_obs_hrs';
mean_weather.total_obs_hrs = total_obs_hrs';

%% make EWS ampm report

k = 1;
for i = 1:1:N_days(1)
    TF = EWS.date == days(i);
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
        date(k) = days(i);
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
        date(k) = days(i);
        ampm(k) = 1;
        freq=[];
        k = k + 1;
    else
        disp('No afternoon EWS records found!');
    end

end
mean_weather_ampm = table;
mean_weather_ampm.ews_date_ampm = date';
mean_weather_ampm.EWS_ampm = ampm';
mean_weather_ampm.mean_datetime_ampm = mean_datetime_ampm';
mean_weather_ampm.weather_median_ampm = weather_median_ampm';
mean_weather_ampm.median_obs_hrs_ampm = median_obs_hrs_ampm';
mean_weather_ampm.total_obs_hrs_ampm = total_obs_hrs_ampm';

%% merge VCD table with EWS daily report
new_table = table;
j = 1;
N_data = size(data);
TF_ref_time_exist = strcmp('ref_sza_utc1',data.Properties.VariableNames);
if sum(TF_ref_time_exist) == 0
    disp('No reference spec time was found! Will not search weather record for ref spec.');
end
for i=1:1:N_data
    TF = data(i,:).date == mean_weather.ews_date;
    TF_ampm = (data(i,:).date == mean_weather_ampm.ews_date_ampm) & (data(i,:).ampm == mean_weather_ampm.EWS_ampm);
    
    
    if (sum(TF_ref_time_exist) > 0) && strcmp(input_table.instrument,'GBS')
        ref_time1 = datetime(datevec(data.ref_sza_utc1));
        ref_time2 = datetime(datevec(data.ref_sza_utc2));
        TF_ref1_isnat = isnat(ref_time1);
        ref_time1(TF_ref1_isnat,:) = ref_time2(TF_ref1_isnat,:);

        UTC_offset = 5;
        EWS_ref_weather = table;
        TF_ref = (data(i,:).date == EWS.date) & (ref_time1(i,:).Hour -UTC_offset == EWS.datetime.Hour);
        if sum(TF_ref) == 1
            EWS_ref_weather.ref_weather = EWS.Weather(TF_ref,:);
            EWS_ref_weather.ref_datetime_in_EWS =  EWS.datetime(TF_ref,:);
        elseif sum(TF_ref) == 0
            disp('Warning: no weather record was found match with time of ref sepc');
            disp(['Year:' num2str(data(i,:).year) '; Day:' num2str(data(i,:).day)]);
            disp('Warning: will try to find record in +1 hr');
            TF_ref = (data(i,:).date == EWS.date) & (ref_time1(i,:).Hour -UTC_offset +1 == EWS.datetime.Hour);
            if sum(TF_ref) == 1
                EWS_ref_weather.ref_weather = EWS.Weather(TF_ref,:);
                EWS_ref_weather.ref_datetime_in_EWS =  EWS.datetime(TF_ref,:);
            elseif sum(TF_ref) == 0
                disp('Warning: no weather record was found match with time of ref sepc + 1hr');
                disp('Warning: will try to find record in -1 hr');
                TF_ref = (data(i,:).date == EWS.date) & (ref_time1(i,:).Hour -UTC_offset -1 == EWS.datetime.Hour);
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

        
    if (sum(TF_ref_time_exist) > 0) && strcmp(input_table.instrument,'GBS')
        if (sum(TF) > 0) & (sum(TF_ampm) > 0) & (sum(TF_ref) > 0) 
            new_table(j,:) = [data(i,:),mean_weather(TF,:),mean_weather_ampm(TF_ampm,:), EWS_ref_weather];
            j = j + 1;
        end
    else
        if (sum(TF) > 0) & (sum(TF_ampm) > 0)
            new_table(j,:) = [data(i,:),mean_weather(TF,:),mean_weather_ampm(TF_ampm,:)];
            j = j + 1;
        end
    end
end

    
mean_delta_o3_table = plot_get_all_weather_impact(input_table.instrument,new_table,save_fig,'_daily')
mean_delta_o3_table_ampm = plot_get_all_weather_impact(input_table.instrument,new_table,save_fig,'_ampm')

