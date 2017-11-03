function UTGBS_CI = pair_CI_weather_v2(input_data,EWS)
%load('H:\work\Eureka\Eureka_weather_station\UTGBS_CI_2010.mat');% load GBS calculated CIs
%load('H:\work\Eureka\Eureka_weather_station\UTGBS_CI_2011.mat');% load GBS calculated CIs

UTGBS_CI = input_data;

% load EWS weather data (with simple lables)
% 0 = clear/mainly clear, 1 = cloudy/mostly cloudy, 2 = other condition
%load('H:\work\Eureka\Eureka_weather_station\EWS_weather2010.mat');
%load('H:\work\Eureka\Eureka_weather_station\EWS_weather2011.mat');
UTGBS_CI.UTC = datenum(UTGBS_CI.DateDDMMYYYY) + UTGBS_CI.Fractionaltime./24;
EWS.UTC = datenum(char(EWS.DateTime),'yyyy-mm-dd HH:MM') + 5/24; % convert from LT to UTC
N = size(UTGBS_CI);
for i = 1:1:N(1)
    [M,I] = min(abs(UTGBS_CI.UTC(i) - EWS.UTC));
    if M <= 1/24;
        UTGBS_CI.Weather_simple_clearL1(i) = EWS.Weather_simple_clearL1(I);
        UTGBS_CI.Weather_simple_clearL2(i) = EWS.Weather_simple_clearL2(I);
        UTGBS_CI.Weather_simple_clearL3(i) = EWS.Weather_simple_clearL3(I);
        UTGBS_CI.Weather_simple_cloudyL1(i) = EWS.Weather_simple_cloudyL1(I);
        UTGBS_CI.Weather_simple_cloudyL2(i) = EWS.Weather_simple_cloudyL2(I);
    else
        warning = datestr(UTGBS_CI.DateDDMMYYYY(i),'yyyy-mm-dd');
        disp(['find one unpaired GBS measurement on ' warning]);
        UTGBS_CI.Weather_simple_clearL1(i) = 3;
        UTGBS_CI.Weather_simple_clearL2(i) = 3;
        UTGBS_CI.Weather_simple_clearL3(i) = 3;
        UTGBS_CI.Weather_simple_cloudyL1(i) = 3;
        UTGBS_CI.Weather_simple_cloudyL2(i) = 3;
    end
end

N = size(UTGBS_CI);
TF = UTGBS_CI.Weather_simple_clearL1 == 3;
missing_weather_records = sum(TF);
p_missing_weather_records = missing_weather_records./N(1)*100;
disp([num2str(missing_weather_records) ' GBS measurements are not paired with weather records']);
disp([num2str(p_missing_weather_records) ' % of records']);
UTGBS_CI(TF,:) = [];

