function scatter_plots_in_differen_weather()
DU = 2.6870e+16;
save_fig = 1;
load('GBS_VCD_2010_2017_MERRA2_2010_2015_EWS.mat');% load GBS/MERRA2/EWS paired data table
data = GBS_VCD_MERRA2_EWS;

data(isnan(data.MERRA2_Ozone),:) = []; % delete days without MERRA2 data
data(isnan(data.mean_vcd),:) = []; % delete days without MERRA2 data


% fig 1

x = data.mean_vcd./DU;
y = data.MERRA2_Ozone;
%dscatter(x,y);
linear_fits(x,y);
xlim([200 600]);
ylim([200 600]);
xlabel('GBS TCO [DU]');
ylabel('MERRA-2 TCO [DU]');
title('All weather types');
print_setting(1/4,save_fig,['MERRA-2_vs_GBS_all_weather']);

% fig 2

TF = (strcmp(data.weather_median_ampm, 'Clear' )) | (strcmp(data.weather_median_ampm, 'Mainly Clear' ));
x = data.mean_vcd(TF,:)./DU;
y = data.MERRA2_Ozone(TF,:);
%dscatter(x,y);
linear_fits(x,y);
xlim([200 600]);
ylim([200 600]);
xlabel('GBS TCO [DU]');
ylabel('MERRA-2 TCO [DU]');
title('Clear & Mainly Clear');
print_setting(1/4,save_fig,['MERRA-2_vs_GBS_Clear_Mainly Clear']);



% fig 3

TF = (strcmp(data.weather_median_ampm, 'Clear' )) ;
x = data.mean_vcd(TF,:)./DU;
y = data.MERRA2_Ozone(TF,:);
%dscatter(x,y);
linear_fits(x,y);
xlim([200 600]);
ylim([200 600]);
xlabel('GBS TCO [DU]');
ylabel('MERRA-2 TCO [DU]');
title('Clear');
print_setting(1/4,save_fig,['MERRA-2_vs_GBS_Clear']);

% fig 4

TF = (strcmp(data.weather_median_ampm, 'Cloudy' )) | (strcmp(data.weather_median_ampm, 'Mostly Cloudy' ));
x = data.mean_vcd(TF,:)./DU;
y = data.MERRA2_Ozone(TF,:);
%dscatter(x,y);
linear_fits(x,y);
xlim([200 600]);
ylim([200 600]);
xlabel('GBS TCO [DU]');
ylabel('MERRA-2 TCO [DU]');
title('Cloudy & Mostly Cloudy');
print_setting(1/4,save_fig,['MERRA-2_vs_GBS_Cloudy_Mostly_Cloudy']);

% fig 5

TF = (strcmp(data.weather_median_ampm, 'Cloudy' )) ;
x = data.mean_vcd(TF,:)./DU;
y = data.MERRA2_Ozone(TF,:);
%dscatter(x,y);
linear_fits(x,y);
xlim([200 600]);
ylim([200 600]);
xlabel('GBS TCO [DU]');
ylabel('MERRA-2 TCO [DU]');
title('Cloudy');
print_setting(1/4,save_fig,['MERRA-2_vs_GBS_Cloudy']);
