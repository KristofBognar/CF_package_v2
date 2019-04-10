function output = reformat_SAOZ_V3()
% this function read in 'raw' SAOZ V3 data, and then reformat it to GBS VCD
% type table, which can be used in function such as "pair_MERRA2_GBS" and
% "weather_impact_MERRA2"

%load('E:\H\work\Eureka\GBS\CI\archive\gbs_saoz_brewer_merra2_ews');
load('C:\Projects\SAOZ\SAOZ_V3_NDACC_SZA.mat');

%data = SAOZ_V3;
data = O3all;
output = table;
am_offset = 14;
pm_offset = 22;
j = 1;
for i = 1:height(data)
    output.UTC_str(j) = datetime(data.Year(i),data.Month(i),data.Day(i),am_offset,0,0);
    output.UTC(j) = datenum(datetime(data.Year(i),data.Month(i),data.Day(i),am_offset,0,0));
    output.DoY(j) = data.DoY(i);
    output.ampm(j) = 0;
    output.mean_vcd(j) = data.O3sr(i);
    output.std_vcd(j) = data.dO3sr(i);
    j = j+1;
    output.UTC_str(j) = datetime(data.Year(i),data.Month(i),data.Day(i),pm_offset,0,0);
    output.UTC(j) = datenum(datetime(data.Year(i),data.Month(i),data.Day(i),pm_offset,0,0));    
    output.DoY(j) = data.DoY(i);
    output.ampm(j) = 1;
    output.mean_vcd(j) = data.O3ss(i);
    output.std_vcd(j) = data.dO3ss(i);
    j = j+1;
end