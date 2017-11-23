function  VCD_merged  = merge_GBS_multi_year_VCDs()
% this function can merge all GBS VCDs from multi-years into a single VCD
% table, please note this merged data table only has GBS VCD info (no Brewer or EWS record)

general_path = 'E:/H/work/Eureka/GBS/CI/';

years = 2010:2017;
VCD_merged = table;
for i = 1:numel(years)
    data_path = [general_path num2str(years(i)) '/CF_newLangely/'];
    VCD = [];
    load([data_path 'temp.mat']);

    if ~isempty(VCD)
        if isempty(VCD_merged)
            VCD_merged = VCD;
        else
            VCD_merged = [VCD_merged;VCD];
        end
    else
        disp(['No VCD table found for year ' num2str(years(i))]);
    end
end

% add UTC serieal time
for i=1:1:height(VCD_merged)
    year = str2num(VCD_merged.year(i,:));
    [day, month] = Julian2Date(year,VCD_merged.day(i));
    hour = (VCD_merged.fd(i) - fix(VCD_merged.fd(i)))*24;
    minute = (hour - fix(hour))*60;
    second = (minute - fix(minute))*60;
    VCD_merged.UTC_str(i,:) = datetime(year,month,day,fix(hour),fix(minute),fix(second));
    VCD_merged.UTC(i,:) = datenum(datevec(VCD_merged.UTC_str(i,:)));
    
end