function EWS = merge_all_EWS()
% this function could merge all EWS data into a single table
years = 1999:2017;
EWS = table;
for i = 1:numel(years)
    
    EWS_1year = [];
    EWS_raw_data_folder = ['E:\H\work\Eureka\Eureka_weather_station\' num2str(years(i)) '\'];
    try
        EWS_1year = make_EWS_weather_table_v2(EWS_raw_data_folder);
    catch
        disp(['Warning: one year is missing: ' num2str(years(i))]);
    end
    if isempty(EWS) & ~isempty(EWS_1year)
        EWS = EWS_1year;
    elseif ~isempty(EWS) & ~isempty(EWS_1year)
        EWS = [EWS;EWS_1year];
    end
        
end