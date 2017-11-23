function data_output = pair_MERRA2_GBS()
% this function will append MERRA2 data to GBS data table.
% if no MERRA2 data avilable, then just fill the MERRA2 columns with place
% holders. 

%load('E:\H\work\Eureka\GBS\CI\weather_impact_v2\weather_impact.mat'); % load GBS data
%GBS = final_table_concat; % this is the GBS-Brewer table, change to any other GBS table, if needed. But need make sure the GBS table has LTC column (Matlab serieal time)
load('E:\H\work\Eureka\GBS\CI\GBS_VCD_2010_2017.mat'); % load GBS VCD data
GBS = VCD_merged; % this is the GBS table, change to any other GBS table, if needed. But need make sure the GBS table has LTC column (Matlab serieal time)
GBS.UTC_str = datetime(datestr(GBS.UTC));
GBS.LTC_str = datetime(datestr(GBS.UTC -5/24));


MERRA2 = prepare_MERRA2();
% prepare a dummy table for the instance that MERRA2 has no data
varnames = MERRA2.Properties.VariableNames;
MERRA2_dummy = array2table(NaN(1,width(MERRA2)),'VariableNames',varnames);
MERRA2_dummy.MERRA2_UTC_str = NaT(1);
MERRA2_dummy.MERRA2_LTC_str = NaT(1);

% loop over GBS data
data_output = table;
for i = 1:1:height(GBS)
    [delta_t,idx] = min(abs(GBS.LTC_str(i) - MERRA2.MERRA2_LTC_str)); 
    MERRA2_temporal_res = duration(3,0,0); % note MERRA-2 has temporal res 3hr
    if delta_t <= MERRA2_temporal_res
        data_output(i,:) = [GBS(i,:),MERRA2(idx,:)];
    elseif delta_t < 2*MERRA2_temporal_res
        disp('Warning: only a MERRA-2 data found within 2times of MERRA2 tempo res!');
        data_output(i,:) = [GBS(i,:),MERRA2(idx,:)];
    else
        disp('Warning: no MERRA-2 data were found could be paried with a GBS measurement! Dummy value will be filled in.');
        data_output(i,:) = [GBS(i,:),MERRA2_dummy];
    end
end

%%
function data = prepare_MERRA2()
load('E:\H\work\MERRA\MERRA2_from_Sophie\MERRA2_table_PV_thermal_TCO_PWV_with_partial_columns_with_surface_o3.mat'); % load MERRA-2 data

data = table;
% only use necessary columns from MERRA-2 table
data.MERRA2_UTC = MERRA2.UTC; 
data.MERRA2_UTC_str = datetime(datestr(MERRA2.UTC)); 
data.MERRA2_LTC_str = datetime(datestr(MERRA2.UTC -5/24));
data.MERRA2_WMO_Tropopauses = MERRA2.WMO_Tropopauses;
data.MERRA2_Dyn_Tropopauses_PV3_5 = MERRA2.Dyn_Tropopauses_PV3_5;
data.MERRA2_Ozone = MERRA2.MERRA_Ozone;










