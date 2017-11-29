function  [VCD_merged,VCD_CF_merged, gbs_brewer_merged, gbscf_brewer_merged]  = merge_GBS_multi_year_VCDs()
% this function can 
% 1) merge multi-year GBS(or SAOZ) VCD (and VCD_CF) tables from multi-years into a single VCD
% table
% 2) merge multi-year GBS-Brewer(or SAOZ-Brewer) (and GBSCF-Brewer) paired tables into a single GBS-Brewer
% paired table

instrument == 'GBS';
%instrument == 'SAOZ';
if strcmp(instrument,'GBS')
    general_path = 'E:/H/work/Eureka/GBS/CI/';
elseif strcmp(instrument,'SAOZ')
    general_path = 'E:/H/work/Eureka/SAOZ/';
end

years = 2010:2017;% the range of years will be merged
VCD_merged = table;
VCD_CF_merged = table;
gbs_brewer_merged = table;
gbscf_brewer_merged = table;
for i = 1:numel(years)
    if strcmp(instrument,'GBS')
        data_path = [general_path num2str(years(i)) '/CF_newLangely/'];
    elseif strcmp(instrument,'SAOZ')
        data_path = [general_path num2str(years(i)) '/CF_newLangely_fixref/'];
    end
    VCD = []; % these lines make sure we won't insert one table from a year more than one time!
    VCD_CF = [];
    gbs_brewer = [];
    gbscf_brewer = [];
    load([data_path 'temp.mat']);

    if ~isempty(VCD) && istable(VCD)
        if isempty(VCD_merged)
            VCD_merged = VCD; % this will merge all normal VCD
            VCD_CF_merged = VCD_CF; % this will merge all CF VCD
        else
            VCD_merged = [VCD_merged;VCD];
            VCD_CF_merged = [VCD_CF_merged;VCD_CF];
        end
    else
        disp(['No VCD table found for year ' num2str(years(i))]);
    end
    
    if (~isempty(gbs_brewer)) && istable(gbs_brewer)
        if isempty(gbs_brewer_merged)
            gbs_brewer_merged = gbs_brewer; % this will merge all GBS/Brewer table
            gbscf_brewer_merged = gbscf_brewer;% this will merge all GBS-CF/Brewer table
        else
            if width(gbs_brewer_merged) ~= width(gbs_brewer) % check if the gbs_brewer table has same columns as previous years'
                [gbs_brewer_merged,gbs_brewer] = shrink_table_for_merge(gbs_brewer_merged,gbs_brewer); % if gbs_brewer tables from different years have different columns, we will just merge the common columns!
                disp('Warning: loaded "gbs_brewer" tables have different width! I will just shrink them to the same size!');
            end
            if width(gbscf_brewer_merged) ~= width(gbscf_brewer)% check if the gbscf_brewer table has same columns as previous years'
                [gbscf_brewer_merged,gbs_brewer] = shrink_table_for_merge(gbscf_brewer_merged,gbs_brewer); % if gbscf_brewer tables from different years have different columns, we will just merge the common columns!
                disp('Warning: loaded "gbscf_brewer" tables have different width! I will just shrink them to the same size!');
            end
            gbs_brewer_merged = [gbs_brewer_merged;gbs_brewer];
            gbscf_brewer_merged = [gbscf_brewer_merged;gbscf_brewer];   
        end
    else
        disp(['No GBS-Brewer table found for year ' num2str(years(i))]);
    end
end

% add UTC serieal time for merged VCD table
for i=1:1:height(VCD_merged)
    year = str2num(VCD_merged.year(i,:));
    [day, month] = Julian2Date(year,VCD_merged.day(i));
    hour = (VCD_merged.fd(i) - fix(VCD_merged.fd(i)))*24;
    minute = (hour - fix(hour))*60;
    second = (minute - fix(minute))*60;
    VCD_merged.UTC_str(i,:) = datetime(year,month,day,fix(hour),fix(minute),fix(second));
    VCD_merged.UTC(i,:) = datenum(datevec(VCD_merged.UTC_str(i,:)));  
end

% add UTC serieal time for merged VCD_CF table
for i=1:1:height(VCD_CF_merged)
    year = str2num(VCD_CF_merged.year(i,:));
    [day, month] = Julian2Date(year,VCD_CF_merged.day(i));
    hour = (VCD_CF_merged.fd(i) - fix(VCD_CF_merged.fd(i)))*24;
    minute = (hour - fix(hour))*60;
    second = (minute - fix(minute))*60;
    VCD_CF_merged.UTC_str(i,:) = datetime(year,month,day,fix(hour),fix(minute),fix(second));
    VCD_CF_merged.UTC(i,:) = datenum(datevec(VCD_CF_merged.UTC_str(i,:)));   
end