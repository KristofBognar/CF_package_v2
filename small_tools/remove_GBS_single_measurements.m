function data_input = remove_GBS_single_measurements()
load('C:\Users\ZhaoX\Documents\paper\CF\fig\gbs_saoz_brewer_merra2_ews_2017_high_quality_dmp_nogbs2017spring.mat');
data_input = GBS_CF;
for year = 2010:2017
    TF_year = str2num(data_input.year) == year;
    data = data_input(TF_year,:);
    for doy = 1:365
        TF_single_day = data.day == doy;
        if sum(TF_single_day) == 1
            data(TF_single_day,:)=[];
        end
    end
    data_input(TF_year,:)=[];
    data_input = [data_input;data];
end
    
