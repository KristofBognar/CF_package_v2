function CF_quality_control()
% this function is a test, which can perform some QC for CF datasets
load('E:\H\work\Eureka\GBS\CI\archive\gbs_saoz_brewer_merra2_ews_2017');

datasets = {'GBS_Brewer_EWS','GBS_CF_Brewer_EWS','SAOZ_Brewer_EWS','SAOZ_CF_Brewer_EWS', ...
    'GBS_MERRA2_EWS','GBS_CF_MERRA2_EWS','SAOZ_MERRA2_EWS','SAOZ_CF_MERRA2_EWS',};

for j = 1:numel(datasets)
    dataset = cell2mat(datasets(j));
    eval(['data =' dataset ';']);
    N_initial = height(data);
    dates = unique(data.date);
    for i=numel(dates):-1:1
        TF = dates(i) == data.date;
        if sum(TF) == 1
            data(TF,:) = [];
        end
    end
    eval([dataset '= data ; ']);
    N_final = height(data);
    N_deleted =  N_initial - N_final;
    pN_deleted =  (N_initial - N_final)/N_initial*100;
    disp([dataset ': ' num2str(N_deleted) ' measurements (single ones) have been deleted, ' num2str(pN_deleted) ' %'])
end

save('gbs_saoz_brewer_merra2_ews_2017_high_quality.mat');