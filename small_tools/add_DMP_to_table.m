function data = add_DMP_to_table()
% this function can add DMP information (height @ Theta=490K, and sPV @
% Theta=490K) to our typical data table, such as GBS_Brewer,
% by Xiaoyi 2017-12-19

load('E:\H\work\Eureka\GBS\CI\archive\gbs_saoz_brewer_merra2_ews_2017_high_quality.mat');
% note that we have DMP information within MERRA2 table, this new MERRA-2
% table can be produced with "make_MERRA2_table_TCO_PWV_column_multi_years"

lists_of_tables = {'GBS_Brewer_EWS','SAOZ_Brewer_EWS','GBS_MERRA2_EWS','SAOZ_MERRA2_EWS','GBS_CF_MERRA2_EWS','SAOZ_CF_MERRA2_EWS', 'SAOZ_V3_reformat_MERRA2_EWS'};
%lists_of_tables = {'GBS_MERRA2_EWS'};

for i =1:numel(lists_of_tables)
    eval(['data = ' cell2mat(lists_of_tables(i)) ';']);
    
    for j = 1:height(data)
       [delta_UTC, idx] = min(abs(data.UTC(j) - MERRA2.UTC));
       if delta_UTC < 0.5 % the successful matching need have delta time less than half day
           data.MERRA2_Theta490K_heigt(j,1) = MERRA2.Theta490K_heigt(idx,1);
           data.MERRA2_sPV_at_Theta490(j,1) = MERRA2.sPV_at_Theta490(idx,1);
       else
           disp('Warning: not matching DMP record was found for a timestamp in the given table');
           data.MERRA2_Theta490K_heigt(j,1) = NaN;
           data.MERRA2_sPV_at_Theta490(j,1) = NaN;
       end
        
    end
    TF = data.MERRA2_sPV_at_Theta490 > 1.6e-4;
    meas_in_vortex = sum(TF)/height(data)*100;
    disp([cell2mat(lists_of_tables(i)) ' has ' num2str(meas_in_vortex) ' % measurements made within the vortex']);
    eval([ cell2mat(lists_of_tables(i)) ' = data;']);
    
end

save('E:\H\work\Eureka\GBS\CI\archive\gbs_saoz_brewer_merra2_ews_2017_high_quality_dmp.mat');





