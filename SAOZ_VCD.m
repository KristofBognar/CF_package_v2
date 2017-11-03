function [vcd_out,vcd_out_vec] = SAOZ_VCD(year,CF_ind,temp_file_nm,fixed_rcd)
%% input %%
%year = '2011';
%temp_file_nm = 'H:\work\Eureka\SAOZ\2011\CF_450_550_minCI\VCD\temp.mat';
load(temp_file_nm);
rcd = fixed_rcd;
min_nbr_pts = 8;
%ideal_sza_range = [86,90];
ideal_sza_range = [86,91];
min_sza = 70;

if CF_ind == true
    dscd_S = dscd_SCF;
else
    dscd_S = dscd_S;
end

addpath('E:\F\Work\VCD\VCD_code\VCD_code_2017');
%% VCD calculation %%
time_pivot = datenum([year '-01-01 00:00:00']) -1;
dscd_S.utc = dscd_S.fd + time_pivot;
dscd_S.ltc = dscd_S.fd + time_pivot - 5/24 ;
dscd_S.ltc_str = datestr(dscd_S.ltc);

dscd_S.vcd_fixref = (dscd_S.mol_dscd + rcd)./dscd_S.amf;


vcd_out_vec = [];
for day = 1:366,
    for ampm = 0:1,
        ind = find(dscd_S.day == day & dscd_S.ampm == ampm);
        if isempty(ind); continue; end
        if length(ind) < min_nbr_pts; continue; end
        
        % pull out vectors
        sza_twi = dscd_S.sza(ind);
        vcd_twi = dscd_S.vcd(ind);
        vcd_fixref_twi = dscd_S.vcd_fixref(ind);
        saa_twi = dscd_S.saa(ind);
        fd_twi = dscd_S.fd(ind);
        sigma_twi = dscd_S.vcd_err(ind);
        dscd_err_twi = dscd_S.err(ind);
        amf_twi = dscd_S.amf(ind);
        
        % pull out relevant indices for range of measurements to be
        % included in the data
        [sza_range, sza_range_j] = get_SZA_indices(sza_twi, ideal_sza_range);
        if min(sza_range) < min_sza; continue; end
        min_j = min(sza_range_j);
        max_j = max(sza_range_j);
        fd_range = [fd_twi(min_j), fd_twi(max_j)];
        saa_range = [saa_twi(min_j), saa_twi(max_j)];
        
        vcd = vcd_twi(min_j:max_j);
        vcd_fixref = vcd_fixref_twi(min_j:max_j);
        dscd_err = dscd_err_twi(min_j:max_j);
        amf = amf_twi(min_j:max_j);
        sigma = sigma_twi(min_j:max_j);
        
        % calculate the weigthed mean
        [w_mean_vcd, w_sigma] = get_weighted_mean(vcd, dscd_err ./ amf);
        [w_mean_vcd_fixref, w_sigma] = get_weighted_mean(vcd_fixref, dscd_err ./ amf);
        w_mean_fd = get_weighted_mean(fd_twi(min_j:max_j), dscd_err ./ amf);
        w_mean_sza = get_weighted_mean(sza_twi(min_j:max_j), dscd_err ./ amf);
        w_mean_saa = get_weighted_mean(saa_twi(min_j:max_j), dscd_err ./ amf);
        %mean_sigma = sqrt(sum(sigma .^2)) ./ length(sigma);
        mean_sigma = mean(sigma);
        
        vcd_out_vec = [vcd_out_vec; day ampm w_mean_fd fd_range...
            w_mean_sza sza_range w_mean_saa saa_range w_mean_vcd...
            w_mean_vcd_fixref w_sigma mean_sigma std(vcd) std(vcd_fixref)];
    end
end
vcd_out.day = vcd_out_vec(:,1);
vcd_out.ampm = vcd_out_vec(:,2);
vcd_out.fd = vcd_out_vec(:,3);
vcd_out.fd_min = vcd_out_vec(:,4);
vcd_out.fd_max = vcd_out_vec(:,5);
vcd_out.sza = vcd_out_vec(:,6);
vcd_out.sza_min = vcd_out_vec(:,7);
vcd_out.sza_max = vcd_out_vec(:,8);
vcd_out.saa = vcd_out_vec(:,9);
vcd_out.saa_min = vcd_out_vec(:,10);
vcd_out.saa_max = vcd_out_vec(:,11);
vcd_out.mean_vcd = vcd_out_vec(:,12); % VCD calculated use daily RCD
vcd_out.mean_vcd_fixref = vcd_out_vec(:,13);% VCD calculated use averaged RCD
vcd_out.sigma_w_vcd = vcd_out_vec(:,14);
vcd_out.sigma_mean_vcd = vcd_out_vec(:,15);
vcd_out.std_vcd = vcd_out_vec(:,16);
vcd_out.std_vcd_fixref = vcd_out_vec(:,17);


    