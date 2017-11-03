function data = ci_o4_testing_plot_v2(input_data,year,plot_path,save_fig)
%% input file %%%

load(input_data);
%load('H:\work\Eureka\GBS\CI\QDOAS_output\u1_2010.mat');
%load('H:\work\Eureka\GBS\CI\2011\UTGBS\u1_2011.mat');
%load('H:\work\Eureka\GBS\CI\2014\UTGBS\u1_2014.mat');
%year = '2014';
%plot_path = 'H:\work\Eureka\GBS\CI\2014\UTGBS\plots\RTM_CI550';

mkdir(plot_path);
cd(plot_path);
size_fig = 1/2;
%save_fig = 0;

plot_CI = 1; % 1 = yes, 0 = no
plot_O4 = 0;
plot_flux = 0;

use_offset_CI = 0; % following Wagner AMT 2016, CI_cal = CI_meas * beta_CI
beta_CI = 0.9;

%% SZA filters %%%
%TF = (data.SZA > 86) & (data.SZA < 91);
%TF = data.SZA < 80;
%data(~TF,:) = [];

%% 1) plot flux
if plot_flux == 1
    nm = data.Properties.VariableNames;
    TF_flux = strfind(nm(1,:),'Fluxes');
    N_fluxes = sum(cell2mat(TF_flux));
    N_QDOASfile = size(data);
    for i = 1:1:N_QDOASfile(2)
        if cell2mat(TF_flux(i)) == 1
            TF_flux_index(i) = 1;
        else
            TF_flux_index(i) = 0;
        end
    end
    TF_flux_index = logical(TF_flux_index);
    fluxes = data(:,TF_flux_index);
    
    %% flux time series
    figure; hold all;
    for i =1:1:N_fluxes
        plot(data.Fractionalday,fluxes{:,i},'.');
    end
    nm(:,~TF_flux_index)=[];
    legend(nm);
    xlabel(['Day of year ' year]);
    ylabel('Flux');
    print_setting(size_fig,save_fig,['flux time series']);
    
    %% flux vs SZA
    figure;hold all;
    for i =1:1:N_fluxes
        plot(data.SZA,fluxes{:,i},'.');
    end
    legend(nm);
    xlabel('SZA');
    ylabel('Flux');
    print_setting(size_fig,save_fig,['flux vs SZA']);
    
    %% CI histogram
    figure;hold all;
    for i =1:1:N_fluxes
        h_hist = histogram(fluxes{:,i});
        bin = 1e3;
        h_hist.Normalization = 'probability';
        h_hist.BinWidth = bin;
    end
    legend(nm);
    xlabel('Flux');
    ylabel('f');
    print_setting(size_fig,save_fig,['flux histogram']);
      
end

%% 2) plot CI
if plot_CI == 1
        
    %CI_1 = data.Fluxes360./data.Fluxes385;
    %CI_1 = data.Fluxes550./data.Fluxes440;
    %CI_1 = data.Fluxes440./data.Fluxes550;
    CI_1 = data.Fluxes450./data.Fluxes550;
    CI_2 = data.Fluxes360./data.Fluxes550;
    CI_3 = data.Fluxes405./data.Fluxes550;
    %CI_4 = data.Fluxes425./data.Fluxes490;
    CI_4 = data.Fluxes490./data.Fluxes550;   
    if use_offset_CI == 1
        CI_1 = CI_1.*beta_CI;
        CI_2 = CI_2.*beta_CI;
        CI_3 = CI_3.*beta_CI;
        CI_4 = CI_4.*beta_CI;
    end
    

    
    %% CI time series
    figure;hold all;
    plot(data.Fractionalday,CI_1,'.');
    plot(data.Fractionalday,CI_2,'.');
    plot(data.Fractionalday,CI_3,'.');
    plot(data.Fractionalday,CI_4,'.');
    %legend('360/385 nm','360/550 nm','405/550 nm','425/490 nm');
    %legend('440/550 nm','360/550 nm','405/550 nm','490/550 nm');
    legend('450/550 nm','360/550 nm','405/550 nm','490/550 nm');
    xlabel(['Day of year ' year]);
    ylabel('Colour Index');
    ylim([0 4]);
    print_setting(size_fig,save_fig,['CI time series']);
    
    %% CI histogram
    figure;hold all;
    h1 = histogram(CI_1);
    h2 = histogram(CI_2);
    h3 = histogram(CI_3);
    h4 = histogram(CI_4);
    bin = 0.05;
    h1.Normalization = 'probability';
    h1.BinWidth = bin;
    h2.Normalization = 'probability';
    h2.BinWidth = bin;
    h3.Normalization = 'probability';
    h3.BinWidth = bin;
    h4.Normalization = 'probability';
    h4.BinWidth = bin;
    %legend('360/385 nm','360/550 nm','405/550 nm','425/490 nm');
    %legend('440/550 nm','360/550 nm','405/550 nm','490/550 nm');
    legend('450/550 nm','360/550 nm','405/550 nm','490/550 nm');
    xlabel('Colour Index');
    ylabel('f');
    xlim([0 4]);
    print_setting(size_fig,save_fig,['CI histogram']);
    
    
    %% CI cloudy screen
    figure;hold all;
    LT = -5;
    LTC = data.Fractionaltime + LT;
    TF = LTC <=0;
    LTC(TF) = LTC(TF) +24;
    plot(LTC,CI_1,'.');
    plot(LTC,CI_2,'.');
    plot(LTC,CI_3,'.');
    plot(LTC,CI_4,'.');
    %legend('360/385 nm','360/550 nm','405/550 nm','425/490 nm');
    %legend('440/550 nm','360/550 nm','405/550 nm','490/550 nm');
    legend('450/550 nm','360/550 nm','405/550 nm','490/550 nm');
    xlabel('Fractional time (local time)');
    ylabel('Colour Index');
    ylim([0 4]);
    print_setting(size_fig,save_fig,['CI cloudy screen']);
    
    %% CI cloudy screen SZA
    figure;hold all;
    plot(data.SZA,CI_1,'.');
    plot(data.SZA,CI_2,'.');
    plot(data.SZA,CI_3,'.');
    plot(data.SZA,CI_4,'.');
    %legend('360/385 nm','360/550 nm','405/550 nm','425/490 nm');
    %legend('440/550 nm','360/550 nm','405/550 nm','490/550 nm');
    legend('450/550 nm','360/550 nm','405/550 nm','490/550 nm');
    xlabel('SZA');
    ylabel('Colour Index');
    ylim([0 4]);
    print_setting(size_fig,save_fig,['CI cloudy screen SZA']);
    
    %% CI cloudy screen SZA density plot
    figure;hold all;
    dscatter(data.SZA,CI_1);
    dscatter(data.SZA,CI_2);
    dscatter(data.SZA,CI_3);
    dscatter(data.SZA,CI_4);
    %legend('360/385 nm','360/550 nm','405/550 nm','425/490 nm');
    %legend('440/550 nm','360/550 nm','405/550 nm','490/550 nm');
    legend('450/550 nm','360/550 nm','405/550 nm','490/550 nm');
    xlabel('SZA');
    ylabel('Colour Index');
    ylim([0 4]);
    print_setting(size_fig,save_fig,['CI cloudy screen SZA density']);
    
    %% CI cloudy screen SZA date plot
    figure;hold all;
    sz = 10;
    month = str2num(datestr(datenum(data.DateDDMMYYYY),'mm'));
    scatter(data.SZA,CI_1,sz,month,'filled');
    scatter(data.SZA,CI_2,sz,month,'filled');
    scatter(data.SZA,CI_3,sz,month,'filled');
    scatter(data.SZA,CI_4,sz,month,'filled');
    %gscatter(data.SZA,CI_1,month);
    %gscatter(data.SZA,CI_2,month);
    %gscatter(data.SZA,CI_3,month);
    %gscatter(data.SZA,CI_4,month);

    %legend('360/385 nm','360/550 nm','405/550 nm','425/490 nm');
    %legend('440/550 nm','360/550 nm','405/550 nm','490/550 nm');
    xlabel('SZA');
    ylabel('Colour Index');
    ylim([0 4]);
    print_setting(size_fig,save_fig,['CI cloudy screen SZA Fractionalday']);
    
    data.CI360550 = CI_2;
    data.CI405550 = CI_3;
    %data.CI440550 = CI_1;
    data.CI450550 = CI_1;
    data.CI490550 = CI_4;
end
%% polynomial fitting
N = size(data);
do_poly = 0; % check the poly of the O4 fitting, 0 = no, 1 = yes
if do_poly == 1
    figure;hold all;
    for i =1:1:N
        x1 = data.O4_VIS_293_a203SlColx1(i);
        x2 = data.O4_VIS_293_a203SlColx2(i);
        x3 = data.O4_VIS_293_a203SlColx3(i);
        x4 = data.O4_VIS_293_a203SlColx4(i);
        x5 = data.O4_VIS_293_a203SlColx5(i);
        x0 = data.O4_VIS_293_a203SlColx0(i);
        p = [x5, x4, x3, x2, x1, x0];
        
        x_test = 425:1:490;
        y_poly = polyval(p,x_test);
        
        plot(x_test,y_poly);
        CI_poly(i) = y_poly(1)/y_poly(66);
    end
    figure;plot(data.Fractionalday,CI_poly,'.');
    hist(CI_poly);
end


%% 3) plot O4
if plot_O4 == 1
    %% O4 time series
    figure; hold all;
    plot(data.Fractionalday,data.O4_VIS_293_a203SlColo4,'.');
    plot(data.Fractionalday,data.O4_VIS_203_a293SlColo4,'.');
    plot(data.Fractionalday,data.O4_VIS_293SlColo4,'.');
    plot(data.Fractionalday,data.O4_VIS_203SlColo4,'.');
    legend('O_4 293K a203K','O_4 203K a293K','O_4 293K','O_4 203K');
    xlabel(['Day of year ' year]);
    ylabel('O_4 dSCD [molec^2/cm^5]');
    print_setting(size_fig,save_fig,['O4 time series']);
    
    %% O4 histogram
    figure; hold all;
    h1 = histogram(data.O4_VIS_293_a203SlColo4);
    h2 = histogram(data.O4_VIS_203_a293SlColo4);
    h3 = histogram(data.O4_VIS_293SlColo4);
    h4 = histogram(data.O4_VIS_203SlColo4);
    bin = 1e42;
    h1.Normalization = 'probability';
    h1.BinWidth = bin;
    h2.Normalization = 'probability';
    h2.BinWidth = bin;
    h3.Normalization = 'probability';
    h3.BinWidth = bin;
    h4.Normalization = 'probability';
    h4.BinWidth = bin;
    legend('O_4 293K a203K','O_4 203K a293K','O_4 293K','O_4 203K');
    xlabel('O_4 dSCD [molec^2/cm^5]');
    ylabel('f');
    print_setting(size_fig,save_fig,['O4 histogram']);
    
    %% O4 difference time series
    figure;hold all;
    for i = 1:1:N(1)
        %mean_o4(i) = (data.O4_VIS_293SlColo4(i) + data.O4_VIS_203SlColo4(i))/2;
        mean_o4(i) = data.O4_VIS_293SlColo4(i);
        delta_1(i) = (data.O4_VIS_293_a203SlColo4(i) - mean_o4(i))/mean_o4(i)*100;
        delta_2(i) = (data.O4_VIS_203_a293SlColo4(i) - mean_o4(i))/mean_o4(i)*100;
        delta_3(i) = (data.O4_VIS_293SlColo4(i) - mean_o4(i))/mean_o4(i)*100;
        delta_4(i) = (data.O4_VIS_203SlColo4(i) - mean_o4(i))/mean_o4(i)*100;
    end
    plot(data.Fractionalday,delta_1,'.');
    plot(data.Fractionalday,delta_2,'.');
    plot(data.Fractionalday,delta_3,'.');
    plot(data.Fractionalday,delta_4,'.');
    ylim([-50 50]);
    legend('O_4 293K a203K','O_4 203K a293K','O_4 293K','O_4 203K');
    xlabel('O_4 dSCD [molec^2/cm^5]');
    ylabel('% difference (vs. O_4 293K)');
    print_setting(size_fig,save_fig,['O4 difference time series']);
    
    %% O4 difference histogram
    figure; hold all;
    TF = abs(delta_1) > 100;delta_1(TF) = [];
    TF = abs(delta_2) > 100;delta_2(TF) = [];
    TF = abs(delta_3) > 100;delta_3(TF) = [];
    TF = abs(delta_4) > 100;delta_4(TF) = [];
    
    h1 = histogram(delta_1);
    h2 = histogram(delta_2);
    h3 = histogram(delta_3);
    h4 = histogram(delta_4);
    bin = 2;
    h1.Normalization = 'probability';
    h1.BinWidth = bin;
    h2.Normalization = 'probability';
    h2.BinWidth = bin;
    h3.Normalization = 'probability';
    h3.BinWidth = bin;
    h4.Normalization = 'probability';
    h4.BinWidth = bin;
    xlim([-50 50]);
    legend('O_4 293K a203K','O_4 203K a293K','O_4 293K','O_4 203K');
    xlabel('% difference (vs. O_4 293K)');
    ylabel('f');
    print_setting(size_fig,save_fig,['O4 difference histogram']);
end