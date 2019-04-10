function combin_box_plots()
% this function can combin all subplots made by "plot_all_boxplots"
%figure_path = 'C:\Users\ZhaoX\Documents\paper\CF\fig\box_plots\spring\';
%figure_path = 'C:\Users\ZhaoX\Documents\paper\CF\fig\box_plots\spring_remake\';
figure_path = 'C:\Users\ZhaoX\Documents\paper\CF\fig\box_plots\no_2017_spring\';

uiopen([ figure_path 'Brewer_vs_GBSandSAOZ_boxplots.fig'],1);
h1 = gca;
uiopen([ figure_path 'MERRA-2_vs_GBSandSAOZ_boxplots.fig'],1);
h2 = gca;
uiopen([ figure_path 'Brewer_vs_GBSandSAOZ_coincidentnumber.fig'],1);
h3 = gca;
uiopen([ figure_path 'MERRA-2_vs_GBSandSAOZ_coincidentnumber.fig'],1);
h4 = gca;

figure;hold all;
h11 = subplot(2,2,1);
copyobj(allchild(h1),h11);
xlim([0.5 6]);
ylim([-10 10]);
ylabel(['TCO % difference (X - Brewer)/Brewer [%]']);
xticklabels({'','','','',''});

h12 = subplot(2,2,2);
copyobj(allchild(h2),h12);
xlim([0.5 6]);
ylim([-10 10]);
ylabel(['TCO % difference (X - MERRA-2)/MERRA-2 [%]']);
xticklabels({'','','','',''});

h21 = subplot(2,2,3);
copyobj(allchild(h3),h21);
xlim([0.5 6]);
ylim([0 600]);
xticklabels({'Clear','Cloudy','Ice Crystals','Mainly Clear','Mostly Cloudy'});
rotateXLabels(gca,45);
legend('UT-GBS','UT-GBS_C_S','SAOZ','SAOZ_C_S','SAOZ_V_3');
ylabel(['Number of coincident measurements with Brewer']);

h22 = subplot(2,2,4);
copyobj(allchild(h4),h22);
xlim([0.5 6]);
ylim([0 600]);
ylabel(['Number of coincident measurements with MERRA-2']);
xticklabels({'Clear','Cloudy','Ice Crystals','Mainly Clear','Mostly Cloudy'});
rotateXLabels(gca,45);