function plot_timeserise(VCD,VCD_CF,brewer_all,plot_path,save_fig)
DU = 2.6870e+16;
brewer_all.datetime = datetime(datevec(brewer_all.UTC));
brewer_all.dayofyear = day(brewer_all.datetime,'dayofyear');
brewer_all.fd = brewer_all.dayofyear + brewer_all.datetime.Hour/24 +brewer_all.datetime.Minute/24/60;

year = str2num(VCD.year(1,:));

TF_year = brewer_all.datetime.Year == year;
brewer = brewer_all(TF_year,:);

TF_DS = strcmp(brewer.ObsCode,'DS');
TF_ZS = strcmp(brewer.ObsCode,'ZS');

brewer_ds = brewer(TF_DS,:);
brewer_zs = brewer(TF_ZS,:);

figure;hold all;
plot(brewer_ds.fd,brewer_ds.ColumnO3,'s');
plot(brewer_zs.fd,brewer_zs.ColumnO3,'s');
plot(VCD.fd,VCD.mean_vcd./DU,'.-');
plot(VCD_CF.fd,VCD_CF.mean_vcd./DU,'.');
plot(VCD.fd,VCD.langley_vcd./DU,'x--');
plot(VCD_CF.fd,VCD_CF.langley_vcd./DU,'x');
legend('Brewer DS','Brewer ZS','GBS', 'GBS-CF', 'GBS (Langley)', 'GBS-CF (Langley)');
ylabel('Ozone VCD [DU]');
xlabel('day of the year');
cd(plot_path);
print_setting(1/2,save_fig,['Timeserise_all']);
