function plot_all_SCIATRAN2_CI_simulations()
addpath('E:\F\Work\MatlabCode\');
database_path = 'E:\H\work\Eureka\GBS\CI\matlab\CF_package_v2\database\';

simulations = ls(database_path);
simulations(1:2,:)=[];
simulations(end,:)=[];
% simulation_types = {'clear_sky','aerosold_lowtran_50km_albedo_0dot06',...
%     'aerosold_lowtran_23km_albedo_0dot06',...
%     ''};

%N = size(simulation_types);
N = size(simulations);
h1 = figure;hold all;
for i = 1:N(1)
    %simulation_type = string(simulation_types(i));
    simulation_type = strtrim(simulations(i,:));
    %data_path = char([database_path + simulation_type + '\']);
    data_path = [database_path simulation_type  '\'];
    plot_path = 'E:\H\work\Eureka\GBS\CI\matlab\plots\'
    try
        model = SCIATRAN2_CI_450_550(data_path,simulation_type,plot_path);
        h2 = gcf;
        close(h2);
        figure(h1);
        plot(model.SZA,model.CI);
        model = [];
    catch
        test;
    end
end
figure(h1);
legend(simulations)