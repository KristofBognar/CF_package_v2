function input_table = read_input_file_smalltools()

input_table = table;

cd ..
cd('CF_package_local');
fid = fopen('CF_input_file_smalltools.txt','r');
%tline = fgetl(fid);

while ~feof(fid)
    tline = fgetl(fid);
    try
        eval(tline);
    catch
        disp('Warning: one required input might not found! See the following line:')
        disp(tline);
    end
end

try
    mkdir(input_table.plot_path);
    status = copyfile('CF_input_file_smalltools.txt', [input_table.concat_data_plot_path 'CF_input_file_smalltools_archive.txt'], 'f');
    if status == 1
        disp('input file has been archived');
    else
        disp('Warning: input file not archived, pls check "read_input_file_smalltools.m" !');
    end
catch
    disp('Warning: input file not archived, pls check "read_input_file_smalltools.m" !');
end