function input_table = read_input_file()

input_table = table;

cd ..
cd('CF_package_local');
fid = fopen('CF_input_file.txt','r');
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
    status = copyfile('CF_input_file.txt', [input_table.plot_path 'CF_input_file_archive.txt'], 'f');
    if status == 1
        disp('input file has been archived');
    else
        disp('Warning: input file not archived, pls check "read_input_file.m" !');
    end
catch
    disp('Warning: input file not archived, pls check "read_input_file.m" !');
end