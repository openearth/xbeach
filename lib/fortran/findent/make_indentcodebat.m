function make_indentcodebat

%% collect all F90 files in src directories
names = {};
dirs = dir('..\..\..\src\makeincludes\*.f90');
for i=1:length(dirs)
    names{end+1} = ['..' filesep '..' filesep '..' filesep 'src' filesep 'makeincludes' filesep dirs(i).name];
end
dirs = dir('..\..\..\src\xbeach\*.f90');
for i=1:length(dirs)
    names{end+1} = ['..' filesep '..' filesep '..' filesep 'src' filesep 'xbeach' filesep dirs(i).name];
end
dirs = dir('..\..\..\src\xbeachlibrary\*.f90');
for i=1:length(dirs)
    names{end+1} = ['..' filesep '..' filesep '..' filesep 'src' filesep 'xbeachlibrary' filesep dirs(i).name];
end

%% generate batch script
fid = fopen('indentcode.bat','wt');
for i=1:length(names)
    fprintf(fid,'%s\n',['findent -o3 <' names{i} '> ' [names{i} '.new']]);
    fprintf(fid,'%s\n',['move /Y  ' [names{i} '.new'] ' '  names{i} ]);
end
fclose(fid);