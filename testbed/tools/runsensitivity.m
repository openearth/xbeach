addpath(pwd)
fi=fopen('sensitivity.txt')
while 1
    l = fgetl(fi)
    if isempty(l), break, end
    runid=deblank(l);
    tag=''
    cd ..
    testbankdir=pwd;
    if filesep=='/'
        release=''
    else
        release='Release\'
    end
    exe=[pwd,filesep,'tools',filesep,'run.bat']
    cases ={ ...
        ['Boers_1C']                     ... %1
        ['Delilah_199010131000'],        ... %2
        ['Deltaflume2006_T01'],          ... %3
        ['Duck_base'],                   ... %4
        ['Ducktest_coarse'],             ... %5
        ['Hard_structures_flat_bottom'], ... %6
        ['Hard_structures_rif'],         ... %7
        ['Humptest_basis'],              ... %8
        ['Lip11D_2E'],                   ... %9
        ['Sbeach_tests_oceancity'],      ... %10
        ['Yu&Slinn_WCI'],                ... %11
        ['Zelt_Case1'],                  ... %12
        ['Zwin_T01'],                    ... %13
        }
    todo= [13]
    for j=1:length(todo)
        i=todo(j);
        rundir=[testbankdir,filesep,runid,filesep,cases{i}];
        mkdir (rundir);
        cd (rundir);
        copyfile (['..',filesep,'..',filesep,'input', ...
            filesep,cases{i},filesep,'*.*'],rundir,'f');
        %% Change params.txt
        l=fgetl(fi);
        numpar=str2num(l);
        for i=1:numpar;
            l=fgetl(fi);
            key=deblank(strtok(l,'='));
            value=deblank(fliplr(strtok(fliplr(l),'=')));
            insertkey('params.txt',key,value);
        end
        %
        eval(['!',exe])
        %plotdata
        %     post
        %     showbank
    end
    cd ([testbankdir,filesep,'tools']);
end
fclose(fi)