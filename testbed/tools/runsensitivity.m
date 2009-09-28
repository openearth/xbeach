testbeddir='d:\data\dano\xbeach\testbed'
exe='D:\data\dano\xbeach\VS2008\Release\xbeach.exe>out'
cd([testbeddir '\tools']);
addpath(pwd)
fi=fopen('sensitivity.txt')
while 1
    l = fgetl(fi)
    if isempty(l)|strcmp(l,'stop'), break, end
    runid=deblank(l);
        %% Read parameters to change
        l=fgetl(fi);
        numpar=str2num(l);
        for ip=1:numpar;
            l=fgetl(fi);
            key{ip}=deblank(strtok(l,'='));
            value{ip}=deblank(fliplr(strtok(fliplr(l),'=')));
        end
    tag=''
    cd ..
    testbankdir=pwd;
    if filesep=='/'
        release=''
    else
        release='Release\'
    end
    cases ={ ...
        ['Boers_1C']                     ... %1
        ['CarrierGreenspan'],            ... %2
        ['Delilah_199010131000'],        ... %3
        ['Deltaflume_T3'],               ... %4
        ['Deltaflume2006_T01'],          ... %5
        ['Deltaflume2006_T04'],          ... %6
        ['Ducktest_coarse'],             ... %7
        ['Hard_structures_flat_bottom'], ... %8
        ['Hard_structures_rif'],         ... %9
        ['Humptest_basis'],              ... %10
        ['Lip11D_2E'],                   ... %11
        ['Sbeach_tests_oceancity'],      ... %12
        ['Yu&Slinn_WCI'],                ... %13
        ['Zelt_Case1'],                  ... %14
        ['Zwin_T01']                     ... %15
        }
    todo= [1:15]
    for j=1:length(todo)
        i=todo(j);
        rundir=[testbankdir,filesep,runid,filesep,cases{i}];
        mkdir (rundir);
        cd (rundir);
        copyfile (['..',filesep,'..',filesep,'input', ...
            filesep,cases{i},filesep,'*.*'],rundir,'f');
        %% Change params.txt
        for ip=1:numpar;
            insertkey('params.txt',key{ip},value{ip});
        end
        %
        eval(['!',exe])
        % plotdata
        %     post
        %     showbank
    end
    cd ([testbankdir,filesep,'tools']);
end
fclose(fi)