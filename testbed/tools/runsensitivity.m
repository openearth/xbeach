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
        ['Delilah_199010131000'],        ... %2
        ['Deltaflume2006_T01'],          ... %3
        ['Deltaflume2006_T04'],          ... %4
        ['Ducktest_coarse'],             ... %5
        ['Hard_structures_flat_bottom'], ... %6
        ['Hard_structures_rif'],         ... %7
        ['Humptest_basis'],              ... %8
        ['Lip11D_2E'],                   ... %9
        ['Sbeach_tests_oceancity'],      ... %10
        ['Yu&Slinn_WCI'],                ... %11
        ['Zelt_Case1'],                  ... %12
        ['Zwin_T01'],                    ... %13
        ['CarrierGreenspan']             ... %14
        }
    todo= [1:14]
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
        plotdata
        %     post
        %     showbank
    end
    cd ([testbankdir,filesep,'tools']);
end
fclose(fi)