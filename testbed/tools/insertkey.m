function insertkey(fname,key_mod,value_mod)
% replace 
fi=fopen(fname)
l=[];i=0;found=0;
while 1
    l = fgetl(fi);
    if ~ischar(l), break, end
    key=deblank(strtok(l,'='));
    %    value=deblank(fliplr(strtok(fliplr(l),'=')))
    if strcmp(key,key_mod);
        l=[key,'   = ',num2str(value_mod)];
        found=1;
    end
    i=i+1;
    s{i}=l;
end
fclose(fi);
fo=fopen(fname,'wt');
for i=1:length(s)
    fprintf(fo,[s{i},'\n']);
end
if ~found
    l=[key_mod,'   = ',num2str(value_mod)]
    fprintf(fo,[l,'\n']);
end
fclose(fo)
