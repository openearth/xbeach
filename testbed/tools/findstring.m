function [lines]=findstring(filename,str)
choice=2;
%=========================================================================
%[FileName,PathName] = uigetfile('*.txt','Select any text file');
%y= [PathName,FileName];
y=filename;
%=========================================================================
%str=input('ENTER YOUR STRING(in single inverted commas)=');
len_str=length(str);
%=========================================================================
%choice=input('TYPE OF SEARCH:CASE SENSITIVE(Enter 1) or INSENSITIVE(Enter 2)=');
%=========================================================================
fid=fopen(y);
%=========================================================================
count=0;
got1=[];
got2=[];
while 1
    tline = fgetl(fid);
      if ~ischar(tline), break, end
        count=count+1;
%=========================================================================
count5=0;
count3=[];
for(t=1:length(tline)-length(str)+1)
count4=0;
for(count2=1:length(str))
%=========================================================================  
    if(choice==2)
    if(str(count2)==tline(t+count2-1)||str(count2)==upper(tline(t+count2-1))||str(count2)==lower(tline(t+count2-1)))
        count4=count4+1;
    end
end
%=========================================================================  
    if(choice==1)
    if(str(count2)==tline(t+count2-1))
        count4=count4+1;
    end
end
%=========================================================================  
end         
count3=[count3 count4];
if(count4==len_str)
    count5=count5+1;
end
end
if(count5>0)
    got1=[got1 count];
    got2=[got2 count5];
end
%=========================================================================          
end
%=========================================================================  
% if(length(got1)==0)
%     disp('SORRY,NOT FOUND');
% end
% % %=========================================================================  
% if(length(got1)~=0)
% disp('LINE NUMBER    NO. OF OCCURENCES')
% max_len=length(int2str(max(got1)));
% for(i=1:length(got1))
%     char1=int2str(got1(i));
%     add_char=max_len-length(char1);
%     if(add_char>0)
%     for(j=1:add_char)
%         char1=[char1 ' '];
%     end
%     end
%     char2=int2str(got2(i));
%     char3=[char1 '                   ' char2];
%     disp(char3)
% end
% end
% % %========================================================================
lines=got1;