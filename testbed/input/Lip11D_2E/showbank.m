% showbank.m
%
% Routine to visualize numerical simulations of LIP experiments
% against actual flume measurements
%
% Starting point: file showbank.m is located in directory worktools
%                 showbank needs file named 'plotdef.inp' as input
%

function showbank ()

[strs,plots] = readPlotDef('plotdef.inp');

% Define some strings for plotting purposes and locate data directory

PN = fliplr(pwd);
[runid,R] = strtok(PN,'\'); runid = fliplr(runid);
[testid,R] = strtok(R,'\'); testid = fliplr(testid);
%[dataid,R] = strtok(R,'\'); dataid = fliplr(dataid);
datadir=['..\..\data\',runid,'\']

% Visualize individual results 
% (per figure & than per parameter)

for idp = 1:length(plots)
   figure(idp)
   clusterData = [];
   count = 0;
   for i = 1:size(plots(idp).params)
      switch lower(plots(idp).type)
         case {'fx','fz','ft'}
            varname=deblank(plots(idp).params(i,:));
            fname = [varname '.tek']
            comp = readTekFirstBlock(fname)
            meas = readTekFirstBlock([datadir fname])
            if ( ~isempty(comp) & ~isempty(meas) )
              [r2,slope,eps,bss]=compskill3(comp,meas,runid,testid,varname);
              command=['!type ',varname,'.err>>..\..\report\',varname,'.err']
               eval(command)
               count = count+1;
               subplot (plots(idp).size(1),plots(idp).size(2),count);
               switch deblank(plots(idp).type)
                  case 'FX', PO = [1 2];
                  case 'FZ', PO = [2 1];
               end
               ph = plot(comp(:,PO(1)),comp(:,PO(2)),'b',meas(:,PO(1)),meas(:,PO(2)),'ro');
               set(gca,'fontsize',6);
               if plots(idp).type == 'FX'
                  set(gca,'xlim',[min(comp(:,PO(1))) max(comp(:,PO(1)))] );
               end
               set(ph,'markersize',3,'linewidth',3);  % Default dot size to 3
               if plots(idp).type == 'FT'
                  set(ph(2),'markersize',1,'linewidth',1);
               end
               xl=get(gca,'xlim');xtext=xl(1)+.1*(xl(2)-xl(1))
               yl=get(gca,'ylim');ytext=yl(1)+.1*(yl(2)-yl(1))
               text(xtext,ytext,['r^2=',round2(r2,2),' m=',round2(slope,2),' \epsilon=',round2(eps,2),' bss=',round2(bss,2)],'fontsize',8)
            end
            
         case 'fs'
            fname = [deblank(plots(idp).params(i,:)) '.sct'];            
            data = readTekFirstBlock(fname);
            if ~isempty(data)
               clusterData = [clusterData; [data i*ones(size(data,1),1)]];
               count = count+1;
               if isempty(plots(idp).cluster)
                  subplot (plots(idp).size(1),plots(idp).size(2),count);
                  mx = max(max(data)); mn = min(min(data));
                  ph = plot([mn mx],[mn mx],'b',data(:,1),data(:,2),'ro');
                  pm = 0.05*(mx-mn);
                  set(gca,'xlim',[mn-pm mx+pm],'ylim',[mn-pm mx+pm]);
                  lh = xlabel('measured'); set(lh,'fontsize',6)
                  lh = ylabel('computed'); set(lh,'fontsize',6)
                  set(gca,'fontsize',6), axis('square')
               end
            end
               
         end
         
         % Set ranges
         if ~isempty(plots(idp).xrange)
            set(gca,'xlim',plots(idp).xrange);
         end
         if ~isempty(plots(idp).yrange)
            set(gca,'ylim',plots(idp).yrange);
         end
         if ~isempty(plots(idp).markersize)
            set(ph(2),'linewidth',plots(idp).markersize);
         end
         set(gca,'ygrid','on') 
         th = title(deblank(plots(idp).params(i,:))); set(th,'fontsize',8);
         
   end
   
   if ~isempty(plots(idp).cluster)
      
      for i = 1:size(plots(idp).cluster,2)
         cluslimits = [0 cumsum(plots(idp).cluster)];
         id = find(clusterData(:,3)>cluslimits(i) & clusterData(:,3)<=cluslimits(i+1));
         subplot(plots(idp).size(1),plots(idp).size(2),i);
         mx = max(max(clusterData(id,[1 2]))); mn = min(min(clusterData(id,[1 2])));
         plot([mn mx],[mn mx],'b',clusterData(id,1),clusterData(id,2),'ro');
         pm = 0.05*(mx-mn);
         set(gca,'xlim',[mn-pm mx+pm],'ylim',[mn-pm mx+pm]);
         lh = xlabel('measured'); set(lh,'fontsize',6)
         lh = ylabel('computed'); set(lh,'fontsize',6)       
         set(gca,'fontsize',6), axis('square'), grid on
         tistr = ['Clustering ' plots(idp).params(cluslimits(i)+1,:) ... 
                  ' to ' plots(idp).params(cluslimits(i+1),:) ];
            th = title(tistr);
            set (th, 'fontsize',8);
            
      end
         
   end
   
   orient tall
   set(gcf,'papertype','a4')
end   

S1a = ['Test ' testid];
S2 = strs.STR2;
S3 = ['Run ' runid];
S4 = strs.STR4;
S5 = strs.STR5;
S6 = strs.STR6;

for idp = 1:length(plots)
%   figure(idp)
%   S1b = plots(idp).mainstr;
%   S6 = plots(idp).figstr;
%   md_paper('portrait',strvcat(S1a,S1b),S2,S3,S4,S5,S6);
%   shh=get(0,'showhiddenhandles');
%   set(0,'showhiddenhandles','on');
   pname = ['..\..\report\',testid '_' runid '_fig' num2str(idp) '.jpg'];
   eval(['print -djpeg ' pname]);
end

% winmngr
%exit


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function readPlotDef                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [strs, plots] = readPlotDef(fname)

% Initialize variables
strs = struct('STR1',[],'STR2',[],'STR3',[], ...
              'STR4',[],'STR5',[],'STR6',[]);
           
% Read plot settings from file 'plotdef.inp'
fid = fopen(fname,'r'); 
str = fgetl(fid);
str = fgetl(fid);
while (strcmp(deblank(str),'END-STRINGS') == 0)
   strs = setfield(strs,str(1:4),str(6:end));
   str = fgetl(fid);   
end

% plots = [];
idp = 1;
while ~feof(fid)
   plots(idp) = struct('type',[],'size',[],'mainstr',[], 'figstr',[], ...
              'params',[],'xrange',[],'yrange',[],'markersize',[],'cluster',[]);
   [dum,str] = strtok(fgetl(fid));
   [plots(idp).type,str] = strtok(str);  % isolate type
   [sz1,str] = strtok(str);
   [sz2,str] = strtok(str);
   [plots(idp).size] = str2num([sz1 ' ' sz2]);
   while ~isempty(deblank(str))
      [fieldname,str] = strtok(str);
      switch lower(fieldname)
      case {'xrange','yrange'}
         [sz1,str] = strtok(str);
         [sz2,str] = strtok(str);
         eval(['plots(idp).' lower(fieldname) ' = [str2num(sz1) str2num(sz2)];']) ;
      case 'markersize'
         [sz1,str] = strtok(str);
         eval(['plots(idp).' lower(fieldname) ' =  str2num(sz1);'])
      case 'cluster'
         [sz1,str] = strtok(str);
         eval(['plots(idp).' lower(fieldname) ' = [plots(idp).' lower(fieldname) ' str2num(sz1)];'])
      end
   end
   [dummy,plots(idp).mainstr] = strtok(fgetl(fid));
   [dummy,plots(idp).figstr] = strtok(fgetl(fid));
   str = fgetl(fid);
   plots(idp).params = [];
   while (strcmp(deblank(str),'END-PLOT') == 0)
      plots(idp).params = strvcat(plots(idp).params,str);
      str = fgetl(fid);   
   end
   idp = idp + 1;   
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function readTekFirstBlock             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data] = readTekFirstBlock(fname)

% Routine to read first block of data from Tekal file

if exist(fname)
   fid = tekal('open',fname);
   data = tekal('read',fid,1);
%    fid = fopen(fname);
%    dummy = fgetl(fid);
%    BLsize = fscanf(fid,'%i',[1 2]);
%    % ad-hoc solution to account for lines with NaN's
%    dummy = fgetl(fid);
%    data = '';
%    for p = 1:BLsize(1)
%      str = fgetl(fid);
%      data = strvcat(data, str);
%    end   
%    data = str2num(data);   
      
   % data = fscanf(fid, '%s',fliplr(BLsize));
   % data = data';
   % data(data(:,2)==0)=nan;
   % fclose(fid);
else
   data = [];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function md_paper                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hBorder=md_paper(cmd,varargin);
% MD_PAPER adds border to plot
%
%   md_paper(PaperType);
%   PaperType can be either 'portrait' or 'landscape'
%   adds "WL | Delft Hydraulics (date and time)" to a page.
%
%   md_paper(PaperType,'String');
%   adds "String (date and time)" to a page.
%
%   md_paper(PaperType,'String1','String2',...)
%   adds a WL | Delft Hydraulics border (currently only A4)
%   right click on the text to edit the texts.
%
%   hBorder=paper('no edit',PaperType,'String1','String2',...)
%   right click editing disabled; for gui environment to be used with
%   paper('edit',hBorder) for editing.
%
%   Use following routines in association with md_paper.m:
%       md_print: to print figures to file or printer
%       winmngr: to close windows. Notice argument option winmngr hvis
%                 for handle visibility

if nargin==0,
  cmd='portrait';
end;

NoEdit=0;
if strcmp(cmd,'no edit'),
  NoEdit=1;
  cmd=varargin{1};
  varargin=varargin(2:end);
end;

switch cmd,
case {'portrait','landscape'},
  hTempBorder=Local_createborder(NoEdit,cmd,varargin{:});
  if nargout>0,
    hBorder=hTempBorder;
  end;
case {'apply','done'},
  Fig=gcbf;
  gcba=get(Fig,'userdata');
  if ~ishandle(gcba),
    delete(gcbf);
    return;
  end;
  hpts=findobj(gcba,'type','text');
  for i=1:7,
    hplottext=findobj(hpts,'flat','tag',sprintf('plottext%i',i));
    hedittext=findobj(Fig,'tag',sprintf('Text%i',i));
    set(hplottext,'string',get(hedittext,'string'));
  end;
  leftpage=findobj(Fig,'tag','LeftPage');
  switch get(get(gcba,'parent'),'paperorientation'),
  case 'portrait',
    if get(leftpage,'value'), %=1
      set(gcba,'xdir','reverse');
    else, %=0
      set(gcba,'xdir','normal');
    end;
  case 'landscape',
    if get(leftpage,'value'), %=1
      set(gcba,'ydir','reverse');
    else, %=0
      set(gcba,'ydir','normal');
    end;
  end;
  if strcmp(cmd,'done'),
    delete(gcbf);
  end;
case 'edit',
  if (nargin==1) & (isempty(gcbf) | ~strcmp(get(gcbf,'selectiontype'),'alt')),
    return;
  elseif nargin==2,
    gcba=varargin{1};
  else,
    gcba=get(gcbo,'parent');
  end;
  
  %first check whether it already exists ...
  Fig=findobj(allchild(0),'type','figure','tag','Border manager for Matlab (c)');
  if ~isempty(Fig),
    for f=transpose(Fig),
      if isequal(get(f,'userdata'),gcba),
        set(f,'visible','on');
        return;
      end;
    end;
  end;
  Fig=ui_paper;

  fig=get(gcba,'parent');
  HandleStr=[' ' num2str(fig)];
  StringStr=get(fig,'name');
  if strcmp(get(fig,'numbertitle'),'on'),
    if isempty(StringStr),
      StringStr=['Figure No.' HandleStr];
    else,
      StringStr=['Figure No.' HandleStr ':' StringStr];
    end;
  end;
  set(Fig,'name',[get(Fig,'name') ' for ' StringStr]);
  
  hpts=findobj(gcba,'type','text');
  for i=1:7,
    hplottext=findobj(hpts,'flat','tag',sprintf('plottext%i',i));
    hedittext=findobj(Fig,'tag',sprintf('Text%i',i));
    set(hedittext,'string',get(hplottext,'string'));
  end;
  set(Fig,'userdata',gcba);
  leftpage=findobj(Fig,'tag','LeftPage');
  switch get(get(gcba,'parent'),'paperorientation'),
  case 'portrait',
    set(leftpage,'value',strcmp(get(gcba,'xdir'),'reverse'));
  case 'landscape',
    set(leftpage,'value',strcmp(get(gcba,'ydir'),'reverse'));
  end;

otherwise,
  error('* Paper orientation should be ''portrait'' or ''landscape''.');
end;


function hBorder=Local_createborder(NoEdit,Orientation,varargin);

BorderColor='k';
border=1;

switch length(varargin),
case {0,1},
  border=0;
  if length(varargin), %=1
    if isempty(varargin{1}),
      plottext=DateStr;
    else,
      plottext=[varargin{1},' (',DateStr,')'];
    end;
  else,
    plottext=['WL | Delft Hydraulics (',DateStr,')'];
  end;
case 7,
  plottext=varargin;
case {2,3,4,5,6},
  plottext=cell(1,7);
  plottext(:)={''};
  plottext(1:length(varargin))=varargin;
  plottext{7}='WL | DELFT HYDRAULICS';
otherwise,
  error('* Too many parameters.');
end;

if NoEdit,
  BDFunction='';
else,
  BDFunction='md_paper edit';
end;

ax=gca;
fg=gcf;
allchld=allchild(fg);
allax=findobj(allchld,'type','axes');
set(fg,'paperunits','centimeter', ...
       'papertype','a4letter', ...
       'paperorientation',Orientation);
xmax=get(fg,'papersize');
ymax=xmax(2);
xmax=xmax(1);
set(fg,'paperposition',[0 0 xmax ymax]);
hBorder=axes('units','normalized', ...
     'position',[0 0 1 1], ...
     'tag','border', ...
     'xlimmode','manual', ...
     'ylimmode','manual', ...
     'xlim',[0 xmax], ...
     'ylim',[0 ymax], ...
     'visible','off');

extramargin2=0.5; % added for in1djt2
extramargin=0.5; % added for hk5djt1

if border,
  switch Orientation,
  case 'portrait',
    line([1 xmax-1 xmax-1 1 1],[1+extramargin2 1+extramargin2 ymax-1-extramargin ymax-1-extramargin 1+extramargin2],'parent',hBorder,'color',BorderColor,'linewidth',1.5);
    b1=0.68*(xmax-2)+1;
    line([1 xmax-1 xmax-1 1 NaN b1 b1 b1 xmax-1],[1.9+extramargin2 1.9+extramargin2 3.7+extramargin2 3.7+extramargin2 NaN 1+extramargin2 3.7+extramargin2 2.8+extramargin2 2.8+extramargin2],'parent',hBorder,'color',BorderColor,'linewidth',1.5);
    b2=(b1+xmax-1)/2;
    line([b2 b2 NaN b2 b2],[1+extramargin2 1.9+extramargin2 NaN 2.8+extramargin2 3.7+extramargin2],'parent',hBorder,'color',BorderColor,'linewidth',1.5);
    text((b1+1)/2,2.8+extramargin2,plottext(1), ...
          'horizontalalignment','center', ...
          'verticalalignment','middle', ...
          'fontname','helvetica', ...
          'tag','plottext1', ...
          'buttondownfcn',BDFunction, ...
          'color',BorderColor);
    text((b2+b1)/2,3.25+extramargin2,plottext(2), ...
          'horizontalalignment','center', ...
          'verticalalignment','middle', ...
          'fontname','helvetica', ...
          'tag','plottext2', ...
          'buttondownfcn',BDFunction, ...
          'color',BorderColor);
    text((b2+xmax-1)/2,3.25+extramargin2,plottext(3), ...
          'horizontalalignment','center', ...
          'verticalalignment','middle', ...
          'fontname','helvetica', ...
          'tag','plottext3', ...
          'buttondownfcn',BDFunction, ...
          'color',BorderColor);
    text((b1+xmax-1)/2,2.35+extramargin2,plottext(4), ...
          'horizontalalignment','center', ...
          'verticalalignment','middle', ...
          'fontname','helvetica', ...
          'tag','plottext4', ...
          'buttondownfcn',BDFunction, ...
          'color',BorderColor);
    text((b2+b1)/2,1.45+extramargin2,plottext(5), ...
          'horizontalalignment','center', ...
          'verticalalignment','middle', ...
          'fontname','helvetica', ...
          'tag','plottext5', ...
          'buttondownfcn',BDFunction, ...
          'color',BorderColor);
    text((b2+xmax-1)/2,1.45+extramargin2,plottext(6), ...
          'horizontalalignment','center', ...
          'verticalalignment','middle', ...
          'fontname','helvetica', ...
          'tag','plottext6', ...
          'buttondownfcn',BDFunction, ...
          'color',BorderColor);
    text((b1+1)/2,1.45+extramargin2,plottext(7), ...
          'horizontalalignment','center', ...
          'verticalalignment','middle', ...
          'fontname','helvetica', ...
          'fontweight','bold', ...
          'tag','plottext7', ...
          'buttondownfcn',BDFunction, ...
          'color',BorderColor);
    for i=1:length(allax),
      set(allax(i),'units','normalized');
      pos_i=get(allax(i),'position');
      n_pos_i(1)=0.055+pos_i(1)*0.89;
      n_pos_i(2)=0.130+pos_i(2)*0.83;
      n_pos_i(3)=pos_i(3)*0.89;
      n_pos_i(4)=pos_i(4)*0.83;
      set(allax(i),'position',n_pos_i);
    end;
  case 'landscape',
    line([1+extramargin2 xmax-1-extramargin xmax-1-extramargin 1+extramargin2 1+extramargin2],[1 1 ymax-1 ymax-1 1],'parent',hBorder,'color',BorderColor,'linewidth',1.5);
    b1=0.32*(ymax-2)+1;
    line([1.9+extramargin2 1.9+extramargin2 3.7+extramargin2 3.7+extramargin2 NaN 1+extramargin2 3.7+extramargin2 2.8+extramargin2 2.8+extramargin2],[1 ymax-1 ymax-1 1 NaN b1 b1 b1 1],'parent',hBorder,'color',BorderColor,'linewidth',1.5);
    b2=(b1+1)/2;
    line([1+extramargin2 1.9+extramargin2 NaN 2.8+extramargin2 3.7+extramargin2],[b2 b2 NaN b2 b2],'parent',hBorder,'color',BorderColor,'linewidth',1.5);
    text(2.8+extramargin2,(ymax-1+b1)/2,plottext(1), ...
          'horizontalalignment','center', ...
          'verticalalignment','middle', ...
          'fontname','helvetica', ...
          'tag','plottext1', ...
          'buttondownfcn',BDFunction, ...
          'color',BorderColor, ...
          'rotation',270);
    text(3.25+extramargin2,(b2+b1)/2,plottext(2), ...
          'horizontalalignment','center', ...
          'verticalalignment','middle', ...
          'fontname','helvetica', ...
          'tag','plottext2', ...
          'buttondownfcn',BDFunction, ...
          'color',BorderColor, ...
          'rotation',270);
    text(3.25+extramargin2,(b2+1)/2,plottext(3), ...
          'horizontalalignment','center', ...
          'verticalalignment','middle', ...
          'fontname','helvetica', ...
          'tag','plottext3', ...
          'buttondownfcn',BDFunction, ...
          'color',BorderColor, ...
          'rotation',270);
    text(2.35+extramargin2,(b1+1)/2,plottext(4), ...
          'horizontalalignment','center', ...
          'verticalalignment','middle', ...
          'fontname','helvetica', ...
          'tag','plottext4', ...
          'buttondownfcn',BDFunction, ...
          'color',BorderColor, ...
          'rotation',270);
    text(1.45+extramargin2,(b2+b1)/2,plottext(5), ...
          'horizontalalignment','center', ...
          'verticalalignment','middle', ...
          'fontname','helvetica', ...
          'tag','plottext5', ...
          'buttondownfcn',BDFunction, ...
          'color',BorderColor, ...
          'rotation',270);
    text(1.45+extramargin2,(b2+1)/2,plottext(6), ...
          'horizontalalignment','center', ...
          'verticalalignment','middle', ...
          'fontname','helvetica', ...
          'tag','plottext6', ...
          'buttondownfcn',BDFunction, ...
          'color',BorderColor, ...
          'rotation',270);
    text(1.45+extramargin2,(ymax-1+b1)/2,plottext(7), ...
          'horizontalalignment','center', ...
          'verticalalignment','middle', ...
          'fontname','helvetica', ...
          'fontweight','bold', ...
          'tag','plottext7', ...
          'buttondownfcn',BDFunction, ...
          'color',BorderColor, ...
          'rotation',270);
    for i=1:length(allax),
      set(allax(i),'units','normalized');
      pos_i=get(allax(i),'position');
      n_pos_i(1)=0.130+pos_i(1)*0.83;
      n_pos_i(2)=0.055+pos_i(2)*0.89;
      n_pos_i(3)=pos_i(3)*0.83;
      n_pos_i(4)=pos_i(4)*0.89;
      set(allax(i),'position',n_pos_i);
    end;
  end;
else, % no border, just text
  text(0.98*xmax,0.02*ymax, ...
       plottext, ...
      'parent',hBorder', ...
      'fontsize',5, ...
      'horizontalalignment','right', ...
      'verticalalignment','bottom');
end;

axes(ax);
units0=get(0,'units');
set(0,'units','centimeters');
maxdim=get(0,'screensize');
maxdim=maxdim(3:4);
if strcmp(Orientation,'landscape'),
  pos1=round([29.6774 20.984]*min(maxdim./[29.6774 20.984]));
  pos2=round([29.6774 20.984]*min(fliplr(maxdim)./[29.6774 20.984]));
  pos=min(pos1,pos2);
else, % 'portrait'
  pos1=round([20.984 29.6774]*min(fliplr(maxdim)./[20.984 29.6774]));
  pos2=round([20.984 29.6774]*min(maxdim./[20.984 29.6774]));
  pos=min(pos1,pos2);
end;   
pos=pos*0.85;
pos=[(maxdim(1)-pos(1))/2 (maxdim(2)-pos(2))/2 pos];
set(fg, ...
 'units','centimeters', ...
   'handlevisibility','off', ...
   'position',pos);
set(fg,'units','pixels');
set(0,'units',units0);
set(fg,'children',[allchld;hBorder]);


function Str=DateStr;
t=[datestr(now,13) ' on ' datestr(now,8) ' '];
x=clock;
if x(3)>3,
  t=[t num2str(x(3)) 'th'];
elseif x(3)==1,
  t=[t num2str(x(3)) 'st'];
elseif x(3)==2,
  t=[t num2str(x(3)) 'nd'];
elseif x(3)==3,
  t=[t num2str(x(3)) 'rd'];
end;
Str=[t ' ' datestr(now,3) ' ' datestr(now,10)];


