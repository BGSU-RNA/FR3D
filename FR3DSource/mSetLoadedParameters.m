%mSetLoadedParameters.m fills in entries in the GUI from a saved search file

Query=Search.Query;

s  = get(handles.SearchPDBs,'String');
ss = get(handles.QueryPDB,'String');

v=[];

for i=1:length(Search.Filenames)
    ff = find(strcmp(upper(s),upper(Search.Filenames{i})));
    if ~isempty(ff),
      v=[v ff(1)];
    end
end
set(handles.SearchPDBs,'Value',v);

if Search.Query.Geometric == 1
    set(handles.Geometric,'Value',1);
    set(handles.NonGeometric,'Value',0);
    set(handles.ViewQuery,'Visible','on');
    set(handles.ReadQuery,'Visible','on');
    set(handles.GenerateMatrix,'Visible','off');

    v=find(strcmp(lower(ss),lower(Search.Query.Filename)));

if isempty(v),
  v = 1;
end

    set(handles.QueryPDB,'Value',v);
    set(handles.QueryPDB,'Visible','on');
    set(handles.QueryPDBTitle,'Visible','on');

% %     mGetPDBfilenames %defines s
% %     set(handles.SearchPDBs,'String',s);
% %     set(handles.SearchPDBs,'Min',1);
% %     set(handles.SearchPDBs,'Max',length(s)+1);
% % 
% %     set(handles.QueryPDB,'String',s);


    NT=Search.Query.NTList;
    A=NT{1};
    for i=2:length(NT)
        A=[A ',' NT{i}];
    end
    set(handles.QueryNTs,'String',A);
    set(handles.QueryNTs,'Visible','on');
    set(handles.QueryNTsTitle,'Visible','on');

    
    
% % %     for i=1:25
% % %         h=findobj('Tag',strcat('ChainPopup',num2str(i)));
% % %         delete(h);
% % %     end
h = findobj('-regexp','Tag','ChainPopup[0-9]');
delete(h);
    
    if length(NT)<13
        x=.054;
        y=.04;
    else
        x=.684/length(NT);
        y=.54/length(NT);
    end

    if isfield(Search.Query,'ChainList'),    
      for i=1:length(NT),
        handles.ChainPopup(i) = uicontrol('Tag',strcat('ChainPopup',num2str(i)),'Style','popupmenu','Units','normalized','Position',[(0.305+x*(i-1)) (0.795) x-.003 .04],'String',Search.Query.ChainList{i},'Background',[1 1 1]);
      end
    end

    set(handles.QueryChains,'Visible','on')%this the text just to the left of the popups
%     mCreateChains
%     for i=1:length(NT)
%         h=findobj('Tag',strcat('ChainPopup',num2str(i)));
%         s=get(h,'String');
%         v=find(strcmp(s,Search.Query.ChainList{i}));
%         set(h,'Value',v);
%     end
  set(handles.NumberOfNTsTitle,'Visible','off'); 
  set(handles.NumberOfNTs,'Visible','off');
  set(handles.GuarCutoff,'Visible','on');
  if isfield(Search.Query,'DiscCutoff')
      set(handles.GuarCutoff,'String',num2str(Search.Query.DiscCutoff));
  end
  set(handles.RelCutoff,'Visible','on');
  if isfield(Search.Query,'RelCutoff')
      set(handles.RelCutoff,'String',num2str(Search.Query.RelCutoff));
  end


else
  set(handles.GuarCutoff,'Visible','off');
  set(handles.RelCutoff,'Visible','off');
    set(handles.Geometric,'Value',0);
    set(handles.NonGeometric,'Value',1);
    set(handles.ViewQuery,'Visible','off');
    set(handles.ReadQuery,'Visible','off');
    set(handles.GenerateMatrix,'Visible','on');
    set(handles.NumberOfNTs,'String',num2str(Search.Query.NumNT));
    set(handles.NumberOfNTs,'Visible','on');
    set(handles.NumberOfNTsTitle,'Visible','on');
    set(handles.QueryPDBTitle,'Visible','off');
    set(handles.QueryNTsTitle,'Visible','off');
    set(handles.QueryChains,'Visible','off');
    set(handles.QueryPDB,'Visible','off');
    set(handles.QueryNTs,'Visible','off');
    for i=1:Search.Query.NumNT
        NT{i}=num2str(i);
    end
    % % %     for i=1:25
    % % %         h=findobj('Tag',strcat('ChainPopup',num2str(i)));
    % % %         delete(h);
    % % %     end
    h = findobj('-regexp','Tag','ChainPopup[0-9]');
    delete(h);
end

set(handles.SearchNameText,'Visible','on');
set(handles.SearchDescriptionText,'Visible','on');
set(handles.GuarCutoffText,'Visible','on');
set(handles.RelCutoffText,'Visible','on');
set(handles.SearchName,'Visible','on');
if isfield(Search.Query,'Number')
    set(handles.SearchName,'String',num2str(Search.Query.Name));%num2str works whether this is character or number
end
set(handles.SearchDescription,'Visible','on');
if isfield(Search.Query,'Description')
    set(handles.SearchDescription,'String',Search.Query.Description);
end
set(handles.Overlap,'Visible','on');
if isfield(Search.Query,'ExcludeOverlap')
    v=Search.Query.ExcludeOverlap;
    if v==0
        v=2;
    end
    set(handles.Overlap,'Value',v);
end
set(handles.DisplayCandidates,'Visible','on');
set(handles.ListCandidates,'Visible','on');

set(handles.RunSearch,'Visible','off');
set(handles.Status,'String','Loaded saved search results. You can examine these results using "Display Candidates" or "List Candidates", or you can repeat the analysis by starting from "Read Query" '); 

%%The matrix
mCreateMatrix_Loaded
