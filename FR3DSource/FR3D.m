% FR3D is Find RNA 3D

% 33(9),36(9),37:39(9),43:47(9)

%<div class="moz-text-flowed" style="font-family: -moz-fixed">
function FR3D(varargin)   %By Ali Mokdad - March 2006

%fprintf('\nThis is FR3D.\n');

%fprintf('Current directory is %s\n', pwd);

%fprintf('Current path is %s\n', path);

%zDisplayNT('1s72',51:60);

%load('C:\FR3D\SearchSaveFiles\2007-08-15_15_23_24-Sarcin_5-nucleotide_in_1s72');
%xDisplayCandidates([],Search);

%figure(2)
%plot(cumsum(rand(1,100)));

% FR3D M-file for FR3D.fig
%      FR3D, by itself, creates a new FR3D or raises the existing
%      singleton*.
%
%      H = FR3D returns the handle to a new FR3D or the handle to
%      the existing singleton*.
%
%      FR3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FR3D.M with the given input arguments.
%
%      FR3D('Property','Value',...) creates a new FR3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FR3D_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FR3D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to runsearch (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FR3D

% Last Modified by GUIDE v2.5 17-Jul-2006 19:31:53

%fprintf('Current path is %s\n', path);

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @FR3D_OpeningFcn, ...
    'gui_OutputFcn',  @FR3D_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FR3D is made visible.
function FR3D_OpeningFcn(hObject, eventdata, handles, varargin)

if ~(exist('FR3DSource') == 7),        % if directory doesn't yet exist
  mkdir('FR3DSource');
end
path(path,[pwd filesep 'FR3DSource']);

if ~(exist('PDBFiles') == 7),        % if directory doesn't yet exist
  mkdir('PDBFiles');
end
path(path,[pwd filesep 'PDBFiles']);

if ~(exist('PrecomputedData') == 7),        % if directory doesn't yet exist
  mkdir('PrecomputedData');
end
path(path,[pwd filesep 'PrecomputedData']);

if ~(exist('SearchSaveFiles') == 7),        % if directory doesn't yet exist
  mkdir('SearchSaveFiles');
end
path(path,[pwd filesep 'SearchSaveFiles']);

savedf=dir(['SearchSaveFiles' filesep '*.mat']);
if length(savedf) == 0,
  savedff=' ';
else
  for i = 1:length(savedf),
    savedff{i,1} = savedf(i).name;
  end
end
set(handles.LOAD,'String',savedff);

[s,snolist] = mGetPDBFilenames;                          % defines PDB lists

if length(s) == 0,
    set(handles.Status,'String','ALERT: there are no PDB files or saved PDB data in the folders "PDBFiles" and "PrecomputedData". Please put some PDB files in these locations and try again!');
    set(handles.ReadQuery,'Visible','Off'); %just to prevent the user from clicking it anyway!
else
  set(handles.SearchPDBs,'String',s);
  set(handles.SearchPDBs,'Min',1);
  set(handles.SearchPDBs,'Max',length(s)+1);
  
  set(handles.SearchPDBs,'Value',1);

  set(handles.QueryPDB,'String',snolist);

  p = find(ismember(snolist,'1s72'));
  if ~isempty(p),
    set(handles.QueryNTs,'String','2692:2694, 2701:2702');
    set(handles.QueryPDB,'Value',p);
  else
    set(handles.QueryNTs,'String','List query nucleotides and chains here');
  end

end


handles.output = hObject;
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function FR3D_OutputFcn(hObject, eventdata, handles)
%varargout{1} = handles.output;
%File = handles.File;


% ----------------------------------------- Loads saved search

function LOAD_Callback(hObject, eventdata, handles)

savedf = dir(['SearchSaveFiles' filesep '*.mat']);
if length(savedf) == 0,
  savedff=' ';
else
  for i = 1:length(savedf),
    savedff{i,1} = savedf(i).name;
  end
end

set(handles.LOAD,'String',savedff);  % Update list without having to reopen GUI

v=get(handles.LOAD,'Value');
f=savedff{v,1};
l= ['SearchSaveFiles' filesep f];      % prepend directory name
load(l)                              % Load search data

if (~isfield(Search,'File')) && (length(Search.Candidates(:,1)) > 0),
  Search = xAddFiletoSearch([],Search);  % add candidate filenames

  if isfield(handles,'File')
    File = handles.File;
    [File,FIndex] = zAddNTData(Search.CandidateFilenames,2,File);
  else
    [File,FIndex] = zAddNTData(Search.CandidateFilenames,2);
  end
  Search = xAddFiletoSearch(File(FIndex),Search);
  save(['SearchSaveFiles' filesep Search.SaveName], 'Search');
  fprintf('Saved updated version of this search\n');

  handles.File = File;
end

mSetLoadedParameters

handles.Search=Search;
guidata(hObject, handles);

% ----------------------------------------- Change query PDB

function QueryPDBedit_Callback(hObject, eventdata, handles)
function QueryPDBedit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ----------------------------------------- Change query nucleotides

function QueryNTs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function QueryNTs_Callback(hObject, eventdata, handles)


% ------------------------------- Executes on button press in ReadQuery

function ReadQuery_Callback(hObject, eventdata, handles)

  %%%Read Query PDB File or just use the one in memory if it is present and same as the new Query PDB
  x=get(handles.QueryPDB,'Value');
  s=get(handles.QueryPDB,'String');
  Query.Filename = s{x};
  set(handles.Status,'String','Reading query PDB, please wait ...');
  drawnow
  %%handles.File.Filename 

  if isfield(handles,'File')
    File = handles.File;
    [File,QIndex]=zAddNTData(Query.Filename,2,File);
  else
    [File,QIndex]=zAddNTData(Query.Filename,2);
  end

  %%%DetermineQuery.NTList %%%%This must be done before running mCreateDynamicGUI

  [Indices,Ch] = zIndexLookup(File(QIndex),get(handles.QueryNTs,'String'));
%   if length(Indices) > 25,
%     set(handles.Status,'String','Error: It is not possible to use this GUI for a motif larger than 25 NTs.');
%     set(handles.GenerateMatrix,'Visible','off');
%   end

  %%%

  %     mCreateDynamicGUI 
  % Creates all dynamic popum menus and edit boxes and hides extra ones 
  % from previous searches

  % Create popupmenus (and delete extra ones) for determining Query.ChainPopup

% % %   for i=1:25
% % %     h=findobj('Tag',strcat('ChainPopup',num2str(i)));
% % %     delete(h);
% % %   end
h = findobj('-regexp','Tag','ChainPopup[0-9]');
delete(h);


  if length(Indices)<13
      x=.054;
      y=.04;
  else
      x=.684/min(25,length(Indices));
      y=.54 /min(25,length(Indices));
  end

  for i=1:min(25,length(Indices)),               % truncate at 25 nucleotides
    Query.NTList{i} = File(QIndex).NT(Indices(i)).Number;
    handles.ChainPopup(i) = uicontrol('Tag',strcat('ChainPopup',num2str(i)),'Style','popupmenu','Units','normalized','Position',[(0.305+x*(i-1)) (0.795) x-.003 .04],'String',Ch{i},'Background',[1 1 1]);
  end

  set(handles.QueryChains,'Visible','on')
                             % this is the text just to the left of the popups

    %Pass some variables created here to be used by other functions
    handles.Indices = Indices;
    handles.NT=Query.NTList;
    handles.Filename = Query.Filename;
    handles.File = File;
    handles.QIndex=QIndex;               % index of File for query file
    handles.NTList=Query.NTList;
    guidata(hObject, handles);

    if get(handles.ViewQuery,'Value')==1
        figure(3)
        clf
        VP.Sugar = 1;
        zDisplayNT(File(QIndex),Indices,VP);
        grid off
        rotate3d on
        zShowInteractionTable(File(QIndex),Indices)
    end
    set(handles.Status,'String','Choose "Query Chains" and click "Generate Interaction Matrix"');
    set(handles.GenerateMatrix,'Visible','on');

set(handles.RunSearch,'Visible','off');
set(handles.ListCandidates,'Visible','off');
set(handles.DisplayCandidates,'Visible','off');

% ------------------------------------------------ Generate interaction matrix

function GenerateMatrix_Callback(hObject, eventdata, handles)

if get(handles.Geometric,'Value') == 1
    NT=handles.NT;
    File=handles.File;
    QIndex=handles.QIndex;

    %%%Read Query.ChainList from GUI
    for i=1:length(NT)
        h=findobj('Tag',strcat('ChainPopup',num2str(i)));
        s=get(h,'String');
        v=get(h,'Value');
        Query.ChainList{i}=s(v);
    end
    %%%End Read Query.ChainList from GUI
    ChainList=Query.ChainList;
    handles.ChainList=ChainList;
    handles.NTlen=length(NT);
    set(handles.Status,'String','Ready to Search');
else
    NTlen=str2num(get(handles.NumberOfNTs,'String'));
    if NTlen>25
        set(handles.Status,'String','25 is the maximum number of nucleotides that can be addressed in this GUI. Ready to "Search"');
    else
        set(handles.Status,'String','Ready to "Search"');
    end
    NTlen=min(25,NTlen);
    for i=1:NTlen
      NT{i}=num2str(i);
    end
    handles.NTlen=NTlen;
end

mCreateMatrix                                  % generate interaction matrix

set(handles.GuarCutoffText,'visible','on');
set(handles.GuarCutoff,'visible','on');
set(handles.RelCutoffText,'visible','on');
set(handles.RelCutoff,'visible','on');
set(handles.Overlap,'visible','on');
set(handles.RunSearch,'visible','on');
set(handles.SearchNameText,'visible','on');
set(handles.SearchName,'visible','on');
set(handles.SearchDescriptionText,'visible','on');
set(handles.SearchDescription,'visible','on');


guidata(hObject, handles);


% ----------------------------- Executes on button press in RunSearch
function RunSearch_Callback(hObject, eventdata, handles)
%if ~isdeployed,clc,end

% if get(handles.Geometric,'Value') == 1

if isfield(handles,'File')
    File=handles.File;
end
if get(handles.Geometric,'Value') == 1
    Query.Filename=handles.Filename;
    Query.ChainList=handles.ChainList;
    Query.NTList=handles.NTList;
end
NTlen=handles.NTlen;
mSpecifyQueryforGUI;

%%Now determine PDB "Filenames" to search in:
x=get(handles.SearchPDBs,'Value');
s=get(handles.SearchPDBs,'String');
for i=1:length(x)
    Filenames{i,1}=s{x(i)};
end
%%

GUIactive = 1;

set(handles.Status,'String','Searching ... Please wait (press Ctrl+C to interrupt)');
drawnow

Verbose = 1;                               % limit notifications
xFR3DSearch

handles.File=File;
%handles.SIndex=SIndex;

if ~isempty(Candidates),
    set(handles.ListCandidates,'Visible','on');
    set(handles.DisplayCandidates,'Visible','on');
    if length(Discrepancy)>1
        s='s';
    else s='';
    end
    str=strcat(num2str(length(Discrepancy)),' candidate',s,' found');
    set(handles.Status,'String',str);

    handles.Search=Search;

else
    set(handles.Status,'String','No candidates found');
    set(handles.ListCandidates,'Visible','off');
    set(handles.DisplayCandidates,'Visible','off');
end
guidata(hObject, handles);


function SearchName_Callback(hObject, eventdata, handles)
function SearchName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SearchPDBs_Callback(hObject, eventdata, handles)
function SearchPDBs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function QueryPDB_Callback(hObject, eventdata, handles)
function QueryPDB_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ViewQuery_Callback(hObject, eventdata, handles)


function NumberOfNTs_Callback(hObject, eventdata, handles)
function NumberOfNTs_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% ------------------ Executes on button press in Geometric.
function Geometric_Callback(hObject, eventdata, handles)
set(handles.Geometric,'Value',1);
set(handles.NonGeometric,'Value',0);

set(handles.NumberOfNTsTitle,'Visible','off');
set(handles.NumberOfNTs,'Visible','off');
set(handles.GenerateMatrix,'Visible','off');

set(handles.QueryPDBTitle,'Visible','on');
set(handles.QueryPDB,'Visible','on');
set(handles.ViewQuery,'Visible','on');
set(handles.QueryNTsTitle,'Visible','on');
set(handles.QueryNTs,'Visible','on');
set(handles.ReadQuery,'Visible','on');

set(handles.Status,'String','Hint: Choose "Query PDB" and "Query NTs", then click "Read Query" - Repeat anytime to restart');
Query.Geometric = 1;


% -----------------  Executes on button press in NonGeometric.
function NonGeometric_Callback(hObject, eventdata, handles)
set(handles.NonGeometric,'Value',1);
set(handles.Geometric,'Value',0)

set(handles.QueryPDBTitle,'Visible','off');
set(handles.ViewQuery,'Visible','off');
set(handles.QueryNTsTitle,'Visible','off');
set(handles.ReadQuery,'Visible','off');

set(handles.QueryPDB,'Visible','off');
set(handles.QueryNTs,'Visible','off');
% h=findobj('Tag','QueryPDB');delete(h);
% h=findobj('Tag','QueryNTs');delete(h);
% % % for i=1:25
% % %     h=findobj('Tag',strcat('ChainPopup',num2str(i)));
% % %     delete(h);
% % % end
h = findobj('-regexp','Tag','ChainPopup[0-9]');
delete(h);

set(handles.QueryChains,'Visible','off')

set(handles.NumberOfNTsTitle,'Visible','on');
set(handles.NumberOfNTs,'Visible','on');
set(handles.GenerateMatrix,'Visible','on');

set(handles.Status,'String','Hint: Input "Number of NTs" (positive integer) and click "Generate Interaction Matrix" - Repeat anytime to restart');
Query.Geometric = 0;


% --------------- 
function SearchDescription_Callback(hObject, eventdata, handles)
function SearchDescription_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------- Executes on button press in Display Candidates
function DisplayCandidates_Callback(hObject, eventdata, handles)

Search=handles.Search;

if isfield(handles,'File')
  File = handles.File;
  [Search,File] = xDisplayCandidates(File,Search);
else
  [Search,File] = xDisplayCandidates([],Search);
end

handles.File=File;
handles.Search = Search;
save(['SearchSaveFiles' filesep Search.SaveName], 'Search');
guidata(hObject, handles);

% --------------- Executes on button press in List Candidates
function ListCandidates_Callback(hObject, eventdata, handles)

Search=handles.Search;
xListCandidates(Search,Inf);

% --- Executes on selection change in Overlap.
function Overlap_Callback(hObject, eventdata, handles)
% hObject    handle to Overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Overlap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Overlap


% --- Executes during object creation, after setting all properties.
function Overlap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Overlap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function GuarCutoff_Callback(hObject, eventdata, handles)
% hObject    handle to GuarCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GuarCutoff as text
%        str2double(get(hObject,'String')) returns contents of GuarCutoff as a double


% --- Executes during object creation, after setting all properties.
function GuarCutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GuarCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%</div>



function RelCutoff_Callback(hObject, eventdata, handles)
% hObject    handle to RelCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RelCutoff as text
%        str2double(get(hObject,'String')) returns contents of RelCutoff as a double


% --- Executes during object creation, after setting all properties.
function RelCutoff_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RelCutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in InteractionText.
function InteractionText_Callback(hObject, eventdata, handles)
xDescribeInteraction

% --- Executes on button press in MaxDistText.
function MaxDistText_Callback(hObject, eventdata, handles)
xDescribeDistance

% --- Executes on button press in NTmaskText.
function NTmaskText_Callback(hObject, eventdata, handles)
xDescribeMask

