% mCreateMatrix is a script that sets up the interaction matrix for FR3D_GUI

% modify the nucleotide index lookup to be aware of the given chain

%%%%%%Create the basepair matrix (and delete extra ones) for determining
%%%%%%Query.Diff, Query.ReqInter, Query.Diagonal
%%%%%%Query.Config added
if length(NT)<13
    x=.054;
    y=.038;
else
    x=.684/length(NT);
    y=.54/length(NT);
end


PreviousConfig=ones(1,25);
for i=1:25
    h=findobj('Tag',strcat('Config',num2str(i)));
    try,PreviousConfig(i)=get(h,'Value');end
    delete(h);
end
% % % h = findobj('-regexp','Tag','Config[0-9]');
% % % delete(h);


for i=1:length(NT)
    str={'','anti','syn'};
    handles.Config(i) = uicontrol('Tag',strcat('Config',num2str(i)),'Style','popupmenu','Units','normalized','Position',[(0.305+x*(i-1)) (0.8-y) x-.003 y-.005],'Background',[1 1 1],'String',str,'Value',PreviousConfig(i));
end

% % % for i=1:25
% % %     hh=findobj('Tag',strcat('BPtexth',num2str(i)));
% % %     delete(hh);
% % %     hv=findobj('Tag',strcat('BPtextv',num2str(i)));
% % %     delete(hv);
% % % end
hh = findobj('-regexp','Tag','BPtexth[0-9]');
hv = findobj('-regexp','Tag','BPtextv[0-9]');
delete(hh,hv)

for i=1:length(NT)
    if get(handles.Geometric,'Value') == 1
        ind = zIndexLookup(File(QIndex),NT{i},ChainList{i});
        Bas = File(QIndex).NT(ind).Base;
    else
        Bas='NT';
    end
    %%%Horizontal title line:
%     handles.BPtexth(i) = uicontrol('Tag',strcat('BPtexth',num2str(i)),'Style','text','Units','normalized','Position',[(0.25+0.057*i) (0.73) .054 .04],'String',strcat(File(QIndex).NT(str2num(NT{i})).Base,NT(i)));
    handles.BPtexth(i) = uicontrol('Tag',strcat('BPtexth',num2str(i)),'Style','text','Units','normalized','Position',[(0.305+x*(i-1)) (0.75-y) x-.003 y-.005],'String',strcat(Bas,NT(i)));
    %%%Vertical title line:
    handles.BPtextv(i) = uicontrol('Tag',strcat('BPtextv',num2str(i)),'Style','text','Units','normalized','Position',[(0.26) (0.725-y*i) .056 y-.005],'String',strcat(Bas,NT(i)));
end

%%%Now the matrix:
% x           = ones(1,25);
% y           = num2str(x);
% y           = regexprep(y,' ','');
% Sequences   = regexprep(y,'1','N');

for i=1:25
    hd=findobj('Tag',strcat('Diagonal',num2str(i)));
    try,PreviousDiagonal{i}=get(hd,'String');end
    delete(hd);
    for j=1:25
        hh=findobj('Tag',strcat('BPType',num2str(i),num2str(j)));
        try,PreviousBPType{i,j}=get(hh,'String');end
        delete(hh);
        hv=findobj('Tag',strcat('Diff',num2str(i),num2str(j)));
        try,PreviousDiff{i,j}=get(hv,'String');end
        delete(hv);
    end
end

% % % hd = findobj('-regexp','Tag','Diagonal[0-9]');
% % % hbp = findobj('-regexp','Tag','BPType[0-9]');
% % % hdif = findobj('-regexp','Tag','Diff[0-9]');
% % % delete(hh,hv,hd,hbp,hdif)


% % % if ~isempty(PreviousMasks)
% % %     PreviousMasks(1)
% % %     PreviousBPType{1,2}
% % %     PreviousDiff{2,1}
% % % end

%MaskList={'N (ACGU)','A','C','G','U','R (AG)','Y (CU)','M (AC)','W (AU)','S (GC)','K (GU)','V (ACG)','H (ACU)','D (AGU)','B (CGU)'};
for i=1:length(NT)
    for j=1:length(NT)
        if j==i
            if isempty(PreviousDiagonal{i}),PreviousDiagonal{i}='N';end
            handles.Diagonal(i) = uicontrol('Tag',strcat('Diagonal',num2str(i)),'Style','edit','Units','normalized','Position',[(0.305+x*(j-1)) (0.73-y*i) x-.003 y-.005],'Background',[1 1 1],'String',PreviousDiagonal{i});
        end
        if i<j
            if isempty(PreviousBPType{i,j}),PreviousBPType{i,j}='';end
            handles.BPType(i,j) = uicontrol('Tag',strcat('BPType',num2str(i),num2str(j)),'Style','edit','Units','normalized','Position',[(0.305+x*(j-1)) (0.73-y*i) x-.003 y-.005],'Background',[1 1 0],'String',PreviousBPType{i,j});
        elseif i>j
            if isempty(PreviousDiff{i,j}),PreviousDiff{i,j}='';end
            handles.Diff(i,j) = uicontrol('Tag',strcat('Diff',num2str(i),num2str(j)),'Style','edit','Units','normalized','Position',[(0.305+x*(j-1)) (0.73-y*i) x-.003 y-.005],'Background',[0 1 1],'String',PreviousDiff{i,j});
        end
    end
end

set(handles.ConfigText,'Visible','on');
set(handles.NTmaskText,'Visible','on');
set(handles.MaxDistText,'Visible','on');
set(handles.InteractionText,'Visible','on');
