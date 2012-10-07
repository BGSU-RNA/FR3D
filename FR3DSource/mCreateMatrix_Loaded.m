%mCreateMatrix_Loaded

%%%%%%Create the basepair matrix (and delete extra ones) for determining
%%%%%%Query.Diff, Query.Edges, Query.Diagonal
%%%%%%Query.Config added
if length(NT)<13
    x=.054;
    y=.04;
else
    x=.684/length(NT);
    y=.54/length(NT);
end

% % % for i=1:25
% % %     h=findobj('Tag',strcat('Config',num2str(i)));
% % %     delete(h);
% % % end
h = findobj('-regexp','Tag','Config[0-9]');
delete(h);

for i=1:length(NT)
    str={'','anti','syn'};
    if isfield(Search.Query,'Config'),
      switch Search.Query.Config{i}
        case 'anti' 
          PreviousConfig(i) = 2;
        case 'syn'  
          PreviousConfig(i) = 3;
        otherwise   
          PreviousConfig(i) = 1;
      end
    else
      PreviousConfig(i) = 1;
    end
    handles.Config(i) = uicontrol('Tag',strcat('Config',num2str(i)),'Style','popupmenu','Units','normalized','Position',[(0.305+x*(i-1)) (0.8-y) x-.003 y-.005],'Background',[1 1 1],'String',str,'Value',PreviousConfig(i));
end

% % % for i=1:25
% % %     hh=findobj('Tag',strcat('BPtexth',num2str(i)));
% % %     delete(hh);
% % %     hv=findobj('Tag',strcat('BPtextv',num2str(i)));
% % %     delete(hv);
% % %     hd=findobj('Tag',strcat('Diagonal',num2str(i)));
% % %     delete(hd);
% % %     for j=1:25
% % %         hh=findobj('Tag',strcat('BPType',num2str(i),num2str(j)));
% % %         delete(hh);
% % %         hv=findobj('Tag',strcat('Diff',num2str(i),num2str(j)));
% % %         delete(hv);
% % %     end
% % % end
hh = findobj('-regexp','Tag','BPtexth[0-9]');
hv = findobj('-regexp','Tag','BPtextv[0-9]');
hd = findobj('-regexp','Tag','Diagonal[0-9]');
hbp = findobj('-regexp','Tag','BPType[0-9]');
hdif = findobj('-regexp','Tag','Diff[0-9]');
delete(hh,hv,hd,hbp,hdif)



for i=1:length(NT)
    if isfield(Search.Query,'NT')
        Bas=Search.Query.NT(i).Base;
        Num=Search.Query.NT(i).Number;
        str=strcat(Bas,Num);
    else
        Bas='NT';
        str=strcat(Bas,num2str(i));
    end
    handles.BPtexth(i) = uicontrol('Tag',strcat('BPtexth',num2str(i)),'Style','text','Units','normalized','Position',[(0.305+x*(i-1)) (0.75-y) x-.003 y-.005],'String',str);
    handles.BPtextv(i) = uicontrol('Tag',strcat('BPtextv',num2str(i)),'Style','text','Units','normalized','Position',[(0.26) (0.725-y*i) .056 y-.005],'String',str);
end

%%%Now the matrix:
%MaskList={'N (ACGU)','A','C','G','U','R (AG)','Y (CU)','M (AC)','W (AU)','S (GC)','K (GU)','V (ACG)','H (ACU)','D (AGU)','B (CGU)'};
PreviousMasks=ones(1,25); %declared
% PreviousDiff={''};
% PreviousBPType={''};

for i=1:length(NT)
    for j=1:length(NT)
        if j==i
            if isfield(Search.Query,'Diagonal')
                PreviousDiagonal{i}=Search.Query.Diagonal{i};
            else PreviousDiagonal{i}='N';
            end
            handles.Diagonal(i) = uicontrol('Tag',strcat('Diagonal',num2str(i)),'Style','edit','Units','normalized','Position',[(0.305+x*(j-1)) (0.73-y*i) x-.003 y-.005],'Background',[1 1 1],'String',PreviousDiagonal{i});
        end
        if i<j
            if isfield(Search.Query,'Edges')
                if ~isempty(Search.Query.Edges{i,j}),
                    PreviousBPType{i,j}=Search.Query.Edges{i,j};
                else
                    PreviousBPType{i,j}='';
                end
            else PreviousBPType{i,j}='';
            end
            handles.BPType(i,j) = uicontrol('Tag',strcat('BPType',num2str(i),num2str(j)),'Style','edit','Units','normalized','Position',[(0.305+x*(j-1)) (0.73-y*i) x-.003 y-.005],'String',PreviousBPType{i,j},'Background',[1 1 0]);
        elseif i>j
            if isfield(Search.Query,'Diff')
                if ~isempty(Search.Query.Diff(i,j)),
                    PreviousDiff{i,j}=Search.Query.Diff{i,j};
                else
                    PreviousDiff{i,j}='';
                end
            else PreviousDiff{i,j}='';
            end
            handles.Diff(i,j) = uicontrol('Tag',strcat('Diff',num2str(i),num2str(j)),'Style','edit','Units','normalized','Position',[(0.305+x*(j-1)) (0.73-y*i) x-.003 y-.005],'String',PreviousDiff{i,j},'Background',[0 1 1]);
        end
    end
end

set(handles.ConfigText,'Visible','on');
set(handles.NTmaskText,'Visible','on');
set(handles.MaxDistText,'Visible','on');
set(handles.InteractionText,'Visible','on');

% Search.Query.Edges
% PreviousBPType
% Search.Query.Diff
% PreviousDiff

