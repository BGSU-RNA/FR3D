% mEditbox(text,title,button) creates a text window
%Ali Mokdad, Aug 2006

function editbox(text,title,fontsize,button) 

if nargin < 3,
  fontsize = 10;
end

% title may be missing, so assign it second;
%if there is a fourth variable (no matter what it is) a pushbutton will be created to close the window

if ~exist('title','var')
    title='';
end

fig = figure('NumberTitle','off','MenuBar','none','Name',title,'Visible','off'); 

%BUG SOLUTION: Visible set to off then to on after putting the text in the 
%edit box so that the slide bar will appear at the top when figure appears

if exist('button','var')
  e=.05;      %this value pushes the edit box up to allow enough space for the OK button
  buttonh=uicontrol('Style','pushbutton','Units','normalized','Position',[.4 0 .2 e],'String','OK','Callback','close');
else 
  e=0;
end

edith = uicontrol('Style','edit','Units','normalized','Position',[0 e 1 1-e],'BackgroundColor',[1 1 1],...
    'Max',length(text),'String',text);

set(edith, 'FontName', 'FixedWidth')
set(edith, 'FontSize', fontsize)
set(edith, 'HorizontalAlignment','Left')

%Max insures multiple lines when needed, so that scroll bar is created
% set(edith,'Max',length(text),'String',text); 

drawnow %drawnow is only necessary on R14SP2
set(fig,'Visible','on')

