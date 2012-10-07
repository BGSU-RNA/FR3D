% mEditbox(text,title,button) creates a text window
%Ali Mokdad, Aug 2006

function editbox(text,title,button) 

% title may be missing, so assign it second;
%if there is a third variable (no matter what it is) a pushbutton will be created to close the window

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

edith = uicontrol('Style','edit','Units','normalized','Position',[0 e 1 1-e],'BackgroundColor',[1 1 1]);

set(edith, 'FontSize', 10);
set(edith, 'HorizontalAlignment','Left');

%[outstring,newpos] = textwrap(edith,text);
%outstring

set(edith,'Max',length(text));
set(edith,'String',text);

drawnow                       %drawnow is only necessary on R14SP2
set(fig,'Visible','on')

