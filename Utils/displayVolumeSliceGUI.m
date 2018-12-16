function displayVolumeSliceGUI(varargin)
%{
Usage : displayVolumeSliceGUI(X)  , where X is the 3D volumetric array
credit: Ohad Menashe
%}

close all

fh = figure('name','VolumeSliceGUI','NumberTitle','off','menubar','none','toolBar','figure');
h = guidata(fh);

if(nargin==0)
    vars = evalin('base','who');
    a=cellfun(@(x) evalin('base',['length(size(' x '))']),vars);
    vars = vars(a==3);
    
    
     if(length(vars)~=1)
         v=listdlg('PromptString','Select volume','SelectionMode','single','ListString',vars);
         
     else
         v = 1;
     end
         h.vol = evalin('base',vars{v});
else
    if(ischar(varargin{1}))
        h.vol = loadModel(varargin{1});
    elseif(length(size(varargin{1}))==3)
        h.vol=varargin{1};
    else
        error('unknonwn input');
    end
end
h.minmax = [min(h.vol(:)) eps+max(h.vol(:))];


guidata(fh,h);
%draw


h.a = axes('parent',fh);
guidata(fh,h);
h.xyz = zeros(4,1);
h.xyz(4) = uibuttongroup('parent',fh,'BorderType','none','SelectionChangeFcn',@dimChange);
h.xyz(3)=uicontrol('Style','radio','parent',h.xyz(4),'string','Z','units','normalized','position',[ 2/3 0 1/3 1]);
h.xyz(2)=uicontrol('Style','radio','parent',h.xyz(4),'string','X','units','normalized','position',[ 1/3 0 1/3 1]);
h.xyz(1)=uicontrol('Style','radio','parent',h.xyz(4),'string','Y','units','normalized','position',[ 0/3 0 1/3 1]);
h.slider = uicontrol('style','slider','parent',fh,'callback',@updateAxes);