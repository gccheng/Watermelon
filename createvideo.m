function varargout = createvideo(varargin)
% CREATEVIDEO MATLAB code for createvideo.fig
%      CREATEVIDEO, by itself, creates a new CREATEVIDEO or raises the existing
%      singleton*.
%
%      H = CREATEVIDEO returns the handle to a new CREATEVIDEO or the handle to
%      the existing singleton*.
%
%      CREATEVIDEO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CREATEVIDEO.M with the given input arguments.
%
%      CREATEVIDEO('Property','Value',...) creates a new CREATEVIDEO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before createvideo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to createvideo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help createvideo

% Last Modified by GUIDE v2.5 17-Feb-2012 16:01:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @createvideo_OpeningFcn, ...
                   'gui_OutputFcn',  @createvideo_OutputFcn, ...
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


% --- Executes just before createvideo is made visible.
function createvideo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to createvideo (see VARARGIN)

% Choose default command line output for createvideo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes createvideo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = createvideo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function txt_picpath_Callback(hObject, eventdata, handles)
% hObject    handle to txt_picpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_picpath as text
%        str2double(get(hObject,'String')) returns contents of txt_picpath as a double


% --- Executes during object creation, after setting all properties.
function txt_picpath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_picpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_vidpath_Callback(hObject, eventdata, handles)
% hObject    handle to txt_vidpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_vidpath as text
%        str2double(get(hObject,'String')) returns contents of txt_vidpath as a double


% --- Executes during object creation, after setting all properties.
function txt_vidpath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_vidpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_vidpath.
function btn_vidpath_Callback(hObject, eventdata, handles)
% hObject    handle to btn_vidpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uiputfile('test.avi','Select video file name');

if 0 ~= file
    set(handles.txt_vidpath, 'String', [path file]);
    set(hObject, 'Value', 1.0);
     
else
    set(handles.btn_vidpath, 'String', '');
    set(hObject, 'Value', 0.0);
end


% --- Executes on button press in btn_picpath.
function btn_picpath_Callback(hObject, eventdata, handles)
% hObject    handle to btn_picpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
test_folder = uigetdir('');

if 0 ~= test_folder
    set(handles.txt_picpath, 'String', test_folder);
    set(hObject, 'Value', 1.0);
     
else
    set(handles.btn_picpath, 'String', '');
    set(hObject, 'Value', 0.0);
end


% --- Executes on button press in btn_start.
function btn_start_Callback(hObject, eventdata, handles)
% hObject    handle to btn_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% The maximum directory depth is 2

vidname = get(handles.txt_vidpath, 'String');
folder = get(handles.txt_picpath, 'String');

listing = dir(folder);
count = size(listing, 1);

picturelist = cell(0);

% get all valid file names
for i=1:count
    name = listing(i).name;
    isdir = listing(i).isdir;
    if ((~isdir) && (name(1)~='.'))
        picturelist = [picturelist; name];
    elseif (isdir && (~strcmp('.',name)) && (~strcmp('..',name)))  
        sub_files = files_in_folder(name);
        picturelist = [picturelist; sub_files];
    end;
end;

aviobj = avifile(vidname, 'fps', 10, 'compression', 'None', 'quality', 100);

%aviObj = VideoWriter(vidname, 'Uncompressed AVI');
%aviObj.FrameRate = 10;
%open(aviObj);

cmap = gray(256);

szFiles = size(picturelist, 1);
for i=1:szFiles
    F_name = [folder '/' picturelist{i}];
    F = imread(F_name);
    M = im2frame(F,cmap);
    aviobj = addframe(aviobj, M);
    %writeVideo(aviObj,F);
end;

aviobj = close(aviobj);
%close(aviObj);


function [files] = files_in_folder(folder)
files = cell(0);
listing = dir(folder);
count = size(listing, 1);

for i=1:count
    name = listing(i).name;
    isdir = listing(i).isdir;
    if ((~isdir) && (name(1)~='.'))
        files = [files; name];
    end;
end;
