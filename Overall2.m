function varargout = Overall2(varargin)
% OVERALL2 MATLAB code for Overall2.fig
%      OVERALL2, by itself, creates a new OVERALL2 or raises the existing
%      singleton*.
%
%      H = OVERALL2 returns the handle to a new OVERALL2 or the handle to
%      the existing singleton*.
%
%      OVERALL2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OVERALL2.M with the given input arguments.
%
%      OVERALL2('Property','Value',...) creates a new OVERALL2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Overall2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Overall2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Overall2

% Last Modified by GUIDE v2.5 23-Aug-2011 16:17:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Overall2_OpeningFcn, ...
                   'gui_OutputFcn',  @Overall2_OutputFcn, ...
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


% --- Executes just before Overall2 is made visible.
function Overall2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Overall2 (see VARARGIN)

% Choose default command line output for Overall2
handles.output = hObject;

rstStAnalyse.dist = [];
rstStAnalyse.x = [];
rstStAnalyse.n = [];
set(handles.hist, 'UserData', rstStAnalyse);

confidence = [];
set(handles.conf, 'UserData', confidence);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Overall2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Overall2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

 
% --- Executes on button press in acvi.
function acvi_Callback(hObject, eventdata, handles)
% hObject    handle to acvi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = ...
     uigetfile({'*.avi';'*.mpeg';'*.mpg';'*.*'},'Open a video file');

if filename
    vid = mmreader([pathname filename]); 
    for i=1:5
        fBuf= read(vid, i);
    end
    fImag = fBuf;
    posImag = get(handles.imag, 'Position');
    axisHeight = posImag(4); axisWidth = posImag(3);
    
    newImag = imresize(fImag, [axisHeight axisWidth]);
    axes(handles.imag);
    imshow(newImag);    
    
    % Save filename and pathname
    set(handles.acvi, 'UserData', [pathname filename]);
    set(handles.imag, 'UserData', fImag);
    
    % Reset histogram data
    rstStAnalyse.dist = [];
    rstStAnalyse.x = [];
    rstStAnalyse.n = [];
    set(handles.hist, 'UserData', rstStAnalyse); 
    
    % Reset confidence data
    confidence = [];
    set(handles.conf, 'UserData', confidence);
    
    % Clear axes
    cla(handles.hist);
    cla(handles.conf);
end;

guidata(hObject, handles);



% --- Executes on button press in stan.
function stan_Callback(hObject, eventdata, handles)
% hObject    handle to stan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Read parameter setting
videofile = get(handles.acvi,'UserData');
duration = uint8(get(handles.dura,'Value'));
filtersize = uint8(get(handles.fisi,'Value'));

% Retrieve saved data and run
conf = get(handles.conf, 'UserData');
rstPrev = get(handles.hist, 'UserData');
[ distance xbin nbin confidence] = st_covmat_analyse(videofile, duration,...
                filtersize, handles, conf, rstPrev.dist, rstPrev.x, rstPrev.n);
disp('Analysis stopped.');

% Save updated data
rstStAnalyse.dist = distance;
rstStAnalyse.x = xbin;
rstStAnalyse.n = nbin;
set(handles.hist, 'UserData', rstStAnalyse);
set(handles.conf, 'UserData', confidence);

%save('distance_bagleftbehind', 'distance');


% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pushed_down = get(hObject, 'Value');
if pushed_down
    disp('Stopping analysis.');
else
    % Read parameter setting
    videofile = get(handles.acvi,'UserData');
    duration = uint8(get(handles.dura,'Value'));
    filtersize = uint8(get(handles.fisi,'Value'));

    % Retrieve saved data and run
    conf = get(handles.conf, 'UserData');
    rstPrev = get(handles.hist, 'UserData');
    [ distance xbin nbin confidence] = st_covmat_analyse(videofile, duration,...
                    filtersize, handles, rstPrev.dist, rstPrev.x, rstPrev.n, conf);
    % Save updated data
    rstStAnalyse.dist = distance;
    rstStAnalyse.x = xbin;
    rstStAnalyse.n = nbin;
    set(handles.hist, 'UserData', rstStAnalyse);
    set(handles.conf, 'UserData', confidence);
end

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.stop, 'Value', 1);
pause(2);
close


function fisi_Callback(hObject, eventdata, handles)
% hObject    handle to fisi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fisi as text
%        str2double(get(hObject,'String')) returns contents of fisi as a double


% --- Executes during object creation, after setting all properties.
function fisi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fisi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sigm_Callback(hObject, eventdata, handles)
% hObject    handle to sigm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sigm as text
%        str2double(get(hObject,'String')) returns contents of sigm as a double


% --- Executes during object creation, after setting all properties.
function sigm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dura_Callback(hObject, eventdata, handles)
% hObject    handle to dura (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dura as text
%        str2double(get(hObject,'String')) returns contents of dura as a double


% --- Executes during object creation, after setting all properties.
function dura_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dura (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function posi_Callback(hObject, eventdata, handles)
% hObject    handle to posi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of posi as text
%        str2double(get(hObject,'String')) returns contents of posi as a double


% --- Executes during object creation, after setting all properties.
function posi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to posi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in piup.
function piup_Callback(hObject, eventdata, handles)
% hObject    handle to piup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fImage = get(handles.imag, 'UserData');
height = size(fImage,1); width = size(fImage,2);
hFigure = figure('Name', 'Pick up a point to observe'); 
set(hFigure,'Position',[400,400,width,height]);
imshow(fImage,'Border', 'tight');
set(hFigure, 'WindowButtonDownFcn',{@my_pickup_Callback, handles, width, height});


% --- Executes on button press in trai.
function trai_Callback(hObject, eventdata, handles)
% hObject    handle to trai (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trai


% --- Executes on button press in dete.
function dete_Callback(hObject, eventdata, handles)
% hObject    handle to dete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dete
