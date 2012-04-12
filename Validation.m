function varargout = Validation(varargin)
% VALIDATION MATLAB code for Validation.fig
%      VALIDATION, by itself, creates a new VALIDATION or raises the existing
%      singleton*.
%
%      H = VALIDATION returns the handle to a new VALIDATION or the handle to
%      the existing singleton*.
%
%      VALIDATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VALIDATION.M with the given input arguments.
%
%      VALIDATION('Property','Value',...) creates a new VALIDATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Validation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Validation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Validation

% Last Modified by GUIDE v2.5 27-Feb-2012 13:45:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Validation_OpeningFcn, ...
                   'gui_OutputFcn',  @Validation_OutputFcn, ...
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


% --- Executes just before Validation is made visible.
function Validation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Validation (see VARARGIN)

% Choose default command line output for Validation
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Validation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Validation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function eb_testvideo_Callback(hObject, eventdata, handles)
% hObject    handle to eb_testvideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eb_testvideo as text
%        str2double(get(hObject,'String')) returns contents of eb_testvideo as a double


% --- Executes during object creation, after setting all properties.
function eb_testvideo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eb_testvideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_testvideo.
function btn_testvideo_Callback(hObject, eventdata, handles)
% hObject    handle to btn_testvideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

test_folder = uigetdir('');

if 0 ~= test_folder
    set(handles.eb_testvideo, 'String', test_folder);
    set(hObject, 'Value', 1.0);
    
    listing = dir(test_folder);
    
    firstFrame = imread([test_folder '/' listing(3).name]);
    axes(handles.imag); imshow(firstFrame);
    set(handles.imag, 'UserData', firstFrame);
    
else
    set(handles.eb_testvideo, 'String', '');
    set(hObject, 'Value', 0.0);
end


% --- Executes on button press in btn_groundtruth.
function btn_groundtruth_Callback(hObject, eventdata, handles)
% hObject    handle to btn_groundtruth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

test_folder = uigetdir('');

if 0 ~= test_folder
    set(handles.eb_groundtruth, 'String', test_folder);
    set(hObject, 'Value', 1.0);
else
    set(handles.eb_groundtruth, 'String', '');
    set(hObject, 'Value', 0.0);
end

function eb_groundtruth_Callback(hObject, eventdata, handles)
% hObject    handle to eb_groundtruth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eb_groundtruth as text
%        str2double(get(hObject,'String')) returns contents of eb_groundtruth as a double


% --- Executes during object creation, after setting all properties.
function eb_groundtruth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eb_groundtruth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_onesite.
function btn_onesite_Callback(hObject, eventdata, handles)
% hObject    handle to btn_onesite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of btn_onesite

checked = get(hObject, 'Value');
fImage = get(handles.imag, 'UserData');

if (1==checked) && (~isempty(fImage))
    height = size(fImage,1); width = size(fImage,2);
    hFigure = figure('Name', 'Pick up a point to observe');
    set(hFigure,'Position',[200,200,width,height]);
    imshow(fImage,'Border', 'tight');
    set(hFigure, 'WindowButtonDownFcn',{@my_pickup_Callback, handles, width, height});
    guidata(hObject, handles);
elseif (1==checked)
    warndlg('Please choose test video first.','warning', 'modal');
    set(hObject, 'Value', 0.0);
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


% --- Executes on button press in btn_test.
function btn_test_Callback(hObject, eventdata, handles)
% hObject    handle to btn_test (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

videodir = get(handles.eb_testvideo, 'String');
filtersize = 3;

if exist('Histograms.mat','file') && exist('mcov.mat','file') && exist('mmu.mat','file')
    HistTrained = 1;
else
    HistTrained = 0;
    warndlg('Model not trained!', 'warning', 'modal');
    return;
end

[ distance xbin nbin confidence] = st_multiple_validate(videodir, ...
        filtersize, handles, HistTrained);
    


% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(hObject, 'Value', 1.0);



function eb_trainvideo_Callback(hObject, eventdata, handles)
% hObject    handle to eb_trainvideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eb_trainvideo as text
%        str2double(get(hObject,'String')) returns contents of eb_trainvideo as a double


% --- Executes during object creation, after setting all properties.
function eb_trainvideo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eb_trainvideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_trainvideo.
function btn_trainvideo_Callback(hObject, eventdata, handles)
% hObject    handle to btn_trainvideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

train_folder = uigetdir('');

if 0 ~= train_folder
    set(handles.eb_trainvideo, 'String', train_folder);
    set(hObject, 'Value', 1.0);
else
    set(handles.eb_trainvideo, 'String', '');
    set(hObject, 'Value', 0.0);
end

% --- Executes on button press in btn_train.
function btn_train_Callback(hObject, eventdata, handles)
% hObject    handle to btn_train (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

videodir = get(handles.eb_trainvideo, 'String');
filtersize = 3;

st_multiple_validate_training(videodir, filtersize, handles, 0);
