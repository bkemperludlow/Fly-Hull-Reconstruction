function varargout = correctionGUI_sam(varargin)
%CORRECTIONGUI_SAM MATLAB code file for correctionGUI_sam.fig
%      CORRECTIONGUI_SAM, by itself, creates a new CORRECTIONGUI_SAM or raises the existing
%      singleton*.
%
%      H = CORRECTIONGUI_SAM returns the handle to a new CORRECTIONGUI_SAM or the handle to
%      the existing singleton*.
%
%      CORRECTIONGUI_SAM('Property','Value',...) creates a new CORRECTIONGUI_SAM using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to correctionGUI_sam_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CORRECTIONGUI_SAM('CALLBACK') and CORRECTIONGUI_SAM('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CORRECTIONGUI_SAM.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help correctionGUI_sam

% Last Modified by GUIDE v2.5 14-Jun-2019 15:27:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @correctionGUI_sam_OpeningFcn, ...
    'gui_OutputFcn',  @correctionGUI_sam_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
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


% --- Executes just before correctionGUI_sam is made visible.
function correctionGUI_sam_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
tmsCurr_default = 17.375 ; %-10 ;
saveEvery_default = 200000 ; %200
startdir_default = 'D:\Fly Data\VNC Sensory Lines\08_09042018\Analysis\Unsorted\Expr_8_mov_008\' ;%'D:\Fly Data\VNC Motor Lines\' ;

% assign default values if no inputs provided
switch nargin
    case 3
        handles.tmsCurr = tmsCurr_default ;
        handles.saveEvery = saveEvery_default ;
        handles.startdir = startdir_default ;
        
        handles.tmsCurr_default = tmsCurr_default ;
    case 4
        handles.tmsCurr = varargin{1} ;
        handles.saveEvery = saveEvery_default ;
        handles.startdir = startdir_default ;
        
        handles.tmsCurr_default = varargin{1} ;
        
    case 5
        handles.tmsCurr = varargin{1} ;
        handles.saveEvery = varargin{2} ;
        handles.startdir = startdir_default ;
        
        handles.tmsCurr_default = varargin{1} ;
    case 6
        handles.tmsCurr = varargin{1} ;
        handles.saveEvery = varargin{2} ;
        handles.startdir = varargin{3} ;
        
        handles.tmsCurr_default = varargin{1} ;
    otherwise
        disp('Wrong number of input arguments')
        keyboard
end
% Choose default command line output for correctionGUI_sam
handles.output = hObject;

% initialize plot handles
handles = initPlots(handles) ;

% try to deal with zoom/keypress conflict
handles.hManager = uigetmodemanager(handles.figure1);
try
    set(handles.hManager.WindowListenerHandles, 'Enable', 'off');  % HG1
catch
    [handles.hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2
end
set(handles.figure1, 'WindowKeyPressFcn',@figure1_WindowKeyPressFcn);
set(handles.figure1, 'KeyPressFcn', []);

% flag to let us know whether or not data is loaded
handles.dataFlag = false ;
handles.datapath_curr = [] ;
handles.dirlist = [] ;

% reference structures to make callbacks easier
adjustTagStruct = struct() ;
slider_handles = findall(0,'Style','Slider') ;
for i = 1:length(slider_handles)
    adjustTagStruct.(slider_handles(i).Tag) = i + 1 ;
end
handles.adjustTagStruct = adjustTagStruct ;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes correctionGUI_sam wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = correctionGUI_sam_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function figure1_WindowKeyPressFcn(hObject, eventdata)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
if handles.dataFlag
    switch eventdata.Key
        case 'v'
            uicontrol(handles.bback);
            bback_Callback(handles.bback,[],handles);
        case 'b'
            uicontrol(handles.back);
            back_Callback(handles.back,[],handles);
        case 'return'
            uicontrol(handles.fwd);
            fwd_Callback(handles.fwd,[],handles);
        case 'f'
            uicontrol(handles.ffwd);
            ffwd_Callback(handles.ffwd,[],handles);
    end
end

% ----------------------------------
% initialize plot handles
function handles = initPlots(handles)
% initialize handle fields for things to plot

scale = 4 ; % 2; % 4
handles.scale = scale ;
handles.Lstub = 3.0*scale ;
handles.Lbar  = 10*scale ;

handles.ignoreFrames = [] ;

% flags to mark whether or not things need to be updated
handles.bodyChangeFlag = false ;
handles.rightWingChangeFlag = false ;
handles.leftWingChangeFlag = false ;

handles.AhatChangeFlag = false ;
handles.stublineChangeFlag = false ;
handles.nrChangeFlag = false ;
handles.phChangeFlag = false ;

handles.spanRightChangeFlag = false ;
handles.chordRightChangeFlag = false ;
handles.spanLeftChangeFlag = false ;
handles.chordLeftChangeFlag = false ;

% intialize plot elements

hold(handles.main_axes, 'on')
handles.hBody = plot3(handles.main_axes, NaN, NaN, NaN, 'g.','MarkerSize',3) ;
handles.hRightWing = plot3(handles.main_axes, NaN, NaN, NaN, 'r.','MarkerSize',3) ;
handles.hLeftWing = plot3(handles.main_axes, NaN, NaN, NaN, 'b.','MarkerSize',3) ;

handles.hAhat = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
    'Color','r','LineWidth',8);
handles.hStubline = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
    'Color','k','LineWidth',8);
handles.hNr = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
    'Color','k','LineWidth',8);
handles.hPh = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
    'Color','r','LineWidth',8);

handles.hSpanRight = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
    'Color','k','LineWidth',4);
handles.hChordRight = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
    'Color','b','LineWidth',4);
handles.hSpanLeft = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
    'Color','k','LineWidth',4);
handles.hChordLeft = line(handles.main_axes,[NaN, NaN],[NaN, NaN],[NaN, NaN],...
    'Color','r','LineWidth',4);

hold(handles.main_axes, 'off')
box(handles.main_axes, 'on')
grid(handles.main_axes, 'on')
axis(handles.main_axes, 'equal')
rotate3d on

%initialize view
azview = -43 ;
elview =  26 ;
view(azview,elview) ;

% --------------------------------
% update plot to reflect changes
function updateDisplay(hObject)
handles = guidata(hObject);

%============================
% determine current frame
%============================
try
    frameCurr = handles.frameCurr ;
catch
    tmsCurr = handles.tmsCurr ;
    tvec = handles.tvec ;
    [~, frameCurr] = min(abs(tvec - tmsCurr)) ;
    handles.frameCurr = frameCurr ;
    handles.tmsCurr = tvec(frameCurr) ;
end

%
%============================
% get voxel info
%============================
row1 = handles.frameStartInd(frameCurr) ;
row2 = handles.frameEndInd(frameCurr) ;
coords = handles.res(row1:row2,2:4) ; % xyz
IDX = handles.RESIDX(row1:row2,:) ;  % (isbody, isrightwing, isleftwing)
bodyRows      = (IDX(:,handles.bodyInd)==1) ;
rightWingRows = (IDX(:,handles.rightWingInd)==1) ;
leftWingRows  = (IDX(:,handles.leftWingInd)==1) ;

scale = handles.scale ; % constant
Lbar = handles.Lbar ; % constant

%============================
% get body vector info
%============================
cb = handles.bodyCM(frameCurr,:) ; % body cm
rollHat = handles.rollHats(frameCurr,:) ; % body axis unit vector

checksum = sum(handles.normRolls(frameCurr,:));
if ((~isfinite(checksum) || norm(handles.normRolls(frameCurr,:))==0) && (frameCurr>1))
    handles.normRolls(frameCurr,:) = handles.normRolls(frameCurr-1,:) ;
    
    %disp('Using roll vector from previous frame') ;
end

% if still zero, start with zero roll
if (norm(handles.normRolls(frameCurr,:))==0)
    handles.normRolls(frameCurr,:) = - cross( rollHat, [0 0 1]) ;
    %disp('starting from rho=0') ;
end
normRoll = handles.normRolls(frameCurr,:) ;
% make sure the roll vector is perpendicular to AHat  and renormalize
normRoll = normRoll - rollHat * dot(normRoll, rollHat) ;
normRoll = normRoll / norm(normRoll) ;
handles.normRolls(frameCurr,:) = normRoll ;

stubvec  = handles.Lstub * RotatePoint(normRoll,[0 0 0],rollHat,90) ;
psiHat = handles.psiHats(frameCurr,:) ;

%============================
% get wing vector info
%============================
cr = handles.rightWingCM(frameCurr,:) ;
cl = handles.leftWingCM(frameCurr,:) ;

span1Hat = handles.span1Hats(frameCurr,:) ;
span2Hat = handles.span2Hats(frameCurr,:) ;
chord1Hat = handles.chord1Hats(frameCurr,:) ;
chord2Hat = handles.chord2Hats(frameCurr,:) ;
%--------------------------------------------------------------------------

%============================
% voxel plotting
%============================
if handles.bodyChangeFlag
    set(handles.hBody, 'XData', coords(bodyRows,1),...
        'YData', coords(bodyRows,2),...
        'ZData', coords(bodyRows,3));
end
if handles.rightWingChangeFlag
    set(handles.hRightWing, 'XData', coords(rightWingRows,1),...
        'YData', coords(rightWingRows,2),...
        'ZData', coords(rightWingRows,3));
end
if handles.leftWingChangeFlag
    set(handles.hLeftWing, 'XData', coords(leftWingRows,1),...
        'YData', coords(leftWingRows,2),...
        'ZData', coords(leftWingRows,3));
end

%============================
% body vector plotting
%============================
if handles.AhatChangeFlag
    x_data = [cb(1), cb(1)+ scale*15*rollHat(1)] ;
    y_data = [cb(2), cb(2)+ scale*15*rollHat(2)] ;
    z_data = [cb(3), cb(3)+ scale*15*rollHat(3)] ;
    
    set(handles.hAhat, 'XData', x_data,...
        'YData', y_data,...
        'ZData', z_data);
end
if handles.stublineChangeFlag
    
    x_data = [cb(1), cb(1)+ stubvec(1)] ;
    y_data = [cb(2), cb(2)+ stubvec(2)] ;
    z_data = [cb(3), cb(3)+ stubvec(3)] ;
    set(handles.hStubline, 'XData', x_data,...
        'YData', y_data,...
        'ZData', z_data);
end
if handles.nrChangeFlag
    x_data = [cb(1)+stubvec(1)-Lbar*normRoll(1), cb(1)+stubvec(1)+Lbar*normRoll(1)] ;
    y_data = [cb(2)+stubvec(2)-Lbar*normRoll(2), cb(2)+stubvec(2)+Lbar*normRoll(2)] ;
    z_data = [cb(3)+stubvec(3)-Lbar*normRoll(3), cb(3)+stubvec(3)+Lbar*normRoll(3)] ;
    set(handles.hNr, 'XData', x_data,...
        'YData', y_data,...
        'ZData', z_data);
end
if handles.phChangeFlag
    x_data = [cb(1),cb(1)+scale*7*psiHat(1)] ;
    y_data = [cb(2),cb(2)+scale*7*psiHat(2)] ;
    z_data = [cb(3),cb(3)+scale*7*psiHat(3)] ;
    set(handles.hPh, 'XData', x_data,...
        'YData', y_data,...
        'ZData', z_data);
end

%============================
% wing vector plotting
%============================
if handles.spanRightChangeFlag
    x_data = [cr(1), cr(1)+ scale*10*span1Hat(1)] ;
    y_data = [cr(2), cr(2)+ scale*10*span1Hat(2)] ;
    z_data = [cr(3), cr(3)+ scale*10*span1Hat(3)] ;
    
    set(handles.hSpanRight, 'XData', x_data,...
        'YData', y_data,...
        'ZData', z_data);
end
if handles.chordRightChangeFlag
    x_data = [cr(1), cr(1)+ scale*8*chord1Hat(1)] ;
    y_data = [cr(2), cr(2)+ scale*8*chord1Hat(2)] ;
    z_data = [cr(3), cr(3)+ scale*8*chord1Hat(3)] ;
    set(handles.hChordRight, 'XData', x_data,...
        'YData', y_data,...
        'ZData', z_data);
end
if handles.spanLeftChangeFlag
    x_data = [cl(1), cl(1)+ scale*10*span2Hat(1)] ;
    y_data = [cl(2), cl(2)+ scale*10*span2Hat(2)] ;
    z_data = [cl(3), cl(3)+ scale*10*span2Hat(3)] ;
    set(handles.hSpanLeft, 'XData', x_data,...
        'YData', y_data,...
        'ZData', z_data);
end
if handles.chordLeftChangeFlag
    x_data = [cl(1), cl(1)+ scale*8*chord2Hat(1)] ;
    y_data = [cl(2), cl(2)+ scale*8*chord2Hat(2)] ;
    z_data = [cl(3), cl(3)+ scale*8*chord2Hat(3)] ;
    set(handles.hChordLeft, 'XData', x_data,...
        'YData', y_data,...
        'ZData', z_data);
end

set(handles.frame_info,'String',...
    ['t(ms)=' num2str(handles.tmsCurr) '  Frame # ' ...
    num2str(handles.frameCurr) ' of ' num2str(handles.Nimages)])

guidata(hObject, handles);

% -------------------------------------
% remove current fly plot from handles
function clearDisplay(hObject)
handles = guidata(hObject);

x_data = [NaN, NaN] ;
y_data = [NaN, NaN] ;
z_data = [NaN, NaN] ;

% voxels
set(handles.hBody, 'XData', NaN, 'YData', NaN,'ZData', NaN);
set(handles.hRightWing, 'XData', NaN,'YData', NaN,'ZData', NaN);
set(handles.hLeftWing, 'XData', NaN,'YData', NaN, 'ZData', NaN);

% body vectors
set(handles.hAhat, 'XData', x_data,'YData', y_data,'ZData', z_data);
set(handles.hStubline, 'XData', x_data,'YData', y_data,'ZData', z_data);
set(handles.hNr, 'XData', x_data,'YData', y_data,'ZData', z_data);
set(handles.hPh, 'XData', x_data,'YData', y_data,'ZData', z_data);

% wing vectors
set(handles.hSpanRight, 'XData', x_data,'YData', y_data,'ZData', z_data);
set(handles.hChordRight, 'XData', x_data,'YData', y_data,'ZData', z_data);
set(handles.hSpanLeft, 'XData', x_data,'YData', y_data,'ZData', z_data);
set(handles.hChordLeft, 'XData', x_data,'YData', y_data, 'ZData', z_data);

set(handles.frame_info,'String','Frame Info')

guidata(hObject, handles);


function frame_info_Callback(hObject, eventdata, handles)
% hObject    handle to frame_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_info as text
%        str2double(get(hObject,'String')) returns contents of frame_info as a double


% --- Executes during object creation, after setting all properties.
function frame_info_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_info (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in data_dir.
function data_dir_Callback(hObject, eventdata, handles)
% hObject    handle to data_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns data_dir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from data_dir
index_selected = get(hObject,'Value');
contents = get(hObject,'String');
item_selected = contents{index_selected};

handles.datapath_curr = item_selected ;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function data_dir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in open_data_dir.
function open_data_dir_Callback(hObject, eventdata, handles)
% hObject    handle to open_data_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% -----------------------------------------------
% see if there's a set directory to start search
try
    startdir = handles.startdir;
catch
    startdir = pwd;
end
% ----------------------------------------------------------------------
% if there's not already data in the directory, just automatically load
if isfield(handles, 'dirlist')
    if ~isempty(handles.dirlist)
        loadFlag = false ;
    else
        loadFlag = true ;
    end
else
    loadFlag = true ;
end
% ---------------------------------------------
% use file picker to get files
tempdirlist = uipickfiles('filter',startdir, 'output', 'struct');
if ~isempty(tempdirlist)
    mat_file_ind = arrayfun(@(x) contains(x.name,'.mat') & ...
        contains(x.name,'Expr'),tempdirlist) ;
    tempdirlist = tempdirlist(mat_file_ind) ;
    tempdirlist = vertcat(handles.dirlist, tempdirlist) ;
    handles.dirlist = tempdirlist;
    set(handles.data_dir,'String',{tempdirlist(~[tempdirlist(:).isdir]).name})
    handles.datapath_curr = tempdirlist(1).name ;
end
if loadFlag
    handles = loadData(hObject, handles) ;
else
    guidata(hObject, handles);
end

% --- Executes on button press in load_data.
function load_data_Callback(hObject, eventdata, handles)
% hObject    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.dataFlag
    clearDisplay(hObject) ;
end
handles = loadData(hObject, handles) ;

% --- Executes on button press in clear_data.
function clear_data_Callback(hObject, eventdata, handles)
% hObject    handle to clear_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clearDisplay(hObject) ;  % reset plot window
handles.dataFlag = false ;      % signal that there's no more data
handles.datapath_curr = [] ;    % empty current data path

guidata(hObject, handles);

% ---------------------------------------------------------------------
% function to load in a fly analysis data file
function handles = loadData(hObject, handles)
if ~isempty(handles.datapath_curr)
    
    % load data file
    data = importdata(handles.datapath_curr) ;
    if isfield(data,'data')
        data = data.data ;
    end
    
    %================================================
    % add fields to handles for relevant plotting data
    %================================================
    % vectors for body/wing orientation
    handles.chord1Hats = data.rightChordHats ;
    handles.chord2Hats = data.leftChordHats ;
    handles.span1Hats  = data.rightSpanHats ;
    handles.span2Hats  = data.leftSpanHats ;
    handles.rollHats   = data.AHat ; % long body axis
    
    handles.Nimages     = data.Nimages ;
    handles.wingLength  = data.wingLength ;
    psiVec  = atan2(data.AHat(:,2),data.AHat(:,1)); % body angle with respect to x axis
    handles.psiHats  = [-sin(psiVec)  cos(psiVec)  zeros(data.Nimages,1)];
    
    phi1Vec   = atan2(data.rightSpanHats(:,2),data.rightSpanHats(:,1));
    handles.phi1Hats = [-sin(phi1Vec) cos(phi1Vec) zeros(data.Nimages,1)];
    
    phi2Vec   = atan2(data.leftSpanHats(:,2),data.leftSpanHats(:,1));
    handles.phi2Hats = [-sin(phi2Vec) cos(phi2Vec) zeros(data.Nimages,1)];
    
    % time information
    df = diff(data.res(:,1)) ;
    frameStartInd = [1 ; find(df==1)+1] ;
    handles.frameStartInd = frameStartInd ;
    handles.frameEndInd   = [frameStartInd(2:end)-1 ; size(data.res,1)] ;
    handles.noffset       = data.params.startTrackingTime - data.params.firstTrackableFrame  ;
    
    tvec = (0:data.Nimages-1) + data.params.startTrackingTime ;
    handles.tvec = tvec / 8000 * 1000 ; % in ms
    
    handles.tmsCurr = handles.tmsCurr_default ;
    
    % voxel info
    handles.res = data.res ;
    handles.RESIDX = data.RESIDX ;
    
    % cm info
    handles.bodyCM = data.bodyCM ;
    handles.rightWingCM = data.rightWingCM ;
    handles.leftWingCM = data.leftWingCM ;
    
    % wing tip info
    handles.rightWingTips = data.rightWingTips ;
    handles.leftWingTips = data.leftWingTips ;
    
    % indexing for body/wings
    handles.bodyInd = data.bodyInd ;
    handles.rightWingInd = data.rightWingInd ;
    handles.leftWingInd = data.leftWingInd ;
    
    %================================================
    % mark that all plot objects need to be updated
    %================================================
    handles = updateChangeFlags(handles,'all') ;
    
    % mark that data is currently loaded
    handles.dataFlag = true ;
    datapath_curr = handles.datapath_curr ;
    [filepath, name, ext] = fileparts(datapath_curr) ;
    if contains(name,'manually_corrected')
        handles.output_path = datapath_curr ;
    else
        name_split = strsplit(name, '_') ;
        ExprNumStr = name_split{2} ;
        MovNumStr = name_split{4} ;
        name_new = ['Expr' ExprNumStr 'mov' MovNumStr ...
            '_Data_manually_corrected' ext] ;
        handles.output_path = fullfile(filepath, name_new) ;
    end
    
    % initialize matrices for adjustments/swaps to data
    handles.adjust = zeros(data.Nimages,19);
    handles.rhoFlag   = false(data.Nimages,1);
    handles.pitchflipFlag = false(data.Nimages,1) ;
    handles.wingSwapFlag = false(data.Nimages,1) ;
    try
        handles.rhoTimes = data.rhoTimes ;
        handles.normRolls   = data.rollVectors ;
    catch
        handles.rhoTimes = [] ;
        handles.normRolls =  zeros(data.Nimages,3) ;
    end
    % enable sliders so that angles can be adjusted
    enableButtons()
else
    disp('No data file selected')
end
guidata(hObject, handles);
updateDisplay(hObject);

% ---------------------------------------------------------------------
% function to store whether or not a given fly objet was changed
function handles = updateChangeFlags(handles,updateType)
% flag all variables for update display (useful for frame shifts)
switch updateType
    case 'all'
        handles.bodyChangeFlag = true ;
        handles.rightWingChangeFlag = true ;
        handles.leftWingChangeFlag = true ;
        
        handles.AhatChangeFlag = true ;
        handles.stublineChangeFlag = true ;
        handles.nrChangeFlag = true ;
        handles.phChangeFlag = true ;
        
        handles.spanRightChangeFlag = true ;
        handles.chordRightChangeFlag = true ;
        handles.spanLeftChangeFlag = true ;
        handles.chordLeftChangeFlag = true ;
    case 'bodyvec'
        handles.AhatChangeFlag = true ;
        handles.stublineChangeFlag = true ;
        handles.nrChangeFlag = true ;
        handles.phChangeFlag = true ;
    case 'rvec'
        handles.spanRightChangeFlag = true ;
        handles.chordRightChangeFlag = true ;
    case 'lvec'
        handles.spanLeftChangeFlag = true ;
        handles.chordLeftChangeFlag = true ;
    case 'rchord'
        handles.chordRightChangeFlag = true ;
    case 'lchord'
        handles.chordLeftChangeFlag = true ;
end
%guidata(hObject, handles);

% ----------------------------------------------
% turn buttons on
function enableButtons()
%handles = guidata(hObject) ;
slider_handles = findall(0,'Style','Slider') ;
set(slider_handles,'Enable','on')
button_handles = findall(0,'Style','pushbutton') ;
set(button_handles,'Enable','on')
toggle_button_handles = findall(0,'Style','togglebutton') ;
set(toggle_button_handles,'Enable','on')

% ----------------------------------------
% return toggle buttons to default status
function resetToggleButtons()
%handles = guidata(hObject) ;
toggle_button_handles = findall(0,'Style','togglebutton') ;
set(toggle_button_handles,'Value',0) ;
defaultColor = get(0,'defaultUicontrolBackgroundColor') ;
set(toggle_button_handles,'BackgroundColor',defaultColor) ;

% -----------------------------------
% return sliders to default position
function resetSliders()
%handles = guidata(hObject) ;
slider_handles = findall(0,'Style','Slider') ;
set(slider_handles,'Value',0) ;

% ---------------------------------
% save *_manually_corrected file
function saveData(hObject)
% function to save manual correction results
handles = guidata(hObject) ;
f = waitbar(0,'') ;
f.Children.Title.Interpreter = 'none' ;
waitbar(0,f, ['Saving data to ' handles.output_path]) ;

% load data file
data = importdata(handles.datapath_curr) ;
if isfield(data,'data')
    data = data.data ;
end

waitbar(0.33)
% replace values in data structure with adjusted values
data.bodyCM = handles.bodyCM ;
data.rightWingCM = handles.rightWingCM ;
data.leftWingCM  = handles.leftWingCM ;

% vectors
data.rightChordHats = handles.chord1Hats ;
data.leftChordHats  = handles.chord2Hats ;
data.rightSpanHats  = handles.span1Hats ;
data.leftSpanHats   = handles.span2Hats ;
data.AHat           = handles.rollHats ;

data.rhoTimes       = handles.rhoTimes ;
data.rollVectors    = handles.normRolls ;
data.ignoreFrames   = handles.ignoreFrames ;

data.RESIDX         = handles.RESIDX ;

if isfield(data,'rightWingTips')
    temp1 = data.rightWingTips ;
    data.rightWingTips([handles.wingSwapFlag],:) = ...
        data.leftWingTips([handles.wingSwapFlag],:) ;
    data.leftWingTips([handles.wingSwapFlag],:) = ...
        temp1([handles.wingSwapFlag],:) ;
end

if isfield(handles, 'manualCorrRangeMS')
    data.manualCorrRangeMS = handles.manualCorrRangeMS ;
end
waitbar(0.66)

save(handles.output_path,'data')

waitbar(1)
close(f)


function [handles, adjust_val] = sliderAdjust(hObject, handles)
% generic function for data adjustment sliders
tag_curr = get(hObject,'Tag') ;
slider_val_curr = get(hObject,'Value') ;
slider_val_prev = ...
    handles.adjust(handles.frameCurr,handles.adjustTagStruct.(tag_curr)) ;
adjust_val = slider_val_curr - slider_val_prev ;
handles.adjust(handles.frameCurr,handles.adjustTagStruct.(tag_curr)) = ...
    slider_val_curr ;


% --- Executes on slider movement.
function body_x_Callback(hObject, eventdata,handles)
% hObject    handle to body_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.bodyCM(frameCurr,1) = handles.bodyCM(frameCurr,1) + adjust_val ;
handles = updateChangeFlags(handles,'bodyvec') ;
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes during object creation, after setting all properties.
function body_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to body_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -80);
set(hObject, 'Max', 80);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/160 , 10/160 ]);
set(hObject, 'Enable', 'off') ;

% --- Executes on slider movement.
function body_y_Callback(hObject, eventdata, handles)
% hObject    handle to body_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.bodyCM(frameCurr,2) = handles.bodyCM(frameCurr,2) + adjust_val ;
handles = updateChangeFlags(handles,'bodyvec') ;
guidata(hObject, handles);
updateDisplay(hObject);

function body_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider 23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -80);
set(hObject, 'Max', 80);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/160 , 10/160 ]);
set(hObject, 'Enable', 'off') ;

% --- Executes on slider movement.
function body_z_Callback(hObject, eventdata, handles)
% hObject    handle to body_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.bodyCM(frameCurr,3) = handles.bodyCM(frameCurr,3) + adjust_val ;
handles = updateChangeFlags(handles,'bodyvec') ;
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes during object creation, after setting all properties.
function body_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to body_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -80);
set(hObject, 'Max', 80);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/160 , 10/160 ]);
set(hObject, 'Enable', 'off') ;

% --- Executes on slider movement.
function body_phi_Callback(hObject, eventdata, handles)
% hObject    handle to body_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.rollHats(frameCurr,:) = ...
    RotatePoint(handles.rollHats(frameCurr,:),[0 0 0],[0 0 1],adjust_val);
handles.normRolls(frameCurr,:) = ...
    RotatePoint(handles.normRolls(frameCurr,:),[0 0 0],[0 0 1],adjust_val);
handles.psiHats(frameCurr,:) = ...
    RotatePoint(handles.psiHats(frameCurr,:),[0 0 0],[0 0 1],adjust_val);
handles = updateChangeFlags(handles,'bodyvec') ;
guidata(hObject, handles);
updateDisplay(hObject);


% --- Executes during object creation, after setting all properties.
function body_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to body_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -180);
set(hObject, 'Max', 180);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/360 , 10/360 ]);
set(hObject, 'Enable', 'off') ;

% --- Executes on slider movement.
function body_pitch_Callback(hObject, eventdata, handles)
% hObject    handle to body_pitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.rollHats(frameCurr,:) = ...
    RotatePoint2(handles.rollHats(frameCurr,:),[0 0 0],...
    handles.psiHats(frameCurr,:),adjust_val);
handles.normRolls(frameCurr,:) = ...
    RotatePoint2(handles.normRolls(frameCurr,:),[0 0 0],...
    handles.psiHats(frameCurr,:),adjust_val);
handles = updateChangeFlags(handles,'bodyvec') ;
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes during object creation, after setting all properties.
function body_pitch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to body_pitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -180);
set(hObject, 'Max', 180);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/360 , 10/360 ]);
set(hObject, 'Enable', 'off') ;

% --- Executes on slider movement.
function body_roll_Callback(hObject, eventdata, handles)
% hObject    handle to body_roll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.normRolls(frameCurr,:) = ...
    RotatePoint(handles.normRolls(frameCurr,:),[0 0 0],...
    handles.rollHats(frameCurr,:),-1*adjust_val);
handles = updateChangeFlags(handles,'bodyvec') ;
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes during object creation, after setting all properties.
function body_roll_CreateFcn(hObject, eventdata, handles)
% hObject    handle to body_roll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -180);
set(hObject, 'Max', 180);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/360 , 10/360 ]);
set(hObject, 'Enable', 'off') ;

% --- Executes on slider movement.
function rx_Callback(hObject, eventdata, handles)
% hObject    handle to rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.rightWingCM(frameCurr,1) = handles.rightWingCM(frameCurr,1) + adjust_val ;
handles = updateChangeFlags(handles,'rvec') ;
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes during object creation, after setting all properties.
function rx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -80);
set(hObject, 'Max', 80);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/160 , 10/160 ]);
set(hObject, 'Enable', 'off') ;

% --- Executes on slider movement.
function r_y_Callback(hObject, eventdata, handles)
% hObject    handle to r_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.rightWingCM(frameCurr,2) = handles.rightWingCM(frameCurr,2) + adjust_val ;
handles = updateChangeFlags(handles,'rvec') ;
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes during object creation, after setting all properties.
function r_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -80);
set(hObject, 'Max', 80);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/160 , 10/160 ]);
set(hObject, 'Enable', 'off') ;

% --- Executes on slider movement.
function r_z_Callback(hObject, eventdata, handles)
% hObject    handle to r_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.rightWingCM(frameCurr,3) = handles.rightWingCM(frameCurr,3) + adjust_val ;
handles = updateChangeFlags(handles,'rvec') ;
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes during object creation, after setting all properties.
function r_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -80);
set(hObject, 'Max', 80);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/160 , 10/160 ]);
set(hObject, 'Enable', 'off') ;

% --- Executes on slider movement.
function r_phi_Callback(hObject, eventdata, handles)
% hObject    handle to r_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;

handles.span1Hats(frameCurr,:) = ...
    RotatePoint(handles.span1Hats(frameCurr,:),[0 0 0],[0 0 1],adjust_val);
handles.chord1Hats(frameCurr,:) = ...
    RotatePoint(handles.chord1Hats(frameCurr,:),[0 0 0],[0 0 1],adjust_val);
%chord1Hat2s(i,:) = RotatePoint(chord1Hat2s(i,:),[0 0 0],[0 0 1],-10);
handles.phi1Hats(frameCurr,:) = ...
    RotatePoint(handles.phi1Hats(frameCurr,:),[0 0 0],[0 0 1],adjust_val);
handles = updateChangeFlags(handles,'rvec') ;
guidata(hObject, handles);
updateDisplay(hObject);


% --- Executes during object creation, after setting all properties.
function r_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -90);
set(hObject, 'Max', 90);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/180 , 10/180 ]);
set(hObject, 'Enable', 'off') ;

% --- Executes on slider movement.
function r_theta_Callback(hObject, eventdata, handles)
% hObject    handle to r_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.span1Hats(frameCurr,:) = ...
    RotatePoint2(handles.span1Hats(frameCurr,:),[0 0 0],...
    handles.phi1Hats(frameCurr,:),-1*adjust_val);
handles.chord1Hats(frameCurr,:) = ...
    RotatePoint2(handles.chord1Hats(frameCurr,:),[0 0 0],...
    handles.phi1Hats(frameCurr,:),-1*adjust_val);
%chord1Hat2s(i,:) = RotatePoint2(chord1Hat2s(i,:),[0 0 0],phi1Hats(i,:),10);
handles = updateChangeFlags(handles,'rvec') ;
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes during object creation, after setting all properties.
function r_theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -90);
set(hObject, 'Max', 90);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/180 , 10/180 ]);
set(hObject, 'Enable', 'off') ;

% --- Executes on slider movement.
function r_eta_Callback(hObject, eventdata, handles)
% hObject    handle to r_eta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.chord1Hats(frameCurr,:) = ...
    RotatePoint(handles.chord1Hats(frameCurr,:),[0 0 0],...
    handles.span1Hats(frameCurr,:),-1*adjust_val);
handles = updateChangeFlags(handles,'rchord') ;
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes during object creation, after setting all properties.
function r_eta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_eta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -180);
set(hObject, 'Max', 180);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/360 , 10/360 ]);
set(hObject, 'Enable', 'off') ;

% --- Executes on slider movement.
function l_x_Callback(hObject, eventdata, handles)
% hObject    handle to l_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.leftWingCM(frameCurr,1) = handles.leftWingCM(frameCurr,1) + adjust_val ;
handles = updateChangeFlags(handles,'lvec') ;
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes during object creation, after setting all properties.
function l_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to l_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -80);
set(hObject, 'Max', 80);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/160 , 10/160 ]);
set(hObject, 'Enable', 'off') ;


% --- Executes on slider movement.
function l_y_Callback(hObject, eventdata, handles)
% hObject    handle to l_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.leftWingCM(frameCurr,2) = handles.leftWingCM(frameCurr,2) + adjust_val ;
handles = updateChangeFlags(handles,'lvec') ;
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes during object creation, after setting all properties.
function l_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to l_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -80);
set(hObject, 'Max', 80);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/160 , 10/160 ]);
set(hObject, 'Enable', 'off') ;


% --- Executes on slider movement.
function l_z_Callback(hObject, eventdata, handles)
% hObject    handle to l_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.leftWingCM(frameCurr,3) = handles.leftWingCM(frameCurr,3) + adjust_val ;
handles = updateChangeFlags(handles,'lvec') ;
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes during object creation, after setting all properties.
function l_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to l_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -80);
set(hObject, 'Max', 80);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/160 , 10/160 ]);
set(hObject, 'Enable', 'off') ;


% --- Executes on slider movement.
function l_phi_Callback(hObject, eventdata, handles)
% hObject    handle to l_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;

handles.span2Hats(frameCurr,:) = ...
    RotatePoint(handles.span2Hats(frameCurr,:),[0 0 0],[0 0 1],-1*adjust_val);
handles.chord2Hats(frameCurr,:) = ...
    RotatePoint(handles.chord2Hats(frameCurr,:),[0 0 0],[0 0 1],-1*adjust_val);
%chord2Hat2s(i,:) = RotatePoint(chord2Hat2s(i,:),[0 0 0],[0 0 1],-10);
handles.phi2Hats(frameCurr,:) = ...
    RotatePoint(handles.phi2Hats(frameCurr,:),[0 0 0],[0 0 1],-1*adjust_val);
handles = updateChangeFlags(handles,'lvec') ;
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes during object creation, after setting all properties.
function l_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to l_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -90);
set(hObject, 'Max', 90);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/180 , 10/180 ]);
set(hObject, 'Enable', 'off') ;

% --- Executes on slider movement.
function l_theta_Callback(hObject, eventdata, handles)
% hObject    handle to l_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.span2Hats(frameCurr,:) = ...
    RotatePoint2(handles.span2Hats(frameCurr,:),[0 0 0],...
    handles.phi2Hats(frameCurr,:),-1*adjust_val);
handles.chord2Hats(frameCurr,:) = ...
    RotatePoint2(handles.chord2Hats(frameCurr,:),[0 0 0],...
    handles.phi2Hats(frameCurr,:),-1*adjust_val);
%chord2Hat2s(i,:) = RotatePoint2(chord2Hat2s(i,:),[0 0 0],phi2Hats(i,:),10);
handles = updateChangeFlags(handles,'lvec') ;
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes during object creation, after setting all properties.
function l_theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to l_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -90);
set(hObject, 'Max', 90);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/180 , 10/180 ]);
set(hObject, 'Enable', 'off') ;

% --- Executes on slider movement.
function l_eta_Callback(hObject, eventdata, handles)
% hObject    handle to l_eta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
frameCurr = handles.frameCurr ;
[handles, adjust_val] = sliderAdjust(hObject, handles) ;
handles.chord2Hats(frameCurr,:) = ...
    RotatePoint(handles.chord2Hats(frameCurr,:),[0 0 0],...
    handles.span2Hats(frameCurr,:),-1*adjust_val);
handles = updateChangeFlags(handles,'lchord') ;
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes during object creation, after setting all properties.
function l_eta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to l_eta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', -180);
set(hObject, 'Max', 180);
set(hObject, 'Value', 0);
set(hObject, 'SliderStep', [1/360 , 10/360 ]);
set(hObject, 'Enable', 'off') ;

% --------------------------------------------------------------------
function ui_save_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
saveData(hObject)

function main_axes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to main_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate main_axes


% --- Executes on button press in bback.
function bback_Callback(hObject, eventdata, handles)
% hObject    handle to bback (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.frameCurr = handles.frameCurr - 10 ;
handles.tmsCurr = handles.tvec(handles.frameCurr) ;
handles = updateChangeFlags(handles,'all') ;
resetToggleButtons()
resetSliders()
if mod(handles.frameCurr,handles.saveEvery) == 0
    saveData(hObject)
end
guidata(hObject, handles);
updateDisplay(hObject);
%resetToggleButtons(h)

% --- Executes on button press in back.
function back_Callback(hObject, eventdata, handles)
% hObject    handle to back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.frameCurr = handles.frameCurr - 1 ;
handles.tmsCurr = handles.tvec(handles.frameCurr) ;
handles = updateChangeFlags(handles,'all') ;
resetToggleButtons()
resetSliders()
if mod(handles.frameCurr,handles.saveEvery) == 0
    saveData(hObject)
end
guidata(hObject, handles);
updateDisplay(hObject);

% --- Executes on button press in fwd.
function fwd_Callback(hObject, eventdata, handles)
% hObject    handle to fwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.frameCurr >= handles.Nimages
    saveData(hObject)
    disp('Movie completed!')
    return
end
handles.frameCurr = handles.frameCurr + 1 ;
handles.tmsCurr = handles.tvec(handles.frameCurr) ;
handles = updateChangeFlags(handles,'all') ;
resetToggleButtons()
resetSliders()
if mod(handles.frameCurr,handles.saveEvery) == 0
    saveData(hObject)
end
guidata(hObject, handles);
updateDisplay(hObject);


% --- Executes on button press in ffwd.
function ffwd_Callback(hObject, eventdata, handles)
% hObject    handle to ffwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.frameCurr + 10) >= handles.Nimages
    saveData(hObject)
    disp('Not enough remaining frames to skip forward')
    return
end
handles.frameCurr = handles.frameCurr + 10 ;
handles.tmsCurr = handles.tvec(handles.frameCurr) ;
handles = updateChangeFlags(handles,'all') ;
resetToggleButtons()
resetSliders()
if mod(handles.frameCurr,handles.saveEvery) == 0
    saveData(hObject)
end
guidata(hObject, handles);
updateDisplay(hObject);


% --- Executes on button press in ignoreFrame.
function ignoreFrame_Callback(hObject, eventdata, handles)
% hObject    handle to ignoreFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button_state = get(hObject,'Value');
frameCurr = handles.frameCurr ;

if button_state == get(hObject,'Max')
    handles.ignoreFrames = sort(unique([handles.ignoreFrames, frameCurr])) ;
    set(hObject,'BackgroundColor','red')
elseif button_state == get(hObject,'Min')
    handles.ignoreFrames = ...
        sort(unique(handles.ignoreFrames(handles.ignoreFrames ~= frameCurr))) ;
    defaultColor = get(0,'defaultUicontrolBackgroundColor') ;
    set(hObject,'BackgroundColor',defaultColor)
end
guidata(hObject, handles);


% --- Executes on button press in save_roll.
function save_roll_Callback(hObject, eventdata, handles)
% hObject    handle to save_roll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
button_state = get(hObject,'Value');
frameCurr = handles.frameCurr ;

if button_state == get(hObject,'Max')
    handles.rhoFlag(frameCurr) = true ;
    handles.rhoTimes = sort(unique([handles.rhoTimes, frameCurr])) ;
    set(hObject,'BackgroundColor','green')
elseif button_state == get(hObject,'Min')
    handles.rhoFlag(frameCurr) = false ;
    handles.rhoTimes = ...
        sort(unique(handles.rhoTimes(handles.rhoTimes ~= frameCurr))) ;
    defaultColor = get(0,'defaultUicontrolBackgroundColor') ;
    set(hObject,'BackgroundColor',defaultColor)
end
guidata(hObject, handles);


% --- Executes on button press in swap_wings.
function swap_wings_Callback(hObject, eventdata, handles)
% hObject    handle to swap_wings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
frameCurr = handles.frameCurr ;

% swap voxels
row1 = handles.frameStartInd(frameCurr) ;
row2 = handles.frameEndInd(frameCurr) ;
IDX = handles.RESIDX(row1:row2,:) ;  % (isbody, isrightwing, isleftwing)
rightWingRows = (IDX(:,handles.rightWingInd)==1) ;
leftWingRows  = (IDX(:,handles.leftWingInd)==1) ;
allRows = row1:row2 ;

handles.RESIDX(allRows(rightWingRows),handles.rightWingInd) = 0 ;
handles.RESIDX(allRows(rightWingRows),handles.leftWingInd) = 1 ;
handles.RESIDX(allRows(leftWingRows),handles.rightWingInd) = 1 ;
handles.RESIDX(allRows(leftWingRows),handles.leftWingInd) = 0 ;

% swap vectors
span1Hats = handles.span1Hats ;
span2Hats = handles.span2Hats ;
chord1Hats = handles.chord1Hats ;
chord2Hats = handles.chord2Hats ;
phi1Hats = handles.phi1Hats ;
phi2Hats = handles.phi2Hats ;

handles.span1Hats(frameCurr,:) = span2Hats(frameCurr,:) ;
handles.chord1Hats(frameCurr,:) = chord2Hats(frameCurr,:) ;
handles.phi1Hats(frameCurr,:) = phi2Hats(frameCurr,:) ;

handles.span2Hats(frameCurr,:) = span1Hats(frameCurr,:) ;
handles.chord2Hats(frameCurr,:) = chord1Hats(frameCurr,:) ;
handles.phi2Hats(frameCurr,:) = phi1Hats(frameCurr,:) ;

% swap cm
rightWingCM = handles.rightWingCM ;
leftWingCM = handles.leftWingCM ;

handles.rightWingCM(frameCurr,:) = leftWingCM(frameCurr,:) ;
handles.leftWingCM(frameCurr,:) = rightWingCM(frameCurr,:) ;

% update structure and plot
handles.wingSwapFlag = ~handles.wingSwapFlag ;
handles = updateChangeFlags(handles,'all') ;
guidata(hObject, handles);
updateDisplay(hObject);


% --- Executes on button press in cluster_wings.
function cluster_wings_Callback(hObject, eventdata, handles)
% hObject    handle to cluster_wings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% -----------------------------------
% read in initial data
frameCurr = handles.frameCurr ;
N_vox_min = 200 ;

row1 = handles.frameStartInd(frameCurr) ;
row2 = handles.frameEndInd(frameCurr) ;
IDX = handles.RESIDX(row1:row2,:) ;  % (isbody, isrightwing, isleftwing)
rightWingRows = (IDX(:,handles.rightWingInd)==1) ;
leftWingRows  = (IDX(:,handles.leftWingInd)==1) ;

N_vox_R = sum(rightWingRows) ; % number of voxels in each wing
N_vox_L = sum(leftWingRows) ;

% --------------------------------------------------------
% smooth/interpolate wing CMs and tips
% 
[rightWingCM_interp, ~, rightWingTips_interp] = ...
    interpolateWingCM(handles,'right') ;
[leftWingCM_interp, ~, leftWingTips_interp] = ...
    interpolateWingCM(handles,'left') ;
% --------------------------------------------------------
% cluster the voxels for whichever wing has more voxels
if (N_vox_R > N_vox_L) && (N_vox_R > N_vox_min)
    wing_str = 'right' ;
elseif (N_vox_R < N_vox_L) && (N_vox_L > N_vox_min)
    wing_str = 'left' ;
else
    disp('Uncertain which wing to cluster')
    return
end
% -------------------------------------------------------------------
% perform clustering
[wingVox, wingRows, label_idx, centroids, badClusterFlag] = ...
    clusterWings(handles, frameCurr, wing_str) ;

% ------------------------------------------------------------------
% did clustering work?
if ~badClusterFlag
    % if we did successfully cluster, we now need to figure out
    % which is left vs right. we'll do this by checking against the
    % interpolated center of mass
    right_centroid_dist = myNorm(centroids -  ...
        repmat(rightWingCM_interp(frameCurr,:),2,1)) ;
    left_centroid_dist = myNorm(centroids -  ...
        repmat(leftWingCM_interp(frameCurr,:),2,1)) ;
    [~, right_idx] = min(right_centroid_dist) ;
    [~, left_idx] = min(left_centroid_dist) ;
    
    % in case both blobs think they should be the same wing
    if (right_idx == left_idx) 
        disp('could not determine left vs right cluster')
    else
        % otherwise, assign the new wing data to our output struct
        wingVoxR = wingVox(label_idx == right_idx,:) ;
        rightCM = centroids(right_idx,:) ;
        wingRowsR = wingRows(:,right_idx) ;
        wingVoxL = wingVox(label_idx == left_idx,:) ;
        leftCM = centroids(left_idx,:) ;
        wingRowsL = wingRows(:,left_idx) ;
        
        % now calculate wing vectors
        rightRefVecs = [handles.bodyCM(frameCurr,:); ...
            handles.bodyCM(frameCurr-1,:); ...
            rightWingTips_interp(frameCurr-1,:)] ;
        [spanHatR, chordHatR, chordAltHatR, ~, wingTipR] = ...
            estimate_wing_vecs(wingVoxR, rightRefVecs, handles.wingLength, ...
            [], [], rightCM) ;
        
        leftRefVecs = [handles.bodyCM(frameCurr,:) ; ...
            handles.bodyCM(frameCurr-1,:); ...
            leftWingTips_interp(frameCurr-1,:)] ;
        [spanHatL, chordHatL, chordAltHatL, ~, wingTipL] = ...
            estimate_wing_vecs(wingVoxL, leftRefVecs, handles.wingLength, ...
            [], [], leftCM) ;
        
        % add to storage arrays
        handles.rightWingCM(frameCurr,:) = rightCM ;
        handles.span1Hats(frameCurr,:) = spanHatR ;
        handles.chord1Hats(frameCurr,:) = chordHatR ;
        %handles.rightChordAltHats(frameCurr,:) = chordAltHatR ;
        handles.rightWingTips(frameCurr,:) = wingTipR ;
        
        handles.leftWingCM(frameCurr,:) = leftCM ;
        handles.span2Hats(frameCurr,:) = spanHatL ;
        handles.chord2Hats(frameCurr,:) = chordHatL ;
        %handles.leftChordAltHats(frameCurr,:) = chordAltHatL ;
        handles.leftWingTips(frameCurr,:) = wingTipL ;
        
        %...including fucking voxels, blurg
        handles.RESIDX(row1:row2,2:3) = [wingRowsR, wingRowsL] ;
        
        % now update display and save changes to data
        handles = updateChangeFlags(handles,'all') ;
        guidata(hObject, handles);
        updateDisplay(hObject);
    end
else
    disp('clustering failed')
end



% --- Executes on button press in save_manualCorrRangeMS.
function save_manualCorrRangeMS_Callback(hObject, eventdata, handles)
% hObject    handle to save_manualCorrRangeMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prompt = {'Enter manual correction start time (ms):',...
    'Enter manual correction end time (ms):'};
dlgtitle = 'Save Manual correction range';
dims = [1 40];
definput = {'-10','30'};
userInput = inputdlg(prompt,dlgtitle,dims,definput) ;

if ~isempty(userInput)
    manualCorrRangeMS = cellfun(@(y) str2double(y), userInput) ;
    if size(manualCorrRangeMS,1) > size(manualCorrRangeMS,2)
        manualCorrRangeMS = manualCorrRangeMS' ;
    end
    handles.manualCorrRangeMS = manualCorrRangeMS ;
else
    handles.manualCorrRangeMS = [] ;
end
guidata(hObject, handles);