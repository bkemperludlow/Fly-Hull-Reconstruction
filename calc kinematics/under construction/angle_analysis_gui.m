function varargout = angle_analysis_gui(varargin)
% PP_GUI MATLAB code for pp_gui.fig
%      PP_GUI, by itself, creates a new PP_GUI or raises the existing
%      singleton*.
%
%      H = PP_GUI returns the handle to a new PP_GUI or the handle to
%      the existing singleton*.
%
%      PP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PP_GUI.M with the given input arguments.
%
%      PP_GUI('Property','Value',...) creates a new PP_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pp_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pp_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pp_gui

% Last Modified by GUIDE v2.5 28-Jun-2016 10:47:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @angle_analysis_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @angle_analysis_gui_OutputFcn, ...
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

%---------------------------------------------------------------------------

function pp_gui_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pp_gui (see VARARGIN)

% Choose default command line output for pp_gui
handles.output = hObject;

axes(handles.axes1); zoom on
axes(handles.axes2); zoom on
axes(handles.axes3); zoom on
axes(handles.axes4); zoom on
axes(handles.axes5); zoom on
axes(handles.axes6); zoom on

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pp_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end
%---------------------------------------------------------------------------

% --- Outputs from this function are returned to the command line.
function varargout = angle_analysis_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

%---------------------------------------------------------------------------

% --- Executes on button press in selectdays.
function selectdata_Callback(hObject, eventdata, handles)
% hObject    handle to selectdays (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    startdir = handles.startdir;
catch
    startdir = 'Z:\';
end
%fileType = { '*.m',   'M-files'   ; '*.mat', 'MAT-files' }
tempdirlist = uipickfiles('filter',startdir, 'output', 'struct','NumFiles',1);
if ~isempty(tempdirlist)
    handles.dirlist = tempdirlist;
try
    %update label of days being plotted;
    labeltxt = '';
    dirlist = handles.dirlist;
    for i = 1:length(dirlist)
        namestr = dirlist(i).name;
        namestr = strsplit(namestr, '\');
        newname = [namestr{end-2},' - ',namestr{end}];
        labeltxt = [labeltxt, newname, ','];
    end
    labeltxt = labeltxt(1:end-1);
    set(handles.daystoplotlabel, 'String', labeltxt);
    startdir = dirlist(end).name; startdir = strsplit(startdir, '\'); 
    startdir = startdir(1:end-1); startdir = strjoin(startdir, '\');
    handles.startdir = startdir;
   
    data = importdata(dirlist) ; 
    %[statslist, dates] = load_stats(dirlist, 0, to_stop,'pellet_count', ...
    %    'srate', 'numtraj');
    %disp(length(statslist));
    %pellets = 0; trialnum = 0;
    %text = {};
    %{
    for i = 1:length(statslist);
        stats = statslist(i);
        try
        pc = [dates{i},' pellets: ', num2str(stats.pellet_count)];
        sr = [dates{i},' success rate: ', num2str(stats.srate.total)];
        sr_l = [dates{i},' success rate_l: ', num2str(stats.srate.laser_succ)];
        sr_nl = [dates{i},' success rate_nl: ', num2str(stats.srate.catch_succ)];
        nt = [dates{i},' num trials: ', num2str(stats.numtraj)];
        tmptext = {pc; sr; sr_l; sr_nl; nt};
        text = [text; tmptext];
        pellets = pellets + stats.pellet_count;
        trialnum = trialnum + stats.trialnum;
        end
    end
    
    handles = update_console(handles, text);
    if length(statslist) > 1
        pc = ['Total pellets: ', num2str(pellets)];
        sr = ['Overall success rate: ', num2str(pellets/trialnum)];
        text = {pc; sr};
        handles = update_console(handles, text);
    end
    %}
catch e
    disp(getReport(e));
end
guidata(hObject, handles);
end
end