% -------------------------------------------------------------------------
% function to get cine file movie info from xml file
%
% NB: this is really just for convenience. parseXML.m will give you all the
% info, but the field structure is complicated
% -------------------------------------------------------------------------
function info_struct = getInfoXML(xml_filename)
% ------------------------------
% initialize output
info_struct = struct ; 
info_struct.filename = xml_filename ; 

% -------------------------------
% convert xml into mat struct
xml_struct = parseXML(xml_filename) ; 

% ---------------------------------------------------------------------
% grab info from struct. right now doing trigger time, frame rate, and
% first im, but should probably expand upon this
% TRIGGER TIME
info_struct.triggerTime = ...
    datenum(xml_struct.chd.CineFileHeader.TriggerTime.Time.Text(1:12),...
        'HH:MM:SS.FFF') ; 
    
% FRAME RATE
if isfield(xml_struct.chd.CameraSetup, 'FrameRate')
    info_struct.frameRate = ...
        str2double(xml_struct.chd.CameraSetup.FrameRate.Text)  ;
elseif isfield(xml_struct.chd.CameraSetup, 'FrameRateDouble')
    info_struct.frameRate = ...
        str2double(xml_struct.chd.CameraSetup.FrameRateDouble.Text)  ;
else
    fprintf('Could not locate frame rate in xml info \n')
    keyboard
end

% FIRST IMAGE
info_struct.firstIm = ...
    str2double(xml_struct.chd.CineFileHeader.FirstMovieImage.Text) ;

end