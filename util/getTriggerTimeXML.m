% -------------------------------------------------------------------------
% quick/dirty function to get movie trigger time from xml file
% 
% NB: this assumes a consistent xml file structure. use "parseXML.m" for a 
% more general (but slower) file read
% -------------------------------------------------------------------------
function triggerTime = getTriggerTimeXML(filename)
fprintf('Under construction! \n')
keyboard

% load xml file
rootNode = xmlread(filename) ; 

% progress through child nodes
chdNode = rootNode.getChildNodes.item(1) ;
cineFileHeaderNode = chdNode.getChildNodes.item(0) ; 
triggerTimeNode = cineFileHeaderNode.getChildNodes.item(7) ; 
triggerTimeAttributes = triggerTimeNode.getAttributes ; 

% childNodes = rootNode.getChildNodes ; 
%chdNode = childNodes.item(1) ; 

end