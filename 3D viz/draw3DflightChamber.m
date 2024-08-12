% -------------------------------------------------------------------------
%  function to draw 3D version of flight chamber for creating figure
%  schematics. NB: was going to put opto stuff here, but i think that's
%  best in different function
%{
fig = figure('Position',[2239, 346, 560, 420])  ;
hold on
ax = gca ; 
parent = hgtransform('Parent',ax);
%}
% -------------------------------------------------------------------------
function grp = draw3DflightChamber(ax, parent, sideLength, crossSectionLength,...
    lineScale, center, rgbFrame, rgbWall, lineStyleStr, objAlpha, resolution, ...
    coilFlag)
% -----------------------
%% params and inputs
if ~exist('center', 'var') || isempty(center)
   center = [0, 0, 0] ;  
end
if ~exist('rgbFrame', 'var') || isempty(rgbFrame)
   rgbFrame = 0.7*[1 1 1] ;  
end
if ~exist('rgbWall', 'var') || isempty(rgbWall)
   rgbWall = [0.6196, 0.7922, 0.8824] ; % light blue 
end
if ~exist('lineStyleStr', 'var') || isempty(lineStyleStr)
   lineStyleStr = '-' ;  
end
if ~exist('objAlpha', 'var') || isempty(objAlpha)
   objAlpha = 1.0 ;  
end
if ~exist('resolution', 'var') || isempty(resolution)
   resolution = 100 ;  
end
if ~exist('coilFlag', 'var') || isempty(coilFlag)
   coilFlag = true ;  
end

lineWidth = lineScale * 0.5 ; % 0.5 is default line thickness, i think
% -------------------------
% wall properties:
wallAlpha = 0.1 ; 

% -------------------------
% coil properties:
coilOuterRadius = sideLength/2 - crossSectionLength ; 
coilInnerRadius = sideLength/2 - 2*crossSectionLength ; 
cylinderRadius = 0.95*coilOuterRadius ; 
rgbCoilMetal = 0.8*[1 1 1] ; 
rgbCoilWire = [184, 115, 51]/255 ;  % copper
coilBorderLW = lineWidth ; 

% -----------------------
% properties of input axis
axes(ax) ;
wasOnHold = ishold ;
hold on ;

% ---------------------------------------------------
%% create box (walls). centered at (0,0,0)
Vertices = sideLength*[0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1] - ...
    (sideLength/2)*[1, 1, 1] ;
Faces = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];

hWall =  patch('Vertices',Vertices,'Faces',Faces,...
    'FaceVertexCData',repmat(rgbWall,size(Faces,1), 1),'FaceColor','flat',...
    'lineStyle',lineStyleStr,'LineWidth',lineWidth) ;
set(hWall,'FaceAlpha',wallAlpha,'EdgeAlpha',wallAlpha)

grp = hgtransform('Parent',parent) ;
wall_grp = hgtransform('Parent',grp);
set(hWall, 'Parent',wall_grp) ;

% change lighting/material properties on walls
material(hWall,'shiny') ; 
% ---------------------------------------------------
%% create frame around walls
% first make square framed got top and bottom. accomplished via loop

translation_list = ...
    (sideLength/2).*[ 1, 0, -1 ; 0, 1, -1 ; -1, 0, -1 ; 0, -1, -1 ; ... % bottom pieces
     1, 0, 1 ; 0, 1, 1 ; -1, 0, 1 ; 0, -1, 1 ; ...                      % top pieces
     1, 1, 0 ; 1, -1, 0 ; -1, -1, 0; -1, 1, 0] ;                        % vertical pieces
N_frame_pieces = 12 ; 
h_frame_array = gobjects(1, N_frame_pieces) ; 

% loop through and make each piece of frame
for i = 1:N_frame_pieces
    if (mod(i,2) == 1) && (i <= 8)
        len_curr = crossSectionLength ;
        wid_curr = sideLength - crossSectionLength ;
        dep_curr = crossSectionLength ;
    elseif (mod(i,2) == 0) && (i <= 8)
        len_curr = sideLength + crossSectionLength ;
        wid_curr = crossSectionLength ;
        dep_curr = crossSectionLength ;
    else
        len_curr = crossSectionLength ; 
        wid_curr = crossSectionLength ; 
        dep_curr = sideLength - crossSectionLength ; 
    end
    
     vert_curr = [0 0 0 ; len_curr 0 0 ; len_curr wid_curr 0 ; 0 wid_curr 0;...
         0 0 dep_curr ; len_curr 0 dep_curr ; len_curr wid_curr dep_curr; ....
         0 wid_curr dep_curr] - 0.5*[len_curr, wid_curr, dep_curr] + ...
         translation_list(i,:); 
     
     h_frame_array(i) = patch('Vertices',vert_curr,'Faces',Faces,...
        'FaceVertexCData',repmat(rgbFrame,size(Faces,1), 1),...
        'FaceColor','flat', 'lineStyle',lineStyleStr,...
        'LineWidth',lineWidth) ;
     set(h_frame_array(i),'FaceAlpha',objAlpha,'EdgeAlpha',objAlpha)
end

% make all frame pieces into one group
frame_grp = hgtransform('Parent',grp);
set(h_frame_array(:), 'Parent',frame_grp) ;

% change lighting/material properties on frame
material(h_frame_array(:),'metal')
set(h_frame_array(:), 'ambientStrength',.9);  
% -----------------------------------------------------------------------
%% add in helmholtz coils, if desired
% each coil consists of two annuli stacked vertically, with a patch in
% between to represent the wire. so four annuli (8 circles) total for the 
% full coil setup
if coilFlag
    % --------------------------------------------------
    % first draw top annulus on bottom face of chamber
    t = linspace(0,2*pi, resolution) ;
    z_list = [-(sideLength/2 + crossSectionLength) ; ...
        (sideLength/2 + crossSectionLength)] ; 
    dz_list = [-1*crossSectionLength ; -crossSectionLength] ; 
    h_annulus_array = gobjects(1,4) ; 
    h_cyl_array = gobjects(1,2) ; 
    h_line_array = gobjects(1,8) ; 
    
    for j = 1:2
        % outer ring
        Xo1 = coilOuterRadius * cos(t) ;
        Yo1 = coilOuterRadius * sin(t) ;
        Zo1 = z_list(j)*ones(size(Xo1))  ;
        %Xo1 = [Xo1, Xo1(1)] ; Yo1 = [Yo1, Yo1(1)] ; Zo1 = [Zo1, Zo1(1)] ;
        % inner ring
        Xi1 = coilInnerRadius * cos(t) ;
        Yi1 = coilInnerRadius * sin(t) ;
        Zi1 = z_list(j)*ones(size(Xi1))  ;
        %Xi1 = [Xi1, Xi1(1)] ; Yi1 = [Yi1, Yi1(1)] ; Zi1 = [Zi1, Zi1(1)] ;
        % color
        C = zeros([size(Xo1,1), 2*size(Xo1,2),  3]) ; % color should be the same for all
        C(:,:,1) = rgbCoilMetal(1) ;
        C(:,:,2) = rgbCoilMetal(2) ;
        C(:,:,3) = rgbCoilMetal(3) ;
        
        % make patch filling space between rings
        h1 = patch([Xi1, Xo1], [Yi1, Yo1], [Zi1, Zo1], C, ...
            'lineStyle','none') ;
        set(h1,'FaceAlpha',objAlpha,'EdgeAlpha',objAlpha)
        h_annulus_array(2*j-1) = h1 ;
        % draw circle borders around annulus (if we use the patch edge we get
        % an annoying line where the circles close)
        h_line_array(4*j-3) = plot3(Xo1, Yo1, Zo1, 'k-',...
            'LineWidth',coilBorderLW) ;
        h_line_array(4*j-2) = plot3(Xi1, Yi1, Zi1, 'k-',...
            'LineWidth',coilBorderLW) ;
        
        % --------------------------------------------------
        % next draw bottom annulus on bottom face of chamber
        Xo2 = coilOuterRadius * cos(t) ;
        Yo2 = coilOuterRadius * sin(t) ;
        Zo2 = (z_list(j) + dz_list(j))*ones(size(Xo2))  ;
        %Xo2 = [Xo2, Xo2(1)] ; Yo2 = [Yo2, Yo2(1)] ; Zo2 = [Zo2, Zo2(1)] ;
        % inner ring
        Xi2 = coilInnerRadius * cos(t) ;
        Yi2 = coilInnerRadius * sin(t) ;
        Zi2 = (z_list(j) + dz_list(j))*ones(size(Xo2))  ;
        %Xi2 = [Xi2, Xi2(1)] ; Yi2 = [Yi2, Yi2(1)] ; Zi2 = [Zi2, Zi2(1)] ;
        % make patch and circles
        h2 = patch([Xi2, Xo2], [Yi2, Yo2], [Zi2, Zo2], C, ...
            'lineStyle','none') ;
        set(h2,'FaceAlpha',objAlpha,'EdgeAlpha',objAlpha)
        h_annulus_array(2*j) = h2 ;
        h_line_array(4*j-1) = plot3(Xo2, Yo2, Zo2, 'k-',...
            'LineWidth',coilBorderLW) ;
        h_line_array(4*j) = plot3(Xi2, Yi2, Zi2, 'k-',...
            'LineWidth',coilBorderLW) ;
        
        % finally, draw cylindrical patch that connects the two
        [Xc, Yc, Zc ] = cylinder(cylinderRadius, resolution) ;
        Zc = Zc * abs(mean(Zi1) - mean(Zi2)) + mean(Zi2) ;
        C = zeros([size(Xc) 3]) ;
        C(:,:,1) = rgbCoilWire(1) ;
        C(:,:,2) = rgbCoilWire(2) ;
        C(:,:,3) = rgbCoilWire(3) ;
        hCoilWire = surf(Xc, Yc, Zc) ;
        set(hCoilWire,'lineStyle','none') ;
        set(hCoilWire,'CData',C) ;
        set(hCoilWire,'FaceAlpha',objAlpha,'EdgeAlpha',objAlpha)
        h_cyl_array(j) = hCoilWire ;
    end
    
    % group coil objects together
    coil_grp = hgtransform('Parent',grp);
    set(h_annulus_array(:), 'Parent',coil_grp) ;
    set(h_line_array(:), 'Parent',coil_grp) ;
    set(h_cyl_array(:), 'Parent',coil_grp) ;
    
    %set material/lighting properties
    material(h_line_array(:) ,'dull')
    material(h_annulus_array(:) ,'metal')
    material(h_cyl_array(:) ,'metal')
    set(h_annulus_array(:), 'ambientStrength',.9);  
    set(h_cyl_array(:), 'ambientStrength',.9);  
end

% -----------------------------------------------------------------------
%% wrap up
% perform global translation, if necessary
T = makehgtform('translate',center) ; 
set(grp, 'Matrix', T) ; 

% final lighting properties
lighting flat

% reset figure's hold properties
if (~wasOnHold)
    hold off ;
end
end