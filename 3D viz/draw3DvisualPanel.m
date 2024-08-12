% -------------------------------------------------------------------------
%  function to draw 3D version of visual panel around flight chamber for
% creating figure schematics.
%
%{
fig = figure('Position',[1, 346, 560, 420])  ;
hold on
ax = gca ;
parent = hgtransform('Parent',ax);

panelRadius = 20 ; 
panelHeight = 30 ; 
angleRange = [pi/2, pi + pi/6] ; 
lineScale = 2 ; 

grp = draw3DvisualPanel(ax, parent, panelRadius,...
    panelHeight, angleRange, lineScale) ; 
%}
% -------------------------------------------------------------------------
function grp = draw3DvisualPanel(ax, parent, panelRadius, ...
    panelHeight, angleRange, lineScale, center, lineStyleStr, ...
    objAlpha, resolution)
% -----------------------
%% params and inputs
if ~exist('ax', 'var') || isempty(ax)
    ax = gca ;
end
if ~exist('center', 'var') || isempty(center)
    center = [0, 0, 0] ;
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

% -------------------------
% some other misc properties of drawing
lineWidth = lineScale * 0.5 ; % 0.5 is default line thickness, i think
rgbPanel = [0.6196, 0.7922, 0.8824] ; % light blue

rgbStripe = [0,0,0] ; 
stripeArcRadius = pi/12 ; 
stripeAlpha = objAlpha ; 

% -----------------------
% properties of input axis
axes(ax) ;
wasOnHold = ishold ;
hold on ;

% initialize total output group
grp = hgtransform('Parent',parent) ;

% -----------------------------------------------------------------------
%% draw in semi cylinder for visual panel
% visual panel (which in reality would be a full cyliner) should be
% displayed as a partial cylinder with stripes (or some kind of pattern) on
% the walls

% --------------------------------------------------
% get coordinates for partial cylinder
t = linspace(angleRange(1),angleRange(2), resolution) ; % angle parameter
z_list = (panelHeight/2).*[-1, 1] ; % bottom and top z position of partial cylinder

% get coordinates for the top and bottom arcs of the partial cylinder
Xbot = panelRadius * cos(t) ;
Ybot = panelRadius * sin(t) ;
Zbot = z_list(1)*ones(size(Xbot))  ;

Xtop = Xbot ;
Ytop = Ybot ;
Ztop = z_list(2)*ones(size(Xtop))  ;

% -------------------------------------------------
%% draw outline of panel
% initialize storage for line objects
h_line_array = gobjects() ;

% make curves for top and bottom arcs
h_line_array(1) = line(ax, Xtop, Ytop, Ztop, 'Color', 'k',...
    'LineWidth',lineWidth, 'LineStyle', lineStyleStr) ;

h_line_array(2) = line(ax, Xbot, Ybot, Zbot, 'Color', 'k',...
    'LineWidth',lineWidth, 'LineStyle', lineStyleStr) ;

% also draw curves connecting the top and bottom arcs
if diff(angleRange) < 2*pi
    h_line_array(3) = line(ax, [Xtop(1), Xtop(1)], [Ytop(1), Ytop(1)], ...
        [Zbot(1), Ztop(1)], 'Color', 'k', 'LineWidth',lineWidth, ...
        'LineStyle', lineStyleStr) ;
    
    h_line_array(4) = line(ax, [Xtop(end), Xtop(end)], [Ytop(end), Ytop(end)], ...
        [Zbot(1), Ztop(1)], 'Color', 'k', 'LineWidth',lineWidth, ...
        'LineStyle', lineStyleStr) ;
end

% group line objects together
line_grp = hgtransform('Parent',grp);
set(h_line_array(:), 'Parent',line_grp) ;

% -------------------------------------------------
%% draw patch to fill in between outline
% NB: if it's a semi-cylinder, use patch. otherwise, just draw cylinder
if diff(angleRange)< 2*pi
    % -------------------------
    % semi cylinder case
    % -------------------------
    % get XYZ coordinates for panel patch
    panelX = [Xbot, fliplr(Xtop)] ;
    panelY = [Ybot, fliplr(Ytop)] ;
    panelZ = [Zbot, fliplr(Ztop)] ;
    
    % get CData for panel
    C = zeros([size(panelX,1), size(panelX,2),  3]) ;
    C(:,:,1) = rgbPanel(1) ;
    C(:,:,2) = rgbPanel(2) ;
    C(:,:,3) = rgbPanel(3) ;
    
    % draw panel
    h_panel = patch(panelX, panelY, panelZ, C, 'lineStyle','none') ;
    set(h_panel,'FaceAlpha',objAlpha,'EdgeAlpha',objAlpha)
    
else
    % -------------------------
    % full cylinder case
    % -------------------------
    % get XYZ coordinates for cylinder
    [Xc, Yc, Zc ] = cylinder(panelRadius, resolution) ;
    Zc = Zc * abs(mean(Ztop) - mean(Zbot)) + mean(Zbot) ;
    
    % set CData for cylinder
    C = zeros([size(Xc) 3]) ;
    C(:,:,1) = rgbPanel(1) ;
    C(:,:,2) = rgbPanel(2) ;
    C(:,:,3) = rgbPanel(3) ;
    
    % draw cylinder as surface plot
    h_panel = surf(Xc, Yc, Zc) ;
    set(h_panel,'lineStyle','none') ;
    set(h_panel,'CData',C) ;
    set(h_panel,'FaceAlpha',objAlpha,'EdgeAlpha',objAlpha)
end



% add the panel to a group
panel_grp = hgtransform('Parent',grp);
set(h_panel, 'Parent', panel_grp) ; 

% ------------------------------------------------------
%% draw panel stripes
% first get indices for stripes from arc parameter (t)
t_stripe = t(1) : stripeArcRadius : t(end) ; 

% initialize storage for stripe patches
h_stripe_array = gobjects() ; 

% initialize counter 
cc = 1 ;

% loop over these to make stripe patches
for k = 1:2:(length(t_stripe)-1)
    % get start and end points (in circle param terms) for current stripe
    t1 = t_stripe(k) ;
    t2 = t_stripe(k+1) ;
    
    % find index that corresponds to this range
    idx = (t >= t1) & (t <= t2) ; 
    
    % get coordinates for current stripe. NB: we're going to make two: one
    % for inside and one for outside (so shrink/grow radius a bit on each)
    stripeRadii = [0.99*panelRadius, (1/0.99)*panelRadius] ;
    for m = 1:length(stripeRadii)
        
        X_stripe_top = stripeRadii(m) * cos(t(idx)) ;
        X_stripe_bot = X_stripe_top ;
        
        Y_stripe_top = stripeRadii(m) * sin(t(idx)) ;
        Y_stripe_bot = Y_stripe_top ;
        
        Z_stripe_top = Ztop(idx) ;
        Z_stripe_bot = Zbot(idx) ;
        
        % set CData for stripe
        C = zeros([size(X_stripe_bot,1), 2*size(X_stripe_bot,2),  3]) ; 
        C(:,:,1) = rgbStripe(1) ;
        C(:,:,2) = rgbStripe(2) ;
        C(:,:,3) = rgbStripe(3) ;
        
        % draw stripe patch
        h_stripe_array(cc) = patch([X_stripe_bot, fliplr(X_stripe_top)],...
            [Y_stripe_bot, fliplr(Y_stripe_top)],...
            [Z_stripe_bot, fliplr(Z_stripe_top)],...
            C, 'lineStyle','none') ;
        
        set(h_stripe_array(cc),'FaceAlpha',stripeAlpha,...
            'EdgeAlpha',stripeAlpha)
        
        % increment counter
        cc = cc + 1 ;
    end
end

% add stripes to group
set(h_stripe_array(:), 'Parent', panel_grp) ; 

% -------------------------------------------
%% set material/lighting properties
material(h_line_array(:) ,'dull')
material(h_panel ,'shiny')
material(h_stripe_array(:) ,'shiny')

set(h_panel, 'ambientStrength',0.9);
set(h_stripe_array(:), 'ambientStrength',0.9);

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

% axis equal
% rotate3d on
end