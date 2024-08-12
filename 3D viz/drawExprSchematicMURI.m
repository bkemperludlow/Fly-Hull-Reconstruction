% -------------------------------------------------------------------------
% script to generate a matlab drawing of the camera/chamber/coil setup
% -------------------------------------------------------------------------
% ------------------------------
%% params
saveFlag = true ; % save output?
insetFlag = false ; % make zoomed-in image of fly?
myaaFlag = false ; % use the "my anti aliasing" tool from FEX? (doesn't work well when remote)

optoFlag = true ;  % add drawing of optogenetic light?
windTubeFlag = false ;  % add drawing of tube that sort of looks like it shoots air?
visualPanelFlag = false ; % add drawing of visual panel?
savePrefix = 'just_opto_silence' ; % any additional string to affix to filenames
pinType = 0 ;    % type of pin on fly. 0 = no pin , 1 = roll pin, 2 = pitch

savePath = 'D:\Dropbox\Dissertation\talk\B EXAM TALK\setup_schematics\' ; % ALTER AS NEEDED

% ----------------------
% overall scales
flyScale = 4.0 ;
scale = 5*flyScale ;
lineScale = flyScale/2.0 ;

% ------------------------------------
% set camera positons and directions
dist_scale = scale*25 ;
center_yz = dist_scale*[1, 0, 0] ;
center_xz = dist_scale*[0, 1, 0] ;
center_xy = dist_scale*[0, 0, 1] ;

% unit vector camera directions
direction_yz = [-1, 0, 0] ;
direction_xz = [0, -1, 0] ;
direction_xy = [0, 0, -1] ;

% camera proportions
camLen = scale*12 ;
camWid = scale*6 ;
camDep = scale*6 ;

% -----------------------------------------
% flight chamber setup
chamberSideLength = scale*15 ;
frameCrossSection = scale*1 ;

% -----------------------------------------
% visual panel setup
panelRadius = 0.45*chamberSideLength ;
panelHeight = chamberSideLength/4 ;
% angleRange = [pi/2 + pi/24, 3*pi/2 - pi/6] + pi/6 ;
angleRange = [0, 2*pi] ;
panelCenter = [0,0, 0.3*chamberSideLength] ; 
panelAlpha = 0.5 ; 

% -------------------------------------------
% wind manifold (to provide gust) setup
% NB: goint to just try to use LED code
tubeLength = scale*10 ;
tubeRadius = scale*0.8 ;
taperLength = scale*6 ;
tubeColor = 0.6*[1,1,1] ;

tube_center = scale*20.*[0.25, -1, 0] ;
tube_direction = [0,0,0] - tube_center ;

% -----------------------------------------
% fly setup -- for multiple flies, create matrix (1 row per fly)
%pinType = 2 ; % 0==no pin, 1==roll pin, 2==pitch
thetab0 = 45 * pi / 180 ;
bodyPos = [[0, 0, 0] ; [-50, -75, -55] ; [80, +100, -30] ; ...
    [-100, +100, +55] ; [+90, -100, +45] ; [+40, -50, +70]] ;
bodyYPR = [[pi/2, pi/4, 0] ; [pi/6, pi/6, -pi/8] ; [pi, pi/3, pi/8] ; ...
    [-pi/4, pi/3, 0] ; [-pi/4, pi/4, 0] ; [-pi/2, pi/6, pi/3]] ;
rightYPR = (pi/180)*[[ 140, 22, 55] ; [120, -15, 75] ; [20, 10, 55] ; ...
    [70, 0, 55] ; [110, 0, 135] ; [95, 0, 30] ] ;
leftYPR = (pi/180)*[[ 140, 22, 55] ; [120, -15, 75] ;  [20, 10, 55] ; ...
    [70, 0, 55] ; [110, 0, 135] ; [95, 0, 30]] ;
N_flies = size(bodyPos,1) ;
flyGridFlag = false ;

% ----------------------------------
% opto led setup
led_center = 100*[-2, 2, 5] ;
led_direction = [0,0,0] - led_center ;
led_direction = led_direction./norm(led_direction) ;
beamColorMap = 'Greens' ; % 'Reds'
ledLength = scale*8 ;
ledRadius = scale*2 ;
beamLength = scale*35 ;

% -----------------------
% misc params
resolution = 100 ;
az = 120 ; % 112 ;
el = 15 ;
lightPos = [0.2638    0.2293   -0.2018] ;

% =========================================================================
%% create figure
fig = figure('PaperPositionMode','auto','Position',[1, 41, 1920, 963],...
    'RendererMode','auto')  ;
hold on
ax = gca ;
parent = hgtransform('Parent',ax);

% draw cameras
cam_yx = draw3Dcamera(ax, parent, camLen, camWid, camDep, center_yz,...
    direction_yz, lineScale) ;
cam_xz = draw3Dcamera(ax, parent, camLen, camWid, camDep, center_xz, ...
    direction_xz, lineScale) ;
cam_xy = draw3Dcamera(ax, parent, camLen, camWid, camDep, center_xy, ...
    direction_xy, lineScale) ;

% draw chamber
flight_chamber = draw3DflightChamber(ax, parent, chamberSideLength, ...
    frameCrossSection, lineScale) ;

% draw visual panel?
if visualPanelFlag
    visual_panel = draw3DvisualPanel(ax, parent, panelRadius, ...
        panelHeight, angleRange, lineScale, panelCenter, [], panelAlpha) ;
end

% TEMP for debugging
axis equal 
rotate3d on 

% put some flies in the chamber
flyGrpList = gobjects(N_flies) ;
for ind = 1:N_flies
    [flyGrp, ~, rightWingGrp, leftWingGrp, dL, ~, ~] = draw3Dfly(ax, ...
        flyScale, resolution, pinType, thetab0, flyGridFlag);
    setFlyDOF(flyGrp, rightWingGrp, leftWingGrp, bodyPos(ind,:), ...
        bodyYPR(ind,:), rightYPR(ind,:) , leftYPR(ind,:), thetab0, dL )
    
    material(flyGrp, 'shiny')
    c = findobj(flyGrp,'Type','surface');
    set(c,'ambientstrength',0.5); % note that this line comes after "material shiny"
    
    flyGrpList(ind) = flyGrp ;
end
% try out the opto LED thing
if optoFlag
    opto_led = draw3DoptoLED(ax, parent, ledLength, ledRadius, ...
        beamColorMap, beamLength, led_center, led_direction) ;
end

if windTubeFlag
    wind_tube = draw3DwindTube(ax, parent, tubeLength, tubeRadius, ...
        taperLength, tube_center, tube_direction, tubeColor) ;
end

% adjust axis properties
axis equal
box on
%grid on

xlabel('X')
ylabel('Y')
zlabel('Z')

view(az, el)

% lighting?
%lighting flat
hlight = light ;
%set(hlight,'position',lightPos) ;
%hlight.Style = 'local' ;

% to save
set(fig,'color','w');
axis off

if saveFlag
    % get optimal (?) resolution
    screen_DPI = get(0,'ScreenPixelsPerInch');
    K = 4 ;
    res = screen_DPI*K ; 
    
    % turn off 'InvertHardCopy' option (tbh, can't remember why...)
    set(fig,'InvertHardcopy','off');
    
   
    % save to png using export graphics
    exportgraphics(fig, fullfile(savePath, ...
        sprintf('%s_setup_schematic_im1.png', savePrefix)), ...
        'Resolution', res)
     
    % try anti-aliasing wity myaa? doesn't work well on remote desktop
    if myaaFlag
        myaa ;
        fprintf('doing weird workaround with myaa')
        F = getframe ;
        imwrite(F.cdata,fullfile(savePath, ...
            sprintf('%s_setup_schematic_im3.png', savePrefix)))
    end
    
    
    print(fig,['-r',num2str(res)], '-dpng', ...
        fullfile(savePath, sprintf('%s_setup_schematic_im2.png', ...
        savePrefix)));
    export_fig(fig, fullfile(savePath, ...
        sprintf('%s_setup_schematic_im4.png', savePrefix)), ...
        '-dpng','-r600')
    %print(gcf, fullfile(savePath, 'setup_schematic.png'),'-dpng','-r600')
    %imwrite(gcf,fullfile(savePath, 'setup_schematic_im.png'))
end

% =========================================================================
%% make inset fig of just fly to show pert
if insetFlag
    h_justFly = figure('PaperPositionMode','auto','Position',[1921, 41, 1920, 963],...
        'RendererMode','auto')  ;
    %     h_justFly = figure('PaperPositionMode','auto','Position', [2143 , 492, 327, 212],...
    %         'RendererMode','auto')  ;
    hold on
    ax = gca ;
    parent = hgtransform('Parent',ax);
    
    flyScaleIncrease = 5 ;
    
    %  draw fly
    [flyGrp, ~, rightWingGrp, leftWingGrp, dL, ~, ~] = draw3Dfly(ax, ...
        flyScaleIncrease*flyScale, resolution, pinType, thetab0, true);
    setFlyDOF(flyGrp, rightWingGrp, leftWingGrp, bodyPos(1,:), bodyYPR(1,:), ...
        rightYPR(1,:) , leftYPR(1,:), thetab0, dL )
    
    material(flyGrp, 'shiny')
    c = findobj(flyGrp,'Type','surface');
    set(c,'ambientstrength',0.9); % note that this line comes after "material shiny"
    
    axis equal
    box on
    %grid on
    
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    
    view(az, el)
    
    % lighting?
    %lighting flat
    hlight = light ;
    lighting gouraud ;
    %set(hlight,'position',lightPos) ;
    %hlight.Style = 'local' ;
    % to save
    set(h_justFly,'color','w');
    axis off
    
    if saveFlag
        screen_DPI = get(0,'ScreenPixelsPerInch');
        K = 4 ;
        set(h_justFly,'InvertHardcopy','off');
        print(h_justFly,['-r',num2str(screen_DPI*K)], '-dpng', ...
            fullfile(savePath, 'just_fly_im2.png'));
        
        if myaaFlag
            myaa ;
            F = getframe ;
            imwrite(F.cdata,fullfile(savePath, 'just_fly_im.png'))
        end
    end
end