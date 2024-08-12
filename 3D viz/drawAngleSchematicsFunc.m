% -------------------------------------------------------------------------
% function to generate 3D fly models that can be used to illustrate body
% and wing kinematic angles
%
% INPUTS:
%   - schematicType: (string) which type of schematic to draw. options are
%       'body', 'wing', and 'wingFrame'
%   - pinType: which type of pin to draw on fly. 0 = no pin; 1 = roll pin;
%       2 = pitch/yaw pin
%   - labelFlag: (boolean) add text labels or not
%   - figSize: 2 element vector giving figure width and height (defaults to
%       full screen)
%   - saveFlag: (boolean) save output or not
%   - savePath: where to save images
%
% -------------------------------------------------------------------------
function [h_main, ax] = drawAngleSchematicsFunc(schematicType, pinType, ...
    labelFlag, figSize, saveFlag, savePath, axisFlag)
% ------------------------
%% inputs/params
if ~exist('schematicType','var') || isempty(schematicType)
    schematicType = 'body' ; % 'body' | 'wing' | 'wingFrame'
end
if ~exist('pinType','var') || isempty(pinType)
    pinType = 0 ;
end
if ~exist('labelFlag','var') || isempty(labelFlag)
    labelFlag = true ;  % add text labels to angles?
end
if ~exist('figSize','var') || isempty(figSize)
    figSize = [20.1667, 11.000] ; % full screen size in inches
end
if ~exist('saveFlag','var') || isempty(saveFlag)
    saveFlag = false ;
end
if ~exist('savePath','var') || isempty(savePath)
    savePath = pwd ; % TO BE ALTERED
end
if ~exist('axisFlag','var') || isempty(axisFlag)
    axisFlag = false ; % TO BE ALTERED
end
% ----------------------
%% overall scales
flyScale = 4.0 ;
scale = 5*flyScale ;
lineScale = flyScale/2.0 ;

arrowLength = 5*scale ;
arrowTipDiameter = 5 ;
arrowTipHeight = 10;
arrowRodDiameter = 2 ;

AHatArrowScale = 6/5 ;
labelTipScale = 1.05 ; % scale arrow tip location to make arrow fit

% colors
arrowColor = [0, 0, 0] ;
curvedColor = [203,24,29]/255 ;
ahatColor = [33,113,181]/255 ;

arrowResolution = 100 ;
arrowLineStyle = 'none' ;% '-' ;  %

arcRadius = 0.65*scale ; %0.75* scale ;

% stroke plane parameters
strokePlaneRadius = 6*scale ;
strokePlaneColor  = [158,202,225]/255 ;  %0.8*[1, 1, 1] ;
strokePlaneAlpha  = 0.3 ; %0.6 ;

% ------------------------
%% misc
flyMaterialType = 'metal' ;
arrowMaterialType = 'dull' ;
flyGridFlag = false ;
flyColorScheme = 'metal' ;

resolution = 100 ;
az = 51 ; % 124 ;%149 ; %104 ; %140 ; %124 ;
el = 16 ; %17 ; %9 ; %20 ; %26 ; %17 ;
lightPos = [0.2638    0.2293   -0.2018] ;

% xlim = [-284.3572  284.4176] ;
% ylim = [-150   150] ;
% zlim = [-60.4879   95.8579] ;

figUnits = 'inches' ;
figPosition = [-0.0833, 0.3333, figSize(1), figSize(1)] ; %[1921, 41, 1920, 963]

% --------------------
%% label info
% font preferences
fontSizeSmall = 6 ; % 12
fontSizeLarge = 6; % 16
fontWeight = 'bold' ;
fontName = 'arial' ;

% lab frame axis labels
xyzLabels = {'x_{lab}', 'y_{lab}', 'z_{lab}'} ;
% body axis labels
xyzBodyLabels = {'x_{body}', 'y_{body}', 'z_{body}'} ;
% body rotation labels
bodyRotLabels = {'roll \rho_{b}', 'pitch \theta_{b}', 'yaw \phi_{b}'} ;
% wing rotation labels
wingRotLabels = {'rotation \eta', 'deviation \theta', 'stroke \phi'} ;

% -----------------------------------------
%% fly DOF
%pinType = 2 ; % 0==no pin, 1==roll pin, 2==pitch
thetab0 = 45 * pi / 180 ;
bodyPos = [0, 0, 0] ;
bodyYPR = [0, pi/4, 0] ;
switch schematicType
    case 'body'
        rightYPR = (pi/180)*[ 140, 22, 55] ;
        leftYPR = (pi/180)*[ 140, 22, 55] ;
    case 'wing'
        rightYPR = (pi/180)*[ 90, 22, 55] ;
        leftYPR = (pi/180)*[ 90, 22, 55] ;
    case 'wingFrame'
        rightYPR = (pi/180)*[ 90, 0, 90] ;
        leftYPR = (pi/180)*[ 90, 0, 90] ;
    otherwise
        fprintf('Invalid schematic type: %s \n', schematicType)
        keyboard
end

% -------------------------------------------------------------------------
%% make figure
h_main = figure('PaperPositionMode','auto','MenuBar','none',...
    'ToolBar','none','DockControls','off','Units',figUnits,...
    'OuterPosition',figPosition,'RendererMode','auto' ) ;
% h_main = figure('PaperPositionMode','auto','Position',figPosition,...
%     'RendererMode','auto')  ;
%     h_justFly = figure('PaperPositionMode','auto','Position', [2143 , 492, 327, 212],...
%         'RendererMode','auto')  ;
hold on
ax = gca ;
parent = hgtransform('Parent',ax);

flyScaleIncrease = 5 ;

% --------------------------------------------------------------------------
%%  draw fly
[flyGrp, bodyGrp, rightWingGrp, leftWingGrp, dL, ~, ~] = draw3Dfly(ax, ...
    flyScaleIncrease*flyScale, resolution, pinType, thetab0, flyGridFlag,...
    flyColorScheme);
setFlyDOF(flyGrp, rightWingGrp, leftWingGrp, bodyPos, bodyYPR, rightYPR , ...
    leftYPR, thetab0, dL )
thoraxRadius = dL/(flyScaleIncrease*flyScale) ;

material(flyGrp, flyMaterialType)
c = findobj(flyGrp,'Type','surface');
set(c,'ambientstrength',0.9); % note that this line comes after "material shiny"

% -------------------------------------------------------------------------
%% get base point of fly
MB = get(flyGrp,'Matrix') ;
bp = [-10 ; 0 ; 0] / 10 * scale ;
v  = [bp ; 1 ] ; % base point on the RAW fly frame (x is body axis)
v2 = MB * v ;
basePoint = v2(1:3) ;

% =========================================================================
%% draw different stuff depending on schematic type
% =========================================================================
switch schematicType
    case 'body'
        % ------------------------------------------------------------------------
        %% lab frame arrows
        if axisFlag
            % define lab axes
            xdirYPRdeg = [ 0  0 0 ] ;
            ydirYPRdeg = [ 90 0 0 ] ;
            zdirYPRdeg = [ 0 90 0 ] ;
            % combine lab frame vectors into matrix
            xyzDirYPRdeg = [xdirYPRdeg; ydirYPRdeg; zdirYPRdeg] ;
            
            % initialize array for tip (to be used for labels)
            arrowTip = [arrowLength + arrowTipHeight, 0, 0] ;
            % scale up a little to avoid overlap
            arrowTip = labelTipScale*arrowTip ;
            xyzArrowTips = repmat(arrowTip,3,1) ;
            
            % initialize parent group for lab frame arrows
            xyzGroup = hgtransform('Parent',ax);
            labArrows = gobjects(size(xyzDirYPRdeg,1),1) ;
            
            % generate arrows and get tip locations
            for dim = 1:size(xyzDirYPRdeg,1)
                % draw arrow
                labArrows(dim) = myArrow(ax, xyzGroup, arrowLength, ...
                    arrowTipDiameter, arrowTipHeight, arrowRodDiameter , ...
                    basePoint, xyzDirYPRdeg(dim,:), arrowColor, ...
                    arrowResolution, arrowLineStyle) ;
                
                % get tip position
                ypr_rad = (pi/180).*xyzDirYPRdeg(dim,:) ;
                rotM=eulerRotationMatrix(ypr_rad(1),ypr_rad(2),ypr_rad(3));
                xyzArrowTips(dim,:) = rotM'*xyzArrowTips(dim,:)' ...
                    + basePoint ;
            end
            %{
        xArrowGrp = myArrow(ax,xyzGroup, arrowLength, arrowTipDiameter, ...
            arrowTipHeight, arrowRodDiameter , basePoint, xdirYPRdeg, arrowColor, ...
            arrowResolution,arrowLineStyle) ;
        yArrowGrp = myArrow(ax,xyzGroup, arrowLength, arrowTipDiameter, ...
            arrowTipHeight, arrowRodDiameter , basePoint, ydirYPRdeg, arrowColor, ...
            arrowResolution,arrowLineStyle) ;
        zArrowGrp = myArrow(ax,xyzGroup, arrowLength, arrowTipDiameter, ...
            arrowTipHeight, arrowRodDiameter , basePoint, zdirYPRdeg, arrowColor, ...
            arrowResolution,arrowLineStyle) ;
            %}
            material(xyzGroup, arrowMaterialType)
            
            % write labels on arrows?
            if labelFlag
                xyzLabelObjs = gobjects(size(xyzDirYPRdeg,1),1) ;
                
                for dim = 1:3
                    xyzLabelObjs(dim) = text(xyzArrowTips(dim,1), ...
                        xyzArrowTips(dim,2),  xyzArrowTips(dim,3), ...
                        xyzLabels{dim},'FontName',fontName,...
                        'FontSize',fontSizeSmall,'FontWeight',fontWeight) ;
                end
            end
        end
        % -------------------------------------------------------------------------
        %% draw  body axis arrow
        aHatArrowGrp = hgtransform('Parent',ax);
        myArrow(ax, aHatArrowGrp, AHatArrowScale*arrowLength, ...
            arrowTipDiameter, arrowTipHeight, arrowRodDiameter , ...
            basePoint, (180/pi).*bodyYPR, ahatColor, arrowResolution, ...
            arrowLineStyle) ;
        AHat = [cos(bodyYPR(1))*cos(bodyYPR(2)), sin(bodyYPR(1))*cos(bodyYPR(2)), ...
            sin(bodyYPR(2))]' ;
        
        material(aHatArrowGrp, arrowMaterialType)
        
        % body axis tip + label
        if labelFlag
            rotM = eulerRotationMatrix(bodyYPR(1), bodyYPR(2), bodyYPR(3)) ;
            AHatTip = [AHatArrowScale*arrowLength + arrowTipHeight, 0, 0] ;
            % scale up a bit to avoid overlap
            AHatTip = labelTipScale.*AHatTip ;
            % rotate and add base point
            AHatTip = (rotM'*AHatTip' + basePoint)' ;
            
            % make label
            AHatLabel = text(AHatTip(1), AHatTip(2), AHatTip(3), ...
                xyzBodyLabels{1}, 'FontName',fontName,'FontSize',fontSizeSmall,...
                'FontWeight',fontWeight, 'Color', ahatColor) ;
        end
        % -------------------------------------------------------------------------
        %% curved arrows
        
        curvedArrowsGroup = hgtransform('Parent',ax);
        %curvedColor = [203,24,29]/255 ;
        
        % ---------------------------
        % BODY ROLL ARROW
        % ---------------------------
        startAngle = 20 ;%-35 ; % 30
        endAngle   = 340 ; %-10 + 295 ; % 30 + 285 % 90+270 ;
        %YPRdeg    = 0*[ bodyYPRdeg(1) -bodyYPRdeg(2) 0 ] ;
        YPRdeg    = (180/pi).*[ bodyYPR(1), -bodyYPR(2), 0 ] ;  %[0 -90 0] ;
        rhoCenterPoint = basePoint + 5*scale*AHat ;
        rhoArrowGroup = myCurvedArrow(ax, curvedArrowsGroup, arcRadius, startAngle, ...
            endAngle,  arrowTipDiameter, 0.5*arrowTipHeight, arrowRodDiameter, ...
            rhoCenterPoint, YPRdeg, curvedColor, arrowResolution,arrowLineStyle) ;
        
        % label for roll arrow
        if labelFlag
            % get center point for this curved arrow
            rhoLabelPoint =  labelTipScale.*(rhoCenterPoint - basePoint) + ...
                basePoint ;
            % displace to the side by arcRadius + a little more
            dispVec = 1.15*labelTipScale*arcRadius*[0; 1.25; -1] ;
            rhoLabelPoint = rhoLabelPoint + dispVec ;
            
            % create label
            rhoLabel = text(rhoLabelPoint(1), rhoLabelPoint(2), ...
                rhoLabelPoint(3), bodyRotLabels{1}, 'FontName',fontName,...
                'FontSize',fontSizeLarge, 'FontWeight','normal', ...
                'Color', curvedColor) ;
        end
        % ---------------------------
        % BODY YAW ARROW
        % ---------------------------
        %body phi
        startAngle = -25 ; %150 ;
        endAngle   = 310 ; %150 + 270 ;% 90+270 ;
        %YPRdeg    = 0*[ bodyYPRdeg(1) -bodyYPRdeg(2) 0 ] ;
        YPRdeg    = [0 0 0] ;
        phiBase = basePoint + [0, 0, 4*scale]' ;
        phiArrowGroup = myCurvedArrow(ax, curvedArrowsGroup, arcRadius, ...
            startAngle, endAngle, arrowTipDiameter, 0.5*arrowTipHeight, ...
            arrowRodDiameter, phiBase, YPRdeg, curvedColor, arrowResolution,...
            arrowLineStyle) ;
        
        % label for yaw arrow
        if labelFlag
            % get center point for this curved arrow
            phiLabelPoint =  labelTipScale.*(phiBase - basePoint) + ...
                basePoint ;
            % displace to the side by arcRadius + a little more
            dispVec =  0.95*labelTipScale*arcRadius*[0; 1; 0.25] ;
            phiLabelPoint = phiLabelPoint + dispVec ;
            
            % create label
            phiLabel = text(phiLabelPoint(1), phiLabelPoint(2), ...
                phiLabelPoint(3), bodyRotLabels{3}, 'FontName',fontName,...
                'FontSize',fontSizeLarge, 'FontWeight','normal', ...
                'Color', curvedColor) ;
        end
        
        % ---------------------------
        % BODY PITCH ARROW
        % ---------------------------
        %body theta
        startAngle = -40 ; %270 ;
        endAngle   = 270 ;  %startAngle + 270 ;% 90+270 ;
        %YPRdeg    = 0*[ bodyYPRdeg(1) -bodyYPRdeg(2) 0 ] ;
        YPRdeg    = [-90 -90 0] ;
        thetaBase = basePoint + [0, 4*scale, 0]' ;
        thetaArrowGroup    = myCurvedArrow(ax, curvedArrowsGroup, arcRadius, ...
            startAngle, endAngle, arrowTipDiameter, 0.5*arrowTipHeight, ...
            arrowRodDiameter, thetaBase, YPRdeg, curvedColor, ...
            arrowResolution,arrowLineStyle) ;
        
        % label for pitch arrow
        if labelFlag
            % get center point for this curved arrow
            thetaLabelPoint =  labelTipScale.*(thetaBase - basePoint) + ...
                basePoint ;
            % displace to the side by arcRadius + a little more
            dispVec =  1.45*labelTipScale*arcRadius*[0; 0.25; -1] ;
            thetaLabelPoint = thetaLabelPoint + dispVec ;
            
            % create label
            thetaLabel = text(thetaLabelPoint(1), thetaLabelPoint(2), ...
                thetaLabelPoint(3), bodyRotLabels{2}, 'FontName',fontName,...
                'FontSize',fontSizeLarge, 'FontWeight','normal', ...
                'Color', curvedColor) ;
        end
        
        % --------------------------------------------------
        % curved arrow properties
        material(curvedArrowsGroup, arrowMaterialType)
        c = findobj(curvedArrowsGroup,'Type','surface');
        set(c,'ambientstrength',0.85); % note that this line comes after "material shiny"
        
    case 'wing'
        % -----------------------------------------------------------------
        %% draw stroke plane
        strokePlaneGroup = hgtransform('Parent',ax);
        
        t = linspace(0,rightYPR(1), resolution) ;
        X = [strokePlaneRadius * cos(t),0] ;  %+ thoraxRadius*scale ;
        Y = [-1*strokePlaneRadius * sin(t),0] - thoraxRadius*scale ;
        Z = zeros(size(X))  ;
        C = zeros([size(X) 3]) ;
        C(:,:,1) = strokePlaneColor(1) ;
        C(:,:,2) = strokePlaneColor(2) ;
        C(:,:,3) = strokePlaneColor(3) ;
        strokePlane = patch(X,Y,Z,C,'lineStyle','none') ;
        set(strokePlane,'FaceAlpha',strokePlaneAlpha,...
            'EdgeAlpha',strokePlaneAlpha)
        set(strokePlane,'parent',strokePlaneGroup) ;
        
        % -----------------------------------------------------------------
        %% draw wing vein arrow
        veinDir = [90, 180, 0] ; %(180/pi)*leftYPR ;
        veinArrowLength = strokePlaneRadius  ;
        wingBase = [0, -1, 0]' ;
        myArrow(ax, rightWingGrp, veinArrowLength, arrowTipDiameter, arrowTipHeight, ...
            arrowRodDiameter , wingBase, veinDir, arrowColor, arrowResolution, ...
            arrowLineStyle) ;
        
        % -----------------------------------------------------------------
        %% body coordinates arrows
        if axisFlag
            xyzBodyArrowGrp = hgtransform('Parent',ax);
            M1 = eulerRotationMatrix(bodyYPR(1), bodyYPR(2), bodyYPR(3)) ;
            %M2 = eulerRotationMatrix(0, -pi/4, 0) ;
            hatVecs = eye(3) ;
            bodyArrows = gobjects(size(hatVecs,1),1) ;
            
            for dim = 1:size(hatVecs)
                hatVecCurr = hatVecs(dim,:) ;
                hatVecBody = M1'*hatVecCurr' ;
                
                yawCurr = atan2(hatVecBody(2), hatVecBody(1)) ;
                pitchCurr = asin(hatVecBody(3)) ;
                YPRdeg = (180/pi)*[yawCurr, pitchCurr, 0] ;
                
                % scale length of body axis arrow
                if (dim == 3)
                    bodyArrowScale = 0.85 ;
                else
                    bodyArrowScale = 1.15 ;
                end
                
                % draw arrow
                bodyArrows(dim) = myArrow(ax, xyzBodyArrowGrp, ...
                    bodyArrowScale*arrowLength, arrowTipDiameter, ...
                    arrowTipHeight, arrowRodDiameter , basePoint, YPRdeg, ...
                    ahatColor, arrowResolution, arrowLineStyle) ;
                
                % make label?
                if labelFlag
                    scaleCurr = (bodyArrowScale + 0.05)*arrowLength + ...
                        arrowTipHeight ;
                    bodyArrowTip = scaleCurr*hatVecBody + basePoint ;
                    text(bodyArrowTip(1), bodyArrowTip(2), bodyArrowTip(3),...
                        xyzBodyLabels{dim}, 'FontName',fontName,...
                        'FontSize',fontSizeSmall, 'FontWeight',fontWeight, ...
                        'Color', ahatColor) ;
                end
            end
            %AHat = [cos(bodyYPR(1))*cos(bodyYPR(2)), sin(bodyYPR(1))*cos(bodyYPR(2)), ...
            %    sin(bodyYPR(2))]' ;
            
            material(xyzBodyArrowGrp, arrowMaterialType)
        end
        % -----------------------------------------------------------------
        %% curved arrows for wing angles
        
        curvedArrowsGroup = hgtransform('Parent',strokePlaneGroup);
        %curvedColor = [203,24,29]/255 ;
        
        % ---------------------------
        % WING PITCH ARROW
        % ---------------------------
        spanVec = [cos(rightYPR(1))*cos(rightYPR(2)), ...
            -sin(rightYPR(1))*cos(rightYPR(2)), ...
            sin(rightYPR(2))]' ;
        startAngle = 180 ;
        endAngle   = 420 ; %25 + 270 ;% 90+270 ;
        %YPRdeg    = 0*[ bodyYPRdeg(1) -bodyYPRdeg(2) 0 ] ;
        YPRdeg    = (180/pi).*[rightYPR(1), 2.5*rightYPR(2), 0 ] ;  %[0 -90 0] ;
        etaCenterPoint = 5.0*scale*spanVec + [ 0  -thoraxRadius*scale 0 ]' ;
        etaArrowGroup = myCurvedArrow(ax, curvedArrowsGroup, arcRadius,...
            startAngle, endAngle,  arrowTipDiameter, 0.5*arrowTipHeight, ...
            arrowRodDiameter, etaCenterPoint, YPRdeg, curvedColor, ...
            arrowResolution,arrowLineStyle) ;
        
        % wing pitch arrow label
        if labelFlag
            % get center point for this curved arrow
            etaLabelPoint =  labelTipScale.*(etaCenterPoint - basePoint) + ...
                basePoint ;
            % displace to the side by arcRadius + a little more
            dispVec = 1.15*labelTipScale*arcRadius*[0; 1.25; +1] ;
            etaLabelPoint = etaLabelPoint + dispVec ;
            
            % create label
            etaLabel = text(etaLabelPoint(1), etaLabelPoint(2), ...
                etaLabelPoint(3), wingRotLabels{1}, 'FontName',fontName,...
                'FontSize',fontSizeLarge, 'FontWeight','normal', ...
                'Color', curvedColor) ;
            
        end
        
        % ---------------------------
        % WING STROKE ANGLE
        % ---------------------------
        startAngle  = 0 ;
        endAngle    = (180/pi)*rightYPR(1) ;
        
        alph = endAngle - 90 ;
        dx   = scale*thoraxRadius ; %*tan(alph*pi/180) ;
        centerPoint = [0, -dx, 0] ; %[ dx  0 0 ] ;
        phiArcRadius   = strokePlaneRadius ;
        YPRdeg      = [0 0 180 ] ;
        phiArrow    = myCurvedArrow(ax, curvedArrowsGroup, phiArcRadius, ...
            startAngle, endAngle, arrowTipDiameter, 0.5*arrowTipHeight, ...
            arrowRodDiameter, centerPoint, YPRdeg, curvedColor, ...
            resolution,arrowLineStyle) ;
        
        % wing stroke arrow label
        if labelFlag
            % get center point for this curved arrow
            phiLabelPoint =  labelTipScale.*centerPoint ;
            % displace to the side by arcRadius + a little more
            dispVec = 0.85*labelTipScale*phiArcRadius*[0.5; -1.25; -0.075] ;
            phiLabelPoint = phiLabelPoint + dispVec ;
            
            % create label
            phiLabel = text(phiLabelPoint(1), phiLabelPoint(2), ...
                phiLabelPoint(3), wingRotLabels{3}, 'FontName',fontName,...
                'FontSize',fontSizeLarge, 'FontWeight','normal', ...
                'Color', curvedColor) ;
            
        end
        % ---------------------------
        % WING DEVIATION ARROW
        % ---------------------------
        startAngle  = 1 ;
        endAngle    = (180/pi)* rightYPR(2) - 1 ;
        YPRdeg      = [-(180/pi)*rightYPR(1),  0,  90 ] ;
        centerPoint =  [ 0  -thoraxRadius*scale 0 ] ;
        dr = sqrt(thoraxRadius^2 + dx^2) ;
        manualFactor = 0.2 ; % 0.75 ;
        thetaArcRadius = phiArcRadius-manualFactor*dr;
        thetaArrow    = myCurvedArrow(ax, curvedArrowsGroup,...
            thetaArcRadius, startAngle, endAngle, ...
            arrowTipDiameter, 0.5*arrowTipHeight, arrowRodDiameter, ...
            centerPoint, YPRdeg, curvedColor, resolution,arrowLineStyle) ;
        
        % wing deviation arrow label
        if labelFlag
            % get center point for this curved arrow
            thetaLabelPoint =  labelTipScale.*centerPoint ;
            % displace to the side by arcRadius + a little more
            dispVec = labelTipScale*thetaArcRadius*[0; -1.3; 0] ;
            thetaLabelPoint = thetaLabelPoint + dispVec ;
            
            % create label
            thetaLabel = text(thetaLabelPoint(1), thetaLabelPoint(2), ...
                thetaLabelPoint(3), wingRotLabels{2}, 'FontName',fontName,...
                'FontSize',fontSizeLarge, 'FontWeight','normal', ...
                'Color', curvedColor,'Rotation',85) ;
            
        end
        
        % ----------------------------------------------
        % set arrow material properties
        material(curvedArrowsGroup, arrowMaterialType)
        c = findobj(curvedArrowsGroup,'Type','surface');
        set(c,'ambientstrength',0.85); % note that this line comes after "material shiny"
        
    case 'wingFrame'
        % remove body and right wing
        delete(bodyGrp)
        delete(rightWingGrp)
        
        % draw left wing vein (x direction)
        veinDir = [-90, 180, 0] ; %(180/pi)*leftYPR ;
        veinArrowLength = 0.85*strokePlaneRadius  ;
        wingBase = [0, 1, 0]' ;
        myArrow(ax, leftWingGrp, veinArrowLength, arrowTipDiameter, ...
            arrowTipHeight, arrowRodDiameter , wingBase, veinDir, ...
            arrowColor, arrowResolution, arrowLineStyle) ;
        
        if labelFlag
            spanVec = [cos(leftYPR(1))*cos(leftYPR(2)), ...
                sin(leftYPR(1))*cos(leftYPR(2)), ...
                sin(leftYPR(2))]' ;
            tipPoint = ...
                0.975*labelTipScale*(veinArrowLength + ...
                arrowTipHeight)*spanVec + [ 0  thoraxRadius*scale 0 ]' ;
            
            text(tipPoint(1), tipPoint(2), tipPoint(3), 'x_{wing}', ...
                'FontName',fontName, 'FontSize',fontSizeSmall, ...
                'FontWeight',fontWeight, 'Color', 'k') ;
        end
        
        % draw left wing y direction
        yDir = [0, 0, 0] ; %(180/pi)*leftYPR ;
        yArrowLength = 0.25*strokePlaneRadius  ;
        wingBase = [0, 1, 0]' ;
        myArrow(ax, leftWingGrp, yArrowLength, arrowTipDiameter, ...
            arrowTipHeight, arrowRodDiameter , wingBase, yDir, ...
            arrowColor, arrowResolution, arrowLineStyle) ;
        
        if labelFlag
            tipPoint = ...
                [0, 0, labelTipScale*(yArrowLength + arrowTipHeight)]' + ...
                [ 0 thoraxRadius*scale  0 ]' ;
            text(tipPoint(1), tipPoint(2), tipPoint(3), 'y_{wing}', ...
                'FontName',fontName, 'FontSize',fontSizeSmall, ...
                'FontWeight',fontWeight, 'Color', 'k') ;
        end
        
        % draw left wing z direction
        zDir = [0, -90, 0] ; %(180/pi)*leftYPR ;
        zArrowLength = 0.25*strokePlaneRadius  ;
        wingBase = [0, 1, 0]' ;
        myArrow(ax, leftWingGrp, zArrowLength, arrowTipDiameter, ...
            arrowTipHeight, arrowRodDiameter , wingBase, zDir, ...
            arrowColor, arrowResolution, arrowLineStyle) ;
        
        if labelFlag
            tipPoint = ...
                [labelTipScale*(zArrowLength + arrowTipHeight), 0, 0]' + ...
                [ 0 thoraxRadius*scale  0 ]' ;
            dispVec = 5*[0, -1, -1]' ;
            tipPoint = tipPoint + dispVec ;
            text(tipPoint(1), tipPoint(2), tipPoint(3), 'z_{wing}', ...
                'FontName',fontName, 'FontSize',fontSizeSmall, ...
                'FontWeight',fontWeight, 'Color', 'k') ;
        end
        
        % update axis limits and view to focus on just left wing
        az = 90 ;
        el = 0 ;
        
        %xlim = [-5, 55] ;
        %ylim = [20, 165] ;
        %zlim = [-65, 80] ;
    otherwise
        fprintf('Invalid schematic type: %s \n', schematicType)
        keyboard
end
% -------------------------------------------------------------------------
%% axis properties and lighting
axis equal
box on
% set(ax,'xlim', xlim, 'ylim',ylim,'zlim',zlim)
%grid on

xlabel('X')
ylabel('Y')
zlabel('Z')

view(az, el)

% lighting?
%lighting flat
%hlight = light ;
%lighting gouraud ;
%set(hlight,'position',lightPos) ;
%hlight.Style = 'local' ;
% to save
set(h_main,'color','w');
axis off

if saveFlag
    % get screen resolution
    screen_DPI = get(0,'ScreenPixelsPerInch');
    K = 4 ;
    
    % basically don't mess with the background color (it's set above)
    set(h_main,'InvertHardcopy','off');
    
    % ----------------------------------------------------------------
    % take shotgun approach to saving, see what gives best appearance
    
    % try the (new?) matlab export graphics function
    exportgraphics(h_main, fullfile(savePath, ['angle_schematic_' ...
        schematicType '_export_graphics.png']),'Resolution',screen_DPI*K)
    
    % regular old print
    print(h_main,['-r',num2str(screen_DPI*K)], '-dpng', ...
        fullfile(savePath, ['angle_schematic_' schematicType '_print.png']));
    
    % export fig
    export_fig(h_main, fullfile(savePath, ['angle_schematic_' ...
        schematicType '_export_fig.png']), '-dpng','-r500', '-a4');
    
    
    myaa ;
    F = getframe ;
    try
        imwrite(F.cdata, fullfile(savePath, ...
            ['angle_schematic_' schematicType '_myaa.png']))
    catch
        fprintf('Failed to save myaa version ...\n')
    end
end
end
