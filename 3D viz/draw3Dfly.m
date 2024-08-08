
function [flyGrp, bodyGrp, rightWingGrp, leftWingGrp, dL, rightWingTip, leftWingTip ] = ....
    draw3Dfly(ax, scale, resolution, pinType, thetab0, gridFlag, colorScheme) %#ok<INUSD>

if ~exist('gridFlag','var') || isempty(gridFlag)
    gridFlag = true ; 
end
if ~exist('colorScheme','var') || isempty(colorScheme)
    colorScheme = 'normal' ; 
end

bodyGridResolution = 20 ;
gridLineWidth = 0.5 ; % 0.5 ; 

headRadius    = 1 ;
thoraxRadius  = 1.2 ;
%abdomenRadius = 1.2 ;
LAMBDA = 1.4 ;
alphaVal = .7 ;


wingSpan      = thoraxRadius * 1.75  ; % 2.2
wingChord     = wingSpan /2   ;
wingThickness = (1/6) ; 

switch colorScheme
    case 'normal'
        pinColor = [.5 .5 1] ;
        eyeColor = [255 85 52] / 255 ; % [0.9 0.0 0.0] ;
        thoraxColor = [194 97 20] / 255 ;
        headColor   = [217 185 136]/255 ;
        abdomenColor = mean([thoraxColor ; headColor]) ;
        rightWingColor = [194 97 20]/255*.6 ; %[1 0 0 ]
        leftWingColor  = [194 97 20]/255*.6  ; % [0.5 0.5 1 ]
    case 'metal'
        pinColor = 0.5*[1 1 1] ;
        eyeColor = 0.5*[1 1 1] ; % [0.9 0.0 0.0] ;
        thoraxColor = 0.6*[1 1 1] ;
        headColor   = 0.7*[1 1 1]  ;
        abdomenColor = 0.7*[1 1 1] ;
        rightWingColor = 0.4*[1 1 1] ; %[1 0 0 ]
        leftWingColor  = 0.4*[1 1 1] ;
    case 'sim'
        pinColor = 0.7*[1 1 1] ;
        eyeColor = [0, 1, 0] ; % [0.9 0.0 0.0] ;
        thoraxColor = [0, 1, 0] ;
        headColor   = [0, 1, 0]  ;
        abdomenColor = [0, 1, 0] ;
        rightWingColor = [1, 0, 0] ; %[1 0 0 ]
        leftWingColor  = [0, 0, 1] ;
    otherwise
        fprintf('Invalid color scheme selection: %s \n', colorScheme)
end
pinLength = thoraxRadius * 2 ;
pinRadius = thoraxRadius * 0.1 ;

veinRadius = thoraxRadius * 0.15 ;

lineStyleStr = 'none' ;
pinLineStyleStr =  'none' ;%'-' ; %
wingLineStyle =  'none' ;%'-' ; %
eyeLineStyle  = 'none' ;

axes(ax) ;
hold on ;

flyGrp       = hgtransform('Parent',ax);
bodyGrp      = hgtransform('Parent',flyGrp);
rightWingGrp = hgtransform('Parent',flyGrp);
leftWingGrp  = hgtransform('Parent',flyGrp);

headGrp     = hgtransform('Parent',bodyGrp);
thoraxGrp   = hgtransform('Parent',bodyGrp);
abdomenGrp  = hgtransform('Parent',bodyGrp);
rightEyeGrp = hgtransform('Parent',headGrp);
leftEyeGrp  = hgtransform('Parent',headGrp);
[Xsph, Ysph, Zsph] = sphere(resolution) ;
[XsphLR, YsphLR, ZsphLR] = sphere(round(resolution/2)) ; % LR=Low Resolution

[XsphCoarse, YsphCoarse, ZsphCoarse] = sphere(bodyGridResolution) ;
[XsphLRCoarse, YsphLRCoarse, ZsphLRCoarse] = sphere(round(bodyGridResolution/2)) ; % LR=Low Resolution

myeps = 1.00 ;

% head
head = surf(XsphLR*headRadius*scale, YsphLR*headRadius*scale, ZsphLR*headRadius*1*scale) ;
C = zeros([size(XsphLR) 3]) ;
C(:,:,1) = headColor(1) ;
C(:,:,2) = headColor(2) ;
C(:,:,3) = headColor(3) ;
set(head,'CData',C) ;
set(head,'linestyle',lineStyleStr) ;

if gridFlag
    % make head grid
    gridDivision = 10 ;
    th = linspace(0,2*pi,200) ;
    xc = cos(th) * scale ;
    yc = xc*0 ;
    zc = sin(th) * scale ;
    hold on ;
    
    for k=0:9
        hcir = plot3(xc,yc,zc,'k-','linewidth',gridLineWidth) ;
        RC   = makehgtform('zrotate',k*2*pi/gridDivision) ;
        tmpGroup = hgtransform('Parent',headGrp);
        set(hcir,'parent',tmpGroup) ;
        set(tmpGroup,'Matrix',RC) ;
    end
    
    dphi = pi / ( gridDivision-1) ;
    phi = (1:(gridDivision-2)) * dphi - pi/2;
    for k=1:length(phi)
        zc = sin(phi(k))* scale + zeros(size(th));
        rad = cos(phi(k))* scale ;
        yc = rad * cos(th) ;
        xc = rad * sin(th) ;
        hcir = plot3(xc,yc,zc,'k-','linewidth',gridLineWidth) ;
        tmpGroup = hgtransform('Parent',headGrp);
        set(hcir,'parent',tmpGroup) ;
    end
    axis equal
end

%{
headGrid = mesh(XsphLR*headRadius*scale*myeps, YsphLR*headRadius*scale*myeps, ZsphLR*headRadius*1*scale*myeps) ;
set(headGrid,'facecolor','none') ;
C = zeros([size(XsphLR) 3]) ;
set(headGrid,'CData',C) ;
set(headGrid,'linestyle','-') ;
%}


% eyes
%----
rightEye = surf(XsphLR*headRadius*scale*0.6, ...
    YsphLR*headRadius*scale*0.6, ... % * 0.85
    ZsphLR*headRadius*2/3*scale*0.6 + headRadius*scale*0.8,'linestyle',lineStyleStr) ;
R1    = makehgtform('xrotate',+pi/2*1.1) ;
C = zeros([size(XsphLR) 3]) ;
C(:,:,1) = eyeColor(1) ;
C(:,:,2) = eyeColor(2) ;
C(:,:,3) = eyeColor(3) ;
set(rightEye,'CData',C) ;
set(rightEye,'Parent',rightEyeGrp) ;
set(rightEyeGrp,'Matrix',R1) ;



leftEye = surf(XsphLR*headRadius*scale*0.6, ...
    YsphLR*headRadius*scale*0.6, ... % *0.85
    ZsphLR*headRadius*2/3*scale*0.6 + headRadius*scale*0.8,'linestyle',lineStyleStr) ;
R2    = makehgtform('xrotate',-pi/2*1.1) ;
set(leftEye,'CData',C) ;

set(leftEye,'Parent',leftEyeGrp) ;
set(leftEyeGrp,'Matrix',R2) ;

S1    = makehgtform('scale',[1 1 2/3 ]) ;
R1    = makehgtform('yrotate',-pi*2/3) ;
T1    = makehgtform('translate',[1.5 0 2/3]*scale) ;
set(head, 'Parent',headGrp) ;
%set(headGrid, 'Parent',headGrp) ;
set(headGrp, 'Matrix',T1*R1*S1) ;

% thorax
thorax = surf(Xsph*thoraxRadius*scale, Ysph*thoraxRadius*scale, Zsph*thoraxRadius*scale,'linestyle',lineStyleStr) ;
C = zeros([size(Xsph) 3]) ;
C(:,:,1) = thoraxColor(1) ;
C(:,:,2) = thoraxColor(2) ;
C(:,:,3) = thoraxColor(3) ;
set(thorax,'CData',C) ;
set(thorax, 'Parent',thoraxGrp) ;

% make thorax grid
if gridFlag
    gridDivision = 20 ;
    th = linspace(0,2*pi,200) ;
    xc = cos(th) * scale * thoraxRadius;
    yc = xc*0 ;
    zc = sin(th) * scale * thoraxRadius;
    hold on ;
    for k=0:9
        hcir = plot3(xc,yc,zc,'k-','linewidth',gridLineWidth) ;
        RC   = makehgtform('zrotate',k*2*pi/gridDivision) ;
        tmpGroup = hgtransform('Parent',thoraxGrp);
        set(hcir,'parent',tmpGroup) ;
        set(tmpGroup,'Matrix',RC) ;
    end
    dphi = pi / ( gridDivision-1) ;
    phi = (1:(gridDivision-2)) * dphi - pi/2;
    for k=1:length(phi)
        zc = sin(phi(k))* scale* thoraxRadius + zeros(size(th));
        rad = cos(phi(k))* scale* thoraxRadius ;
        yc = rad * cos(th) ;
        xc = rad * sin(th) ;
        hcir = plot3(xc,yc,zc,'k-','linewidth',gridLineWidth) ;
        tmpGroup = hgtransform('Parent',thoraxGrp);
        set(hcir,'parent',tmpGroup) ;
    end
end

% abdomen
abdomen = surf(Xsph*headRadius*scale, Ysph*headRadius*scale, Zsph*headRadius*7/3*scale,'linestyle',lineStyleStr) ;
C = zeros([size(Xsph) 3]) ;
C(:,:,1) = abdomenColor(1) ;
C(:,:,2) = abdomenColor(2) ;
C(:,:,3) = abdomenColor(3) ;
set(abdomen,'CData',C) ;
R1 = makehgtform('yrotate',-75*pi/180) ;
T1 = makehgtform('translate',[-headRadius*4/3*scale 0 0]) ;
set(abdomen, 'Parent',abdomenGrp) ;
set(abdomenGrp,'Matrix',T1*R1) ;

% make abdomen grid
if gridFlag 
    gridDivision = 20 ;
    th = linspace(0,2*pi,200) ;
    xc = cos(th) * headRadius*scale ;
    yc = xc*0 ;
    zc = sin(th) * headRadius*7/3*scale;
    hold on ;
    for k=0:9
        hcir = plot3(xc,yc,zc,'k-','linewidth',gridLineWidth) ;
        RC   = makehgtform('zrotate',k*2*pi/gridDivision) ;
        tmpGroup = hgtransform('Parent',abdomenGrp);
        set(hcir,'parent',tmpGroup) ;
        set(tmpGroup,'Matrix',RC) ;
    end
    dphi = pi / ( gridDivision-1) ;
    phi = (1:(gridDivision-2)) * dphi - pi/2;
    for k=1:length(phi)
        zc = sin(phi(k))* headRadius*7/3*scale + zeros(size(th));
        rad = cos(phi(k))* headRadius*scale ;
        yc = rad * cos(th) ;
        xc = rad * sin(th) ;
        hcir = plot3(xc,yc,zc,'k-','linewidth',gridLineWidth) ;
        tmpGroup = hgtransform('Parent',abdomenGrp);
        set(hcir,'parent',tmpGroup) ;
    end
end



% now rotate the whole fly to look real
R1 = makehgtform('xrotate',pi) ;
R2 = eye(4) ; % makehgtform('yrotate',-pi/4) ;
set(bodyGrp,'Matrix',R2*R1) ;

% magnetic pin
switch (pinType)
    case 1 % roll pin
        [Z, X, Y] = cylinder ;
        Y = (Y - 0.5)*scale*pinLength ;
        Z = Z*scale*pinRadius - thoraxRadius*scale;
        X = X*pinRadius*scale ;
        pin = surf(X,Y,Z,'linestyle',pinLineStyleStr) ;
        
        C = zeros([size(X) 3]) ;
        C(:,:,1) = pinColor(1) ;
        C(:,:,2) = pinColor(2) ;
        C(:,:,3) = pinColor(3) ;
        set(pin,'CData',C) ;
        
        % disks that seal the pin cylinder
        % disk 1
        t = linspace(0,2*pi, resolution) ;
        X = scale*pinRadius * cos(t) ;
        Y = zeros(size(X))  + 0.5*scale*pinLength ;
        Z = scale*pinRadius * sin(t) - thoraxRadius*scale;
        C = zeros([size(X) 3]) ;
        C(:,:,1) = pinColor(1) ;
        C(:,:,2) = pinColor(2) ;
        C(:,:,3) = pinColor(3) ;
        hdisk1 = patch(X,Y,Z,C,'linestyle',pinLineStyleStr) ;
        
        % disk2
        Y = Y - scale*pinLength ;
        hdisk2 = patch(X,Y,Z,C,'linestyle',pinLineStyleStr) ;
        
        pinGrp = hgtransform('Parent',thoraxGrp);
        set([pin hdisk1 hdisk2],'Parent',pinGrp) ;
        
        R1 = makehgtform('yrotate',-pi/4) ;
        set(pinGrp,'Matrix',R1) ;
    case 2
        [Z, Y, X] = cylinder ;
        Y = Y*pinRadius*scale ;
        Z = Z*scale*pinRadius ;  
        X = -X*scale*pinLength ;
        pin = surf(X,Y,Z,'linestyle',pinLineStyleStr) ;
        
        C = zeros([size(X) 3]) ;
        C(:,:,1) = pinColor(1) ;
        C(:,:,2) = pinColor(2) ;
        C(:,:,3) = pinColor(3) ;
        set(pin,'CData',C) ;
        
        % disks that seal the pin cylinder
        % disk 1
        t = linspace(0,2*pi, resolution) ;
        Y = scale*pinRadius * cos(t) ;
        X = zeros(size(Y)) ;
        Z = scale*pinRadius * sin(t) ;
        C = zeros([size(X) 3]) ;
        C(:,:,1) = pinColor(1) ;
        C(:,:,2) = pinColor(2) ;
        C(:,:,3) = pinColor(3) ;
        hdisk1 = patch(X,Y,Z,pinColor,'linestyle',pinLineStyleStr) ;
        
        % disk2
        X = X - scale*pinLength ;
        hdisk2 = patch(X,Y,Z,pinColor,'linestyle',pinLineStyleStr) ;
        
        pinGrp = hgtransform ;
        set([pin hdisk1 hdisk2],'Parent',pinGrp) ;
        R1 = makehgtform('yrotate',-pi/4) ;
        T1 = makehgtform('translate',[0 0 -scale*thoraxRadius]) ;
        T1R1 = T1*R1 ;
        set(pinGrp,'Matrix',T1R1) ;
        
        set(pinGrp,'Parent',thoraxGrp)
end
% antennea

% }--- WINGS ---{
% right wing

dy = LAMBDA*wingSpan*scale;
rightWing = surf(Xsph*scale*wingChord-wingChord*scale, ... 
    Ysph*scale*wingSpan - dy, ...
    Zsph*scale*wingThickness );
C = zeros([size(Xsph) 3]) ;
C(:,:,1) = rightWingColor(1) ;
C(:,:,2) = rightWingColor(2) ;
C(:,:,3) = rightWingColor(3) ;
set(rightWing,'CData',C,'linestyle',wingLineStyle) ;            

rightWingTip = [-wingChord*scale ; - scale*wingSpan - dy ; 0];
% right wing vein
[X, Z, Y ] = cylinder(veinRadius, resolution) ;
Y = Y * scale * wingSpan * LAMBDA - dy;
X = X * scale ;
Z = Z * scale ;
rightVein = surf(X,Y,Z) ;
C = zeros([size(X) 3]) ;
C(:,:,1) = rightWingColor(1) ;
C(:,:,2) = rightWingColor(2) ;
C(:,:,3) = rightWingColor(3) ;
set(rightVein,'CData',C,'linestyle',wingLineStyle) ; 

% disk that seals the right wing vein
t = linspace(0,2*pi, resolution) ;
X = scale*veinRadius * cos(t) ;
Y = zeros(size(X))  - dy ;
Z = scale*veinRadius * sin(t) ;
C = zeros([size(X) 3]) ;
C(:,:,1) = rightWingColor(1) ;
C(:,:,2) = rightWingColor(2) ;
C(:,:,3) = rightWingColor(3) ;
hdisk3 = patch(X,Y,Z,rightWingColor,'linestyle',wingLineStyle) ;

Y = Y + scale * wingSpan * LAMBDA ;
hdisk4 = patch(X,Y,Z,rightWingColor,'linestyle',wingLineStyle) ;

alpha([rightWing rightVein hdisk3 hdisk4],alphaVal) ;
set([rightWing rightVein hdisk3 hdisk4], 'Parent',rightWingGrp) ;

%R1 = makehgtform('yrotate',thetab0) ;
%T1 = makehgtform('translate',[0 -thoraxRadius*scale 0  ]) ;
dL = thoraxRadius*scale ;
%set(rightWingGrp,'Matrix',R1) ;


leftWing = surf(Xsph*scale*wingChord-wingChord*scale, ... 
    Ysph*scale*wingSpan + dy, ...
    Zsph*scale*wingThickness );
C = zeros([size(Xsph) 3]) ;
C(:,:,1) = leftWingColor(1) ;
C(:,:,2) = leftWingColor(2) ;
C(:,:,3) = leftWingColor(3) ;
set(leftWing,'CData',C,'linestyle',wingLineStyle) ;            

leftWingTip = [rightWingTip(1) ; -rightWingTip(2) ; rightWingTip(3)] ;
% left wing vein
[X, Z, Y ] = cylinder(veinRadius, resolution) ;
Y = Y * scale * wingSpan * LAMBDA - 0*dy;
X = X * scale ;
Z = Z * scale ;
leftVein = surf(X,Y,Z) ;
C = zeros([size(X) 3]) ;
C(:,:,1) = leftWingColor(1) ;
C(:,:,2) = leftWingColor(2) ;
C(:,:,3) = leftWingColor(3) ;
set(leftVein,'CData',C,'linestyle',wingLineStyle) ; 

% disks that seals the left wing vein
t = linspace(0,2*pi, resolution) ;
X = scale*veinRadius * cos(t) ;
Y = zeros(size(X))  - 0*dy ;
Z = scale*veinRadius * sin(t) ;
C = zeros([size(X) 3]) ;
C(:,:,1) = leftWingColor(1) ;
C(:,:,2) = leftWingColor(2) ;
C(:,:,3) = leftWingColor(3) ;
hdisk5 = patch(X,Y,Z,leftWingColor,'linestyle',wingLineStyle) ;

Y = Y + scale * wingSpan * LAMBDA;
hdisk6 = patch(X,Y,Z,leftWingColor,'linestyle',wingLineStyle) ;

alpha([leftWing leftVein hdisk5 hdisk6],alphaVal) ;
set([leftWing leftVein hdisk5 hdisk6], 'Parent',leftWingGrp) ;

end

