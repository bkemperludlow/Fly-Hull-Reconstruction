%plots chord vectors from idealizedWingStroke. This is a mess. Need to
%clean up before real use. Should be made into a function

%% Initial preferences. Needs to be adjusted for each use. These will be function inputs
wingSide = 'R' ; %adjust this for different wings 
spacing = 1 ; %in units of frame number
flipInd = 18 ; %Pick which halfstroke you want to look at. Indexed by flipTimes. If 0, you just get the first (potentially partial) halfstroke
strokeType = 'fwd' ;
overlay = 0 ; %This is either 1 or 0. If 0, make a new plot. If 1, plot on an extant one

%% Load in data about the particular wingstroke. Should work for real data too
if isfield(data.params,'span')
    span = data.params.span ;
else 
    span = .002 ; %m
end
if isfield(data.params,'chord')
    chord = data.params.chord ;
else
    chord = .0007 ; %m
end
dt = 1/data.params.fps ;
startTrackingTime = data.params.startTrackingTime ;

if wingSide == 'R'
    fwdFlipTimes = data.fwdFlipTimesR ;
    backFlipTimes = data.backFlipTimesR ;
    wingTip = data.rightWingTips ;
    chordHat = data.rightChordHats ;
    spanHat = data.rightSpanHats ;
    F_rot = data.F_rotR ;          
    %F_lift = data.F_liftR ;
    %F_drag = data.F_dragR ;
    F_trans = data.F_transR ;
    %U_t = data.U_tR ;
    %torque = data.torqueR ;
elseif  wingSide == 'L'
    fwdFlipTimes = data.fwdFlipTimesL ;
    backFlipTimes = data.backFlipTimesL ;
    wingTip = data.leftWingTips ;
    chordHat = data.leftChordHats ;
    spanHat = data.leftSpanHats ;
    F_rot = data.F_rotL ;      
    %F_lift = data.F_liftL ;
    %F_drag = data.F_dragL ;
    F_trans = data.F_transL ;
    %U_t = data.U_tL ;
    %torque = data.torqueL ;
else
    disp('Check wingSide')
end

%THIS IS WHAT'S ACTUALLY PLOTTED
bodyCM = data.bodyCM * 50e-6 ; %now in real units
base = wingTip*50e-6 - (span/2)*spanHat ; %wingTip needs to be multiplied by 50e-6 for real data
scaledChord = (chord/2)*chordHat ; 
scaledSpan = (span/2)*spanHat ;
%U_tHat = U_t ./ repmat(myNorm(U_t),1,3) ;

%% Find time interval over which chords will be plotted

%This is a little awkward, but I just want to pick out one half-stroke.
%Below, I pick out the flip time corresponding to t1, then I find the next
%(opposite flip time) later. So if I want to look at e.g. the first forward
%stroke, I take flipInd = 1 and strokeType = 'fwd', find the first backward
%flip, and then find the next forward flip time after that 
if strcmp(strokeType,'fwd')
    if flipInd == 0
        t1 = startTrackingTime ;
    else
        t1 = backFlipTimes(flipInd) ;
    end
    laterFlips = find(fwdFlipTimes > t1) ;
    t2 = fwdFlipTimes(laterFlips(1)) ;
elseif strcmp(strokeType,'back')
    if flipInd == 0
        t1 = startTrackingTime ;
    else
        t1 = fwdFlipTimes(flipInd) ;
    end
    laterFlips = find(backFlipTimes > t1) ;
    t2 = backFlipTimes(laterFlips(1)) ;
else
    disp('Check strokeType variable')
end

startFrame = int16(t1/dt - startTrackingTime + 2) ;
endFrame = int16(t2/dt - startTrackingTime + 2) ;
N = (endFrame - startFrame)/int16(spacing) + 1; 

%% Make the color change depending on which wing we're looking at.
%Unnecessary, but, then again, w/e. Black is early time, bright is later
%To change: different color for backstroke vs forward stroke?
colorShift = linspace(0,1,N) ;
if wingSide == 'R'
    color = [colorShift', zeros(N,2)] ;
elseif  wingSide == 'L'
    color = [zeros(N,2),  colorShift'] ; 
else
    disp('Check wingSide')
end
%% Defines the figure

if overlay == 0
    chordplot = figure ;
elseif overlay == 1
    try 
        figure(chordplot) ;
    catch 
        disp('NO EXTANT CHORD PLOT')
    end
else
    disp('Check overlay variable. Needs to be either 0 or 1')
end

%myPlotVector([0 0 0], .5e-3*data.AHat(1,:),[.5 .5 .5], 3) ; %plots body axis hat. Need to fix if this becomes dynamic
hold on

j = 1; %My indexing is all fucked up, don't want to figure it out now

for i = startFrame:spacing:endFrame ;
    myPlotVector(base(i,:),scaledChord(i,:),color(j,:),2.5) ; %chord vector
    myPlotVector(bodyCM(i,:),scaledSpan(i,:),color(j,:)) ; %span vector
    myPlotVector(base(i,:),50*F_rot(i,:),[0 .5 0],2,'--') ;
    %myPlotVector(base(i,:),50*F_trans(i,:),[0 .5 0],2,'--') ; %total translational force
    %myPlotVector(base(i,:),50*F_lift(i)*[0 0 1],[0 .5 0]) ; %lift
    %myPlotVector(base(i,:),-50*F_drag(i)*U_tHat(i,:),[.5 .5 0]) ; %drag
    %myPlotVector(bodyCM(i,:),.5e5*torque(i,:),color(j,:)) ;  %total torque
    j = j+1 ;
end

%set(gca, 'fontsize', 14) ;
xlabel('X') ; ylabel('Y') ; zlabel('Z') ; 
box on; grid on;
axis equal ;
axis tight ;
    
    
