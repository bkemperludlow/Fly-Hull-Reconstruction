%% Set inputs

defineConstantsScript

ExprNum = 7 ;
MovNum = 8 ; 
PitchType = 'up' ;
wingSide = 'L' ;
smoothFlag = false ;

%% Load data

if MovNum < 10
    zstr = '00';
elseif MovNum < 100
    zstr = '0';
else
    zstr = '';
end

if strcmp(PitchType,'up')
    datapath = strcat('F:\luca\Analysis\pitch up\Expr_',...
        num2str(ExprNum), '_mov_',zstr,num2str(MovNum)) ;
elseif strcmp(PitchType,'down')
    datapath = strcat('F:\luca\Analysis\pitch down\Expr_',...
        num2str(ExprNum), '_mov_',zstr,num2str(MovNum)) ;
else
    disp('check PitchType')
    return ;
end
cd(datapath)

datafilename = strcat(datapath,'\Expr',num2str(ExprNum), ...
    'mov',zstr,num2str(MovNum), '_Data_manually_corrected.mat') ;

load(datafilename) ;

%% Check to make sure 'data' has the right fields
if (isfield(data,'ignoreFrames'))
    ignoreFrames = data.ignoreFrames ;
else
    ignoreFrames = [] ;
end

if (isfield(data,'correctionTime'))
    correctionTime=data.correctionTime;
elseif(isfield(data,'manualCorrRangeMS'))
    correctionTime=data.manualCorrRangeMS;
end

if isfield(data.params,'pulseLengthMS')
    pulseLength = data.params.pulseLengthMS ;
elseif (~isfield(data.params,'pulseLengthMS')) && (ExprNum == 7)
    pulseLength = 5.8 ; %ms
else
    pulseLength = 8 ; %ms
end

%% Get wing times and angles
if isfield(data,'fwdFlipTimesR') && isfield(data,'backFlipTimesR')
    fwdFlipTimesR = data.fwdFlipTimesR ;
    backFlipTimesR = data.backFlipTimesR ;
    fwdFlipTimesL = data.fwdFlipTimesL ;
    backFlipTimesL = data.backFlipTimesL ;
else
    [fwdFlipTimesR, backFlipTimesR, fwdFlipTimesL, backFlipTimesL, ~,~, data]...
        = saveWingFlipsAndAngles(ExprNum,MovNum,PitchType) ;
end

t = data.t ; 

if strcmp(wingSide,'L')
    backFlipTimes = backFlipTimesL ; 
    fwdFlipTimes = fwdFlipTimesL ;
    phi = data.anglesBodyFrame(:,PHIL) ; 
    eta = data.anglesBodyFrame(:,ETAL) ; 
elseif strcmp(wingSide,'R')
    backFlipTimes = backFlipTimesR ; 
    fwdFlipTimes = fwdFlipTimesR ; 
    phi = -data.anglesBodyFrame(:,PHIR) ; 
    eta = data.anglesBodyFrame(:,ETAR) ; 
else
    disp('check wingSide')
    return ;
end

for q = 1:length(eta)
    while eta(q) > 360
        eta(q) = eta(q) - 360 ;
    end
    while eta(q) < 0
        eta(q) = eta(q) + 360 ;
    end
end
%% Find times during maneuver and before perturbation
% I should probably make the time intervals I'm defining more rigorous
i2 = find(t == 0) ;
try 
    i1 = i2 - 120 ;
catch
    i1 = 1 ; 
end

t_maneuver1 = .015 ;
t_maneuver2 = .025 ;
[~,i3] = min(abs(t - t_maneuver1)) ; 
[~,i4] = min(abs(t - t_maneuver2 )) ; 

%% Make figure
h = figure ; 
set(h,'name',['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum) ' etaVSphi'],'numbertitle','off')
plot(phi(i1:i2),eta(i1:i2),'Color',[0 .5 0],'LineWidth',1.5)
hold on
plot(phi(i3:i4),eta(i3:i4),'m','LineWidth',1.5)


axis equal
set(gca,'xlim',[0 180])
set(gca,'ylim',[0 180])
grid on

xlabel('\phi')
ylabel('\eta')
title(['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum) ' Wing ' wingSide])
legend({'Pre-Pert.',[num2str(t_maneuver1) '-' num2str(t_maneuver2) ' ms']},'location','northeast')
    
    
    
    
    
    
    
    
