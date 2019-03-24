%% beginning stuff
mainfig = figure ;
ax = gca ;
hold on

defineConstantsScript

Expr = [7 7 7 7 7 7 ];
Mov = [7 8 9 15 19 54] ;
PitchType = 'up' ;
angleType = 1 ;

%   angleType - integer that indicates which type of angle to plot
%       0 = stroke angle
%       1 = wing pitch
%       2 = deviation angle

dt = 1/8000 ;

runSize = length(Expr) ; 

for i = 1:runSize

    %% Load data
    if Mov(i) < 10
        zstr = '00';
    elseif Mov(i) < 100
        zstr = '0';
    else
        zstr = '';
    end
    
    if strcmp(PitchType,'up')
        datapath = strcat('F:\luca\Analysis\pitch up\Expr_',...
            num2str(Expr(i)), '_mov_',zstr,num2str(Mov(i))) ;
    elseif strcmp(PitchType,'down')
        datapath = strcat('F:\luca\Analysis\pitch down\Expr_',...
            num2str(Expr(i)), '_mov_',zstr,num2str(Mov(i))) ;
    else
        disp('check PitchType')
        return ;
    end
    cd(datapath)
    
    datafilename = strcat(datapath,'\Expr',num2str(Expr(i)), ...
        'mov',zstr,num2str(Mov(i)), '_Data_manually_corrected.mat') ;
    
    load(datafilename) ;
    
    %% Get relevant variables and smooth
    fwdFlipTimesR = data.fwdFlipTimesR ;
    backFlipTimesR = data.backFlipTimesR ;
    if isfield(data,'t')
        t = data.t ;
    else
        t1_temp = data.params.startTrackingTime/8000 ;
        t2_temp = data.params.endTrackingTime/8000 ;
        dt_temp = 1/data.params.fps ;
        t = t1_temp:dt_temp:t2_temp ;
    end
    if (isfield(data,'ignoreFrames'))
        ignoreFrames = data.ignoreFrames ;
    else
        ignoreFrames = [] ;
    end
   
    switch angleType
        case 0
            angleR = -data.anglesBodyFrame(:,PHIR) ;
            angleL = data.anglesBodyFrame(:,PHIL) ;
            EstErr = 1 ;
            labelStr = 'Wing Stroke Angle' ;
            
        case 1
            angleR = (180/pi)*mod(unwrap((pi/180)*data.anglesBodyFrame(:,ETAR)),2*pi) ;
            angleL = (180/pi)*mod(unwrap((pi/180)*data.anglesBodyFrame(:,ETAL)),2*pi) ;
            EstErr = 1 ;
            labelStr = 'Wing Pitch Angle' ;
            
        case 2
            angleR = data.anglesBodyFrame(:,THETAR) ;
            angleL = data.anglesBodyFrame(:,THETAL) ;
            EstErr = 2 ;
            labelStr = 'Wing Deviation Angle' ;
    end
    
    if (~isempty(ignoreFrames))
        angleR(ignoreFrames) = NaN ;
        angleL(ignoreFrames) = NaN ;
    end
    
    indR = find(~isnan(angleR)) ;
    indL = find(~isnan(angleL)) ;
    ind = intersect(indR,indL) ;
    currtvecR = t(ind) ;
    currtvecL = t(ind) ;
    currangleR = angleR(ind) ;
    currangleL = angleL(ind) ;
    
    [sp_angleR, angleR_smooth, ~] = mySplineSmooth(currtvecR,currangleR,EstErr) ;
    [sp_angleL, angleL_smooth, ~] = mySplineSmooth(currtvecL,currangleL,EstErr) ;
    
    %% Fix phase relationship
    firstFlipInd = find(fwdFlipTimesR > -.01, 1, 'first') ;
    t_diff = abs(t - fwdFlipTimesR(firstFlipInd)) ;
    [~,i1] = min(t_diff) ;
    
    angle = .5*(angleR_smooth(i1:end)+angleL_smooth(i1:end)) ;
    Ntemp = length(angle) ;
    tms = -10:(1000*dt):(-10+(Ntemp-1)*1000*dt) ;
    
    hold on
    plot(tms,angle,'-','LineWidth',1,'Color',[1 1 1]*.8)
    

end

