function hmetadata = plotMetaData_mk2(XVarName,YVarName,xlim,...
    fitFlag,textFlag)
%hfront = plotMetaData_mk2('PhiFront','pitchAcceleration',[-23 30],true, false) ;
%-----------------------------------------------------------------------
%Plots metadata from tables.
%
%Potential variable names:
%   -PhiFront
%   -PhiBack
%   -PhiAmp
%   -PhiMid
%   -DeltaPitch
%   -pitchAcceleration
%   -raw versions of the above (e.g. rawPhiFront)
%-----------------------------------------------------------------------
hmetadata = figure ;
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 4 4]);

if (~exist('textFlag','var'))
    textFlag = false ;
end
%fitFlag = true ;

textDistX = 1 ;
textDistY = 5e4 ;

cd('F:\luca\Analysis\metadata\Stroke Angle\')
pitchUp = importdata('pitchUp10302014.mat');
pitchDown = importdata('pitchDown10302014.mat');
T = union(pitchUp, pitchDown) ; 
%T = pitchDown ; 

%XVarName = 'PhiFront' ;
%YVarName = 'pitchAcceleration' ;
XErrVal = 1 ; %Assume degrees. Should fix this to make more general
YErrVarName = 'AccelFitErr' ;

XVarInd = find(strcmp(T.Properties.VariableNames,XVarName)) ;
if isempty(XVarInd) 
    disp('check X variable name')
    T.Properties.VariableNames
    return ;
end
YVarInd = find(strcmp(T.Properties.VariableNames,YVarName)) ;
if isempty(YVarInd)
    disp('check Y variable name')
    T.Properties.VariableNames
    return ;
end
YErrVarInd = find(strcmp(T.Properties.VariableNames,YErrVarName)) ;

if isempty(YErrVarInd)
    disp('check error variable name')
    T.Properties.VariableNames
    return ;
end

X = T{:,XVarInd} ;
Y = T{:,YVarInd} ;
ErrY = T{:,YErrVarInd} ;
ErrX = XErrVal*ones(size(X)) ;

%Do the linear fit if wanted
if fitFlag
    weights = ErrY.^(-1) ;
    A = [X, ones(size(X))] ;
    [fitcoeffs, stdx, mse, S] = lscov(A,Y) ;
    hline = plot([min(X), max(X)], 1e-5*[fitcoeffs(1)*min(X)+fitcoeffs(2), fitcoeffs(1)*max(X)+fitcoeffs(2)],...
        'Color',[1 1 1]*.5,'LineWidth',3) ;
    sse = sum((Y - A*fitcoeffs).^2);
    sst = sum((Y - mean(Y)).^2);
    rsq = 1 - sse/sst;
end

if strcmp(XVarName,'PhiFront')
    PlotColor = [0.91 0.41 0.17] ; %deep carrot orange
elseif strcmp(XVarName,'PhiBack')
    PlotColor = [0.5859 0.4805 0.7109] ; %lavender
else
    PlotColor = [0 0 0] ;
end
hold on
h2 = herrorbar(X,Y/1e5,ErrX) ; 
h1 = errorbar(X,Y/1e5,ErrY/1e5,'ko','MarkerSize',8,'MarkerFaceColor',PlotColor,'LineWidth',2) ;
%lavender =  [0.5859 0.4805 0.7109]
%deep carrot orange = [0.91 0.41 0.17]
set(h2(2),'linestyle','none') ;
set(h2(1),'Color',[0 0 0]) ;
set(h2(1),'LineWidth',2) ;
set(gca,'fontsize',10)

xlabel(XVarName); 
ylabel(YVarName);
%title('Wing Actuation');
%legend([h1 hline],{'Data',['R^2 = ' num2str(rsq,2)]},'location','northeast') 
%legend([h1 hline],{'D','R'},'location','northeast') 
set(gca,'xlim',xlim) ;
%xlabel('\phi_{front} [deg]')
%ylabel('Pitch Acceleration [deg/s^2]')
%title('Pitch Acceleration vs Front Stroke Angle')
if textFlag
    for i = 1:length(X)
        text(X(i)+textDistX,(Y(i)+textDistY)/1e5,[num2str(T.Expr(i)) ...
            '.' num2str(T.Mov(i))])
    end
end





