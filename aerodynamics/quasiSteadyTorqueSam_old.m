%Find torques on the fly body due to quasi-steady forces.
%We assume that forces act at the center of pressure (CP) which is 70% of 
%the length along the line from base of the wing to tip of the wing (Fry+Sayaman+Dickinson, 2005)
%
%By Sam Whitehead
%-------------------------------------------------------------------------------------
%% Load data

plotFlag = 1 ;

%datafilename = 'idealizedWingstrokeData_phasediff.mat' ; 
%load(datafilename) ;

defineConstantsScript
%{
exprNum = 13 ;
movNum  = 7 ;

if (movNum<10)
    zstr = '00' ;
elseif (movNum<100)
    zstr = '0' ;
else
    zstr = '';
end

datapath = ['F:\luca\Analysis\pitch up\Expr_' num2str(exprNum) '_mov_' zstr num2str(movNum) '\' ] ;
cd(datapath)

datafilename = [ datapath ...
    'Expr' num2str(exprNum) 'mov' zstr num2str(movNum) '_Data_manually_corrected.mat' ] ; %

load(datafilename) ;
%}
if (isfield(data,'ignoreFrames'))
    ignoreFrames = data.ignoreFrames ;
else
    ignoreFrames = [] ;
end

%{
if (isfield(data,'correctionTime'))
    correctionTime=data.correctionTime;
elseif(isfield(data,'manualCorrRangeMS'))
    correctionTime=data.manualCorrRangeMS;
end
%}

startTime = data.params.startTrackingTime ;
endTime = data.params.endTrackingTime ;

span = .0025 ; %in meters
chord = .0007 ; %in meters

%% Find forces

[F_liftR, F_dragR, F_rotR, wingTip_smoothR, U_t_projHatR, spanHatR, chordHat_projR, t] = quasiSteadyForceSam_old(data,'R',plotFlag) ; 
[F_liftL, F_dragL, F_rotL, wingTip_smoothL, U_t_projHatL, spanHatL, chordHat_projL, ~] = quasiSteadyForceSam_old(data,'L',plotFlag) ;

%% Find centers of pressure

%bodyCM = data.bodyCM_old*50e-6 ; %Need to use bodyCM_old for Expr7mov65 for some reason
%epsilonVec = data.AHat - repmat([0 0 1],size(data.AHat,1),1).*repmat(dot(data.AHat,repmat([0 0 1],size(data.AHat,1),1),2),1,3) ;
%epsilonHat = epsilonVec./repmat(myNorm(epsilonVec),1,3) ;
%epsilon = epsilonHat*.0001615 - .0001615*repmat([0 0 1],size(epsilonHat,1),1) ;
bodyCM = data.bodyCM*50e-6 ;%+ epsilon(1)*bodyX + epsilon(2)*axisHats ;
CPR = wingTip_smoothR*50e-6 - .3*span*spanHatR - bodyCM ; %+ epsilon ; %Need to watch out for units here
CPL = wingTip_smoothL*50e-6 - .3*span*spanHatL - bodyCM ;%+ epsilon ; 

%% Calculate torques
% We assume that drag is antiparallel to wing tip velocity and lift is
% perpendicular to both drag and span

liftVecR = cross(spanHatR,U_t_projHatR) ; 
indR = find(liftVecR(:,3) < 0) ; 
if (~isempty(indR))
    liftVecR(indR,:) = -1*liftVecR(indR,:) ;
end
FHatR = (liftVecR - U_t_projHatR)./repmat(myNorm(liftVecR - U_t_projHatR),1,3) ;

liftVecL = cross(spanHatL,U_t_projHatL) ; 
indL = find(liftVecL(:,3) < 0) ; 
if (~isempty(indL))
    liftVecL(indL,:) = -1*liftVecL(indL,:) ;
end
FHatL = (liftVecL - U_t_projHatL)./repmat(myNorm(liftVecL - U_t_projHatL),1,3) ;

F_transR =  - repmat(F_dragR,1,3).*U_t_projHatR  + repmat(F_liftR,1,3).*liftVecR ;
F_transL =  - repmat(F_dragL,1,3).*U_t_projHatL  + repmat(F_liftL,1,3).*liftVecL ;
F_rotRVec = repmat(F_rotR,1,3).*FHatR ;
F_rotLVec = repmat(F_rotL,1,3).*FHatL ;

torqueR = cross(CPR, F_transR) ; %added minus sign 13/10/14 to account for error in idealizedwingstroke
torqueL = cross(CPL, F_transL) ; 

torque = torqueR + torqueL ;

%Save things to data

data.F_liftR = F_liftR ;        %right wing lift coefficients
data.F_liftL = F_liftL ;        %left wing lift coefficients
data.F_dragR = F_dragR ;        %right wing drag coefficients
data.F_dragL = F_dragL ;        %left wing drag coefficients
data.F_transR = F_transR ;      %right wing total translational force vector
data.F_transL = F_transL ;      %left wing total translational force vector
data.F_rotR = F_rotR ;          %right wing rotational force magnitude
data.F_rotL = F_rotL ;          %left wing rotational force magnitude
data.torqueR = torqueR ;        %right wing torque
data.torqueL = torqueL ;        %left wing torque
data.torques = torque ;         %total torque
data.t = t' ;                   %time

%save(datafilename, 'data')

%% Make figure
% N = length(t) ;
% htorque = figure ;
% 
% 
% for frameNum = 383:N 
% 
%     Axvec = [-100*data.AHat(frameNum,1) 100*data.AHat(frameNum,1)] ;
%     Ayvec = [-100*data.AHat(frameNum,2) 100*data.AHat(frameNum,2)] ;
%     Azvec = [-100*data.AHat(frameNum,3) +100*data.AHat(frameNum,3)] ;
%     plot3(Axvec,Ayvec,Azvec,'k-','LineWidth',4)
%     box on ; grid on ;
%     hold on ;
%     
%     Rwingvec_x = [0 100*spanHatR(frameNum,1)] ;
%     Rwingvec_y = [0 100*spanHatR(frameNum,2)] ;
%     Rwingvec_z = [0 100*spanHatR(frameNum,3)] ;
%     
%     Lwingvec_x = [0 100*spanHatL(frameNum,1)] ;
%     Lwingvec_y = [0 100*spanHatL(frameNum,2)] ;
%     Lwingvec_z = [0 100*spanHatL(frameNum,3)] ;
%     
%     plot3(Rwingvec_x,Rwingvec_y,Rwingvec_z,'r-','LineWidth',3.5)
%     plot3(Lwingvec_x,Lwingvec_y,Lwingvec_z,'b-','LineWidth',3.5)
%     
%     torquevec_x = [0 torque(frameNum,1)*1e10] ;
%     torquevec_y = [0 torque(frameNum,2)*1e10] ;
%     torquevec_z = [0 torque(frameNum,3)*1e10] ;
%     
%     plot3(torquevec_x,torquevec_y,torquevec_z,'Color',[0 .5 0],'LineWidth',3)
%     xlabel('X [vox]')
%     ylabel('Y [vox]')
%     zlabel('Z [vox]')
%     title(['t = ' num2str(1000*t(frameNum)) ' ms'])
%     axis([-100 100 -100 100 -100 100])
%     pause(1) ;
%     cla
% end