
startFlipInd = 1;
endFlipInd = find(data.fwdFlipTimesR > 0,1,'first') - 1;
strokeNum = startFlipInd - endFlipInd ;

oldTorqueVec = zeros(strokeNum,1) ;
forceMat = zeros(strokeNum,2) ;
%A = [1 0; 0 -1] ; %constraints
%b = [.0005 ; 0] ;%[.0005 ; .0005] ;
lb = [-0.0006 -0.0006] ;
ub = [0.0006 0.0006] ;

F_transL = data.F_transL ;
F_transR = data.F_transR ;

try
    pitchTorque = data.pitchTorque ;
catch
    torques = data.torques ;
    axisHats = data.AHat ;

    bodyY = cross(axisHats,repmat([0 0 1],length(torques),1)) ;
    bodyY = bodyY ./ repmat(myNorm(bodyY),1,3) ;

    pitchTorque = dot(torques,bodyY,2) ;
end
axisHats = data.AHat ;
t = data.t ;
N = length(t) ;
backFlipTimes = data.backFlipTimesR ; %arbitrary
bodyCM = data.bodyCM ;

bodyY = cross(axisHats,repmat([0 0 1],N,1)) ;
bodyY = bodyY ./ repmat(myNorm(bodyY),1,3) ;

bodyX = cross(bodyY,axisHats) ;

F_x = dot((F_transL + F_transR),bodyX,2) ;
F_z = dot((F_transL + F_transR),axisHats,2) ;

for i = startFlipInd:endFlipInd
    [~,i1] = min(abs(t - backFlipTimes(i))) ;
    [~,i2] = min(abs(t - backFlipTimes(i+1))) ;
    
    oldTorqueVec(i) = sum(pitchTorque(i1:i2)) ;
    forceMat(i,1) = -sum(F_z(i1:i2)) ;
    forceMat(i,2) = sum(F_x(i1:i2)) ;
end

[epsilon,resnorm,residual,exitflag,output,lambda] = lsqlin(forceMat,oldTorqueVec,[],[],...
    [],[],lb,ub,[0;0]) ;

%epsilon is in form given by notes, i.e. epsilon(1) is the component in the
%x direction (I call the body axis 'z' and the pitch torque axis 'y') and
%epsilon(2) is in the body axis direction (z)

df = diff(data.res(:,1)) ;
frameStartInd = [1 ; find(df==1)+1] ;
frameEndInd   = [frameStartInd(2:end)-1 ; size(data.res,1)] ;
clear df ;

FrameRef= 200 ;
bodyCM_new = bodyCM(FrameRef,:) + epsilon(1)/(50e-6)*bodyX(FrameRef,:)...
    + epsilon(2)/(50e-6)*axisHats(FrameRef,:) ; 
    

row1       = frameStartInd(FrameRef) ;
row2       = frameEndInd(FrameRef) ;
coords     = data.res(row1:row2,2:4) ; % xyz
IDX        = data.RESIDX(row1:row2,:) ;
bodyRows   = (IDX(:,data.bodyInd)==1) ;
bodyVoxels = coords(bodyRows,:) ;
Nvox       = sum(bodyRows) ;

figure;
plot3(bodyVoxels(:,1),bodyVoxels(:,2),bodyVoxels(:,3),'bo');
axis equal;
hold on;
plot3(bodyCM(FrameRef,1),bodyCM(FrameRef,2),bodyCM(FrameRef,3),'gx',...
    'MarkerSize',20,'LineWidth',6)
plot3(bodyCM_new(1),bodyCM_new(2),bodyCM_new(3),'rx',...
    'MarkerSize',20,'LineWidth',6)

