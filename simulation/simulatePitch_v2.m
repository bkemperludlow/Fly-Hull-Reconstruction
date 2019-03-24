%% Runs an ode solver on bodyPitchODE.m or bodyPitchODE_CTRL.m
%Notes:
%       -
%
%--------------------------------------------------------------------------
%% Initialize
simType = 2 ; % 0 = uncontrolled, 1 = controlled, 2 = controlled with perturbation
plotFlag = true ;
movieFlag = false ; 
saveFlag = true ;

savePath = 'F:\luca\Analysis\simulation\NoControl' ;

%Time range
ti          =   0 ; %seconds
tf          =   1 ;
dt          =   .001 ; 
tspan = ti:dt:tf ;

%Controller model:
deltaT = .006 ; %seconds
K_i = .3 ; %.3 ;
K_p = .007 ; %.007 ; %seconds

%perturbation characteristics
pulseStart = .350 ; %seconds
pulseEnd = .355 ;
pulseStrength = 7e3 ; %rad/s^2 (so it's an acceleration not a torque). Roughly estimated from movies %7e3

%Specify initial values
y_0         =   0 ;
ydot_0      =   0 ;
z_0         =   0 ;
zdot_0      =   0 ;
theta_0     =   pi/4 ;
thetadot_0  =   0 ;
s_0 = [y_0, ydot_0, z_0, zdot_0, theta_0, thetadot_0] ;

%Define color for array
N = length(tspan) ;
colorArray = linspace(0,1,N);
colorVec = [colorArray' zeros(N,2) ] ;

%% Cases for different simulation types
switch simType
    %Uncontrolled
    case 0
        sol = ode45(@bodyPitchODE, [ti tf], s_0);
        titlestr = '(Uncontrolled)' ;
        
    %Controlled
    case 1 
        pertFlag = false ; 
        sol = dde23(@(t,s,Z)bodyPitchODE_CTRL(t,s,Z,K_i,K_p,...
            pertFlag,pulseStart, pulseEnd, pulseStrength),...
            deltaT, s_0, [ti tf] );
        titlestr = '(Controlled -- No Perturbation)' ;
    
    %Controlled with perturbation
    case 2 
        pertFlag = true ; 
        sol = dde23(@(t,s,Z)bodyPitchODE_CTRL(t,s,Z,K_i,K_p,...
            pertFlag,pulseStart, pulseEnd, pulseStrength),...
            deltaT, s_0, [ti tf] );
        titlestr = '(Controlled -- With Perturbation)' ;
end

sint = deval(sol, tspan) ; %evaluate solution at tspan
theta = sint(5,:) ; %radians
y = sint(1,:) ; % meters
z = sint(3,:) ; % meters

if plotFlag
    %Plot
    htrajectory = figure ;
    set(gca,'fontsize',12)
    hold on
    for i = 1:N
        plot(1000*y,1000*z,'ko','MarkerSize',6,'MarkerFaceColor',colorVec(i,:))
    end
    xlabel('y [mm]')
    ylabel('z [mm]')
    title(['CoM Trajectory ' titlestr])
    
    hpitch = figure ;
    set(gca,'fontsize',12)
    plot(tspan*1000,(180/pi)*theta, 'ko','MarkerSize',6,...
        'MarkerFaceColor',[0 0 .8])
    xlabel('t [ms]')
    ylabel('\theta [deg]')
    title(['Pitch Angle vs Time ' titlestr])
end

if movieFlag
   bodyPitchODE_makeMovie(tspan,y,z,theta,savePath) ;
end

if saveFlag
    cd(savePath)
    controlsim = struct ;
    controlsim.tspan = tspan ;
    controlsim.sint = sint ;
    controlsim.pulseStart = pulseStart ;
    controlsim.pulseEnd = pulseEnd ;
    controlsim.pulseStrength = pulseStrength ;
    controlsim.K_i = K_i ;
    controlsim.K_p = K_p ;
    controlsim.deltaT = deltaT ;
    controlsim.s_0 = s_0 ;
    
    save controlsim controlsim
end




