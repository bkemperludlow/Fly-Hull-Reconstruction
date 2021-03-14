% -------------------------------------------------------------------------
% script to run minimization of forces/torques for parameterized kinematics
% -------------------------------------------------------------------------
% -------------------------------
%% angle conversions
DEG2RAD = (pi/180) ; 
RAD2DEG = (180/pi) ; 

% -------------------------------------------
%% define initial values
params = defineQuasiSteadyParams ; 
omega0 = params.omega ; 
phi_f0 = DEG2RAD*(params.phi_f ); 
phi_b0 = DEG2RAD*(params.phi_b) ;
psi_m0 = DEG2RAD*(params.psi_m) ;
psi_00 = DEG2RAD*(params.psi_0) ;
del_psi0 = params.del_psi ;
span0 = params.span ; 
hinge_x0 = params.r_hinge ; 

x0 = [omega0, phi_f0, phi_b0, psi_m0, psi_00, del_psi0, span0, hinge_x0] ; 

% varyFlag = [1, 1, 1, 1, 1, 1, 1, 1] ; 
varyFlag = [1, 1, 1, 1, 1, 1, 0, 0] ; 
% -----------------------------
%% set solver constraints
% angular frequency
if varyFlag(1)
    omega_min = 2*pi*180 ; 
    omega_max = 2*pi*270 ; 
else
    omega_min = omega0 ; 
    omega_max = omega0 ; 
end

% ventral stroke amp (rad)
if varyFlag(2)
    phi_f_min = DEG2RAD*10 ; 
    phi_f_max = DEG2RAD*50 ; 
else
    phi_f_min = phi_f0  ; 
    phi_f_max = phi_f0  ; 
end

% dorsal stroke amp (rad)
if varyFlag(3)
    phi_b_min = DEG2RAD*140 ; 
    phi_b_max = DEG2RAD*190 ;
else
    phi_b_min = phi_b0  ; 
    phi_b_max = phi_b0  ; 
end

% wing rotation amplitude (rad)
if varyFlag(4)
    psi_m_min = DEG2RAD*45 ; 
    psi_m_max = DEG2RAD*70 ;
else
    psi_m_min = psi_m0  ; 
    psi_m_max = psi_m0  ; 
end

% wing rotation midpoint (rad)
if varyFlag(5)
    psi_0_min = DEG2RAD*70 ; 
    psi_0_max = DEG2RAD*110 ;
else
    psi_0_min = psi_00  ; 
    psi_0_max = psi_00  ; 
end

% wing rotation phase offset (rad)
if varyFlag(6)
    del_psi_min = -1*DEG2RAD*110 ; 
    del_psi_max = -1*DEG2RAD*60 ;
else
    del_psi_min = del_psi0  ; 
    del_psi_max = del_psi0  ; 
end

% span length (meters)
if varyFlag(7)
    span_min = 1.0e-3 ; 
    span_max = 3e-3 ; 
else
    span_min = span0 ; 
    span_max = span0 ; 
end

 % length along long body axis from CM to hinge (m)
if varyFlag(8)
    hinge_x_min = 1e-4 ;
    hinge_x_max = 0.5e-3 ; 
else
    hinge_x_min = hinge_x0 ;
    hinge_x_max = hinge_x0 ; 
end

% but into lower/upper bounds vectors
lb = [omega_min, phi_f_min, phi_b_min, psi_m_min, psi_0_min, del_psi_min,...
    span_min, hinge_x_min] ; 
ub = [omega_max, phi_f_max, phi_b_max, psi_m_max, psi_0_max, del_psi_max,...
    span_max, hinge_x_max] ; 


% ------------------------------------------------------------------
%% use linear inequality to keep total stroke amp reasonably sized
% A = [0, 0, 0, 0, 0; 0, 1, -1, 0, 0 ; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; ...
%     0, 0, 0, 0, 0;] ;
% b = [0, -1*DEG2RAD*120, 0, 0, 0]' ;
% --------------------------------------------------
%% set other constraint types as empty
A = [] ;
b = [] ; 
Aeq = [] ; 
beq = [] ; 

% ------------------------------------------------
%% run solver 
options = optimoptions(@fmincon, 'Display', 'iter', ...
            'FiniteDifferenceType','central') ; 
            
[x, fval] = fmincon(@stabilityCostFunc,x0,A,b,Aeq,beq,lb,ub,[],options) ;

% ------------------------------------------------
%% convert phi_f, phi_b to phi_0, phi_m
phi_0 = (x(2) + x(3))./2 ; 
phi_m = (x(3) - x(2))./2 ; 

% ------------------------------------------------
%% print results
fprintf('f = %f \n', x(1)/(2*pi))
fprintf('phi_m = %f \n', RAD2DEG*phi_m)
fprintf('phi_0 = %f \n', RAD2DEG*phi_0)
fprintf('psi_m = %f \n', RAD2DEG*x(4))
fprintf('psi_0 = %f \n', RAD2DEG*x(5))
fprintf('del_psi = %f \n', x(6))
fprintf('span = %f \n', x(7))
fprintf('hinge_x = %f \n', x(8))

fprintf('\n total cost = %f \n', fval)