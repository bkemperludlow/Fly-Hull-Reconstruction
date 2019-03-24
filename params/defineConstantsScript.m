
%addpath .\fdaM

PSI    = 1 ; % body yaw
BETA   = 2 ; % body pitch
PHIR   = 3 ;
THETAR = 4 ;
ETAR   = 5 ;
PHIL   = 6 ;
THETAL = 7 ;
ETAL   = 8 ;
RHO    = 9 ; % body roll

PHIB   = PSI ; % body yaw
THETAB = BETA ; % body pitch
PSIB   = RHO ;  % body roll


rad2deg = 180 / pi ;
deg2rad = pi / 180 ;
newton2mg   = 1e5 ;

%pixels_to_meters = (1/232) * 1e-2 ; % 1cm=232pixels, 1pix=(1/232)cm=(1/23200)m