
%wandPoints_1
M1 = getWandPoints(17,45,'xy_1.cin', 'xz_1.cin', 'yz_1.cin','wandPoints1.csv');
M2 = getWandPoints(17,45,'xy_2.cin', 'xz_2.cin', 'yz_2.cin','wandPoints2.csv');
M3 = getWandPoints(17,45,'xy_3.cin', 'xz_3.cin', 'yz_3.cin','wandPoints3.csv');
%M4 = getWandPoints(17,45,'xy_4.cin', 'xz_4.cin', 'yz_4.cin','wandPoints4.csv');
%M5 = getWandPoints(17,45,'xy_5.cin', 'xz_5.cin', 'yz_5.cin','wandPoints5.csv');
%M6 = getWandPoints(17,45,'xy_6.cin', 'xz_6.cin', 'yz_6.cin','wandPoints6.csv');
%M7 = getWandPoints(17,45,'xy_7.cin', 'xz_7.cin', 'yz_7.cin','wandPoints7.csv');
%M8 = getWandPoints(17,45,'xy_8.cin', 'xz_8.cin', 'yz_8.cin','wandPoints8.csv');

M = [ M1; M2; M3] ; % ;M5];  %; M3 ; M4] ;  % ; M4 ;M5]; %; M3];

dlmwrite('wandPoints.csv', M) ;
easyWand5

