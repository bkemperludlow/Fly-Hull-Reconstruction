% file type
ext = '.cin' ; 

%wandPoints_1
M1 = getWandPoints(17,45,['xy_1' ext], ['xz_1' ext], ['yz_1' ext], ...
    'wandPoints1.csv');
M2 = getWandPoints(17,45,['xy_2' ext], ['xz_2' ext], ['yz_2' ext],...
    'wandPoints2.csv');
M3 = getWandPoints(17,45,['xy_3' ext], ['xz_3' ext], ['yz_3' ext],...
    'wandPoints3.csv');
%M4 = getWandPoints(17,45,'xy_4.cin', 'xz_4.cin', 'yz_4.cin','wandPoints4.csv');
%M5 = getWandPoints(17,45,'xy_5.cin', 'xz_5.cin', 'yz_5.cin','wandPoints5.csv');
%M6 = getWandPoints(17,45,'xy_6.cin', 'xz_6.cin', 'yz_6.cin','wandPoints6.csv');
%M7 = getWandPoints(17,45,'xy_7.cin', 'xz_7.cin', 'yz_7.cin','wandPoints7.csv');
%M8 = getWandPoints(17,45,'xy_8.cin', 'xz_8.cin', 'yz_8.cin','wandPoints8.csv');

% close parallel pool
delete(gcp) 

M = [ M1; M2; M3] ; % ;M5];  %; M3 ; M4] ;  % ; M4 ;M5]; %; M3];

dlmwrite('wandPoints.csv', M) ;
easyWand5

