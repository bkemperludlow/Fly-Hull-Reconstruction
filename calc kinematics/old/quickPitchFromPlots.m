h = findobj(gca,'Type','line') ;
xd = get(h,'XData');
yd = get(h,'YData');

t = xd{1,1} ;
y = yd{1,1} ; 

zeropoint = find(t == 0) ;
baseyd = y(zeropoint) ;
%[maxyd, maxind] = max(y) ;
[minyd, minind] = min(y(1,1:800)) ;

%i1 = find(y <= .1*(maxyd-baseyd)+baseyd) ; 
i1 = find(y >= .1*(minyd-baseyd)+baseyd) ; 
%i2 = find(i1 > maxind ) ; 
i2 = find(i1 > minind ) ; 
i3 = min(i1(i2)) ;


t_c = t(i3) 
%delta_pitch = maxyd - baseyd
delta_pitch = minyd - baseyd
close all