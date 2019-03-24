function s = RotatePoint(r,p1,p2,theta)
% rotates the vector 'r' by 'theta' degrees about the axis defined by p1 and p2
% based on:
% http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html

theta = theta * pi / 180 ; % convert to radians

x = r(1) ;
y = r(2) ;
z = r(3) ;

u = p2(1) - p1(1) ;
v = p2(2) - p1(2) ;
w = p2(3) - p1(3) ;

nrm = sqrt(u^2 + v^2 + w^2) ;
u = u / nrm ;
v = v / nrm ;
w = w / nrm ;

a = 0 ; 
b = 0 ;
c = 0 ;

s = zeros(size(r)) ;

s(1) = ( a*(v^2+w^2) - u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(theta)) + ...
    x*cos(theta) + (-c*v+b*w-w*y+v*z)*sin(theta) ;

s(2) = ( b*(u^2+w^2) - v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(theta)) + ...
    y*cos(theta) + (c*u-a*w+w*x-u*z)*sin(theta) ;

s(3) = ( c*(u^2+v^2) - w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(theta)) + ...
    z*cos(theta) + (-b*u+a*v-v*x+u*y)*sin(theta) ;


end