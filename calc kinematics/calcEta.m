
function eta = calcEta(s, c, leftRight, plotFlag)
% calculates the wing pitch angle eta in degrees
% input parameters:
% s - span unit vector
% c - chord unit vector
% leftRight - a string that indicates whether it's the right or left wing.
% only the lower cae of the first character is considered ('l' or 'r')

if (size(s,1)==1)
    s = s' ;
end
if (size(c,1)==1)
    c = c' ;
end

if (nargin==3)
    plotFlag = false ;
end
delta    = 1000*eps ;
% -----------------------
% determine left or right
% -----------------------

chr = lower(leftRight(1)) ;

switch(chr)
    case 'r'
        wingFlag = 1 ;
    case 'l'
        wingFlag = -1 ;
    otherwise
        disp('illegal value for the leftRight input argument') ;
        eta = NaN ;
        return
end


% -----------------------------
% calculate eta - using vectors
% -----------------------------

phihat = - cross(s, [0 ; 0 ; 1]) * wingFlag ;
phihat = phihat / norm(phihat) ;

eta = acos( dot(c, phihat) ) * 180 / pi ;

% take care of eta's sign
if (eta~=0)
    v3 = cross(phihat, c) ; % v3 should be parallel to s
    v3 = v3 / norm(v3) ;
    % if v3 and s are on the same direction, then sgn=+1
    % if v3 and s are on opposite direction then sgn=-1
    
    sgn = dot(v3, s) ; % sgn is either 1 or -1
    sgn2 = round(sgn) ; % to prevent cases where sgn = +-1 +1 epsilon
    if (abs(sgn-sgn2)>delta)
        disp('ERR in the calculation of eta. Check this.') ;
        %keyboard ;
    end
    
    eta = eta * sgn2 * wingFlag ;
end
% print the result
% disp(['method 1: eta=' num2str(eta) ]) ;

if (eta<0) 
    eta = eta + 360 ;
end

% ----
% plot
% ----
if (plotFlag)
    v1 = - cross(s,[0 ;0 ;1]) * wingFlag ; % projection of c onto the xy plane
    v1 = v1 / norm(v1) ;                   % normalize
    v2 = - cross(v1, s) * wingFlag ;
    origin = [ 0 ; 0 ; 0 ] ;
    
    figure ;
    hold on ;
    myplot(origin, origin,'ks') ;
    myplot(origin, s, 'ro-') ;
    %{
    myplot(s, c, 'bo-') ;
    myplot(s, v1,'ko--') ;
    myplot(s, v2,'ks--') ;
    myplot(s, phihat,'gd-') ;
    %}
    myplot(origin, c, 'bo-') ;
    myplot(origin, v1,'ko--') ;
    myplot(origin, v2,'ks--') ;
    myplot(origin, phihat,'gd-') ;
    hold off ;
    axis equal ;
    grid on ;
    box on ;
    view(3) ;
    xlabel('x') ; ylabel('y') ; zlabel('z') ;
end
end

function myplot(v1,v2,colstr)
v3 = v1+v2 ;
plot3( [v1(1) v3(1)], [v1(2) v3(2)], [v1(3) v3(3)],colstr,'linewidth',2) ;
end