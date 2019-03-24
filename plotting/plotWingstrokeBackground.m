function plotWingstrokeBackground(ax, backFlipTimes, fwdFlipTimes, mycolor, transparencyFlag)
% plots gray backgrouds behind the BACK strokes
% ax - axis to plot on
% miny, maxy - vertical limits of rectrangles to plot (since this function
% will be called before the data is plotted)
% backFlipTimes, fwdFlipTimes - self explanatory

if (~exist('transparencyFlag'))
    transparencyFlag = false ;
end
    
axis(ax);
hold on ;
%mycolor = [1 1 1 ] * 0.8 ;
faceAlpha = 0.5 ;
ylim = get(ax, 'ylim') ;
a1 = ylim(1) ; 
a2 = ylim(2) ;
clear ylim

try
    dt = fwdFlipTimes(2) - fwdFlipTimes(1) ;
catch
    dt = backFlipTimes(2) - backFlipTimes(1) ;
end

% check whether the first flip is back or fwd flip
if (fwdFlipTimes(1) < backFlipTimes(1))
    fwdind = 1 ;
    backind = 1 ;
else    
    % add a fictitions fwd flip
    fwdFlipTimes = [ (fwdFlipTimes(1)-2*dt) fwdFlipTimes] ;
    fwdind = 1 ;
    backind = 1 ;
end

if (fwdFlipTimes(end) > backFlipTimes(end))
    % add a fictitions backflip time
    backFlipTimes = [backFlipTimes (backFlipTimes(end)+2*dt)] ;
    fwdind = 1 ;
    backind = 1 ;
end

Nf = length(fwdFlipTimes) ;
Nb = length(backFlipTimes) ;
% color the BACK stroke in gray backgroud


while(fwdind<=Nf && backind<=Nb)
    
    ts = fwdFlipTimes(fwdind) ;
    tf = backFlipTimes(backind) ;
    hf = fill([ts tf tf ts ts ] , [ a1 a1 a2 a2 a1],mycolor) ;
    set(hf,'linestyle','none') ;
    uistack(hf,'bottom') ; %added by SW
    if (transparencyFlag)
        set(hf,'facealpha',faceAlpha);
    end
    
    fwdind = fwdind + 1 ;
    backind = backind + 1 ;
    
end

end
