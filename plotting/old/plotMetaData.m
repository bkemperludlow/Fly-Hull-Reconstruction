%Plots metadata from tables

testDist = 2e5 ; 
textFlag = true ; 
haccel = figure;

S = readtable('pitchDown10302014.csv');
T = readtable('pitchUp10302014.csv');
%--------------------------------------------------------------------------------
subplot(4,2,1)
plot(T.PhiFront, T.pitchAcceleration,'b.','MarkerSize',20)
hold on
plot(S.PhiFront, S.pitchAcceleration,'r.','MarkerSize',20)
title('Pitch Acceleration vs Phi Front')
grid on;
if textFlag
    for i = 1:length(T.PhiFront)
        text(T.PhiFront(i),T.pitchAcceleration(i)+textDist,[num2str(T.Expr(i)) ...
            '.' num2str(T.Mov(i))], 'horizontalalignment','center')
    end
    for j = 1:length(S.PhiFront)
        text(S.PhiFront(j),S.pitchAcceleration(j)+textDist,[num2str(S.Expr(j)) ...
            '.' num2str(S.Mov(j))], 'horizontalalignment','center')
    end
end
%--------------------------------------------------------------------------------
subplot(4,2,2)
plot(T.RawPhiFront, T.pitchAcceleration,'b.','MarkerSize',20)
hold on
plot(S.RawPhiFront, S.pitchAcceleration,'r.','MarkerSize',20)
title('Pitch Acceleration vs Raw Phi Front')
grid on;
if textFlag
    for i = 1:length(T.RawPhiFront)
        text(T.RawPhiFront(i),T.pitchAcceleration(i)+textDist,[num2str(T.Expr(i)) ...
            '.' num2str(T.Mov(i))], 'horizontalalignment','center')
    end
    for j = 1:length(S.RawPhiFront)
        text(S.RawPhiFront(j),S.pitchAcceleration(j)+textDist,[num2str(S.Expr(j)) ...
            '.' num2str(S.Mov(j))], 'horizontalalignment','center')
    end
end
%--------------------------------------------------------------------------------
subplot(4,2,3)
plot(T.PhiBack, T.pitchAcceleration,'b.','MarkerSize',20)
hold on
plot(S.PhiBack, S.pitchAcceleration,'r.','MarkerSize',20)
title('Pitch Acceleration vs Phi Back')
grid on;
if textFlag
    for i = 1:length(T.PhiBack)
        text(T.PhiBack(i),T.pitchAcceleration(i)+textDist,[num2str(T.Expr(i)) ...
            '.' num2str(T.Mov(i))], 'horizontalalignment','center')
    end
    for j = 1:length(S.PhiBack)
        text(S.PhiBack(j),S.pitchAcceleration(j)+textDist,[num2str(S.Expr(j)) ...
            '.' num2str(S.Mov(j))], 'horizontalalignment','center')
    end
end
%--------------------------------------------------------------------------------
subplot(4,2,4)
plot(T.RawPhiBack, T.pitchAcceleration,'b.','MarkerSize',20)
hold on
plot(S.RawPhiBack, S.pitchAcceleration,'r.','MarkerSize',20)
title('Pitch Acceleration vs Raw Phi Back')
grid on;
if textFlag
    for i = 1:length(T.RawPhiBack)
        text(T.RawPhiBack(i),T.pitchAcceleration(i)+textDist,[num2str(T.Expr(i)) ...
            '.' num2str(T.Mov(i))], 'horizontalalignment','center')
    end
    for j = 1:length(S.RawPhiBack)
        text(S.RawPhiBack(j),S.pitchAcceleration(j)+textDist,[num2str(S.Expr(j)) ...
            '.' num2str(S.Mov(j))], 'horizontalalignment','center')
    end
end
%--------------------------------------------------------------------------------
subplot(4,2,5)
plot(T.PhiMid, T.pitchAcceleration,'b.','MarkerSize',20)
hold on
plot(S.PhiMid, S.pitchAcceleration,'r.','MarkerSize',20)
title('Pitch Acceleration vs Phi Mid')
grid on;
if textFlag
    for i = 1:length(T.PhiMid)
        text(T.PhiMid(i),T.pitchAcceleration(i)+textDist,[num2str(T.Expr(i)) ...
            '.' num2str(T.Mov(i))], 'horizontalalignment','center')
    end
    for j = 1:length(S.PhiMid)
        text(S.PhiMid(j),S.pitchAcceleration(j)+textDist,[num2str(S.Expr(j)) ...
            '.' num2str(S.Mov(j))], 'horizontalalignment','center')
    end
end
%--------------------------------------------------------------------------------
subplot(4,2,6)
plot(T.RawPhiMid, T.pitchAcceleration,'b.','MarkerSize',20)
hold on
plot(S.RawPhiMid, S.pitchAcceleration,'r.','MarkerSize',20)
title('Pitch Acceleration vs Raw Phi Mid')
grid on;
if textFlag
    for i = 1:length(T.RawPhiMid)
        text(T.RawPhiMid(i),T.pitchAcceleration(i)+textDist,[num2str(T.Expr(i)) ...
            '.' num2str(T.Mov(i))], 'horizontalalignment','center')
    end
    for j = 1:length(S.RawPhiMid)
        text(S.RawPhiMid(j),S.pitchAcceleration(j)+textDist,[num2str(S.Expr(j)) ...
            '.' num2str(S.Mov(j))], 'horizontalalignment','center')
    end
end
%--------------------------------------------------------------------------------
subplot(4,2,7)
plot(T.PhiAmp, T.pitchAcceleration,'b.','MarkerSize',20)
hold on
plot(S.PhiAmp, S.pitchAcceleration,'r.','MarkerSize',20)
title('Pitch Acceleration vs Phi Amp')
grid on;
if textFlag
    for i = 1:length(T.PhiAmp)
        text(T.PhiAmp(i),T.pitchAcceleration(i)+textDist,[num2str(T.Expr(i)) ...
            '.' num2str(T.Mov(i))], 'horizontalalignment','center')
    end
    for j = 1:length(S.PhiAmp)
        text(S.PhiAmp(j),S.pitchAcceleration(j)+textDist,[num2str(S.Expr(j)) ...
            '.' num2str(S.Mov(j))], 'horizontalalignment','center')
    end
end
%--------------------------------------------------------------------------------
subplot(4,2,8)
plot(T.RawPhiAmp, T.pitchAcceleration,'b.','MarkerSize',20)
hold on
plot(S.RawPhiAmp, S.pitchAcceleration,'r.','MarkerSize',20)
title('Pitch Acceleration vs RawPhiAmp')
grid on;
if textFlag
    for i = 1:length(T.RawPhiAmp)
        text(T.RawPhiAmp(i),T.pitchAcceleration(i)+textDist,[num2str(T.Expr(i)) ...
            '.' num2str(T.Mov(i))], 'horizontalalignment','center')
    end
    for j = 1:length(S.RawPhiAmp)
        text(S.RawPhiAmp(j),S.pitchAcceleration(j)+textDist,[num2str(S.Expr(j)) ...
            '.' num2str(S.Mov(j))], 'horizontalalignment','center')
    end
end