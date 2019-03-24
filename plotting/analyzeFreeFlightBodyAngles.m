function [h_pitch_hist, h_pitch_acf, h_roll_hist, h_roll_acf] ...
    = analyzeFreeFlightBodyAngles(flyLine) 

dataPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\Free Flight Analysis' ;

cd(dataPath)
load freeFlightBodyAnglesStruct_new

experimentalPitchColor = [.7 0 0] ; %[255 102 0]/255 ;
controlPitchColor = [0 .7 0] ; %[170,238,255]/255 ;
experimentalRollColor = [191 87 0]/255 ;
controlRollColor = [0 119 204]/255 ;
lineWidth = 1.5 ;

MB204BInd = find([freeFlightBodyAnglesStruct(:).ExprNum] == 1  | ...
    [freeFlightBodyAnglesStruct(:).ExprNum] == 2) ; % | ...
    %[freeFlightBodyAnglesStruct(:).ExprNum] == 8) ; 
MB258CInd = find([freeFlightBodyAnglesStruct(:).ExprNum] == 11  | ...
    [freeFlightBodyAnglesStruct(:).ExprNum] == 10 ) ; %| ...
    %[freeFlightBodyAnglesStruct(:).ExprNum] == 12 | ...
    %[freeFlightBodyAnglesStruct(:).ExprNum] == 13) ; 
Gal4_9281Ind = find([freeFlightBodyAnglesStruct(:).ExprNum] == 9 | ...
    [freeFlightBodyAnglesStruct(:).ExprNum] == 10 | ...
    [freeFlightBodyAnglesStruct(:).ExprNum] == 11) ; 
%EmptyVectorInd = find([freeFlightBodyAnglesStruct(:).ExprNum] == 14 | ...
%    [freeFlightBodyAnglesStruct(:).ExprNum] == 15 ) ; 

switch flyLine
    case 'MB204B'
        freeFlightBodyAnglesStruct = freeFlightBodyAnglesStruct(MB204BInd) ; 
        plotLabel = 'MB204B' ;
    case 'MB258C'
        freeFlightBodyAnglesStruct = freeFlightBodyAnglesStruct(MB258CInd) ;
        plotLabel = 'MB258C' ;
    case '9281(2)'
        freeFlightBodyAnglesStruct = freeFlightBodyAnglesStruct(Gal4_9281Ind) ;
        plotLabel = '9281(2)' ;
    otherwise
        disp('All fly lines')
        plotLabel = 'All' ;
end


controlInd = find([freeFlightBodyAnglesStruct(:).flyType] == 2) ;
experimentalInd = find([freeFlightBodyAnglesStruct(:).flyType] == 1) ;
%Nmovies = length(freeFlightBodyAnglesStruct) ;

nLags = 200 ;
xbins_pitch = 0:3:90 ;
xbins_roll = -45:4.5:45 ;

xlim_pitch = [min(xbins_pitch) max(xbins_pitch)] ;
xlim_roll = [min(xbins_roll) max(xbins_roll)] ;

%% get acfs for plots

for i = 1:length(freeFlightBodyAnglesStruct)
    bodyPitch = freeFlightBodyAnglesStruct(i).pitch_smoothed ;   
    bodyRoll = freeFlightBodyAnglesStruct(i).bodyRoll ; 
    
    [acf_pitch, lags_pitch] = autocorr(bodyPitch, min([nLags, length(bodyPitch)-1])) ;
    [acf_roll, lags_roll] = autocorr(bodyRoll, min([nLags, length(bodyRoll)-1])) ;
    
    freeFlightBodyAnglesStruct(i).acf_pitch = acf_pitch ; 
    freeFlightBodyAnglesStruct(i).lags_pitch = lags_pitch ; 
    freeFlightBodyAnglesStruct(i).acf_roll = acf_roll ; 
    freeFlightBodyAnglesStruct(i).lags_roll = lags_roll ; 

end
%% make some histograms
control_pitch_pdf = makePitchHist(freeFlightBodyAnglesStruct, controlInd, xbins_pitch) ;
experimental_pitch_pdf = makePitchHist(freeFlightBodyAnglesStruct, experimentalInd, xbins_pitch) ;
control_roll_pdf = makeRollHist(freeFlightBodyAnglesStruct, controlInd, xbins_roll) ;
experimental_roll_pdf = makeRollHist(freeFlightBodyAnglesStruct, experimentalInd, xbins_roll) ;

% pitch histograms
max_pitch_pdf = max([experimental_pitch_pdf, control_pitch_pdf]) ;
JSD_pitch = jensen_shannon_divergence(control_pitch_pdf , experimental_pitch_pdf) 

h_pitch_hist = figure('Position', [38 705 316 256],'PaperPositionMode', 'auto');
subplot(2,1,1)
bar(xbins_pitch, control_pitch_pdf, ...
    'faceColor',controlPitchColor) ;
set(gca,'xlim',xlim_pitch)
set(gca,'ylim',[0 max_pitch_pdf])
ylabel('PDF')

subplot(2,1,2)
bar(xbins_pitch, experimental_pitch_pdf, ...
    'faceColor',experimentalPitchColor) ; 
set(gca,'xlim',xlim_pitch)
set(gca,'ylim',[0 max_pitch_pdf])
ylabel('PDF')
xlabel('Body Pitch Angle, \theta_B , [deg]')

%------------------------------------------------------------------------------
% roll histograms
max_roll_pdf = max([experimental_roll_pdf, control_roll_pdf]) ;
JSD_roll = jensen_shannon_divergence(control_roll_pdf , experimental_roll_pdf) 

h_roll_hist = figure('Position',[40 371 316 256],'PaperPositionMode', 'auto');
subplot(2,1,1)
bar(xbins_roll, control_roll_pdf, ...
    'faceColor',controlRollColor) ; 
set(gca,'xlim',xlim_roll)
set(gca,'ylim',[0 max_roll_pdf])
ylabel('PDF')


subplot(2,1,2)
bar(xbins_roll, experimental_roll_pdf, ...
    'faceColor',experimentalRollColor) ; 
set(gca,'xlim',xlim_roll)
set(gca,'ylim',[0 max_roll_pdf])
ylabel('PDF')
xlabel('Body Roll Angle, \rho_B , [deg]')

%% make acf plots

h_pitch_acf = figure('Position', [368 754 285 254],'PaperPositionMode', 'auto');
hold on

for j = controlInd
    plot([freeFlightBodyAnglesStruct(j).lags_pitch]/8, ...
        [freeFlightBodyAnglesStruct(j).acf_pitch],'-','LineWidth',lineWidth,...
        'Color',controlPitchColor) 
    %keyboard ;
end
for k = experimentalInd
    plot([freeFlightBodyAnglesStruct(k).lags_pitch]/8, ...
        [freeFlightBodyAnglesStruct(k).acf_pitch],'-','LineWidth',lineWidth,...
        'Color',experimentalPitchColor) 
    %keyboard ;
end
axis tight
xlabel('Lags [ms]')
ylabel('Pitch Autocorrelation')

%------------------------------------------------------------------------------
h_roll_acf = figure('Position', [388 436 285 254],'PaperPositionMode', 'auto');
hold on

for j = controlInd
    plot([freeFlightBodyAnglesStruct(j).lags_roll]/8, ...
        [freeFlightBodyAnglesStruct(j).acf_roll],'-','LineWidth',lineWidth,...
        'Color',controlRollColor) 
    %keyboard ;
end
for k = experimentalInd
    plot([freeFlightBodyAnglesStruct(k).lags_roll]/8, ...
        [freeFlightBodyAnglesStruct(k).acf_roll],'-','LineWidth',lineWidth,...
        'Color',experimentalRollColor) 
    %keyboard ;
end
axis tight
xlabel('Lags [ms]')
ylabel('Roll Autocorrelation')

%{
cd('G:\b1 paper\figures (rough)\free flight plots')
print(h_pitch_hist,'pitch_hist.svg','-dsvg')
print(h_roll_hist,'roll_hist.svg','-dsvg')
print(h_pitch_acf,'pitch_ACF.svg','-dsvg')
print(h_roll_acf,'roll_ACF.svg','-dsvg')
%}
end

%---------------------------------------------------------

%% functions to make plots
%---------------------------------------------------------
function pdf = makePitchHist(struct, ind, xbins)
pitchList = [] ;

for k = ind
    pitch_temp = struct(k).pitch_smoothed ;
    pitchList = [pitchList; pitch_temp] ;
end    

[f_pitch, x_pitch] = hist(pitchList, xbins) ; 
pdf = f_pitch ./ trapz(x_pitch, f_pitch) ; 

end
%---------------------------------------------------------

function pdf = makeRollHist(struct, ind, xbins)
rollList = [] ;

for k = ind
    roll_temp = struct(k).bodyRoll ;
    rollList = [rollList; roll_temp] ;
end    

[f_roll, x_roll] = hist(rollList, xbins) ; 
pdf = f_roll ./ trapz(x_roll, f_roll) ; 
end
%---------------------------------------------------------

    