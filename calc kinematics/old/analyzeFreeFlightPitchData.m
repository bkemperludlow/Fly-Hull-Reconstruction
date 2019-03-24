dataPath = 'G:\Janelia Flies\kir2.1 flies\Analysis\No Perturbation\Compiled Data' ;

cd(dataPath)
load freeFlightPitchStruct

kirColor = [.7 0 0] ; %[255 102 0]/255 ;
whiteEyesColor = [0 .7 0] ; %[170,238,255]/255 ;

kirRollColor = [80 33 0]/255 ;
whiteEyesRollColor = [0 119 204]/255 ;

for j = 3:7
   idx = find([freeFlightPitchStruct(:).ExprNum] == j) ;
   
   figure('Name',['Experiment Number ' num2str(j)])
   hold on
   for m = 1:length(idx)
       %fucking delta
       t = freeFlightPitchStruct(idx(m)).time ;
       pitchAngle = freeFlightPitchStruct(idx(m)).bodyPitchAngle ;
       plot(1000*t, pitchAngle)
      
   end
    axis tight
    grid on
end

%Bad movies: 17, 31

%% make some histograms

pitchCell = cell(1,5) ;
timeCell = cell(1,5) ;
for j = 3:7
   idx = find([freeFlightPitchStruct(:).ExprNum] == j) ;
   pitchArray = [] ;
   timeArray = [] ;
   for m = 1:length(idx)      
       pitchArray = [pitchArray; freeFlightPitchStruct(idx(m)).bodyPitchAngle] ;
       timeArray = [timeArray , freeFlightPitchStruct(idx(m)).time ] ;     
   end
   pitchCell{j-2} = pitchArray ; 
   timeCell{j-2} = timeArray ; 
end

nBins = 50 ;
[f_all, x_all] = hist([pitchCell{1}; pitchCell{2}; pitchCell{4}; pitchCell{5}],nBins) ;
[f_3, ~] = hist(pitchCell{1},x_all) ;
[f_4, ~] = hist(pitchCell{2},x_all) ;
[f_6, ~] = hist(pitchCell{4},x_all) ;
[f_7, ~] = hist(pitchCell{5},x_all) ;


hhist = figure('Position',[140 100 500 700],'PaperPositionMode','auto')   ; 
subplot(4,1,1)
bar(x_all, f_3/trapz(x_all,f_3), 'faceColor',kirColor) 
set(gca,'xlim', [-5 70] )
set(gca,'ylim', [0 0.085] )
set(gca,'fontsize',12)
legend({'204b > kir2.1'},'location','northwest') %, '204b > w1118', '258c > kir2.1', '258c > w1118'})

subplot(4,1,2)
bar(x_all, f_4/trapz(x_all,f_4), 'faceColor',whiteEyesColor) 
set(gca,'xlim', [-5 70] )
set(gca,'ylim', [0 0.085] )
set(gca,'fontsize',12)
legend({'204b > w1118'},'location','northwest')

subplot(4,1,3)
bar(x_all, f_6/trapz(x_all,f_6), 'faceColor',kirColor) 
set(gca,'xlim', [-5 70] )
set(gca,'ylim', [0 0.085] )
set(gca,'fontsize',12)
legend({'258c > kir2.1'},'location','northwest')

subplot(4,1,4)
bar(x_all, f_7/trapz(x_all,f_7), 'faceColor',whiteEyesColor) 
set(gca,'xlim', [-5 70] )
set(gca,'ylim', [0 0.085] )
set(gca,'fontsize',12)
legend({'258c > w1118'},'location','northwest')


%% Look at dat autocorrelation
nLags = 200 ;
h_autocorr = figure ;
hold on
for j = 3:7
   idx = find([freeFlightPitchStruct(:).ExprNum] == j) ;
   if (j == 3) || (j ==6)
       plotColor = kirColor ;
   else
       plotColor = whiteEyesColor ;
   end
   
   for m = 1:length(idx)
       %fucking delta
       
       t = 1000*freeFlightPitchStruct(idx(m)).time ;
       pitchAngle = freeFlightPitchStruct(idx(m)).bodyPitchAngle ;
       [acf, lags] = autocorr(pitchAngle, min([nLags, length(pitchAngle)-1])) ;
       plot(lags, acf, 'Color',plotColor,'LineWidth',2)
       disp(m) ;
   end
    axis tight
    %grid on
end

%}