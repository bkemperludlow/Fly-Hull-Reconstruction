% -------------------------------------------------------------------------
% will be back to improve on this 
% -------------------------------------------------------------------------
%% path info
dataPath = 'D:\Fly Data\Simulation\Expr_96_mov_009' ;
dataFn = 'Expr_96_mov_009_simFit.mat' ; 

saveFlag = true ; 
savePath = fullfile('D:\Dropbox\Cohen Group Meeting 08-28-2020', ...
    'fitSim2DataExample.png') ; 

sim_data_struct = importdata(fullfile(dataPath, dataFn)) ; 

b1_paper_plot_pref 


% ----------------------------------------------
%% read data struct
tms = 1000.*sim_data_struct.t ;  % convert to ms

thetaB_sim = (180/pi).*sim_data_struct.thetaB_sim ; % convert to deg
thetaB_sim_raw = (180/pi).*sim_data_struct.thetaB_sim_raw ;
thetaB_sim_raw = thetaB_sim_raw - thetaB_sim_raw(1) ; 
thetaB_data = (180/pi).*sim_data_struct.thetaB_data ; 

pulseStart = 1000.*sim_data_struct.pulseStart ; % convert to ms
pulseEnd = 1000.*sim_data_struct.pulseEnd ;

xlim = [min(tms), max(tms)] ; 
ymin = min([thetaB_sim' ; thetaB_sim_raw' ; thetaB_data]) - 5; 
ymax = max([thetaB_sim' ; thetaB_sim_raw' ; thetaB_data]) + 5; 
ylim = [ymin, ymax] ; 

raw_ind = 1200 + (1:length(thetaB_sim)) ; % FIX THIS!
% --------------------------------------------
%% make figure
h_main = figure('OuterPosition', [691   474   426   370]) ;  
hold on

% make pert bar
tsfvec = [pulseStart pulseEnd pulseEnd pulseStart pulseStart] ;
avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
hf = fill(tsfvec , avec,'y') ;
set(hf,'facecolor',[255 238 170]/255,'facealpha',1,'edgecolor','none') ;
set(hf,'HandleVisibility','off')

% plot raw sim pitch data
plot(tms, thetaB_sim_raw(raw_ind), '-', 'Color', 0.6*[1, 1, 1],...
    'LineWidth', 1.5)

% plot smoothed sim data
plot(tms, thetaB_sim, 'k-', 'LineWidth', 1.5)

% plot actual data
plot(tms, thetaB_data, 'b-', 'LineWidth', 1.5)

% ----------------------------------
%% axis properties
set(gca, 'xlim', xlim, 'ylim', ylim)
ylabel('Body Pitch Angle, \theta_b , (deg)')
xlabel('Time (ms)') 

legend({'raw sim', 'filt sim', 'data'},'location','southeast')


% -------------------------------------
%% save?
if saveFlag
   print(h_main, savePath, '-dpng', '-r300') 
end