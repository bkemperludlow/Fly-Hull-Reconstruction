function pValue = myFTest(P_controller_struct, PI_controller_struct) 

N_P = length(P_controller_struct.fwdFlipTimes) ; 
N_PI = length(PI_controller_struct.fwdFlipTimes) ; 

df_P = N_P - 3 ; %P controller has 3 parameters 
df_PI = N_PI - 4;  %PI controller has 4 parameters

rmse_P = P_controller_struct.rms ;
RSS_P = N_P*(rmse_P)^2 ;
rmse_PI = PI_controller_struct.rms ;
RSS_PI = N_PI*(rmse_PI)^2 ;

% from http://tinyurl.com/h63kc4d
Fstatistic = ((RSS_P - RSS_PI)/(df_P - df_PI))/ (RSS_PI/df_PI) ;

pValue = fcdf(Fstatistic,df_P,df_PI) ;

end