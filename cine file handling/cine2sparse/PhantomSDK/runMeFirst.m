% these commands should be exectuted in the beginning of the matlab session
% in case the compiler location variable is not set on the computer path, create it temporeraly:
% setenv('MW_MINGW64_LOC','C:\ProgramData\MATLAB\SupportPackages\R2017a\3P.instrset\mingw_492.instrset')
setenv('MW_MINGW64_LOC','C:\ProgramData\MATLAB\SupportPackages\R2020b\3P.instrset\mingw_w64.instrset')

% add SDK to path
folder_path = [pwd,'\SDK 13.4.788.0'];
addpath(pwd,genpath(folder_path))
rmpath(genpath([folder_path,'\bin\Win32']))
LoadPhantomLibraries