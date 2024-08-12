clear
close all
clc


path='G:\Google Drive huji\Cine_sparse_files\Exp_29_11_2020\'
listing = dir(path)
for k=1:1:length(listing)
    movfile=strfind(listing(k).name,'mov');
    file_end=strfind(listing(k).name,'.m');   
    if isempty(movfile)==0 && isempty(file_end)==0
        name_splt=split(listing(k).name,'_');
        mkdir([path,name_splt{1}])
        moveFrom=[path,listing(k).name];
        moveto=[path,name_splt{1}]
        movefile(moveFrom,moveto)
    end
end
    