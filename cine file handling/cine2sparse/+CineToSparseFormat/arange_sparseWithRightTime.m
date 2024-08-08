clear
close all
clc


path='G:\Google Drive huji\Cine_sparse_files\Exp_23_9_2020\'
% listing = dir(path)
% for k=1:1:length(listing)
%     movfile=strfind(listing(k).name,'mov');
%     %        movNum=str2double(listing(k).name(movfile+3:end));
%     movNum = str2double(regexp(listing(k).name, '\d+', 'match'));
%     
%     if isempty(movfile)==0
%         dirName=listing(k).name;
%         listing_file = dir([path,'/',dirName,'/']);
%         
%         for kcam=1:1:length(listing_file)
%             
%             camfile=strfind(listing_file(kcam).name,'cam');
%             if isempty(camfile)==0
%                 load([path,'/',dirName,'/',listing_file(kcam).name]);
%                 numCam=str2num(listing_file(kcam).name(camfile+3));
%                 time{movNum,1}=movNum;
%                 time{movNum,numCam}=metaData.xmlStruct.CineFileHeader.TriggerTime.Time.Text;
%                 
%                 
%             end
%             
%             
%         end
%     end
% end
% save([path,'\','time'],'time')

%%
load([path,'\','time'],'time')

path2='G:\Google Drive huji\Cine_sparse_files\Exp_23_9_2020_cam_order\Exp_23_9_2020_Reorder\'
mkdir(path2);
for k=1:1:size(time,1)
    movdir=sprintf('mov%d\\',time{k,1})
    mkdir([path2,movdir]);
    listing = dir([path,movdir])
%     for klist=1:1:length(listing)
%         NNop_dir=strfind(listing(klist).name,'NNoutput');
%         if NNop_dir==1 
%            copyfile([path,movdir,'\NNoutput'],[path2,movdir,'\NNoutput'],'f')
%     end
%         Hull_op_dir=strfind(listing(klist).name,'Hull_op');
%        if Hull_op_dir==1 
%            copyfile([path,movdir,'\Hull_op'],[path2,movdir,'\Hull_op'],'f')
%     end 
%     end
    
    
    movname=sprintf('mov%d_cam2_sparse.mat',time{k,1})
    copyfile([path,movdir,movname],[path2,movdir],'f')
    
    cam2Time = str2double(regexp(time{k,2}, '\d+', 'match'));
    cam3Time = str2double(regexp(time{k,3}, '\d+', 'match'));
    cam4Time = str2double(regexp(time{k,4}, '\d+', 'match'));
    
    comp2_to3=sum(cam2Time(1:4)==cam3Time(1:4))==4;
    comp2_to4=sum(cam2Time(1:4)==cam4Time(1:4))==4;
    
    if comp2_to3==1
        movname_cam3=sprintf('mov%d_cam3_sparse.mat',time{k,1})
        copyfile([path,movdir,movname_cam3],[path2,movdir],'f')
        
    else
        for mov_ind=1:1:size(time,1)
            if isempty(time{mov_ind,3})==0
            cam3Time = str2double(regexp(time{mov_ind,3}, '\d+', 'match'));
            comp2_to3=sum(cam2Time(1:4)==cam3Time(1:4))==4;
            if comp2_to3==1
                movname_cam3=sprintf('mov%d_cam3_sparse.mat',time{mov_ind,1});
                movdir_tmp=sprintf('mov%d\\',time{mov_ind,1})
                movname_cam3_tmp=sprintf('mov%d_cam3_sparse.mat',time{k,1});
              
                copyfile([path,movdir_tmp,movname_cam3],[path2,movdir,movname_cam3_tmp],'f');
            end
            end
            
        end
    end
    
    
    if comp2_to4==1
        movname_cam4=sprintf('mov%d_cam4_sparse.mat',time{k,1});
        copyfile([path,movdir,movname_cam4],[path2,movdir],'f');
        
    else
        for mov_ind=1:1:size(time,1)
                        if isempty(time{mov_ind,4})==0

            cam4Time = str2double(regexp(time{mov_ind,4}, '\d+', 'match'));
            comp2_to4=sum(cam2Time(1:4)==cam4Time(1:4))==4;
            if comp2_to4==1
                movname_cam4=sprintf('mov%d_cam4_sparse.mat',time{mov_ind,1});
                movdir_tmp=sprintf('mov%d\\',time{mov_ind,1})
                movname_cam4_tmp=sprintf('mov%d_cam4_sparse.mat',time{k,1});
              
                copyfile([path,movdir_tmp,movname_cam4],[path2,movdir,movname_cam4_tmp],'f');
            end
                        end
            
            
        end
    end
    
    
    
    
end


%     moveFrom=[path,listing(k).name];
%         moveto=[path,name_splt{1}]
%         movefile(moveFrom,moveto)