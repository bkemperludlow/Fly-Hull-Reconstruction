function [tip]=find_pin_mk4(pincoords_yz,pincoords_xz,pincoords_xy,easyWandData,dlt)
%Order for image tipe should be yz xz xy
[commonTime, idx_xz,idx_xy] = intersect(pincoords_xz(:,3),pincoords_xy(:,3));
pincoords_yz= nan(size(commonTime,1),2);

imagetip = [pincoords_yz,pincoords_xz(idx_xz,1:2),pincoords_xy(idx_xy,1:2),commonTime];
imagetip(:,2:2:6)=repmat(513,size(imagetip,1),3)-imagetip(:,2:2:6);

Nframes = size(imagetip,1);
tip=zeros(Nframes,3);
[Kxz,Kyz,Kxy,Rxz,Ryz,Rxy] = find_transform2(easyWandData, dlt,[2 1 3]);
voxelSize=0.00005;


for i=1:Nframes
    
    imageTip_yz=imagetip(i,1:2);
    imageTip_xz=imagetip(i,3:4);
    imageTip_xy=imagetip(i,5:6);
    
    if (sum(~isnan(imagetip(i,1:6)))== 0 )
        
        G1= reconstructpoint(Ryz,Kyz,imageTip_yz,Rxy,Kxy,imageTip_xy);
        G2= reconstructpoint(Ryz,Kyz,imageTip_yz,Rxz,Kxz,imageTip_xz);
        G3= reconstructpoint(Rxy,Kxy,imageTip_xy,Rxz,Kxz,imageTip_xz);
        G = [G1';G2';G3'];
        
        tip(i,1:3) = mean(G);
        tip(i,1:3)=round(tip(i,:)/voxelSize);
    else
        if ~isnan(imageTip_yz(1)) && ~isnan(imageTip_xz(1))
            
            
            tip(i,1:3)= reconstructpoint(Ryz,Kyz,imageTip_yz,Rxz,Kxz,imageTip_xz);
            tip(i,1:3)=round(tip(i,:)/voxelSize);
            
        elseif ~isnan(imageTip_yz(1)) && ~isnan(imageTip_xy(1))
            
            tip(i,1:3)= reconstructpoint(Ryz,Kyz,imageTip_yz,Rxy,Kxy,imageTip_xy);
            tip(i,1:3)=round(tip(i,:)/voxelSize);
            
        elseif ~isnan(imageTip_xz(1)) && ~isnan(imageTip_xy(1))
            
            tip(i,1:3)= reconstructpoint(Rxy,Kxy,imageTip_xy,Rxz,Kxz,imageTip_xz);
            tip(i,1:3)=round(tip(i,:)/voxelSize);
            
        else
            
            tip(i,1:3)=[NaN NaN NaN];
        end
        
    end
end

tip = [tip,imagetip(:,7)];

end


function [Kxz,Kyz,Kxy,Rxz,Ryz,Rxy] = find_transform2(easyWandData, dlt,order)

%load(calibration_directory)
%dlt=load(coefs_directory);
%order is a 1x3 vector  such that order(1)= index of camera YZ for calibration
%                                 order(2)= index of camera XZ for calibration
%                                 order(3)= index of camera XY for calibration

dltxy=dlt(:,order(3));
dltxz=dlt(:,order(2));
dltyz=dlt(:,order(1));

pp=easyWandData.ppts;
f=easyWandData.focalLengths';

f=f(:,order);
pp=pp(:,[2*order(1)-1,2*order(1),2*order(2)-1,2*order(2),2*order(3)-1,2*order(3)]);

Kyz=[f(1),0,pp(1);0,f(1),pp(2);0,0,1];
Kxz=[f(2),0,pp(3);0,f(2),pp(4);0,0,1];
Kxy=[f(3),0,pp(5);0,f(3),pp(6);0,0,1];

Axz=[dltxz(1),dltxz(2),dltxz(3),dltxz(4);dltxz(5),dltxz(6),dltxz(7),dltxz(8);dltxz(9),dltxz(10),dltxz(11),1];
Ayz=[dltyz(1),dltyz(2),dltyz(3),dltyz(4);dltyz(5),dltyz(6),dltyz(7),dltyz(8);dltyz(9),dltyz(10),dltyz(11),1];
Axy=[dltxy(1),dltxy(2),dltxy(3),dltxy(4);dltxy(5),dltxy(6),dltxy(7),dltxy(8);dltxy(9),dltxy(10),dltxy(11),1];

Rxz= Kxz \ Axz ; % inv(Kxz)*Axz;
Ryz= Kyz \ Ayz ; % inv(Kyz)*Ayz;
Rxy= Kxy \ Axy ; % inv(Kxy)*Axy;

Rxz=Rxz/max([norm(Rxz(1:3,1)),norm(Rxz(1:3,2)),norm(Rxz(1:3,3))]);
Ryz=Ryz/max([norm(Ryz(1:3,1)),norm(Ryz(1:3,2)),norm(Ryz(1:3,3))]);
Rxy=Rxy/max([norm(Rxy(1:3,1)),norm(Rxy(1:3,2)),norm(Rxy(1:3,3))]);


end


function [G]=reconstructpoint(R1,K1,imageTip_1,R2,K2,imageTip_2)


U1 =   R1(1:3,1:3)' * (K1 \ ([imageTip_1,1]'));
O1 = - R1(1:3,1:3)' * R1(1:3,4);

U2 =   R2(1:3,1:3)' * (K2 \ ([imageTip_2,1]'));
O2 = - R2(1:3,1:3)' * R2(1:3,4);

A1=[-dot(U1,U1),dot(U1,U2);-dot(U1,U2),dot(U2,U2)];
B1=[-dot(O2-O1,U1);-dot(O2-O1,U2)];

lambda1=linsolve(A1,B1);

G= (O1+O2+lambda1(1)*U1+lambda1(2)*U2)'/2;

end