
%load calibration result
load('C:\Users\Fruit Flies\Desktop\temp\calibration_easyWandData.mat')

g=1; %image you want to check in the videos
nb=120; %number of points in the grid
mylength=0.02; %boundaries of the volume in real space

%detect bigger circles in each image

im1=imread('C:\Users\Fruit Flies\Desktop\temp\xy.tif');
[c_xy,radius_xy]=imfindcircles(im1,[35 50],'ObjectPolarity','dark');
c_xy=[c_xy(1),511-c_xy(2)];

im2=imread('C:\Users\Fruit Flies\Desktop\temp\yz.tif');
[c_yz,radius_yz]=imfindcircles(im2,[35 50],'ObjectPolarity','dark');
c_yz=[c_yz(1),511-c_yz(2)];

im3=imread('C:\Users\Fruit Flies\Desktop\temp\xz.tif');
[c_xz,radius_xz]=imfindcircles(im3,[35 50],'ObjectPolarity','dark');
c_xz=[c_xz(1),511-c_xz(2)];

%create volume where the sphere should be

x_min = -mylength;     % in real-space units (m)
x_max = mylength;
x_mesh = nb;

y_min = -mylength;
y_max = +mylength;
y_mesh = nb;

z_min = -mylength;
z_max = +mylength;
z_mesh = nb;

X=linspace(x_min,x_max,x_mesh);
Y=linspace(y_min,y_max,y_mesh);
Z=linspace(z_min,z_max,z_mesh);

dlt=load('C:\Users\Fruit Flies\Desktop\temp\calibration_dltCoefs.csv');
dlt_xz=dlt(:,3);
dlt_yz=dlt(:,2);
dlt_xy=dlt(:,1);
V = false(x_mesh, y_mesh, z_mesh) ;
count = 0;

%for each point in real space, check if its projection in the three cameras
%images is inside each circle

for i=1:x_mesh
    for j=1:y_mesh
        for k=1:z_mesh
            
             pixel_xz=dlt_inverse(dlt_xz,[X(i),Y(j),Z(k)]);
             
            if (norm(pixel_xz-c_xz)>radius_xz)
                continue ;
            end
            
             pixel_yz=dlt_inverse(dlt_yz,[X(i),Y(j),Z(k)]);
             
            if (norm(pixel_yz-c_yz)>radius_yz)
                continue ;
            end
            
             pixel_xy=dlt_inverse(dlt_xy,[X(i),Y(j),Z(k)]);
            
            if   norm(pixel_xy-c_xy)<radius_xy 
                
                V(i,j,k) = true ;
                count = count + 1 ;
            end
            
            
        end
    end
end

% allocate "result"
RES = zeros(count, 3) ;
c=0;
for i=1:x_mesh
    for j=1:y_mesh
        for k=1:z_mesh
            
            if(V(i,j,k))
                c=c+1;
                RES(c,:) = [X(i), Y(j), Z(k)] ;
            end
        end
    end
end

scatter3(RES(:,1),RES(:,2),RES(:,3))
axis equal;

