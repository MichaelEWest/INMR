
%This file written to "simulate" a 2D Cartesian MRI sequence
%Assume cat is a representation of T1s or T2s of a homogeneous image
clear 
clc
cat = imread('CatGrayscale.gif');
birb = imread('BirdGrayscale.gif');
%%
%Gradient strengths-- assume 30 mT/m max
%assume voxel spacing is 1mm center to center (512 ~ .5 meters)
%      __
% RF  |  |
%    _|  |________________________________________________
%     0  1   
%            1    2
%             ____
%GPE (x)     /____\
%    _______/______\______________________________________
%           \ ____ /
%            \____/
%
%                    2           3
%GFE (y)             _________________________________
%    ______         /                                 \___
%          \       /
%           \_____/
%            1   2
%                                ()
%                             ()()()()
%                           ()()    ()()
% SIG___________________()()()        ()()()_______________
%                    2           3             
%
%
%       |          TE            |
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIALIZE CLASSICAL "VECTOR MATRIX"

%512x512x3 matrix, a.k.a. the image size and 3 cartesian vector components
%per voxel
rho0 = zeros(512,512,4);
%Assume each voxel has the same initial magnetization (in z-direction)
for nn = 1:1:512
    rho0(:,:,3) = 1; 
    rho0(:,:,4) = 1;  %(4th element always 1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONSTANTS
G_max = 3; %G/cm

N_PE = 72;   %Number of phase encodes
ImageSize = 50;  %"size" of image, cm
Del_y = 0.1; %Spatial resolution of voxels, cm
gratio = 4248; %Gyromagnetic ratio, Hz/G

T1scale = 5000;
T2scale = 0.05;
T1map = double(birb)/256*T1scale;
T2map = double(cat)/256*T2scale;

Del_ky = 1/(N_PE*Del_y);  %k-space voxel spacing
ky_max = 1/2*(N_PE-1)*Del_ky;  %cm^-1
Tau_PE = ky_max/gratio/G_max;  %Phase encode time 


%RF PULSE
%Assume equal excitation for first pulse (puts every voxel in x-direction)
rho1 = Rot(rho0,Ry(pi/2));


%PHASE ENCODE GRADIENT (X)  (1-2)
GPE = 3;
GFE1 = 3;
rho1 = xGrad(rho1,ImageSize,gratio*GPE*Tau_PE);
rho1 = yGrad(rho1,ImageSize,gratio*GFE1*Tau_PE);
rho2 = Relax(rho1,Tau_PE,T1map,T2map);


%FORM ECHO (2-3)
GFE2=-0.75;
Tau_RO = Tau_PE*-GFE2/GFE1;
%rho2 = yGrad


function val = xGrad(rho,ImageSize,thetam)  
    for mm = 1:1:512
        xdisp(mm) = -ImageSize/2+ImageSize/512*mm;
        for nn = 1:1:512
            val(mm,nn,:) = Rz(thetam*xdisp(mm))*squeeze(rho(mm,nn,:));
        end
    end
end
function val = yGrad(rho,ImageSize,thetam)  
    for mm = 1:1:512
        for nn = 1:1:512
            ydisp(mm) = -ImageSize/2+ImageSize/512*nn;
            val(mm,nn,:) = Rz(thetam*ydisp(mm))*squeeze(rho(mm,nn,:));
        end
    end
end
function val = Relax(rho,tau,T1map,T2map)
    for mm = 1:1:512
        for nn =1:1:512
            T1 = T1map(mm,nn);
            T2 = T2map(mm,nn);
            E1 = exp(-tau/T1);
            E2 = exp(-tau/T2);
            Relmat =  [E2 0 0 0;
                       0 E2 0 0;
                       0 0 E1 1-E1;
                       0 0  0 1];
           val(mm,nn,:) = Relmat*squeeze(rho(mm,nn,:));
        end
    end


end
function val = Rot(rho,rmat)
    for mm = 1:1:512
        for nn = 1:1:512
            val(mm,nn,:) = rmat*squeeze(rho(mm,nn,:));
        end
    end
end
function val = Rx(theta)
    val = [1 0           0          0;
           0 cos(theta) -sin(theta) 0;
           0 sin(theta)  cos(theta) 0;
           0 0           0          1];
end
function val = Ry(theta)
    val = [ cos(theta) 0 sin(theta) 0;
            0          1 0          0;
           -sin(theta) 0 cos(theta) 0;
            0          0 0          1];
end
function val = Rz(theta)
    val = [cos(theta) -sin(theta) 0 0;
           sin(theta)  cos(theta) 0 0;
           0           0          1 0;
           0           0          0 1];
end


