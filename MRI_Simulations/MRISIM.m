
%This file written to "simulate" a 2D Cartesian MRI sequence
%Assume cat is a representation of T1s or T2s of a homogeneous image
%
%
%
% THIS VERSION IS SHIT. It tries to assume 1 magnetization vector per
% voxel. BAD IDEA! So many wonky refocusings happening.
%
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
%CONSTANTS
G_max = 3; %G/cm
GFE = -0.75;
GPEmax = 3;
GFEwind = 3;

N_PE = 32;   %Number of phase encodes
N_RO = 32;   %Number of readout points
ImageSize = 16;  %"size" of image, cm
IMS = 32;      %Number of pixels per side of square image
Del_y = 0.1; %Spatial resolution of voxels, cm
gratio = 4248; %Gyromagnetic ratio, Hz/G

T1scale = 5000;
T2scale = 0.0005;

      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INITIALIZE CLASSICAL "VECTOR MATRIX"

%512x512x3 matrix, a.k.a. the image size and 3 cartesian vector components
%per voxel
rho0 = zeros(IMS,IMS,4);
%Assume each voxel has the same initial magnetization (in z-direction)
for nn = 1:1:IMS
    rho0(:,:,3) = 1; 
    rho0(:,:,4) = 1;  %(4th element always 1)
end

cat = imresize(cat,IMS/512);
birb = imresize(birb,IMS/512);

T1map = double(birb)/256*T1scale;
T2map = double(cat)/256*T2scale;

Del_ky = 1/(N_PE*Del_y);  %k-space voxel spacing
ky_max = 1/2*(N_PE-1)*Del_ky;  %cm^-1
Tau_PE = ky_max/gratio/G_max;  %Phase encode time 

Tau_RO = -2*Tau_PE*GFEwind/GFE;

%RF PULSE
%Assume equal excitation for first pulse (puts every voxel in x-direction)


for tt =1:1:N_PE
    rho1 = Rot(rho0,Ry(pi/2),IMS);
    GPE(tt) = -GPEmax+2*GPEmax*(tt-1)/(N_PE-1);
    %PHASE ENCODE GRADIENT (X)  (1-2)
    rho1 = xGrad(rho1,ImageSize,IMS,gratio*GPE(tt)*Tau_PE);  %Apply phase encode
    rho1 = yGrad(rho1,ImageSize,IMS,gratio*GFEwind*Tau_PE); %Apply freq encode
    rho2 = Relax(rho1,Tau_PE,IMS,T1map,T2map);

    imshow(rho2(:,:,1))
    %FORM ECHO  (2-4)
    %rho2 = yGrad

    for nn =1:1:N_RO
        Taustep(nn) = Tau_RO/N_RO*nn;
        rho3(:,:,:,nn) = yGrad(rho2,ImageSize,IMS,gratio*GFE*Taustep(nn));
        Signalx(nn,tt) = sum(rho3(:,:,1,nn),'all');
        Signaly(nn,tt) = sum(rho3(:,:,2,nn),'all');
        Signalz(nn,tt) = sum(rho3(:,:,3,nn),'all');
    end
end

figure(2)
subplot(2,2,1)
imshow(cat)

subplot(2,2,2)
var =(ifft2(cat));
imshow(fftshift(real(var)))

subplot(2,2,3)
imshow(Signalx)

subplot(2,2,4)
Reconimg =(fft2(Signalx+1i*Signaly));
imshow(real(Reconimg))


%%
function val = xGrad(rho,ImageSize,IMS,thetam)  
    for mm = 1:1:IMS
        xdisp(mm) = -ImageSize/2+ImageSize/(IMS-1)*(mm-1);
        for nn = 1:1:IMS
            val(mm,nn,:) = Rz(thetam*xdisp(mm))*squeeze(rho(mm,nn,:));
        end
    end
end
function val = yGrad(rho,ImageSize,IMS,thetam)  
    for mm = 1:1:IMS
        for nn = 1:1:IMS
            ydisp(mm) = -ImageSize/2+ImageSize/(IMS-1)*(nn-1);
            val(mm,nn,:) = Rz(thetam*ydisp(mm))*squeeze(rho(mm,nn,:));
        end
    end
end
function val = Relax(rho,tau,IMS,T1map,T2map)
    for mm = 1:1:IMS
        for nn =1:1:IMS
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
function val = Rot(rho,rmat,IMS)
    for mm = 1:1:IMS
        for nn = 1:1:IMS
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


