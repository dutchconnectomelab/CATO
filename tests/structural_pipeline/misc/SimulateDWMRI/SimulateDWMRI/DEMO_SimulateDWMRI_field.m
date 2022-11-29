function S=DEMO_SimulateDWMRI_field()
%-fanDTasia ToolBox------------------------------------------------------------------
% This Matlab script is part of the fanDTasia ToolBox: a Matlab library for Diffusion 
% Weighted MRI (DW-MRI) Processing, Diffusion Tensor (DTI) Estimation, High-order 
% Diffusion Tensor Analysis, Tensor ODF estimation, Visualization and more.
%
% A Matlab Tutorial on DW-MRI can be found in:
% http://www.cise.ufl.edu/~abarmpou/lab/fanDTasia/tutorial.php
%
%-CITATION---------------------------------------------------------------------------
% If you use this software please cite the following work:
% A. Barmpoutis et al. "Adaptive kernels for multi-fiber reconstruction", 
% In the Proceedings of IPMI, 2009, pp. 338-349. 
%
%-DESCRIPTION------------------------------------------------------------------------
% This demo script shows how to use SimulateDWMRI.m in order to generate a synthetic
% DW-MRI (Diffusion-Weighted MRI) field of size 32x32.
%
%-USE--------------------------------------------------------------------------------
% S=DEMO_SimulateDWMRI_field;
%
% S: is the simulated DWMRI signal
%
%-DISCLAIMER-------------------------------------------------------------------------
% You can use this source code for non commercial research and educational purposes 
% only without licensing fees and is provided without guarantee or warrantee expressed
% or implied. You cannot repost this file without prior written permission from the 
% authors. If you use this software please cite the following work:
% A. Barmpoutis et al. "Adaptive kernels for multi-fiber reconstruction", 
% In the Proceedings of IPMI, 2009, pp. 338-349. 
%
%-AUTHOR-----------------------------------------------------------------------------
% Angelos Barmpoutis, PhD
% Computer and Information Science and Engineering Department
% University of Florida, Gainesville, FL 32611, USA
% abarmpou at cise dot ufl dot edu
%------------------------------------------------------------------------------------


b_value=1500; %SimulateDWMRI.m generates signal using b=1500s/mm2
S0=1;  %Here we set the S0
%Here we construct a set of 21 unit vectors
GradientOrientations=[0.1639 0.5115 0.8435;0.1176 -0.5388 0.8342;0.5554 0.8278 -0.0797;-0.4804 0.8719 0.0948;0.9251 -0.0442 0.3772;0.7512 -0.0273 -0.6596;0.1655 -0.0161 0.9861;0.6129 -0.3427 0.7120;0.6401 0.2747 0.7175;-0.3724 -0.3007 0.8780;-0.3451 0.3167 0.8835;0.4228 0.7872 0.4489;0.0441 0.9990 0.0089;-0.1860 0.8131 0.5515;0.8702 0.4606 0.1748;-0.7239 0.5285 0.4434;-0.2574 -0.8032 0.5372;0.3515 -0.8292 0.4346;-0.7680 -0.4705 0.4346;0.8261 -0.5384 0.1660;0.9852 -0.0420 -0.1660];


%The next lines implement the simulation of the DW-MRI signal
S=zeros(32,32,size(GradientOrientations,1));
for i=1:size(GradientOrientations,1)
    fprintf(1,'.');
    for x=1:32
        for y=1:32
            
           f1_flag=0;
           f2_flag=0;
           
           %Define the underlying fiber orientation of fiber bundle #1
           if x*x+y*y>16*16 & x*x+y*y<32*32
               v=[y/x -1 0];v=v/sqrt(v*v');
               fiber_orientation1=v;
               f1_flag=1;
           end
           
           %Define the underlying fiber orientation of fiber bundle #2 (i.e. we have 2-fiber crossing) 
           if x<y+10 & x>y-10
               fiber_orientation2=[sqrt(2)/2 sqrt(2)/2 0];
               f2_flag=1;
           end
           
           %Here we fill the rest of the voxels
           if f1_flag==0 & f2_flag==1
               fiber_orientation1=fiber_orientation2;
           elseif f1_flag==1 & f2_flag==0
               fiber_orientation2=fiber_orientation1;
           elseif f1_flag==0 & f2_flag==0
               fiber_orientation1=[0 0 1];
               fiber_orientation2=[0 0 1];
           end
           
           
            S(x,y,i)=S0*(SimulateDWMRI(fiber_orientation1,GradientOrientations(i,:))+SimulateDWMRI(fiber_orientation2,GradientOrientations(i,:)))/2;
        end
    end
end

fprintf(1,'\nDo you want to save the simulated dataset in fanDTasia (FDT) format\n');
fprintf(1,'so that you can open it with fanDTasia the on-line DW-MRI processing tool?\n');
fprintf('For YES type 1, for NO type 0.\n');
a=input('Answer [1,0] :');
if a==1
    for x=1:32
        for y=1:32
            fdt(x,y,1,1)=S0;
            fdt(x,y,1,[2:size(GradientOrientations,1)+1])=S(x,y,:);
        end
    end
    writeFDT('synthetic_dataset.fdt',fdt);
    f=fopen('synthetic_dataset.txt','w');
    fprintf(f,'%.8f %.8f %.8f %.2f\n',1,0,0,10);
    for i=1:size(GradientOrientations,1)
        fprintf(f,'%.8f %.8f %.8f %.2f\n',GradientOrientations(i,1),GradientOrientations(i,2),GradientOrientations(i,3),1500);
    end
    fclose(f);
    
    fprintf(1,'\nFiles saved: synthetic_dataset.fdt, synthetic_dataset.txt\n');
end


fprintf(1,'\nIf you use SimulateDWMRI.m please cite the following work:\n');
fprintf(1,'A. Barmpoutis et al. "Adaptive kernels for multi-fiber reconstruction",\n'); 
fprintf(1,'In the Proceedings of IPMI, 2009, pp. 338-349.\n');