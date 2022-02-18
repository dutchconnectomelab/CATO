function S=DEMO_SimulateDWMRI_voxel()
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
% DW-MRI (Diffusion-Weighted MRI) dataset of 1 voxel.
%
%-USE--------------------------------------------------------------------------------
% S=DEMO_SimulateDWMRI_voxel;
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

%Underlying fiber orientation 1
fiber_orientation1=[cos(20*pi/180) sin(20*pi/180) 0];

%Underlying fiber orientation 2 (i.e. we have a fiber crossing in this voxel)
fiber_orientation2=[cos(100*pi/180) sin(100*pi/180) 0];

%The next lines implement the simulation of the DW-MRI signal
S=zeros(size(GradientOrientations,1),1);
for i=1:size(GradientOrientations,1)
    S(i)=S0*(SimulateDWMRI(fiber_orientation1,GradientOrientations(i,:))+SimulateDWMRI(fiber_orientation2,GradientOrientations(i,:)))/2;
end

fprintf(1,'\nIf you use SimulateDWMRI.m please cite the following work:\n');
fprintf(1,'A. Barmpoutis et al. "Adaptive kernels for multi-fiber reconstruction",\n'); 
fprintf(1,'In the Proceedings of IPMI, 2009, pp. 338-349.\n');