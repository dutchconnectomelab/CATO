function out=SimulateDWMRI(fiber_orientation,dwmri_gradient_orientation)
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
% This method computes an accurate continuous approximation of 
% the DW-MRI signal using adaptive kernels (A. Barmpoutis et al. 
% "Adaptive kernels for multi-fiber reconstruction", In the 
% Proceedings of IPMI, 2009). The approximation is based on Fig. 3
% of the article which shows the plot of the DW-MRI signal attenuation 
% from water molecules whose diffusion is restricted inside a cylindrical 
% fiber of radius rho=5 micrometers and length L=5mm using b-value = 1500s/mm2.
%
%-USE--------------------------------------------------------------------------------
% out=SimulateDWMRI(fiber_orientation,dwmri_gradient_orientation)
%
% Parameters: fiber_orientation is a 3D unit vector
%            dwmri_gradient_orientation is a 3D unit vector
% Example: SimulateDWMRI([1 0 0],[0 1 0])
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

fiber(1,[1:3])=fiber_orientation([1:3]);
gradient(1,[1:3])=dwmri_gradient_orientation([1:3]);
fiber=fiber/sqrt(fiber*fiber');
gradient=gradient/sqrt(gradient*gradient');
cosine=fiber*gradient';
cp=[0.0672 0.1521 0.3091 0.4859 0.6146];
num_of_cp=length(cp);
out=0;
for k=1:num_of_cp
    out=out+cp(k)*basisN(abs(cosine),k,num_of_cp);
end
out=out+cp(num_of_cp)*basisN(abs(cosine),num_of_cp+1,num_of_cp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=basisN(argument,segment_id,total_segments)
%Evaluates the 2nd-order spline basis
%argument must be a real in [0...1]
%segment_id must be an integer in [1...total_segments]
argument=1-argument;
for i=1:length(argument)
	if (argument(i)>= (segment_id-2)/total_segments) & (argument(i)<(segment_id-1)/total_segments)
        t=(argument(i)-(segment_id-2)/total_segments)/((segment_id-1)/total_segments-(segment_id-2)/total_segments);
        out(i)=0.5*t*t;
	elseif (argument(i)>= (segment_id-1)/total_segments) & (argument(i)<segment_id/total_segments)
        t=(argument(i)-(segment_id-1)/total_segments)/((segment_id)/total_segments-(segment_id-1)/total_segments);
        out(i)=-t*t+t+0.5;
	elseif (argument(i)>= (segment_id)/total_segments) & (argument(i)<(segment_id+1)/total_segments)
        t=(argument(i)-(segment_id)/total_segments)/((segment_id+1)/total_segments-(segment_id)/total_segments);
        out(i)=0.5*(1-t)*(1-t);
	else
        out(i)=0;
    end
end