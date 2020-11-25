function [alphaVec, X2, r, exitFlag] = lmsolverDTI(alphaVec,B,S,weight)
% LMSOLVERDTI Non-linear tensor fitting using Levenberg–Marquardt algorithm
%
%   INPUT VARIABLES
%   alphaVec:
%   Initial tensor estimate.
%
%   B:
%   Gradient matrix (B-matrix).
%
%   S:
%   Signal intensities of a voxel.
%
%   weight:
%   Weighting of scan directions in non-linear fitting.
%
%   OUTPUT VARIABLES
%   alphaVec:
%   Updated tensor estimate.
%
%   X2:
%   Chi-square of the model.
%
%   r:
%   Residual or difference between the experimental and estimated values.
%
%   exitFlag:
%   exitFlag = 0 is succes and exitFlag = 4 indicates that the maximum
%   number of iterations was reached without convergence

%   Modified from: http://people.duke.edu/~hpgavin/m-files/lm.m
%   by Henri Gavin, Dept. Civil & Environ. Engineering, Duke Univ. 4 May 2016
%   modified from: http://octave.sourceforge.net/optim/function/leasqr.html
%   using references by
%   Press, et al., Numerical Recipes, Cambridge Univ. Press, 1992, Chapter 15.
%   Sam Roweis       http://www.cs.toronto.edu/~roweis/notes/lm.pdf
%   Manolis Lourakis http://www.ics.forth.gr/~lourakis/levmar/levmar.pdf
%   Hans Nielson     http://www2.imm.dtu.dk/~hbn/publ/TR9905.ps
%   Mathworks        optimization toolbox reference manual
%   K. Madsen, H.B., Nielsen, and O. Tingleff
%   http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
%   
%   LICENSE based on leasqr.m from the optimization package:
%   Copyright (C) 1992-1994 Richard Shrager
%   Copyright (C) 1992-1994 Arthur Jutan <jutan@charon.engga.uwo.ca>
%   Copyright (C) 1992-1994 Ray Muzic <rfm2@ds2.uh.cwru.edu>
%   Copyright (C) 2010-2019 Olaf Till <i7tiol@t-online.de>
%
%   This program is free software; you can redistribute it and/or modify it under
%   the terms of the GNU General Public License as published by the Free Software
%   Foundation; either version 3 of the License, or (at your option) any later
%   version.
%   This program is distributed in the hope that it will be useful, but WITHOUT
%   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
%   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
%   details.
%   You should have received a copy of the GNU General Public License along with
%   this program; if not, see <http://www.gnu.org/licenses/>.


% In comparison with original:
% y = dependent data = S
% p = parameters = alpha
% t = independent data = B
% weight = weight
% func = @(t,p,c) exp(B*p)


epsilon_1 = 1e-3;
epsilon_2 = 1e-3;
epsilon_3 = 1e-1;
epsilon_4 = 1e-1;
lambda_0 = 0.01;
MaxIter = 100; % FASTER
lambda_up = 11;
lambda_down = 9;

exitFlag = 0;

if nargin < 4
    weight = abs(1/(S'*S)) .* ones(length(S), 1);
end

Npar = length(alphaVec);

%%

Shat = exp(B * alphaVec);
J = B .* repmat(S, [1 7]);
JtWJ  = J' * ( J .* ( weight * ones(1,Npar) ) ); % JtWJ = J'*W*J;
r = S - Shat; % residual error between model and data
JtWdy = J' * ( weight .* r );
X2 = Shat' * (Shat .* weight);


if ( max(abs(JtWdy)) < epsilon_1 )
    exitFlag = 1;
    return;
end

counter = 0;
lambda  = lambda_0;

DoF = numel(r);

while true
    counter = counter + 1;
    
    h = ( JtWJ + lambda*diag(diag(JtWJ)) ) \ JtWdy;
    
    alpha_try = alphaVec + h;
    
    % optional: apply constraints
    
    r = S - exp(B * alpha_try);
    
    X2_try = r' * (r .* weight);
    
    rho = (X2 - X2_try) / (h' * (lambda * h + JtWdy));
    
    if ( max(abs(h ./ alphaVec)) < epsilon_2 )
%         fprintf(' **** Convergence in Parameters **** \n')
%         fprintf(' **** epsilon_2 = %e\n', epsilon_2);
        break;
    end    
    
    if rho > epsilon_4
        X2 = X2_try;
        alphaVec = alpha_try;
        lambda = max(lambda / lambda_down, 1e-7);
        
        Shat = exp(B * alphaVec);
        J = B .* repmat(S, [1 7]);
        JtWJ  = J' * ( J .* ( weight * ones(1,Npar) ) );
        r = S - Shat;
        JtWdy = J' * ( weight .* r );        
    else
        lambda = min(lambda * lambda_up, 1e7);
    end
    
    if ( max(abs(JtWdy)) < epsilon_1 ) && ( counter > 2 )
%         fprintf(' **** Convergence in r.h.s. ("JtWdy")  **** \n')
%         fprintf(' **** epsilon_1 = %e\n', epsilon_1);
        break;
    end
    if ( X2/DoF < epsilon_3 &&  counter > 2 )
%         fprintf(' **** Convergence in reduced Chi-square  **** \n')
%         fprintf(' **** epsilon_3 = %e\n', epsilon_3);
        break;
    end
    if ( counter == MaxIter )
%         disp(' !! Maximum Number of Iterations Reached Without Convergence !!')
            exitFlag = 4;
        break;
    end
    
end

