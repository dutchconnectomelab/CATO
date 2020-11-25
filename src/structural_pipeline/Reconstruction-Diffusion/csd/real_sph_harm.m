function real_sh = real_sph_harm(m, n, theta, phi)
% real_sh = REAL_SPH_HARM (m, n, theta, phi) computes real spherical
% harmonics for each m (nparams x 1) and n (nparams x 1) of each location
% theta (npoints x 1) and phi (npoints x 1).
%
% Where the real harmonic $Y^m_n$ is defined to be:
%
% Imag($Y^m_n$) * sqrt(2)     if m > 0
% $Y^0_n$                     if m = 0
% Real($Y^|m|_n$) * sqrt(2)   if m < 0
%
% real_sh is the real harmonic $Y^m_n$ sampled at `theta` and `phi`.

% prepare variables
mabs = abs(m(:));
n = n(:);
phi = phi(:);
theta = theta(:);
theta = cos(theta);

% calculate the spherical harmonics sh at point (phi, theta)
sh = nan(length(n), length(theta));

for i = 1:length(n)
    
    mi= mabs(i);
    ni = n(i);
    
    Lmn = legendre(ni, theta);
    Lmn = Lmn(mi+1, :)';
    
    a1 = (2 * ni + 1) / (4*pi);
    a2 = factorial(ni-mi) / factorial(ni+mi);
    C = sqrt(a1*a2);
    
    sh(i, :) = C * Lmn .* exp(1i * mi * phi);
    
end

% convert spherical harmonics to real spherical harmonics.
real_sh = sh;
real_sh(m > 0, :) = imag(sh(m > 0, :));
real_sh(m <= 0, :) = real(sh(m <= 0, :));
real_sh(m ~= 0, :) = sqrt(2) * real_sh(m ~= 0, :);

end
