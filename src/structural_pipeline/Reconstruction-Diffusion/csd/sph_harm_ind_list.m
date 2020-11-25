function [m_list, n_list] = sph_harm_ind_list(sh_order)
% function [m_list, n_list] = SPH_HARM_IND_LIST(sh_order)
%
% SPH_HARM_IND_LIST provides the order m (n_harmonics x 1) and
% degree (n_harmonics x 1) of all symmetric spherical harmonics of degree
% less or equal to sh_order.

if mod(sh_order, 2) ~= 0 || sh_order<0
    error('sph_harm_ind_list:InvalidInput', ...
        'sh_order must be an even integer >=0');
end

% n_list describes the degree
n_range = 0:2:(sh_order + 1);
n_list = repelem(n_range, n_range * 2 + 1);

% m_list are the associated order
m_list = nan(1, length(n_list));

offset = 0;
for ii = n_range
    m_list((offset:offset + 2 * ii)+1) = -ii:ii;
    offset = offset + 2 * ii + 1;
end
