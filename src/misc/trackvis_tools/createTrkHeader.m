function header = createTrkHeader(varargin)
% CREATETRKVHEADER Create header for TrackVis files (in TRK-format).
%
% HEADER = createTrkHeader(VARARGIN) creates HEADER using default
% values or parameters specified as name-value pairs in VARARGIN.
%
% See: http://www.trackvis.org/docs/?subsect=fileformat

assert(mod(numel(varargin), 2) == 0, ...
    'CATO:createTrkHeader:WrongNumberArgs', ...
    'Optional input arguments must come in pairs.');

% Default values for header:
header.id_string = 'TRACK';
header.dim = [0 0 0];
header.voxel_size = [0 0 0];
header.origin = [0 0 0];
header.n_scalars = 0;
header.scalar_name = '';
header.n_properties = 0;
header.property_name = '';
header.vox_to_ras = zeros(4);
header.reserved = '';
header.voxel_order = '';
header.pad2 = '';
header.image_orientation_patient = [0 0 0 0 0 0];
header.pad1 = '';
header.invert_x = 0;
header.invert_y = 0;
header.invert_z = 0;
header.swap_xy = 0;
header.swap_yz = 0;
header.swap_zx = 0;
header.n_count = 0;
header.version = 2;
header.hdr_size = 1000;

% Update default values with name-value parameter pairs provided as input
% arguments.
while ~isempty(varargin)
    varName = varargin{1};
    varValue = varargin{2};
    
    header.(varName) = varValue;
    
    varargin(1:2) = [];
end