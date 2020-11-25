function header = readTrkHeader(fid)
% READTRKVHEADER Read header from fiber file (in TRK-format).
%
% HEADER = readTrkHeader(FID) reads the HEADER of the TrackVis file that is
% pointed to by the file identifier FID.
%
% See: http://www.trackvis.org/docs/?subsect=fileformat

% Based on:
%   trk_read.m by John Colby
%   https://github.com/johncolby/along-tract-stats/blob/master/trk_read.m
%   which is distributed under the GPL-3.0 License

header.id_string = fread(fid, [1, 6], '*char*1');
header.dim = fread(fid, [1, 3], 'int16');
header.voxel_size = fread(fid, [1, 3], 'float');
header.origin = fread(fid, [1, 3], 'float');
header.n_scalars = fread(fid, 1, 'int16');
header.scalar_name = fread(fid, [20, 10], '*char*1')';
header.n_properties = fread(fid, 1, 'int16');
header.property_name = fread(fid, [20, 10], '*char*1')';
header.vox_to_ras = fread(fid, [4, 4], 'float');
header.reserved = fread(fid, [1, 444], '*char*1');
header.voxel_order = fread(fid, [1, 4], '*char*1');
header.pad2 = fread(fid, [1, 4], '*char*1');
header.image_orientation_patient = fread(fid, [1, 6], 'float');
header.pad1 = fread(fid, [1, 2], '*char*1');
header.invert_x = fread(fid, 1, '*uchar');
header.invert_y = fread(fid, 1, '*uchar');
header.invert_z = fread(fid, 1, '*uchar');
header.swap_xy = fread(fid, 1, '*uchar');
header.swap_yz = fread(fid, 1, '*uchar');
header.swap_zx = fread(fid, 1, '*uchar');
header.n_count = fread(fid, 1, 'int');
header.version = fread(fid, 1, 'int');
header.hdr_size = fread(fid, 1, 'int');