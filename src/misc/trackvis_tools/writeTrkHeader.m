function writeTrkHeader(fid, header)
% WRITETRKHEADER Write TrackVis header to file.
%
% writeTrkHeader(FID, HEADER) writes HEADER to file with identifier FID.
%
% See: http://www.trackvis.org/docs/?subsect=fileformat

% Based on:
%   trk_write.m by John Colby
%   https://github.com/johncolby/along-tract-stats/blob/master/trk_write.m
%   which is distributed under the GPL-3.0 License

% Check orientation
[~, ix] = max(abs(header.image_orientation_patient(1:3)));
[~, iy] = max(abs(header.image_orientation_patient(4:6)));
iz = 1:3;
iz([ix iy]) = [];
assert(isequal([ix iy iz], [1 2 3]), ...
    ['Image orientation of the original image is assumed LPS. Other ', ...
    'orientations are not supported.']);

% Write TrackVis header
fwrite(fid, charpad(header.id_string, 1, 6), 'char*1');
fwrite(fid, header.dim, 'int16');
fwrite(fid, header.voxel_size, 'float');
fwrite(fid, header.origin, 'float');
fwrite(fid, header.n_scalars, 'int16');
fwrite(fid, charpad(header.scalar_name, 10, 20)', 'char*1');
fwrite(fid, header.n_properties, 'int16');
fwrite(fid, charpad(header.property_name, 10, 20)', 'char*1');
fwrite(fid, header.vox_to_ras, 'float');
fwrite(fid, charpad(header.reserved, 1, 444), 'char*1');
fwrite(fid, charpad(header.voxel_order, 1, 4), 'char*1');
fwrite(fid, charpad(header.pad2, 1, 4), 'char*1');
fwrite(fid, header.image_orientation_patient, 'float');
fwrite(fid, charpad(header.pad1, 1, 2), 'char*1');
fwrite(fid, header.invert_x, 'uchar');
fwrite(fid, header.invert_y, 'uchar');
fwrite(fid, header.invert_z, 'uchar');
fwrite(fid, header.swap_xy, 'uchar');
fwrite(fid, header.swap_yz, 'uchar');
fwrite(fid, header.swap_zx, 'uchar');
fwrite(fid, header.n_count, 'int');
fwrite(fid, header.version, 'int');
fwrite(fid, header.hdr_size, 'int');

end

function charsPadded = charpad(chars, n, m)
charsPadded = repmat(char(0), n ,m);
charsPadded(1:size(chars, 1), 1:size(chars,2)) = chars;
end