function hdr = load_nifti_partially(niftifile, voxelMask)
% Load large 4D nifti files partially
%
% Adapted from load_nifti.m:
%    Original Author: Doug Greve
%    Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%    Terms and conditions for use, reproduction, distribution and contribution
%    are found in the 'FreeSurfer Software License Agreement' contained
%    in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%    https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense

    % load header
    hdr = load_nifti_hdr(niftifile);

    % get number of voxels and gradient directions
    nVoxels = prod(hdr.dim(2:4));
    ndirections = hdr.dim(5);

    % check datatype
    switch(hdr.datatype)
        case 4
            datatype = 'int16';
            nbytes = 2;
        case 8
            datatype = 'int32';
            nbytes = 4;
        case 16
            datatype = 'single';
            nbytes = 4;
        case 64
            datatype = 'double';
            nbytes = 8;
        case 512
            datatype = 'uint16';
            nbytes = 2;
        case 768
            datatype = 'uint32';
            nbytes = 4;
        otherwise
            error('Data type not supported.');
    end

    % preallocate memory for output
    hdr.partialvol = zeros(nnz(voxelMask), ndirections, datatype);

    % read data
    fid = fopen(niftifile, 'r', hdr.endian);

    fseek(fid,round(hdr.vox_offset),'bof');
    for i = 1:ndirections
         thisVol = fread(fid, nVoxels, datatype);
         hdr.partialvol(:, i) = thisVol(voxelMask);
    end

    fclose(fid);

    % rescale if necessary
    if(hdr.scl_slope ~= 0)
%         fprintf('Rescaling NIFTI: slope = %g, intercept = %g\n', ...
%             hdr.scl_slope, hdr.scl_inter);
        hdr.partialvol = hdr.partialvol * hdr.scl_slope + hdr.scl_inter;
    end
end
