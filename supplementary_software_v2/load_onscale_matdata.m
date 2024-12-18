function data_out = load_onscale_matdata(dirname, filename, nrows, downsample_t)
% LOAD_ONSCALE_MATDATA - Load data from directory of numbered OnScale *.mat outputs
%
% data_out = load_onscale_matdata(dirname, filename, nrows, downsample_t)
%
% Parameters
% ----------
% dirname : [string] path to folder containing a collection of OnScale results
% filename : [string] file prefix of *.mat files (e.g. filename1.mat, filename2.mat)
% nrows : [int] number of rows to return from each data file array
% downsample_t : [int] load in time points, skipping this number of files between each load
%
% Returns
% -------
% data_out : [struct] contains the arrays read from the *.mat files
%
% Author: John J. Pitre, Jr.
%
% Pitre, JJ, MA Kirby, DS Li, TT Shen, RK Wang, M O'Donnell, and I Pelivanov.
%    Nearly-incompressible transverse isotropy (NITI) of cornea elasticity: 
%    model and experiments with acoustic micro-tapping OCE. Scientific Reports 
%    (2020).
%
% ---
    if nargin < 4
        downsample_t = 1;
    end
    
    % Get all mat files in directory dirname starting with filename
    datafile_info = dir([fullfile(dirname, filename), '*.mat']);
    nfiles = length(datafile_info);
    nt = length(1:downsample_t:nfiles);
    
    % Convert datafile struct into a cell array of filenames
    datafiles = cell(nfiles, 1);
    for i = 1:nfiles
        datafiles{i} = datafile_info(i).name;
    end
    
    % Sort filenames in numerical order
    datanumber = zeros(nfiles, 1);
    for i = 1:nfiles
        datanumber(i) = str2double(datafiles{i}((length(filename)+1):(end-4)));
    end
    [~, sortorder] = sort(datanumber);
    datafiles = datafiles(sortorder);
    
    % Flags for import variables
    flags.itime = 0;
    flags.ftime = 0;
    flags.xcrd = 0;
    flags.ycrd = 0;
    flags.xcrs = 0;
    flags.ycrs = 0;
    flags.xvel = 0;
    flags.yvel = 0;
    flags.xdsp = 0;
    flags.ydsp = 0;
    flags.pldt = 0;
    
    % Load first data file and preallocate fields
    data_in = load(fullfile(dirname, datafiles{1}));
    data_fields = fieldnames(data_in);
    if any(endsWith(data_fields, 'itime'))
        data_out.step = zeros(nt, 1);
        flags.itime = 1;
    end
    if any(endsWith(data_fields, 'ftime'))
        data_out.t = zeros(nt, 1);
        flags.ftime = 1;
    end
    if any(endsWith(data_fields, 'xcrd'))
        tmp_var = getvar(data_in, 'xcrd');
        data_out.x = tmp_var;
        flags.xcrd = 1;
    end
    if any(endsWith(data_fields, 'ycrd'))
        tmp_var = getvar(data_in, 'ycrd');
        if nargin < 3 || isempty(nrows)
            nrows = length(tmp_var);
        end
        data_out.y = tmp_var(1:nrows);
        flags.ycrd = 1;
    end
    if any(endsWith(data_fields, 'xcrs'))
        tmp_var = getvar(data_in, 'xcrs');
        tmp_size = size(tmp_var);
        if nargin < 3 || isempty(nrows)
            nrows = tmp_size(2);
        end
        data_out.xs = tmp_var(:, 1:nrows)';
        flags.xcrs = 1;
    end
    if any(endsWith(data_fields, 'ycrs'))
        tmp_var = getvar(data_in, 'ycrs');
        tmp_size = size(tmp_var);
        if nargin < 3 || isempty(nrows)
            nrows = tmp_size(2);
        end
        data_out.ys = tmp_var(:, 1:nrows)';
        flags.ycrs = 1;
    end
    if any(endsWith(data_fields, 'xvel'))
        tmp_var = getvar(data_in, 'xvel');
        tmp_size = size(tmp_var);
        if nargin < 3 || isempty(nrows)
            nrows = tmp_size(2);
        end
        data_out.xvel = zeros(nrows, tmp_size(1), nt);
        flags.xvel = 1;
    end
    if any(endsWith(data_fields, 'yvel'))
        tmp_var = getvar(data_in, 'yvel');
        tmp_size = size(tmp_var);
        if nargin < 3 || isempty(nrows)
            nrows = tmp_size(2);
        end
        data_out.yvel = zeros(nrows, tmp_size(1), nt);
        flags.yvel = 1;
    end
    if any(endsWith(data_fields, 'xdsp'))
        tmp_var = getvar(data_in, 'xdsp');
        tmp_size = size(tmp_var);
        if nargin < 3 || isempty(nrows)
            nrows = tmp_size(2);
        end
        data_out.xdsp = zeros(nrows, tmp_size(1), nt);
        flags.xdsp = 1;
    end
    if any(endsWith(data_fields, 'ydsp'))
        tmp_var = getvar(data_in, 'ydsp');
        tmp_size = size(tmp_var);
        if nargin < 3 || isempty(nrows)
            nrows = tmp_size(2);
        end
        data_out.ydsp = zeros(nrows, tmp_size(1), nt);
        flags.ydsp = 1;
    end
    if any(endsWith(data_fields, 'pldt'))
        tmp_var = getvar(data_in, 'pldt');
        tmp_size = size(tmp_var);
        data_out.pldt = zeros(tmp_size(2), nt);
        flags.pldt = 1;
    end
    
    counter = 1;
    for i = 1:downsample_t:nfiles
        data_in = load(fullfile(dirname, datafiles{i}));
        if flags.itime
            tmp_var = getvar(data_in, 'itime');
            data_out.step(counter) = tmp_var;
        end
        if flags.ftime
            tmp_var = getvar(data_in, 'ftime');
            data_out.t(counter) = tmp_var;
        end
        if flags.xvel
            tmp_var = getvar(data_in, 'xvel');
            tmp_var = tmp_var';
            tmp_var = tmp_var(1:nrows,:,:);
            data_out.xvel(:,:,counter) = tmp_var;
        end
        if flags.yvel
            tmp_var = getvar(data_in, 'yvel');
            tmp_var = tmp_var';
            tmp_var = tmp_var(1:nrows,:,:);
            data_out.yvel(:,:,counter) = tmp_var;
        end
        if flags.xdsp
            tmp_var = getvar(data_in, 'xdsp');
            tmp_var = tmp_var';
            tmp_var = tmp_var(1:nrows,:,:);
            data_out.xdsp(:,:,counter) = tmp_var;
        end
        if flags.ydsp
            tmp_var = getvar(data_in, 'ydsp');
            tmp_var = tmp_var';
            tmp_var = tmp_var(1:nrows,:,:);
            data_out.ydsp(:,:,counter) = tmp_var;
        end
        if flags.pldt
            tmp_var = getvar(data_in, 'pldt');
            data_out.pldt(:,counter) = tmp_var(2,:)';
        end
        counter = counter + 1;
    end
    
    
function field_var = getvar(datastruct, varname)
    data_fields = fieldnames(datastruct);
    field_name = data_fields(endsWith(data_fields, varname));
    field_var  = getfield(datastruct, field_name{1});