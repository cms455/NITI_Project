function data_out = load_onscale_flexdata(dirname, filename, nrows, downsample_t, startrow)
% LOAD_ONSCALE_FLEXDATA - Load data from directory of numbered OnScale *.flxdato files
%
% data_out = load_onscale_flexdata(dirname, filename, nrows, downsample_t, startrow)
%
% Parameters
% ----------
% dirname : [string] path to folder containing a collection of OnScale results
% filename : [string] file prefix of *.flxdato files (e.g. filename1.flxdato, filename2.flxdato)
% nrows : [int] number of rows to return from each data file array
% downsample_t : [int] load in time points, skipping this number of files between each load
% startrow : [int] index of the first row to return from each data file array
%
% Returns
% -------
% data_out : [struct] contains the arrays read from the *.flxdato files
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
    startrow = 1;
end
if nargin < 5
    startrow = 1;
end

% Get all mat files in directory dirname starting with filename
datafile_info = dir([fullfile(dirname, filename), '*']);
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
    datanumber(i) = str2double(datafiles{i}((length(filename)+1):end));
end
[~, sortorder] = sort(datanumber);
datafiles = datafiles(sortorder);

% Load first data file and preallocate fields
data_in = read_flxdato(fullfile(dirname, datafiles{1}));
nfields = length(data_in.DRecord);
field_names = cell(nfields, 1);
for i = 1:nfields
    field_names{i} = replace(data_in.DRecord(i).A.Name, ' ', '');
end

data_out.t = zeros(nt, 1);
for i = 1:nfields
    switch field_names{i}
    case 'xcrd'
        data_out.x = data_in.DRecord(i).B;
    case 'ycrd'
        data_out.y = data_in.DRecord(i).B(startrow:(startrow+nrows-1));
    case 'xcrs'
        tmp_size = size(data_in.DRecord(i).B);
        if nargin < 3 || isempty(nrows)
            nrows = tmp_size(2);
        end
        if ~isempty(data_in.DRecord(i).B)
            data_out.xs = data_in.DRecord(i).B(:, startrow:(startrow+nrows-1))';
        else
            data_out.xs = [];
        end
    case 'ycrs'
        tmp_size = size(data_in.DRecord(i).B);
        if nargin < 3 || isempty(nrows)
            nrows = tmp_size(2);
        end
        if ~isempty(data_in.DRecord(i).B)
            data_out.ys = data_in.DRecord(i).B(:, startrow:(startrow+nrows-1))';
        else
            data_out.ys = [];
        end
    case 'xvel'
        tmp_size = size(data_in.DRecord(i).B);
        if nargin < 3 || isempty(nrows)
            nrows = tmp_size(2);
        end
        data_out.xvel = zeros(nrows, tmp_size(1), nt);
    case 'yvel'
        tmp_size = size(data_in.DRecord(i).B);
        if nargin < 3 || isempty(nrows)
            nrows = tmp_size(2);
        end
        data_out.yvel = zeros(nrows, tmp_size(1), nt);
    case 'xdsp'
        tmp_size = size(data_in.DRecord(i).B);
        if nargin < 3 || isempty(nrows)
            nrows = tmp_size(2);
        end
        data_out.xdsp = zeros(nrows, tmp_size(1), nt);
    case 'ydsp'
        tmp_size = size(data_in.DRecord(i).B);
        if nargin < 3 || isempty(nrows)
            nrows = tmp_size(2);
        end
        data_out.ydsp = zeros(nrows, tmp_size(1), nt);
    case 'pldt'
        tmp_size = size(data_in.DRecord(i).B);
        data_out.pldt = zeros(tmp_size(2), nt);
    case 'tmpr'
        tmp_size = size(data_in.DRecord(i).B);
        if nargin < 3 || isempty(nrows)
            nrows = tmp_size(2);
        end
        data_out.tmpr = zeros(nrows, tmp_size(1), nt);
    end
end


% Fill fields
counter = 1;
for j = 1:downsample_t:nfiles
    data_in = read_flxdato(fullfile(dirname, datafiles{j}));
    for i = 1:nfields
        switch field_names{i}
        case 'xvel'
            tmp_size = size(data_in.DRecord(i).B);
            if nargin < 3 || isempty(nrows)
                nrows = tmp_size(2);
            end
            data_out.t(counter) = data_in.DRecord(i).A.Time;
            data_out.xvel(:,:,counter) = data_in.DRecord(i).B(:, startrow:(startrow+nrows-1))';
        case 'yvel'
            tmp_size = size(data_in.DRecord(i).B);
            if nargin < 3 || isempty(nrows)
                nrows = tmp_size(2);
            end
            data_out.t(counter) = data_in.DRecord(i).A.Time;
            data_out.yvel(:,:,counter) = data_in.DRecord(i).B(:, startrow:(startrow+nrows-1))';
        case 'xdsp'
            tmp_size = size(data_in.DRecord(i).B);
            if nargin < 3 || isempty(nrows)
                nrows = tmp_size(2);
            end
            data_out.t(counter) = data_in.DRecord(i).A.Time;
            data_out.xdsp(:,:,counter) = data_in.DRecord(i).B(:, startrow:(startrow+nrows-1))';
        case 'ydsp'
            tmp_size = size(data_in.DRecord(i).B);
            if nargin < 3 || isempty(nrows)
                nrows = tmp_size(2);
            end
            data_out.t(counter) = data_in.DRecord(i).A.Time;
            data_out.ydsp(:,:,counter) = data_in.DRecord(i).B(:, startrow:(startrow+nrows-1))';
        case 'pldt'
            tmp_size = size(data_in.DRecord(i).B);
            data_out.t(counter) = data_in.DRecord(i).A.Time;
            data_out.pldt(:,counter) = data_in.DRecord(i).B(2,:)';
        case 'tmpr'
            tmp_size = size(data_in.DRecord(i).B);
            if nargin < 3 || isempty(nrows)
                nrows = tmp_size(2);
            end
            data_out.t(counter) = data_in.DRecord(i).A.Time;
            data_out.tmpr(:,:,counter) = data_in.DRecord(i).B(:, startrow:(startrow+nrows-1))';
        end
    end
    counter = counter + 1;
end