function Data = read_flxdato(FName)
% READ_FLXDATO - Read data from OnScale *.flxdato output file
%
% This function is a modified version of the function provided by PZFlex
% (now OnScale) for reading *.flxdato files, accessed from:
%
% https://www.mathworks.com/matlabcentral/fileexchange/36170-matlab-script-to-read-in-flxdato-files-from-flexlab?s_tid=prof_contriblnk
% Accessed on Febr 27, 2019.
%
% The original function has been modified to account for cases where the 
% data arrays are large (> 8 MB for double arrays, > 4 MB for single arrays).
% Modified by John J. Pitre, Jr.
%
% Data = read_flxdato(FName)
%
% Parameters
% ----------
% FName : [string] *.flxdato file to read
%
% Returns
% -------
% Data : [struct] structure containing data read from file
%
% ---


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Data OUT1 Reader
%   Description:    Reads PZFlex data file into Matlab
%   Author:         Weidlinger
%   Date:           12/09/11
%   Version:        1.0
%
%   History:        1.0         Initial version
%                   2019-04-18  Edit to allow reading large arrays in 8M
%                               byte blocks for float64 data - JJP
%                   2019-04-21  Edit to ignore data block for variables
%                               with data length zero, fixed bug in large
%                               array import - JJP
%                   2020-04-05  Edit to allow reading large arrays in 4M
%                               byte blocks for float32 data - JJP
%
%
% Copyright (c) 2012, PZFlex
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
%   
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set precision of open / close bytes
OCForm = 'int32';
% Open the file
Fid = fopen(FName); % Replace FName with name of flex data file
                    % e.g. 'FILENAME.flxdato'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Headers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read header record 1 open byte
OpenBytes = fread(Fid,1,OCForm);
% Decide precision used in file
if (OpenBytes == 8)
    
    % Single precision
    IntForm = 'int32';
    FloatForm = 'float32';
    
elseif (OpenBytes == 12)
    
    % Double precision
    IntForm = 'int64';
    FloatForm = 'float64';
    
else
    
    % Unknown header format
    disp('Error: Unknown header format.');
    Data = -1;
    return;
    
end;
% Read remainder of header record 1
Header1.Label = fread(Fid,4,'*char')';
Header1.Type = fread(Fid,1,IntForm);
CloseBytes = fread(Fid,1,OCForm);
% Check header and write to main structure
if (OpenBytes ~= CloseBytes)    
    Data = -1;
    return;   
end;
% Save to structure
Data.Header1 = Header1;
% Read header record 2
OpenBytes = fread(Fid,1,OCForm);             
Header2.FType = fread(Fid,20,'*char')';
Header2.FForm = fread(Fid,20,'*char')';
Header2.Date = fread(Fid,80,'*char')';
Header2.Code = fread(Fid,20,'*char')';
Header2.User = fread(Fid,200,'*char')';
Header2.Extra1 = fread(Fid,20,'*char')';
Header2.Extra2 = fread(Fid,20,'*char')';
CloseBytes = fread(Fid,1,OCForm); 
% Check header and write to main structure
if (OpenBytes ~= CloseBytes)    
    Data = -1;
    return;   
end;
% Save to structure
Data.Header2 = Header2;
% Read header record 3
OpenBytes = fread(Fid,1,OCForm);             
Header3.Title = fread(Fid,200,'*char')';
Header3.Tag = fread(Fid,80,'*char')';
CloseBytes = fread(Fid,1,OCForm); 
% Check header and write to main structure
if (OpenBytes ~= CloseBytes)    
    Data = -1;
    return;   
end;
% Save to structure
Data.Header3 = Header3;
% Read header record 4
OpenBytes = fread(Fid,1,OCForm);             
Header4.Version = fread(Fid,1,IntForm);
Header4.NPartitions = fread(Fid,1,IntForm);
Header4.IPartition = fread(Fid,1,IntForm); 
Header4.IntBytes = fread(Fid,1,IntForm);
Header4.RealBytes = fread(Fid,1,IntForm);
Header4.Char1Bytes = fread(Fid,1,IntForm);
Header4.Char2Bytes = fread(Fid,1,IntForm);
Header4.Char3Bytes = fread(Fid,1,IntForm); 
Header4.Extra1 = fread(Fid,1,IntForm); 
Header4.Extra2 = fread(Fid,1,IntForm); 
CloseBytes = fread(Fid,1,OCForm); 
% Check header and write to main structure
if (OpenBytes ~= CloseBytes)    
    Data = -1;
    return;   
end;
% Save to structure
Data.Header4 = Header4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data records
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise index
Index = 1;
% Loop through data records
while (feof(Fid) == 0)

    % Read data record A
    OpenBytes = fread(Fid,1,OCForm);
    
    % If cannot read openbytes the file end has been reached
    if (isempty(OpenBytes))
        
        break;
        
    end;
    
    % Continue record A
    DRecord(Index).A.Name = fread(Fid,20,'*char')';
    DRecord(Index).A.Length = fread(Fid,1,IntForm);
    DRecord(Index).A.NTmStp = fread(Fid,1,IntForm);
    DRecord(Index).A.Time = fread(Fid,1,FloatForm);
    DRecord(Index).A.NDim = fread(Fid,1,IntForm);
    DRecord(Index).A.IRange = fread(Fid,1,IntForm);
    DRecord(Index).A.JRange = fread(Fid,1,IntForm);
    DRecord(Index).A.KRange = fread(Fid,1,IntForm);
    DRecord(Index).A.IType = fread(Fid,1,IntForm);
    DRecord(Index).A.IPart = fread(Fid,1,IntForm);
    DRecord(Index).A.IExt1 = fread(Fid,1,IntForm); 
    DRecord(Index).A.IExt2 = fread(Fid,1,IntForm); 
    CloseBytes = fread(Fid,1,OCForm);
    
    % Check header
    if (OpenBytes ~= CloseBytes)
        Data = -1;
        return;
    end;
    
    
    
    % Interpret file type
    if (DRecord(Index).A.IPart == 0)

        % File is type 0, Record A and B(1)

        % Set data format and size
        if (DRecord(Index).A.IType == 1)

            % Read reals
            DataForm = FloatForm;

            % Set data size
            NRead = DRecord(Index).A.Length;
            MatSize = [DRecord(Index).A.IRange ...
                DRecord(Index).A.JRange ...
                DRecord(Index).A.KRange];

        elseif (DRecord(Index).A.IType == 2)

            % Read integers
            DataForm = IntForm;

            % Set data size
            NRead = DRecord(Index).A.Length;
            MatSize = [DRecord(Index).A.IRange ...
                DRecord(Index).A.JRange ...
                DRecord(Index).A.KRange];

        elseif (DRecord(Index).A.IType < 0)

            % Read characters
            DataForm = '*char';

            % Set data size
            NRead = DRecord(Index).A.Length * ...
                abs(DRecord(Index).A.IType);
            MatSize = [abs(DRecord(Index).A.IType) ...
                DRecord(Index).A.Length ...
                1];

        end;


        % Read Record B(1)
        
        % JJP - If the variable is empty, don't do anything - for example if we
        % export xcrs or ycrs for a non-skewed grid, the header information is
        % available, but there is no data written (and so no open or close
        % bytes).
        if DRecord(Index).A.Length > 0
            
            % JJP - For very large arrays, the binary seems to be written in 8
            % million byte blocks, each with their own open and close bytes.
            % The last blocks may be less than 8 million. The following loop
            % accounts for this. If the array will require > 8M bytes, compute
            % the number of reads required and read in blocks, checking open
            % and close bytes for each. 
            if strcmp(DataForm, 'float64')
                % Preliminary constants, etc
                bytes_per_double = 8;
                bytes_per_block = 8E6;
                doubles_per_block = bytes_per_block/bytes_per_double;
                data_bytes_to_read = NRead*bytes_per_double;
                nblocks = ceil(data_bytes_to_read/8E6);

                if nblocks > 1 % We need to read in multiple blocks
                    RawData = zeros(doubles_per_block, nblocks);

                    % Read in all but the last block - these are all 8M bytes of
                    % doubles along with their open and close bytes
                    for iblock = 1:(nblocks-1)
                        OpenBytes = fread(Fid,1,OCForm);
                        RawData(:, iblock) = fread(Fid, doubles_per_block, DataForm);
                        CloseBytes = fread(Fid,1,OCForm);
                        % Check header
                        if (OpenBytes ~= CloseBytes)
                            Data = -1;
                            return;
                        end;
                    end

                    % The last block is shorter
                    doubles_remaining = NRead - (nblocks-1)*doubles_per_block;
                    OpenBytes = fread(Fid,1,OCForm);
                    RawData(1:doubles_remaining, nblocks) = fread(Fid, doubles_remaining, DataForm);
                    CloseBytes = fread(Fid,1,OCForm);
                    % Check header
                    if (OpenBytes ~= CloseBytes)
                        Data = -1;
                        return;
                    end;

                    % Reshape and store results in data record
                    DRecord(Index).B = reshape(RawData(1:NRead), MatSize);

                else % There is only one block to read, so use original code
                    OpenBytes = fread(Fid,1,OCForm);
                    RawData = fread(Fid,NRead,DataForm)';
                    DRecord(Index).B = reshape(RawData,MatSize);
                    CloseBytes = fread(Fid,1,OCForm);
                    % Check header
                    if (OpenBytes ~= CloseBytes)
                        Data = -1;
                        return;
                    end;
                end

                
            % JJP - For very large arrays, the floats seems to be written in
            % 4 million byte blocks, each with their own open and close bytes.
            % The last blocks may be less than 4 million. The following loop
            % accounts for this. If the array will require > 4M bytes, compute
            % the number of reads required and read in blocks, checking open
            % and close bytes for each. 
            elseif strcmp(DataForm, 'float32') 
                
                % Preliminary constants, etc
                bytes_per_float = 4;
                bytes_per_block = 4E6;
                floats_per_block = bytes_per_block/bytes_per_float;
                data_bytes_to_read = NRead*bytes_per_float;
                nblocks = ceil(data_bytes_to_read/4E6);

                if nblocks > 1 % We need to read in multiple blocks
                    RawData = zeros(floats_per_block, nblocks);

                    % Read in all but the last block - these are all 4M bytes of
                    % singles along with their open and close bytes
                    for iblock = 1:(nblocks-1)
                        OpenBytes = fread(Fid,1,OCForm);
                        RawData(:, iblock) = fread(Fid, floats_per_block, DataForm);
                        CloseBytes = fread(Fid,1,OCForm);
                        % Check header
                        if (OpenBytes ~= CloseBytes)
                            Data = -1;
                            return;
                        end;
                    end

                    % The last block is shorter
                    floats_remaining = NRead - (nblocks-1)*floats_per_block;
                    OpenBytes = fread(Fid,1,OCForm);
                    RawData(1:floats_remaining, nblocks) = fread(Fid, floats_remaining, DataForm);
                    CloseBytes = fread(Fid,1,OCForm);
                    % Check header
                    if (OpenBytes ~= CloseBytes)
                        Data = -1;
                        return;
                    end;

                    % Reshape and store results in data record
                    DRecord(Index).B = reshape(RawData(1:NRead), MatSize);

                else % There is only one block to read, so use original code
                    OpenBytes = fread(Fid,1,OCForm);
                    RawData = fread(Fid,NRead,DataForm)';
                    DRecord(Index).B = reshape(RawData,MatSize);
                    CloseBytes = fread(Fid,1,OCForm);
                    % Check header
                    if (OpenBytes ~= CloseBytes)
                        Data = -1;
                        return;
                    end;
                end
                
            % JJP - I'm not sure if there are similar issues with large
            % arrays of integers. These should default to the original code
            % for now.
            else
                OpenBytes = fread(Fid,1,OCForm);
                RawData = fread(Fid,NRead,DataForm)';
                DRecord(Index).B = reshape(RawData,MatSize);
                CloseBytes = fread(Fid,1,OCForm);
                % Check header
                if (OpenBytes ~= CloseBytes)
                    Data = -1;
                    return;
                end;
            end
        else % the data field is empty
            DRecord(Index).B = [];
        end
        
    else

        % File is type 1, Record A, A1 and B(2)

        % Read Record A1
        OpenBytes = fread(Fid,1,OCForm);
        DRecord(Index).A1.IBegin = fread(Fid,1,IntForm); 
        DRecord(Index).A1.IEnd = fread(Fid,1,IntForm);
        DRecord(Index).A1.JBegin = fread(Fid,1,IntForm); 
        DRecord(Index).A1.JEnd = fread(Fid,1,IntForm);
        DRecord(Index).A1.KBegin = fread(Fid,1,IntForm); 
        DRecord(Index).A1.KEnd = fread(Fid,1,IntForm);
        CloseBytes = fread(Fid,1,OCForm);

        % Check header
        if (OpenBytes ~= CloseBytes)
            Data = -1;
            return;
        end;


        % Set data format and size
        if (DRecord(Index).A.IType == 1)

            % Read reals
            DataForm = FloatForm;

            % Set data size
            MatSize = [(DRecord(Index).A1.IEnd - DRecord(Index).A1.IBegin + 1) ...
                       (DRecord(Index).A1.JEnd - DRecord(Index).A1.JBegin + 1) ...
                       (DRecord(Index).A1.KEnd - DRecord(Index).A1.KBegin + 1)];
            NRead = MatSize(1) * MatSize(2) * MatSize(3);

        elseif (DRecord(Index).A.IType == 2)

            % Read integers
            DataForm = IntForm;

            % Set data size
            MatSize = [(DRecord(Index).A1.IEnd - DRecord(Index).A1.IBegin + 1) ...
                       (DRecord(Index).A1.JEnd - DRecord(Index).A1.JBegin + 1) ...
                       (DRecord(Index).A1.KEnd - DRecord(Index).A1.KBegin + 1)];
            NRead = MatSize(1) * MatSize(2) * MatSize(3);

        elseif (DRecord(Index).A.IType < 0)

            % Read characters
            DataForm = '*char';

            % Set data size
            MatSize = [abs(DRecord(Index).A.IType) ...
                       (DRecord(Index).A1.IEnd - DRecord(Index).A1.IBegin + 1) ...
                       1];
            NRead = MatSize(1) * MatSize(2);

        end;


        % Read Record B(2)
        OpenBytes = fread(Fid,1,OCForm);
        RawData = fread(Fid,NRead,DataForm)';
        DRecord(Index).B = reshape(RawData,MatSize);
        CloseBytes = fread(Fid,1,OCForm);

        % Check header
        if (OpenBytes ~= CloseBytes)
            Data = -1;
            return;
        end;

    end;


    % If data is character then flip contents of record B
    if (DRecord(Index).A.IType < 0)

        DRecord(Index).B = DRecord(Index).B';

    end;
    

    
    % Increment index
    Index = Index + 1;

    
end;
% Save to structure
Data.DRecord = DRecord';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close the file
fclose(Fid);