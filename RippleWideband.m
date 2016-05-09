classdef RippleWideband < WidebandManager
    %% Access FLAC-compressed wideband data from the Ripple NIP Box.
    %
    % Note that the "channel" argument for member functions is actually
    % just a channel index and therefore runs from 1 to channelCount(),
    % even if the actual channel numbers (in self.channel) are not actually
    % continuous. Hopefully this makes coding a lot easier since Matlab
    % doesn't have anything like an iterator.
    
    
    properties(Constant=true)
        REF_GROUP_SIZE = 32; % Size of the re-referencing groups
    end
    
    properties(GetAccess=public, SetAccess=protected)
        dataFiles   % Filenames for individual data files
        scaleFactor % Converts A/D units to physical units (see channels!)
        header      % Struct array containing file metadata
    end
    
    methods
        function self = RippleWideband(basedir ,ignoreHeader)
            %% RippleWideband: Create a Wideband reader for the specified data directory
            % Ideally, the directory contains a mat/txt header (from 
            % rippleToFlac versions after v726909...). 
            %
            % If no mat files are present, it tries to guess
            % metadata for earlier versions without the headers. 
            %
            % If multiple mat files are present, you're on your own!
            
            if(exist('ignoreHeader', 'var') && ignoreHeader)
                self.init_without_header(basedir);
            else
            
                f = dir(fullfile(basedir, '*.mat'));
                if isempty(f)
                    warning('RippleWideband:MissingMatHeader', ...
                        'No header file in basedir, checking for FLAC...');
                
                    self.header = [];
                    self.init_without_header(basedir);
                elseif length(f) == 1
                    self.init_from_header(basedir, f.name)
                else
                    % This is a work-around for the spike-sorter spam that ends
                    % up in the data directories (bletch). Of course, it
                    % doesn't quite work because the channels have been
                    % reordered too.
                    matnames = strrep({f.name}, '.mat', '');
                    txtdir = dir(fullfile(basedir, '*.txt'));
                    txtnames = strrep({txtdir.name}, '.txt', '');

                    intersect_names = intersect(matnames, txtnames);
                    if length(intersect_names) == 1
                        self.init_from_header(basedir, [matnames{1} '.mat']);
                    else
                        error('RippleWideband:TooManyMat', ...
                        'Multiple mat files found--there should only be one!');
                    end
                end
            end
        end         
        
        function data = getRawByIndex(self, chan_index, start, stop)
            %% getRawByIndex Return unscaled (i.e, in AD int16s) data          
            % data = getRawByIndex(..., start, stop) returns raw A/D values             
            % sample# start and sample# stop for the specified channel. 
            %
            % If start and stop are not provided, return all available data.
            
            if length(chan_index) > 1        
                if ~(exist('start', 'var') && exist('stop', 'var'))
                    start = 1;
                    stop = self.nSamples;
                end
                
                data = nan(length(chan_index), stop-start+1);
                for ch=1:length(chan_index)                   
                        data(ch,:) = audioread(self.dataFiles{chan_index(ch)}, [start, stop], 'native');                                                                
                end
            else % This does exactly the same thing, avoiding an extra copy
                if exist('start', 'var') && exist('stop', 'var')
                    data = audioread(self.dataFiles{chan_index}, [start, stop], 'native')';
                else
                    data = audioread(self.dataFiles{chan_index}, 'native')';
                end
            end
        end
        
        
        function data = getDataByIndex(self, chan_index, start, stop)
            %% getDataByIndex Return scaled (i.e, in microvolts) data
            % data = getDataByIndex(..., start, stop) returns data between 
            % sample# start and sample# stop for the specified channel.
            %
            % If start and stop are not provided, return all available data. 
            if ~(exist('start', 'var') && exist('stop', 'var'))
                    start = 1;
                    stop = self.nSamples;
            end
            
            data = bsxfun(@times, ...
                double(self.getRawByIndex(chan_index, start, stop)), ...
                self.scaleFactor(chan_index));
        end
        
        function data = getDataByTime(self, chan_index, t_start, t_stop)
            %% getDataByTime Return scaled (i.e, in microvolts) data
            % data = getDataByIndex(..., start, stop) returns data between 
            % start and stop seconds (Ripple time-base) for the specified 
            % channel.
            %
            % If start and stop are not provided, return all available data. 
            
            if ~exist('t_start', 'var') && ~exist('t_stop', 'var')
                if length(chan_index) > 1
%                     data = nan(length(chan_index), stop-start+1, 'double');
                    for ch=1:length(chan_index)
                        % You don't actually need a cast here since audioread
                        % only updates a proper subset of the data and thus
                        % inherits data's type. That said, leave it
                        % because this seems like something that The Mathworks
                        % may break at any moment.
                        data(ch,:) = double(audioread(self.dataFiles{chan_index(ch)}, 'native'));
                    end
                else
                    data = double(audioread(self.dataFiles{chan_index}, 'native'))';
                end
            else
                start = floor(t_start * self.samplingRate);
                if start < 1
                    start = 1;
                end
                
                stop = floor(t_stop * self.samplingRate);
                
                if length(chan_index) > 1
                    data = nan(length(chan_index), stop-start+1, 'double');
                    for ch=1:length(chan_index)
                        % You don't actually need a cast here since audioread
                        % only updates a proper subset of the data and thus
                        % inherits data's type. That said, leave it
                        % because this seems like something that The Mathworks
                        % may break at any moment.
                        data(ch,:) = double(audioread(self.dataFiles{chan_index(ch)}, [start, stop], 'native'));
                    end
                else
                    % This does exactly the same thing, avoiding an extra copy
                    % (or any copies, if the optimizer is smart). Obviously,
                    % you need a cast here or the multiplication is wrong
                    % below.
                    if exist('start', 'var') && exist('stop', 'var')
                        data = double(audioread(self.dataFiles{chan_index}, [start, stop], 'native'))';
                    end
                end
            end
            
            % You *do* need bsxfun here!
            data =  bsxfun(@times, data, self.scaleFactor(chan_index));
        end
       
        
        function refGroup = getRefGroup(self, chan_index)
            %% getRefGroup Return re-referencing group for the specified channel *index*.
            refGroup = self.refGroups(chan_index);
        end
        
        function n = channelCount(self)
            %% channelCount Returns the total number of channels
            n = length(self.channels);
        end
        
    end
    
    methods(Access=protected)
        
        function verifyFiles(self)
            %% verifyFiles(): Ensure that all files can be opened and have
            % the same length, bit depth, etc. This is mostly paranoia, but
            % hey...
            ai = audioinfo(self.dataFiles{1});
            comment = ai.Comment;
            
            assert(self.channelCount == length(self.dataFiles), ...
                'RippleWideband:MissingFiles', 'Missing FLAC files?');
            
            for ii=1:length(self.dataFiles)
                ai = audioinfo(self.dataFiles{ii}); 
                if ai.NumChannels ~= 1
                    error('RippleWideband:FlacChannels', ...
                        'Only one channel per FLAC file is supported');
                end
                
                if ai.BitsPerSample ~= 16
                    error('RippleWideband:FlacBits', ...
                        'Only 16 bits/sample FLAC files are supported');
                end
                
                if ai.SampleRate ~= self.samplingRate
                    error('RippleWideband:FlacFS', ...
                        'Sampling Rate varies across files: found %d and %d', ...
                        fs, ai.SampleRate);
                end
                
                if ai.TotalSamples ~= self.nSamples
                    error('RippleWideband:FlacLength', ...
                        'Files have different numbers of samples');
                end
                
                if ~isequal(ai.Comment, comment)
                    warning('RippleWideband:FlacComment', ...
                        'Files have different comments/scale factors');
                end
            end
        end
        
        
        function r = findRefGroups(self)
            %% Assign channels to reference groups. This assigns all channels
            % from the same headstage to the same reference group.
            if isnumeric(self.channels)
                channels = self.channels;
            else
                channels = [self.channels.number]';
            end
            
            % This is a little hacky, but the idea is that we want to
            % calculate channel//32, since each headstage has 32 channels.
            % After that, we make the values consecutive (since we may not
            % have 4 headstages plugged into each port).
            r = ceil(channels/self.REF_GROUP_SIZE);
            
            total_counts = tabulate(r);
            n_missing  = cumsum(total_counts(:,2)==0);
            for ii=1:size(total_counts, 1)
                if total_counts(ii,2) && n_missing(ii)
                    index = (r==total_counts(ii,1));
                    r(index) = r(index) - n_missing(ii);
                end
            end
        end
        
        
        function init_from_header(self, basedir, matfile)
            %% Initialize the class using information from the header .mat file
            
            self.header = load(fullfile(basedir, matfile));
            self.samplingRate  = self.header.sampling_frequency;
            
            self.channels = self.header.channels;           
            self.header = rmfield(self.header, 'channels'); %Might as well save 500k
            
            self.dataFiles = cell(self.channelCount(), 1);
            self.scaleFactor = nan(self.channelCount(), 1);
            for ii=1:length(self.dataFiles)
                self.dataFiles{ii} = fullfile(basedir, self.channels(ii).filename);
                self.scaleFactor(ii) = self.channels(ii).d2a_scale_factor;
            end
            
            
            ai = audioinfo(self.dataFiles{1});
            self.nSamples = ai.TotalSamples;
            
            self.verifyFiles();
            self.refGroups = self.findRefGroups();
        end
        
        
        function init_without_header(self, basedir) 
            %% Older data directories are missing the mat/txt header files
            % that describe the recording parameters. This makes a
            % reasonable guess at them. Hopefully we can eventually phase
            % this out as everything gets converted but until then...
            
            
            %% Find data files in the directory
            f = dir(fullfile(basedir, '*_ch*.flac'));
            if isempty(f)
                error('WidebandManager:NoFiles', ...
                    'No FLAC files found in %s', basedir');
            end
                        
            self.dataFiles = fullfile(basedir, sort({f.name}));
            ch = regexpi({f.name}, '_ch(\d+).flac', 'tokens'); 
            for ii=1:length(ch)
                self.channels(ii).number = str2double(ch{ii}{1}{1});
            end
            %%self.channels.number = cellfun(@(x) str2double(x{1}{1}), ch)';
                       
          
            %% Pull out metadata and check it
            ai = audioinfo(self.dataFiles{1});
            self.samplingRate = ai.SampleRate;
            self.nSamples = ai.TotalSamples;
            
            %% Be paranoid and check all the files for length/sampling rate
            % (This causes N different files to be opened, so you may want
            %  to skip this if file access is very slow)
            self.verifyFiles();
            
           
            if self.samplingRate == 30000
                self.scaleFactor = repmat(0.208319345683157, size(self.channels));
            else
                % Need to replace this with actual plexon scale
                self.scaleFactor = repmat(1, size(self.channels));  %#ok<RPMT1>
            end
            
            %% Set up the reference data
            self.refGroups = self.findRefGroups();    
        end
    end    
    
    methods(Hidden=true)
        %% Hide pointless functions inherited from handle
        function addlistener(varargin) ; end
        function ge(varargin); end
        function gt(varargin); end
        function le(varargin); end
        function lt(varargin); end
    end
end



