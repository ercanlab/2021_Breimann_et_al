function varargout = matVis(varargin)
%% matVis(data1, data2, ... , 'PropertyName', PropertyValue, ...)
%**************************************************************************
% GUI for displaying data sets of arbitrary dimension using 2D images and 
% 1D profiles. Key features of matVis include
% - Easy data input: load variables from the workspace or various image file
%   formats including Zeiss lsm-files and .mat-files
% - Easy navigation through data sets of arbitrary dimension - Easy to
% learn user interface with tooltips for each control item - Visualize
% multiple data sets of identical dimensions in parallel - Histogram-guided
% contrast adjustment including gamma correction and
%   algorithmic thresholding
% - Large variety of colormaps
% - Image and plot filter
% - Image projections including image tiling
% - RGB display including color-coded max- and mean-projections ('RGB
%   stretch')
% - Extensive region of interest manager
% - Video generator, which can combine the content of mutiple windows into  
%   a single customizable movie
% - Specify scales, names and units for the data dimensions
% - Use of alpha data including contrast and gamma adjustments
% - 2D histogram for image data / alpha data with histogram guided 2D 
%   contrast adjustment
% - Figure and function handles can be accessed from 'outside', making
%   matVis easy to integrate with other GUIs or Matlab functions
% - Save configuration (settings, window positions, ...) for future use
% - Use 'startPar' argument to selectively overwrite configuration
%   paramaters (see below for list of start parameters)
% - Check for updates from within the main GUI
%  
% Input Arguments
%   -none-                      Select image file(s) (including matrices 
%                               saved as .mat-file, multi-image tif, Zeiss 
%                               lsm-files, RGB files, animated gif and most 
%                               other image formats) via uigetfile. Files 
%                               can be combined into a single matrix or 
%                               opened "in parallel".
% OR
%   data1, data2, ...           Data sets of arbitrary but identical dimension
% OR
%   'pathname\filename.ext'     Path- + filename of desired file (works only for
%                               one file)
%
% Optional Arguments            'PropertyName', PropertyValue
%       'dimScale'              Scale for dimensions, i.e. values of bounds
%                               of 'pixel-locations' (like the x and y paramteres in
%                               imagesc(x,y,d))
%       'dimNames'              Cell array containing dimension names, e.g.
%                               names = {'x';'y';'z';'time';'lambda'};
%                               Make sure that size(names,1) == ndims(data{1}).
%                               If omitted, dimensions will be named
%                               Dim1  , Dim2, Dim3, ...
%       'dimUnits'              Units of dimensions (simply for display purposes)
%       'matNames'              Names of matrices displayed in Figure Name.
%                               Needs to be specified as a cell array of
%                               strings, even if only one matrix is
%                               specified.
%                               Default is the name of the variable in the matlab workspace, but if an
%                               expression is used (such as d(:,:,1)), it will be empty leaving
%                               matVis to display an empty string.
%       'alphaMap'              Alphamap, has to be of same size as data
%                               matrix. Useful for masking images. Background color will
%                               be set to black.
%       'startPar'              List of configuration settings. These
%                               settings 'override' the custom settings saved in a customConfig
%                               file (if available). The list should be a cell array with the
%                               common {'propertyName1'; 'propertyValue1'; 'propertyName2';
%                               'propertyValue2';...} structure. The following properies can be
%                               set:
%                               xySel             x- and y-dim for image display
%                               zoomVal           Matrix indicating zoom setting
%                               plotDim           Dimensions along which plots will be displayed
%                               plotMean          Number indicating status of plot-average-mode: 0: 1x1, 1: 3x3, 2: 5x5, 3: [1x1,3x3,5x5], 4: zoom area, 5: RGB
%                               currPos           Vector indicating starting position
%                               cmMinMax          Matrix of size 2 x number of data sets for colormap limits
%                               cmap              Number of colormap from list of available colormaps:  {'Gray';'Gray (Range)'; 'Jet'; 'HSV'; 'Hot'; 'Cool';
%                                                   'Red 1';'Red 2';'Green 1';'Green 2';'Blue 1';'Blue 2'; 'Rainbow1';'Rainbow2';'Rainbow3';'Rainbow4';
%                                                   'Blue-Gray-Yellow (0 centered)';'Magenta-Gray-Green (0 centered)'}
%                               aspRatio          Two element vector
%                               rgbDim            Number of RGB dimensions
%                               projDim           Number of projection dimensions
%                               projMethod        Number or string: {'None';'Max';'Min';'Mean';'Std';'Var'; 'Tile'}
%                               windowVisibility  Binary vector indicating visibility of [imageWin zoomWin plotWin]
%                               windowPositions   Structure s containing the fields s.gui, s.imageWin, s.zoomWin and s.plotWin. In case there is one data set, 
%                                                   the values of the fields are four element vectors [x0, y0, width, height], in case there are nMat data sets, 
%                                                   the values are matrices of size nMat x 4.
%                                                   WARNING: gui Position includs TooltipDisplay!!!
%
% Output Argument (optional)
%    Structure containing the following fields:
%       - figHandles: figure handles to all figures, including main GUI
%       - fctHandles: function handles of functions useful for 'data
%         exchange' or 'external' update of matVis settings
%       - settings: structure containing information about current settings
%         (currently zoom and zoomXY)
%       - dataName: cell array containing names of the loaded data
%
% Mouse actions
%   Left click and drag                      Draw zoom region
%   Right click                              Unzoom
%   Left click and drag in Zoom Rectangle    Move zoom region (pan)
%   Middle click or Shift + left click       Move Position Lines
%   Double click                             Copy content of current figure
%                                            to clipboard (either as bitmap or vector graphics).
%                                            Works only for Windows OS.
%   Right click in main gui                  Bring all visibile windows on
%                                            top of the screen
%   Scroll Wheel (Matlab 2007a or later)     Zoom in and out (for zoom and
%                                            image windows)
%
% Keyboard
%   Instead of using sliders or text boxes to change the current position values,
%   you can use the number keys on the keyboard instead. Press '1' to increase
%   the current value of Dim1 by one, press 'control'+'1' to decrease its
%   value by one,... This is however only possible if the main gui is the
%   current figure (i.e. the selected window) and none of its controls is active
%   (click in some empty space within the gui if it doesn't work as you
%   expect).
%   Works also only for the first nine dimensions :(
%
% RGB Display
%   RGB display can be toggled with the RGB button. You can use RGB display
%   for each dimension that is not used as either 'x' or 'y' dimension. As
%   you press the RGB button you switch between all possible dimensions.
%   The current image is displayed as the green channel, the preceding as
%   red and the succeeding as the blue channel. If there are only two
%   values in the selected dimension the blue channel will be left empty.
%   If you select the 'Stretch RGB' option, a colormap along the complete
%   dimension will be created such that red corresponds to high intensities
%   in the initial part, green to high intensities in the intermediate part,
%   and blue to high intensities in the late part of the selected
%   dimension. Changing the position of the slider of this dimension does
%   not have any effect.
%
% Note that some functions (e.g. histogram, export data) are so
% far only supported for the first data set! Using them while multiple data
% sets are loaded might lead to errors!
%
% See end of this file for a list of known bugs and planned feature
% implementations.
% See http://www.colors-and-contrasts.com/Documents/matVisGuide.pdf for a
% complete manual (can also be accessed from the main Gui of matVis).
% Manual is not completely up to date. Please contact me if you have 
% questions / bug reports / feature requests!
%
%**************************************************************************
%               Copyright 13.06.2012, S. Junek (sjunek@gwdg.de)
%
versionNumber = 1.102;  % Current version number of matVis
%% Check Matlab version and installed toolboxes
v = version;
% v = num2str(v(1:3));
% if v < 7
%     errordlg('matVis is not supported for MATLAB versions prior v7 (R14). Sorry!','matVis not supported');
%     return;
% end
oldMATLAB = verLessThan('MATLAB','8.4'); % to differenciate for HG2 
mlVer = ver;
withDipimage = any(strfind([mlVer.Name],'DIPimage'));
withImageProcessingTB = any(strfind([mlVer.Name],'Image Processing Toolbox'));
withStatisticsTB = any(strfind([mlVer.Name],'Statistics Toolbox'));  %#ok
% Parameters that might be overwritten by optional arguments
customDimScale = 0;  % Flag to indicate whether dimension scales are provided
dimNames = [];       % Dimensions names
withDimUnits = 0;    % Dimension units
fromFile = 0;        % Flag indicating whether data are loaded from file
os = computer;       % Operating system
macScaleFactor = [1.05 1.05]; % Scaling factor to adjust width and height of GUI and uicontrols for Mac OS-X
%% Read Data ...
% ... from Files
if nargin == 0 || ischar(varargin{1})
  debugMatVis = 0;
  if nargin && strcmp(varargin{1},'debug') && (islogical(varargin{2}) || isscalar(varargin{2}))
    debugMatVis = varargin{2};
  end
    if nargin == 0 || strcmp(varargin{1},'debug')
        %Choose files
        [f,p] = uigetfile( '*.mat;*.tif;*.tiff; *.jpg; *.bmp; *.gif; *.png; *.lsm; *.da; *.dat',...
            'Select one or more files of identical dimension!', 'MultiSelect', 'on');
        if isequal(f, 0)
            varargout{1} = [];
            return;
        end
    else
        f = varargin{1};
        p = '';
    end
    allNames = '';
    fromFile = 1;
    %Convert to array if there is only one data set
    if ~isa(f, 'cell')
        ff{1} = f;
        f = ff;
    end
    if numel(f) > 2
        q = questdlg('Sort files alphabetically?','Sort files','Sort in ascending order','Sort in descending order','Keep selection order','Sort in ascending order');
        switch q
            case 'Sort in ascending order'
                f = sort(f);
            case 'Sort in descending order'
                f = sort(f,1,'descend');
        end
    end
    %Load files
    isCustomTif = 0;
    for i = 1:size(f,2)
        [pn, fn, ext] = fileparts(f{i}); %#ok
        display(['Loading file ', f{i}]);
        if strcmp(ext, '.tif') || strcmp(ext, '.tiff')
            %Try to read as custom tif
            w = imfinfo([p f{i}]);
            if any(w(1).StripOffsets == 832)
                [data{i}, par{i}] = readCustomTif([p,f{i}]); %#ok
                isCustomTif = 1;
                if i==1
                    filePath = p;
                end
                %                 if ~isequal(par{i}.dim, par{1}.dim)
                %                     display('Error: Data sets need to have identical dimensions!');
                %                     return;
                %                 end
                try
                    assignin('base', ['tifPar_',fn], par{i});
                catch    %#ok
                    assignin('base', ['tifPar_',num2str(i)], par{i});
                end
                display('File description: ');
                display(par{i}.description);
                defaultColormap{i} = []; %#ok
                % Data scaling and names, assume it's the same for all files
                if i==1
                    dimScale(1,:) = par{i}.scale.x*[1 size(data{i},1)];
                    dimScale(2,:) = par{i}.scale.y*[1 size(data{i},2)];
                    dimScale(3,:) = abs(par{i}.scale.z)*[1 size(data{i},3)];
                    dimScale(4,:) = (par{i}.config.aq.frt+par{i}.config.aq.frdt)*[1 size(data{i},4)];
%                     if numel(size(data{i})) > 4
                        dimScale(5,:) = [1 size(data{i},5)];
%                     end
                    customDimScale = 1;
                    dimNames{1} = 'x';
                    dimNames{2} = 'y';
                    dimNames{3} = 'z';
                    dimNames{4} = 't';
                    dimNames{5} = 'l';
                    dimUnits{1} = 'µm';
                    dimUnits{2} = 'µm';
                    dimUnits{3} = 'µm';
                    dimUnits{4} = 's';
                    dimUnits{5} = 'channel';
                    withDimUnits = 1;
                    sz = ones(1,5);
                    sz(1:numel(size(data{1}))) = size(data{1});
                    dimScale(sz==1,:) = [];
                    dimUnits(sz==1) = [];
                    dimNames(sz==1) = [];
                    sz(sz==1) = [];
%                     w = size(data{1});
%                     w(w==1) = [];
%                     dimScale = dimScale(1:numel(w),:);
%                     dimUnits = dimUnits(1:numel(w));
                    pxWidth = diff(dimScale,1,2)'./(sz-1);
                end
                data{i} = squeeze(data{i});
            end
        end
        %Read other formats
        if ~isCustomTif
            %(Animated) Gif
            if strcmp(ext, '.mat')
                ww = load([p,f{i}]);
                ww = struct2cell(ww);
                data{i} = ww{1};  %#ok
                clear ww;
                defaultColormap{i} = [];  %#ok
            elseif strcmp(ext, '.gif')
                [data{i},defaultColormap{i}] = imread([p,f{i}], 'frames', 'all');  %#ok
                data{i} = squeeze(data{i});  %#ok
                %Tif (possibly multi-image)
            elseif strcmp(ext, '.tif') || strcmp(ext, '.tiff')
                tic
                % Speed optimized version, not tested on color tifs
                ww = imfinfo([p,f{i}]);  % Tif file information
                nCol = ww(1).SamplesPerPixel;  % Number of color channels (NOT TESTED!)
                nImg = numel(ww);  % Number of frames
                if ww(1).BitsPerSample(1) <= 16
                    out = zeros(ww(1).Height, ww(1).Width,nCol, nImg,'uint16');
                else
                    out = zeros(ww(1).Height, ww(1).Width, nCol, nImg);
                end
                readLength = [ww(1).Width ww(1).Height]; % order has to be reversed for correct reading of data (don't know why)
                gl_fid = fopen ([p,f{i}], 'r', 'l');
                for j=1:nImg
                    for colIdx=1:nCol
                        fseek(gl_fid, ww(j).StripOffsets(colIdx) , 'bof');
                        if ww(j).BitsPerSample(colIdx) <=8
                            out(:,:,colIdx,j) = uint16(fread(gl_fid,readLength, '*uint8'))';
                        elseif ww(j).BitsPerSample(colIdx) <=16
                            out(:,:,colIdx,j) = fread(gl_fid, readLength, '*uint16')';
                        else ww(j).BitsPerSample(colIdx)
                            out(:,:,colIdx,j) = fread(gl_fid,readLength, 'uint32')';
                        end
                    end
                end
                fclose(gl_fid);
                data{i} = squeeze(out);
                clear out;
                toc
                % Old (slow) version:
%                 ww = imfinfo([p,f{i}]);
%                 im1 = imread([p f{i}],1);
%                 switch ww(1).BitDepth
%                     case 8
%                         d{i} = zeros([size(im1),numel(ww)],'uint8');  %#ok
%                     case 16
%                         d{i} = zeros([size(im1),numel(ww)],'uint16'); %#ok
%                 end
%                 if strcmp(ww(1).PhotometricInterpretation, 'RGB')
%                     data{i} = imread([p,f{i}]);  %#ok
%                 else
%                     for j = 1:numel(ww)
%                         data{i}(:,:,j) = imread([p,f{i}],j,'Info',ww);  %#ok
%                     end
%                 end
                defaultColormap{i} = []; %#ok
            elseif strcmp(ext, '.lsm')
                dimNames{1} = 'x';
                dimNames{2} = 'y';
                dimNames{4} = 'l';
                dimNames{3} = 'z';
                dimNames{5} = 't';
                if i==1  % Data size and scaling, assume it's the same for all files
                    if any(which('lsminfo'))
                        ww = lsminfo([p,f{i}]);
                        lsmDim = ww.DIMENSIONS;
                        lsmDim(4) = ww.NUMBER_OF_CHANNELS;
                        lsmDim(5) = ww.TIMESTACKSIZE;
                        lsmDim(lsmDim == 0) = 1;
                        dimScale(1,:) = ww.VoxelSizeX*1e6*[1 lsmDim(1)];
                        dimScale(2,:) = ww.VoxelSizeY*1e6*[1 lsmDim(2)];
                        dimScale(4,:) = [1 ww.NUMBER_OF_CHANNELS];
                        dimScale(3,:) = ww.VoxelSizeZ*1e6*[1 lsmDim(3)];
                        dimScale(5,:) = ww.TIMESTACKSIZE*[1 lsmDim(5)];
                        dimUnits{1} = 'µm';
                        dimUnits{2} = 'µm';
                        dimUnits{4} = 'channel';
                        dimUnits{3} = 'µm';
                        dimUnits{5} = 's';
                        customDimScale = 1;
                        withDimUnits = 1;
                        % remove singular dimensions
                        dimScale(lsmDim==1,:) = [];
                        dimUnits(lsmDim==1) = [];
                        dimNames(lsmDim==1) = [];
                        pxWidth = diff(dimScale,1,2)'./(lsmDim(lsmDim~=1)-1);
                        if iscell(ww(1).ScanInfo.BITS_PER_SAMPLE)
                            bd =  ww(1).ScanInfo.BITS_PER_SAMPLE{1};
                        else
                            bd =  ww(1).ScanInfo.BITS_PER_SAMPLE(1);
                        end
                    else
                        lsmDim = str2double(inputdlg({'x';'y';'z';'z';'t';'BitDepth'},'Enter data size for all 5 dimensions and bit depth'));
                        bd = lsmDim(6);
                        lsmDim(6) = [];
                        dimNames(lsmDim==1) = [];
                    end
                end
                ww = imfinfo([p,f{i}]);
                nCol = ww(1).SamplesPerPixel;
                nImg = floor(numel(ww)/2); % Every second image is a thumbnail of the preceding full-res image
                if ww(1).BitsPerSample(1) <= 16
                    data{i} = zeros(ww(1).Height, ww(1).Width,nCol, nImg,'uint16');
                else
                    data{i} = zeros(ww(1).Height, ww(1).Width, nCol, nImg);
                end
                readLength = [ww(1).Width ww(1).Height]; % order has to be reversed for correct reading of data (don't know why)
                gl_fid = fopen ([p,f{i}], 'r', 'l');
                for j=1:nImg
                    for col=1:nCol
                        fseek(gl_fid, ww(2*j-1).StripOffsets(col) , 'bof');
                        if ww(2*j-1).BitsPerSample(col) <=8
                            data{i}(:,:,col,j) = uint16(fread(gl_fid,readLength, '*uint8'))';
                        elseif ww(2*j-1).BitsPerSample(col) <=16
                            data{i}(:,:,col,j) = fread(gl_fid, readLength, '*uint16')';
                        else ww(2*j-1).BitsPerSample(col)
                            data{i}(:,:,col,j) = fread(gl_fid,readLength, 'uint32')';
                        end
                    end
                end
                fclose(gl_fid);
                data{i} = permute(reshape(permute(data{i},[2 1 3 4 5]), lsmDim([1 2 4 3 5])),[1 2 4 3 5]);
                data{i} = squeeze(data{i}); %  in case  of a single color channel
                defaultColormap{i} = [];             %#ok
            elseif strcmp(ext, '.da')
                data{i} = readNeuroPlex([p,f{i}]); %#ok
                defaultColormap{i} = [];             %#ok
            elseif strcmp(ext, '.dat')
                data{i} = importdata([p,f{i}]); %#ok
                defaultColormap{i} = [];             %#ok
                %Everything else
            else
                [data{i},defaultColormap{i}] = imread([p,f{i}]);  %#ok
            end
            if ~isequal(size(data{1}), size(data{i}))
                display('Error: Data sets need to have identical dimensions!');
                return;
            end
        end
        %Check Input Dimensions
        varName{i} = f{i};        %#ok %Name of file
        allNames = [allNames, f{i}, '; ']; %#ok
    end
    allNames(end-1:end) = [];
    if numel(data) > 1
        q = questdlg('Combine files into single data set or keep as individual sets?','Combine files','Keep separate','Combine into single set','Keep separate');
        if strcmp(q, 'Combine into single set')
            % Reduce name of combined data set using name of first and last
            % file
            if numel(data) > 2
                allNames = [varName{1} ' ... ' varName{end}];
            end
            if customDimScale
               pxWidth = diff(dimScale,1,2)'./(size(data{1})-1);
               dimScale(end+1,:) = [1 size(dimScale,1)+1];
            end
            if withDimUnits
                dimUnits{end+1} = '';
            end
            varName = [];
            varName{1} = allNames;
            % Create size vector and dimension index for reshape data after
            %  converting cell array into matrix. The cells are
            %  concatenated along the 2nd dimension, so the matrix has to
            %  be first reshaped to separate the two dimensions and then
            %  permuted to put the file-number in the last dimension.
            dim = size(data{1});
            if numel(dim) > 2
                dim(4:end+1) = dim(3:end);
                dim(3) = numel(data);
            else
                dim(3) = numel(data);
            end
            ind = 1:numel(dim);
            ind(3) = ind(end);
            ind(end) = 3;
            dd = data;
            data = [];
            dd = permute(reshape(cell2mat(dd), dim),ind);
            data{1} = dd;
            clear dd;
            if ~isempty(dimNames)
               dimNames{end+1} = 'Files'; 
            end
            if customDimScale
                pxWidth(end+1) = 1;
            end
        end
    end
    nMat = numel(data);
    if isempty(dimNames)
        for i = 1:ndims(data{1})
            dimNames{i} = ['Dim', num2str(i)]; %#ok
        end
    end
    withAlpha = 0;
    startPar = [];
    numberOptionalArgs = 0;
    % ... from Arguments
else
    %Check for optional arguments, identified by string
    debugMatVis = 0;
    for i=1:nargin
        optionalArgs(i) = isa(varargin{i}, 'char');  %#ok
    end
    numberOptionalArgs = sum(optionalArgs);
    optionalArgs = find(optionalArgs);
    %Determine number of data matrices
    nMat = nargin - 2*numberOptionalArgs;
    %Read Data and create long string containing all variable names
    allNames = '';
    for i = 1:nMat
        if isa(varargin{i}, 'dip_image')
            data{i} = double(abs(varargin{i}));  %#ok
        else
            data{i} = varargin{i};  %#ok
        end
        varName{i} = inputname(i);   %#ok %Name of input variable
        allNames = [allNames, varName{i}, '; '];  %#ok
    end
    allNames(end-1:end) = [];
    %Set flags/default values of possible arguments
    %AlphaMap
    withAlpha = 0;
    startPar = [];
    %DimensionNames
    for i = 1:ndims(varargin{1})
        dimNames{i} = ['Dim', num2str(i)];  %#ok
    end
    %Read optional arguments
    for i = 1:numberOptionalArgs
        optionalArgIdentifier{i} = varargin{optionalArgs(i)};
        val = varargin{optionalArgs(i)+1};
        switch  optionalArgIdentifier{i}
            case 'alphaMap'
                if nMat == 1
                    alphaMap{1} = squeeze(val);
                elseif isnumeric(val)
                    for j = 1:nMat
                        alphaMap{j} = squeeze(val);  %#ok
                    end
                else
                    for j = 1:nMat
                        alphaMap{j} = squeeze(val{j});  %#ok
                    end
                end
                withAlpha = 1;
                currAlphaMap = [];
            case 'dimNames'
                if length(val) ~= ndims(varargin{1})
                    error(sprintf('Dimension of matrix and number of ''dimNames'' have to be equal!\nUse <a href="matlab:help matVis">help matVis</a> for help.'));
                end
                dimNames = val;
            case 'dimUnits'
                if length(val) ~= ndims(varargin{1})
                    error(sprintf('Dimension of ''dimUnits'' mut fit matrix dimension!\nUse <a href="matlab:help matVis">help matVis</a> for help.'));
                end
                dimUnits = val;
                withDimUnits = 1;
            case 'matNames'
                if length(val) ~= length(data)
                    error(sprintf('Dimension of ''matNames'' mut fit number of provided data sets!\nUse <a href="matlab:help matVis">help matVis</a> for help.'));
                end
                varName = val;
                allNames = '';
                for j=1:length(varName)
                    allNames = [allNames, varName{j}, '; '];  %#ok
                end
                allNames(end-1:end) = [];
            case 'startPar'
                startPar = val;
            case 'dimScale'
                dimScale = val;
                customDimScale = 1;
                pxWidth = diff(dimScale,1,2)'./(size(data{1})-1);
            case 'debug'
                debugMatVis = val;
        end
    end
    %Check Input Dimensions
    defaultColormap{1} = [];
    isCustomTif = 0;
end
% start debug modus recognition
fctCount = [];                           %Function Counter for debug mode
if debugMatVis, debugMatVisFcn(1); end

%Check for 1D data or singular dimensions
if size(size(data{1}),2) == 2 && (size(data{1},1)==1 || size(data{1},2)==1)
    disp('Error: matVis not suitable for (1D) vectors. Ever heard of the function ''plot''?');
    figure; plot(data{1});
    return
elseif  sum(size(data{1}) == 1) > 0
    dimNames(size(data{1}) == 1) = [];
    for i=1:nMat
        data{i} = squeeze(data{i});  %#ok
    end
end
%% Set initial values
dim = size(data{1});                     %Dimensions of data set
nDim = size(dim,2);                      %Number of dimensions
if customDimScale
    pxWidth(pxWidth==0) = 1;  % Avoid error messages
end
if ~withDimUnits
    dimUnits = cell(1,nDim);
end
maxVal = zeros(1,nMat*(1+withAlpha));
minVal = zeros(1,nMat*(1+withAlpha));
maxValInd = zeros(1,nMat*(1+withAlpha));
minValInd = zeros(1,nMat*(1+withAlpha));
cmMinMax = zeros(nMat*(1+withAlpha),2);
for i=1:nMat
    %         d = data{i};
    [maxVal(i) maxValInd(i)] = max(data{i}(:));         %Maximum Data Value
    [minVal(i) minValInd(i)] = min(data{i}(:));         %Minimum Data Value
    sldLimVal(i,:) = [minVal(i) maxVal(i)];
    
    if minVal(i) == maxVal(i)
        minVal(i) = maxVal(i)-1; % to avoid error messages
        if isempty(varName{i})
            display(['Warning: All numbers in the data set #' num2str(i) ' equal ' num2str(minVal(i)) '.']);
        else
            display(['Warning: All numbers in the data set ''' varName{i} ''' equal ' num2str(minVal(i)) '.']);
        end
    end
    cmMinMax(i,:) = [minVal(i) maxVal(i)]'; %Colormap Limits
end
if withAlpha
    for i=1:nMat
        [maxVal(nMat+i) maxValInd(nMat+i)] = max(alphaMap{i}(:));         %Maximum Data Value    % Used to be nanmax
        [minVal(nMat+i) minValInd(nMat+i)] = min(alphaMap{i}(:));         %Minimum Data Value   % Used to be nanmin
        sldLimVal(nMat+i,:) = [minVal(nMat+i) maxVal(nMat+i)];
        if minVal(nMat+i) == maxVal(nMat+i)
            minVal(nMat+i) = maxVal(nMat+i)-1; % to avoid error messages
            if isempty(varName{i})
                display(['Warning: All numbers in alpha map #' num2str(i) ' equal ' num2str(maxVal(nMat+i)) '.']);
            else
                display(['Warning: All numbers in alpha map to data set ''' varName{i} ''' equal ' num2str(maxVal(nMat+i)) '.']);
            end
        end
    end
    if all(maxVal(nMat+1:end)<=1) && all(minVal(nMat+1:end)>=0) % If all alpha values are between 0 and 1, use [0,1] as alpha limits
        cmMinMax(nMat+1:2*nMat,1) = 0; %Colormap Limits
        cmMinMax(nMat+1:2*nMat,2) = 1; %Colormap Limits
    else
        cmMinMax(nMat+i,:) = [minVal(nMat+1:2*nMat); maxVal(nMat+1:2*nMat)]'; %Colormap Limits
    end
end

clear d;
%Set initial plot limits
for i=1:nDim
    plotXLim(i,:) = [1 dim(i)]; %#ok
end
if customDimScale
    plotXLimScale = dimScale;
end
for i=1:nMat
    plotYLim(i,:) = [minVal(i) maxVal(i)]; %#ok
end
%Find position of maximum value (adopted from ind2sub function) and set as
%starting position; if data set is too large to find this position set
%center of data set as starting position instead
try
    maxPosLin = maxValInd(1+withAlpha);
    minPosLin = minValInd(1+withAlpha);
    k = [1 cumprod(dim(1:end-1))];
    for i = nDim:-1:1
        viMax = rem(maxPosLin-1  , k(i)) + 1;
        vjMax = (maxPosLin - viMax)/k(i) + 1;
        maxPos(i) = vjMax; %#ok
        maxPosLin = viMax;
        viMin = rem(minPosLin-1  , k(i)) + 1;
        vjMin = (minPosLin - viMin)/k(i) + 1;
        minPos(i) = vjMin; %#ok
        minPosLin = viMin;
    end
    currPos = maxPos;
catch    %#ok
    currPos = round(dim/2);                 %Starting Position
    maxPos = currPos;
    minPos = currPos;
end
savedPos = [];                           %Matrix for saved positions using small buttons to the right
savedZoom = [];                          %Matrix for saved zoom using small buttons to the right

% peallocation of handles
if oldMATLAB
  currIm = [];                             %Handle of Figure of Current Image
  imHandle = [];                           %Handle of Figure in Image Window
  zoomHandle = [];                         %Handle of Figure in Zoom Window
  tifParFig = [];                          %Handle to Figure displaying CustomTif parameter
  imAx = [];                               %Handle of Axes in Image Window
  zoomAx = [];                             %Handle of Axes in Zoom Window
  zoomReg = [];                            %Handle of Rectangle indicating Zoom Region
  tempWin = [];                            %Temp. window used for user input (aspect ratio, axes limits)
  exportWin = [];                          %Temp. window used for data export
  subPlotHandles = [];                     %Handles to subplot axes in Plots Window
  subPlotPlots = [];                       %Handles to data displayed in Plots Window
  subPlotPlotsAlpha = [];                  %Handles to data displayed in Plots Window
  posLine = [];                            %Lines in plots indicating current value of diplayed dimension
  lineHorZoom = [];                        %Handles to position lines in Image and Zoom Windows
  lineVertZoom = [];
  lineHorIm = [];
  lineVertIm = [];
  rectPlotArea_Zoom = [];
  rectPlotArea_Im = [];
  histPlots = [];                          %Handle for Histogram Plot
  histObj = [];                            %Handle for Histogram Objects (lines)

  roiWin = [];                             %Handle to Roi Manager Window
  roiImage = [];                           %Handle to image displaying current Roi
  roiAxes = [];                            %Handle to axes displaying current Roi
  profileWin = [];                         %Handle to Profile Manager Window
  profileImage = [];                           %Handle to image displaying current profile
  profileTraceWin = [];
  profileAxes = [];                            %Handle to axes displaying current profile
  tempRoiPropWin = [];                         %Handle of Figure for ROI update of ROI settings
else
  % if the array of handles is defined as empty matrix (double) some functions (like colorbar) do not work
  % correct predefinition is only required for arrays of objects!!!
  % so most predefinitions are not required
  
%   currIm = matlab.ui.Figure.empty;                             %Handle of Figure of Current Image
  imHandle = matlab.ui.Figure.empty;                           %Handle of Figure in Image Window
  zoomHandle = matlab.ui.Figure.empty;                         %Handle of Figure in Zoom Window
  tifParFig = matlab.ui.Figure.empty;                          %Handle to Figure displaying CustomTif parameter
  imAx = matlab.graphics.axis.Axes.empty;                               %Handle of Axes in Image Window
  zoomAx = matlab.graphics.axis.Axes.empty;                             %Handle of Axes in Zoom Window
  zoomReg = [];                            %Handle of Rectangle indicating Zoom Region
  tempWin = matlab.ui.Figure.empty;                            %Temp. window used for user input (aspect ratio, axes limits)
  % tempWin = []; is used later in matVis
  exportWin = matlab.ui.Figure.empty;                          %Temp. window used for data export
  subPlotHandles = matlab.graphics.chart.primitive.Line.empty;              %Handles to subplot axes in Plots Window
  subPlotPlots = matlab.graphics.chart.primitive.Line.empty;                %Handles to data displayed in Plots Window
  subPlotPlotsAlpha = matlab.graphics.chart.primitive.Line.empty;           %Handles to data displayed in Plots Window
  posLine = matlab.graphics.chart.primitive.Line.empty;                     %Lines in plots indicating current value of diplayed dimension
  lineHorZoom = matlab.graphics.chart.primitive.Line.empty;                 %Handles to position lines in Image and Zoom Windows
  lineVertZoom = matlab.graphics.chart.primitive.Line.empty;
  lineHorIm = matlab.graphics.chart.primitive.Line.empty;
  lineVertIm = matlab.graphics.chart.primitive.Line.empty;
  rectPlotArea_Zoom = matlab.graphics.primitive.Rectangle.empty;
  rectPlotArea_Im = matlab.graphics.primitive.Rectangle.empty;
  histPlots = matlab.ui.Figure.empty;                                       %Handle for Histogram Plot
  histObj = matlab.graphics.chart.primitive.Line.empty;                     %Handle for Histogram Objects (lines)

  roiWin = matlab.ui.Figure.empty;                                          %Handle to Roi Manager Window
  roiImage = matlab.graphics.primitive.Image;                               %Handle to image displaying current Roi
  roiAxes = matlab.graphics.axis.Axes.empty;                                %Handle to axes displaying current Roi
  profileWin = matlab.ui.Figure.empty;                                      %Handle to Roi Manager Window
  profileImage = matlab.graphics.primitive.Image;                           %Handle to image displaying current profile
  profileTraceWin = matlab.ui.Figure.empty;               
  profileAxes = matlab.graphics.axis.Axes.empty;                            %Handle to axes displaying current profile
  tempRoiPropWin = matlab.ui.Figure.empty;                                  %Handle of Figure for ROI update of ROI settings
end
  
currImVal = [];                          %Values of current image (currIm might contain altered values due to filtering and/or gamma values ~= 1)
currAlphaMapVal = [];                    %Values of current alpha map (currAlphaMap might contain altered values due to filtering and/or  gamma values ~= 1)
currContrastSel = 1;                     %Dataset currently selected for contrast selection / histogram display
linkContrastSettings = 1;                %Manual contrast adjustment: adjust contrast of all data sets together in a linear way
xySel = [1 2];                           %Selected Dimensions for Images
zoomVal = zeros(nDim,2);
for i=1:nDim
    zoomVal(i,1) = 1;
    zoomVal(i,2) = dim(i);
end
zoomValXY = [1,1,dim(xySel(2)),dim(xySel(1))]; %Zoom Value [left down width height]
if customDimScale
    zoomPxXY = [dimScale(2,1) dimScale(1,1) dimScale(2,2)-dimScale(2,1) dimScale(1,2)-dimScale(1,1)];
end
zoomStart = [];                          %Zoom Value for Paning Zoom Window
aspRatio = [1 1];                        %Aspect Ratio for image and zoom display
pStart = [];                             %Start Position for Paning Zoom Window
plotSel = zeros(1  ,nDim);                 %Index of dimensions for selected plots (will be set in gui configuration)
nPlots = 0;                              %Number of plots to be displayed
plotDim = [];                            %Dimensionnumbers of selected plots
rgbDim = [];                             %Dimension used for RGB display
rgbCount = 0;                            %Count for flipping through RGB dimension (zero means no RGB display)
findExtremumFlag = 0;                    % Flag indicating whether at each change of position the cursor (i.e. selected position) should automatically jump to the min or max
findExtrMnMx = '';
findExtrZmIm = '';
%rgbLimits = [];                          %Limits for RGB channels in channel display mode
rgbVal = [];                             %rgb map for stretched RGB mode
rgbValPlots = [];                        %rgb map for stretched RGB mode used for plots
colorcodePlots = [];                     %colorcode for plots: used for plots based on ROIs
rgbStretchSldVal  = [0 1];               %Slider value for RGB stretch mode
cmap = gray(255);                        %Current colormap
plotValues = [];                         %Values used in plots
plotValuesAlpha = [];                    %AlphaValues used in plots
isPlaying = 0;                           %Playing status
playDim = [];                            %Dimension which was selected for play mode
projMethod = 0;                          %Number indicating method for projection.
%0 - no projection, 1 - maximum
%projection, 2 - mean projection,...
projDim = 3;
currWinScale = 100;                      %Window-pixel scale for 100% button
if withAlpha
    currGamma = ones(1,2*nMat);
else
    currGamma = ones(1,nMat);
end
histVal = [];                            %Histogram Values
guiHistVal = [];                         %Histrogram Values for GUI histogram
dragging =  0;           %Force update of GUI hist values during
histValCurrIm = [];                      %Histogram Values of current image
globalHist = [];                         %Values for Global Histogram (over complete (first) data set)
calculatingGlobalHist = 0;               % Flag indicating that global hist is currently calculated (to avoid errors during calculation)
%dataPerc = [];                          %Percentiles of data, filled after updateGlobalHist is called for the first time
%Function not yet implemented
exportCount = 0;                         %Count of exported data subsets
isBusy = 0;                              %Flag for busy-state of matVis
shuttingDown = 0;                        %Flag for shutting matVis down
currRoi = [];                            %Current Roi
roiMin = [];                             %Min of current Roi
roiMax = [];                             %Max of current Roi
roiMean = [];                            %Mean of current Roi
roiSize = [];                            %Number of pixels of Roi
roiList = [];                            %List of Rois
roiLine = [];                            %Handle to lines indicating Rois
roiCenterIndicator = [];                 %Handle to '+' (text) indicating center of current Roi
nRois = 0;                               %Number of Rois
roiText = [];                            %Handle to text (numbers) for Rois
roiName = [];                            %Roi name of currently selected roi above roiAxes
roiListbox = [];                         %Listbox in Roi Gui containing all ROIs; selected ROI is "current Roi"
roiBtRename = [];                        %Button to rename roi
roiBtShowNames = [];                     %Button to show or hide Roi text in figure windows
roiBtReplace = [];                       %Button to exchange existing Roi by new one.
roiPopDataExport = [];                   %Handle to pop menu to choose dimension for Roi based data export
tbRoiShift = [];                         %Handle for toggle button for shifting of Rois
tbRoiRotate = [];                        %Handle for toggle button for rotating Rois
tbRoiScale = [];                         %Handle for toggle button for scaling Rois
tb_newRoi = [];                          %Handle for toggle buttons to creat new rois
bt_copyRoi = [];
bt_deleteRoi = [];
bt_roiExport = [];
bt_exportRoiData = [];
newRoiSize = 5;                          %Size of "one-click" ROIs (radius, squares will have 2*newRoiSize+1 side length)
jump2ROIPos_cb = 1;                      %Handle to checkbox indidacting whether the current position chould be changed when a new ROI is selected (jump to position where ROI was created)
etxt_newRoiName = [];
txt_RoiPropWinTxt = [];                  %Handle Textfield in Figure for ROI update of ROI settings

profilePoints = [];
profileInwork = 0;
profileComplete = 0;
currProfile = [];                            %Current profile
profileMin = [];                             %Min of current profile
profileMax = [];                             %Max of current profile
profileMean = [];                            %Mean of current profile
profileSize = [];                            %Number of pixels of profile
profileStruct = struct;                      %Structure used for creating and updating profiles ("finisehed" profiles are then collected in profileList)
profileList = [];                            %List of profiles
profileLine = [];                            %Handle to lines indicating profiles
nProfiles = 0;                               %Number of profiles
profileText = [];                            %Handle to text (numbers) for profiles
profileName = [];                            %profile name of currently selected profile above profileAxes
profileListbox = [];                         %Listbox in profile Gui containing all profiles; selected profile is "current profile"
profileBtRename = [];                        %Button to rename profile
profileBtShowNames = [];                     %Button to show or hide profile text in figure windows
profileBtReplace = [];                       %Button to exchange existing profile by new one.
profilePopDataExport = [];                   %Handle to pop menu to choose dimension for profile based data export
tbProfileShift = [];                         %Handle for toggle button for shifting of profiles
tbProfileRotate = [];                        %Handle for toggle button for rotating profiles
tbProfileScale = [];                         %Handle for toggle button for scaling profiles
tb_newProfile = [];                          %Handle for toggle buttons to creat new profiles
bt_deleteProfile = [];
bt_profileExport = [];
bt_exportProfileData = [];
isMatVis2DHist = 0;                        % Flag to indicate whether current matVis is called to display 2D hist
calledFrom2DHist = 0;                   % Flag indicating whether matVis update was called as 2D hist (to avoid recrusive callbacks between main matVis and 3D Hist matVis)
matVis2DHist = [];                      % Handle to matVis that will/might be called as 2D hist-matVis
with2DHist = 0;                         % Flag to indicate whether current matVis has associated 2D hist matVis open
x_2DHist = [];                          % x vector of 2D histogram data
y_2DHist = [];                          % y vector of 2D histogram data
% fctLevel = 0;                            %Function level for debug mode
% fctCount = [];                           %Function Counter for debug mode
% debugIndent = 4;                         %Indentation width for debug mode output
pop_profileProjection = [];              % Placeholder for handle for projection selection  in profile GUI (will be newly defined, i.e. no difference between "new" and "old" Matlab)
% Paramteres for movie recording
movdata.parts = [];
movdata.set.movname = 'movie_file_name'; % preset name for movie output
movdata.set.movdir = pwd; % preset dir for movie output
movdata.set.movcodec = 'MPEG-4'; % presetting for movie output
movdata.set.movFrameRate = 10; % presetting for movie FrameRate
movdata.set.movquality = 75; % presetting for movie quality (only for MPEG-4 and Motion JPEG AVI)
movdata.set.movcompress = 10; % presetting for compression rate (only Archival and Motion JPEG 2000)
movdata.set.movsize = [800 380]; % presetting for movie frame size
if ~withAlpha
    movdata.set.movbg = [0.99 0.92 0.8]; % presetting for movie background color
    movdata.set.movPartsbg = [0.73 0.83 0.96]; % presetting for movie parts color
    movdata.set.movTitleColor = [0 0 0]; % presetting for movie title color
    movdata.set.movLabelColor = [0 0 0]; % presetting for movie label color
else
    movdata.set.movbg = [0 0 0]; % presetting for movie background color
    movdata.set.movPartsbg = [0 0 0]; % presetting for movie parts color
    movdata.set.movTitleColor = [0.8 0 0]; % presetting for movie title color
    movdata.set.movLabelColor = [1 1 1 ]; % presetting for movie label color
end
movdata.set.prevscale = 1; % presetting for preview scaling (value of popup; 1:100%, 2:75%, 3:50%)
movdata.set.partunits = 1; % presetting for Movie Part Units (value of popup; 1:normalized, 2:pixels)
movdata.set.title = 'Movie Title'; % presetting for Movie Title
movdata.set.titlepos = [.5 .92]; % presetting for Movie Title Position
movdata.set.titlesize = 12; % presetting fot Movie Title FontSize in pt
movdata.set.titlestatus = 1; % presetting fot Movie Title FontSize in pt
movdata.rec = 0; % status of recordbutton
movdata.fastCaptureMode = true; %false; % use 'hardcopy' instead of 'export_fig'
%% Configuration of Windows
%Custom configuration from config file
%Window Properties

% Ensure root units are pixels and get the size of the screen
set(0,'Units','pixels')
scnSize = get(0,'ScreenSize'); %
monSize = get(0,'MonitorPosition');
monSizeMin = min(monSize(:,[1 2]),[],1);         % min x y position
monSizeMax = [max(sum(monSize(:,[1 3]),2)-1,[],1) max(sum(monSize(:,[2 4]),2)-1,[],1)]; % max x y position

% Create main gui already here to be able to find the correct width of the
% window borders (depend on OS)
% Calculation of MenuBar width does not work. For some reason, the figure
% size is not changed when showing/hiding the menu bar unless in debug mode

gui = figure('Menubar','none');
% borderWidthMenuBars = get(gui,'OuterPosition') - get(gui, 'Position');
set(gui,'Visible','off','Menubar','none');
borderWidth = get(gui,'OuterPosition') - get(gui, 'Position');
winWidthSide = abs(borderWidth(1));
winWidthTop  = borderWidth(4)-winWidthSide;
% winWidthMenuBar = borderWidthMenuBars(4) - borderWidth(4); %Will be used for toggling menu bars (see function toggleMenuBars)
winWidthMenuBar = 48;
dimHeight     = 40; %Height used for control elements of each dimension in the gui (in pixel)
guiSize = [320 nDim*dimHeight+460];
sz = ([...                                             %  add 'min' before bracket to enable square windows
    floor((scnSize(3)-(guiSize(1)+6*winWidthSide))/2)...      % (screen width - GUI) / 2 for Image and Zoom wim
    floor((scnSize(4)-2*(winWidthSide+winWidthTop))/2)]);   % screen hight / 2 first line: Image+zoom, second line: Plot win
winSize = sz;
%Window Visibility
defaultConfig.winVis.imageWin = 1;       %Default: 1
defaultConfig.winVis.zoomWin  = 1;       %Default: 1
defaultConfig.winVis.plotWin  = 1;       %Default: 1

%Window Position
%For one data set
defaultConfig.winPos.one.imageWin = ...
    [guiSize(1) + 3*winWidthSide, scnSize(4) - winSize(2) - winWidthTop,  winSize];       %Default: [337, 545,  450, 450];
defaultConfig.winPos.one.zoomWin  = ...
    [guiSize(1) + 5*winWidthSide + winSize(1), scnSize(4) - winSize(2) - winWidthTop,  winSize];       %Default: [800, 545,  450, 450];
defaultConfig.winPos.one.plotWin  = ...
    [guiSize(1) + 3*winWidthSide, scnSize(4)-2*(winSize(2)+ winWidthTop) - winWidthSide , 2*(winSize(1) + winWidthSide), winSize(2)];       %Default: [5,    10, 1250, 500];
%For multiple data sets
sz = ([...                                             % add 'min' before brackat to enable square windows
    floor((scnSize(3)-(guiSize(1)+8*winWidthSide))/4)...      % (screen width - GUI) / 2 for Image and Zoom wim
    floor((scnSize(4)-nMat*(winWidthSide+winWidthTop))/nMat)]);   % screen hight / 2 first line: Image+zoom, second line: Plot win
winSize = sz;

% winWidth = round(900 / nMat);
% winHeight = 200;
% plotHeight = round((540-nMat*30) / nMat);
for i =  1:nMat
    defaultConfig.winPos.mult.imageWin(i,:) = ... [337+(i-1)*(winWidth+8), 794, winWidth,  winHeight]; %Default:[337+(i-1)*(winWidth+8), 794, winWidth,  winHeight];
        [guiSize(1) + 3*winWidthSide, scnSize(4) - i*(winSize(2) + winWidthTop + winWidthSide) + winWidthSide,  winSize];
    defaultConfig.winPos.mult.zoomWin(i,:) = ... [337+(i-1)*(winWidth+8), 555, winWidth,  winHeight];  %Default:[337+(i-1)*(winWidth+8), 555, winWidth,  winHeight];
        [guiSize(1) + 5*winWidthSide + winSize(1), scnSize(4) - i*(winSize(2) + winWidthTop + winWidthSide) + winWidthSide,  winSize];
    defaultConfig.winPos.mult.plotWin(i,:) =  ... [5, 10+(i-1)*(plotHeight+30),     1250, plotHeight];   %Default:[5, 10+(i-1)*(plotHeight+30),     1250, plotHeight];
        [guiSize(1) + 7*winWidthSide + 2*winSize(1), scnSize(4) - i*(winSize(2) + winWidthTop + winWidthSide) + winWidthSide,  2*(winSize(1)), winSize(2)];
end
% defaultConfig.winPos.mult.gui = [5 815-nDim*20 320 nDim*20+180];       %Default: [5 815-nDim*20 320 nDim*20+180]);
% defaultConfig.winPos.one.gui  = [5 815-nDim*20 320 nDim*20+180];       %Default: [5 815-nDim*20 320 nDim*20+180]);
defaultConfig.winPos.mult.gui = [winWidthSide scnSize(4)-guiSize(2)-winWidthTop guiSize];       %Default: [5 815-nDim*20 320 nDim*20+180]);
defaultConfig.winPos.one.gui  = [winWidthSide scnSize(4)-guiSize(2)-winWidthTop guiSize];       %Default: [5 815-nDim*20 320 nDim*20+180]);

% Image and Zoom Window Options
%Aspect Ratio 1:1
defaultConfig.aspectRatio = 1;           %Default: 1
%Colorbar Display
defaultConfig.colorbar = 0;              %Default: 0
%Colormap
defaultConfig.colormap = 1;              %Default: 1 (gray)
%Gamma
defaultConfig.gamma = ones(1,nMat);
%RGB Mode
defaultConfig.RGB = 0;
%Colormap Mode (Global, Local or Manual)
defaultConfig.colormapMode = 'Global';   %Default: 'Global'
%Lines Visibility
defaultConfig.lineVis = 1;               %Default: 1
%Menu Bar
defaultConfig.menuBarVis = 0;            %Default: 0
%Update of histogram during playback
defaultConfig.playHist = 0;
%Update of histogram during change of position indicators
defaultConfig.moveHist = 0;
% Link figure position/size
defaultConfig.linkFigSize = 1;

% Plot Options
%Plot Dimensions
defaultConfig.plotDim = [1 2];           %Default: [1 2]
%Plot Zoom
defaultConfig.plotZoom = 0;              %Default: 0
%Scale Plot
defaultConfig.plotScale = 0;       %Default: 0
%Marker
defaultConfig.marker = 0;                %Default: 0
%Average
defaultConfig.plotMean = 0;              %Default: 0 (no averaging)
% Tooltips
defaultConfig.tooltips = 1;              %Default: 1 (display tooltips)

if nMat == 1
    defaultConfig.winPos = defaultConfig.winPos.one;
else
    defaultConfig.winPos = defaultConfig.winPos.mult;
end

[matVisPath,f,e] = fileparts(which('matVis.m'));
compName = strtrim(getenv('Computername'));
configFile = [];
if exist(fullfile(matVisPath,['matVisConfig_' compName '.mat']), 'file')
    configFile = ['matVisConfig_' compName];
elseif exist(fullfile(matVisPath,'matVisConfig.mat'), 'file')
    configFile = 'matVisConfig';
end
if ~isempty(configFile)
    load(fullfile(matVisPath,configFile));
    customConfig.plotDim(customConfig.plotDim > nDim) = [];
    gp = customConfig.winPos.gui;
    if size(customConfig.winPos.imageWin,1) ~= nMat
        customConfig.winPos = defaultConfig.winPos;
    end
    customConfig.winPos.gui = gp;
    customConfig.winPos.gui(2) = customConfig.winPos.gui(2)-(guiSize(2)-customConfig.winPos.gui(4));
    customConfig.winPos.gui(4) = guiSize(2);
    % Check for missing fields in the config structure
    configFields = fieldnames(defaultConfig);
    flag = 0;
    for i = 1:numel(configFields)
        if ~isfield(customConfig,configFields{i})
            customConfig.(configFields{i}) = defaultConfig.(configFields{i});
            flag = 1;
        end
    end
    if flag
        save(fullfile(matVisPath,'matVisConfig.mat'), 'customConfig');
    end
else
    customConfig = defaultConfig;
end
% currConfig = customConfig; %Not necessary after removing currConfig
% option

% Start values (if specified)
updateProj = 0;  % Needs to be updated if projection paramter is specified
updateRGB = 0;  % Needs to be updated if projection paramter is specified
if ~isempty(startPar)
    for i=1:length(startPar)/2
        propName = startPar{2*i-1};
        propVal  = startPar{2*i};
        switch propName
            case 'xySel'   % checked
                xySel = propVal;
            case 'zoomVal'   % checked
                zoomVal = propVal;
            case 'plotDim'   % checked
                customConfig.plotDim = propVal;
            case 'plotMean'   % checked
                customConfig.plotMean = propVal;
            case 'currPos'  % checked
                if isa(propVal(1),'char')
                    if strcmp(propVal, 'min')
                        currPos = minPos;
                    elseif strcmp(propVal, 'max')
                        currPos = maxPos;
                    end
                else
                    currPos = propVal;
                end
            case 'cmMinMax'   % checked
                cmMinMax = propVal;
            case 'cmap'       % checked
                customConfig.colormap = propVal;
            case 'aspRatio'   % checked
                customConfig.aspectRatio = propVal;
            case 'rgbDim'     % seems ok
                rgbDim = propVal;
                updateRGB = 1;
            case 'projDim'   % checked
                updateProj = 1;
                projDim = propVal;
            case 'projMethod'   % checked
                if isa(propVal(1),'char')
                    switch propVal
                        case 'max'
                            projMethod = 1;
                        case 'min'
                            projMethod = 2;
                        case 'mean'
                            projMethod = 3;
                        case 'std'
                            projMethod = 4;
                        case 'var'
                            projMethod = 5;
                        case 'tile'
                            projMethod = 6;
                    end
                else
                    projMethod = propVal;
                end
            case 'windowVisibility'
                customConfig.winVis.imageWin = propVal(1);
                customConfig.winVis.zoomWin  = propVal(2);
                customConfig.winVis.plotWin  = propVal(3);
            case 'windowPositions'
                if isfield(propVal,'imageWin')
                    customConfig.winPos.imageWin = propVal.imageWin;
                elseif isfield(propVal,'gui')
                    customConfig.winPos.imageWin(1:2) = customConfig.winPos.imageWin(1:2)+propVal.gui(1:2)-customConfig.winPos.gui(1:2);
                end
                if isfield(propVal,'zoomWin')
                    customConfig.winPos.zoomWin = propVal.zoomWin;
                elseif isfield(propVal,'gui')
                    customConfig.winPos.zoomWin(1:2) = customConfig.winPos.zoomWin(1:2)+propVal.gui(1:2)-customConfig.winPos.gui(1:2);
                end
                if isfield(propVal,'plotWin')
                    customConfig.winPos.plotWin = propVal.plotWin;
                elseif isfield(propVal,'gui')
                    customConfig.winPos.plotWin(1:2) = customConfig.winPos.plotWin(1:2)+propVal.gui(1:2)-customConfig.winPos.gui(1:2);
                end
                if isfield(propVal,'gui')
                    propVal.gui(3:4) = customConfig.winPos.gui(3:4); % Ignore width and height of GUI, since it is fixed by matVis
                    customConfig.winPos.gui = propVal.gui;
                end
                  
                
            case 'objVisibility'
                customConfig.lineVis = propVal;
            case 'calledAs2DHist'
                sendFrom2DHist = propVal{1};  % Function handle used when matVis is called to display 2D histogram
                xVal2DHIst = propVal{2};
                yVal2DHIst = propVal{3};
                isMatVis2DHist = 1;  % Flag to indicate current matVis is called to display 2D hist
            case 'recordFct'
                recordFct = val;
                withRecord = 1;
        end
    end
end
zoomValXY([1,3]) = zoomVal(xySel(2),:);
zoomValXY([2,4]) = zoomVal(xySel(1),:);
%% Icons
icon_matVis = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,50,34,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,101,218,252,34,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,151,252,252,0,0;
    0,67,252,252,34,168,252,235,101,34,168,252,235,67,0,0,0,0,67,168,218,252,252,185,84,0,252,252,252,252,252,151;
    0,118,252,252,235,252,252,252,235,202,252,252,252,218,0,0,0,151,252,252,252,252,252,252,67,50,252,252,252,252,252,118;
    0,151,252,252,151,0,185,252,252,168,17,168,252,252,0,0,134,252,252,118,17,118,252,252,17,0,67,252,252,101,0,0;
    0,185,252,235,17,0,151,252,252,17,0,118,252,218,0,34,252,252,151,0,0,168,252,218,0,0,118,252,252,50,0,0;
    0,252,252,151,0,0,185,252,185,0,0,168,252,185,0,118,252,252,34,0,0,218,252,168,0,0,151,252,252,0,0,0;
    34,252,252,101,0,0,252,252,118,0,0,218,252,118,0,185,252,252,0,0,67,252,252,118,0,0,185,252,185,0,0,0;
    84,252,252,67,0,50,252,252,84,0,0,252,252,101,0,185,252,252,67,17,202,252,252,101,0,0,252,252,218,67,17,0;
    118,252,252,0,0,84,252,252,34,0,67,252,252,50,0,118,252,252,252,252,185,252,252,67,0,0,235,252,252,252,17,0;
    168,252,202,0,0,118,252,252,0,0,101,252,252,0,0,0,134,235,235,134,17,252,252,67,0,0,101,218,252,235,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,29,7,0,0,0,0,0,0,0,0,0,0,0,0,0;
    58,108,108,108,14,0,0,0,0,0,79,108,108,86,0,0,72,108,101,22,0,0,0,0,0,0,0,0,0,0,0,0;
    22,108,108,108,43,0,0,0,0,7,108,108,108,50,0,0,108,108,108,50,0,0,0,0,0,0,0,0,0,0,0,0;
    0,101,108,108,72,0,0,0,0,43,108,108,108,14,0,0,50,108,79,7,0,0,0,0,0,0,0,0,0,0,0,0;
    0,65,108,108,101,0,0,0,0,65,108,108,86,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,36,108,108,108,29,0,0,0,101,108,108,50,0,0,0,108,108,108,29,0,0,0,0,36,79,108,101,79,43,0,0;
    0,0,101,108,108,58,0,0,22,108,108,108,14,0,0,0,108,108,108,29,0,0,0,58,108,108,108,108,108,58,0,0;
    0,0,65,108,108,86,0,0,58,108,108,86,0,0,0,0,108,108,108,29,0,0,7,108,108,108,50,36,58,29,0,0;
    0,0,36,108,108,108,7,0,86,108,108,50,0,0,0,0,108,108,108,29,0,0,29,108,108,108,22,0,0,0,0,0;
    0,0,7,108,108,108,43,7,108,108,108,14,0,0,0,0,108,108,108,29,0,0,7,101,108,108,108,65,14,0,0,0;
    0,0,0,72,108,108,58,36,108,108,86,0,0,0,0,0,108,108,108,29,0,0,0,43,101,108,108,108,108,43,0,0;
    0,0,0,43,108,108,86,58,108,108,50,0,0,0,0,0,108,108,108,29,0,0,0,0,7,50,101,108,108,101,7,0;
    0,0,0,7,108,108,108,86,108,108,14,0,0,0,0,0,108,108,108,29,0,0,0,0,0,0,7,108,108,108,29,0;
    0,0,0,0,79,108,108,108,108,86,0,0,0,0,0,0,108,108,108,29,0,0,7,94,58,36,50,108,108,108,14,0;
    0,0,0,0,50,108,108,108,108,50,0,0,0,0,0,0,108,108,108,29,0,0,36,108,108,108,108,108,108,58,0,0;
    0,0,0,0,14,108,108,108,108,14,0,0,0,0,0,0,108,108,108,29,0,0,14,58,79,108,101,79,43,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
icon_matVis(:,:,2) =    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,33,22,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,66,142,164,22,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,98,164,164,0,0;
    0,44,164,164,22,109,164,153,66,22,109,164,153,44,0,0,0,0,44,109,142,164,164,120,55,0,164,164,164,164,164,98;
    0,77,164,164,153,164,164,164,153,131,164,164,164,142,0,0,0,98,164,164,164,164,164,164,44,33,164,164,164,164,164,77;
    0,98,164,164,98,0,120,164,164,109,11,109,164,164,0,0,87,164,164,77,11,77,164,164,11,0,44,164,164,66,0,0;
    0,120,164,153,11,0,98,164,164,11,0,77,164,142,0,22,164,164,98,0,0,109,164,142,0,0,77,164,164,33,0,0;
    0,164,164,98,0,0,120,164,120,0,0,109,164,120,0,77,164,164,22,0,0,142,164,109,0,0,98,164,164,0,0,0;
    22,164,164,66,0,0,164,164,77,0,0,142,164,77,0,120,164,164,0,0,44,164,164,77,0,0,120,164,120,0,0,0;
    55,164,164,44,0,33,164,164,55,0,0,164,164,66,0,120,164,164,44,11,131,164,164,66,0,0,164,164,142,44,11,0;
    77,164,164,0,0,55,164,164,22,0,44,164,164,33,0,77,164,164,164,164,120,164,164,44,0,0,153,164,164,164,11,0;
    109,164,131,0,0,77,164,164,0,0,66,164,164,0,0,0,87,153,153,87,11,164,164,44,0,0,66,142,164,153,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,53,13,0,0,0,0,0,0,0,0,0,0,0,0,0;
    106,198,198,198,26,0,0,0,0,0,145,198,198,158,0,0,132,198,185,40,0,0,0,0,0,0,0,0,0,0,0,0;
    40,198,198,198,79,0,0,0,0,13,198,198,198,92,0,0,198,198,198,92,0,0,0,0,0,0,0,0,0,0,0,0;
    0,185,198,198,132,0,0,0,0,79,198,198,198,26,0,0,92,198,145,13,0,0,0,0,0,0,0,0,0,0,0,0;
    0,119,198,198,185,0,0,0,0,119,198,198,158,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,66,198,198,198,53,0,0,0,185,198,198,92,0,0,0,198,198,198,53,0,0,0,0,66,145,198,185,145,79,0,0;
    0,0,185,198,198,106,0,0,40,198,198,198,26,0,0,0,198,198,198,53,0,0,0,106,198,198,198,198,198,106,0,0;
    0,0,119,198,198,158,0,0,106,198,198,158,0,0,0,0,198,198,198,53,0,0,13,198,198,198,92,66,106,53,0,0;
    0,0,66,198,198,198,13,0,158,198,198,92,0,0,0,0,198,198,198,53,0,0,53,198,198,198,40,0,0,0,0,0;
    0,0,13,198,198,198,79,13,198,198,198,26,0,0,0,0,198,198,198,53,0,0,13,185,198,198,198,119,26,0,0,0;
    0,0,0,132,198,198,106,66,198,198,158,0,0,0,0,0,198,198,198,53,0,0,0,79,185,198,198,198,198,79,0,0;
    0,0,0,79,198,198,158,106,198,198,92,0,0,0,0,0,198,198,198,53,0,0,0,0,13,92,185,198,198,185,13,0;
    0,0,0,13,198,198,198,158,198,198,26,0,0,0,0,0,198,198,198,53,0,0,0,0,0,0,13,198,198,198,53,0;
    0,0,0,0,145,198,198,198,198,158,0,0,0,0,0,0,198,198,198,53,0,0,13,172,106,66,92,198,198,198,26,0;
    0,0,0,0,92,198,198,198,198,92,0,0,0,0,0,0,198,198,198,53,0,0,66,198,198,198,198,198,198,106,0,0;
    0,0,0,0,26,198,198,198,198,26,0,0,0,0,0,0,198,198,198,53,0,0,26,106,145,198,185,145,79,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
icon_matVis(:,:,3) =    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,10,12,2,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,7,12,12,0,0;
    0,3,12,12,2,8,12,11,5,2,8,12,11,3,0,0,0,0,3,8,10,12,12,9,4,0,12,12,12,12,12,7;
    0,6,12,12,11,12,12,12,11,10,12,12,12,10,0,0,0,7,12,12,12,12,12,12,3,2,12,12,12,12,12,6;
    0,7,12,12,7,0,9,12,12,8,1,8,12,12,0,0,6,12,12,6,1,6,12,12,1,0,3,12,12,5,0,0;
    0,9,12,11,1,0,7,12,12,1,0,6,12,10,0,2,12,12,7,0,0,8,12,10,0,0,6,12,12,2,0,0;
    0,12,12,7,0,0,9,12,9,0,0,8,12,9,0,6,12,12,2,0,0,10,12,8,0,0,7,12,12,0,0,0;
    2,12,12,5,0,0,12,12,6,0,0,10,12,6,0,9,12,12,0,0,3,12,12,6,0,0,9,12,9,0,0,0;
    4,12,12,3,0,2,12,12,4,0,0,12,12,5,0,9,12,12,3,1,10,12,12,5,0,0,12,12,10,3,1,0;
    6,12,12,0,0,4,12,12,2,0,3,12,12,2,0,6,12,12,12,12,9,12,12,3,0,0,11,12,12,12,1,0;
    8,12,10,0,0,6,12,12,0,0,5,12,12,0,0,0,6,11,11,6,1,12,12,3,0,0,5,10,12,11,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,50,13,0,0,0,0,0,0,0,0,0,0,0,0,0;
    100,188,188,188,25,0,0,0,0,0,138,188,188,150,0,0,125,188,175,38,0,0,0,0,0,0,0,0,0,0,0,0;
    38,188,188,188,75,0,0,0,0,13,188,188,188,88,0,0,188,188,188,88,0,0,0,0,0,0,0,0,0,0,0,0;
    0,175,188,188,125,0,0,0,0,75,188,188,188,25,0,0,88,188,138,13,0,0,0,0,0,0,0,0,0,0,0,0;
    0,113,188,188,175,0,0,0,0,113,188,188,150,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,63,188,188,188,50,0,0,0,175,188,188,88,0,0,0,188,188,188,50,0,0,0,0,63,138,188,175,138,75,0,0;
    0,0,175,188,188,100,0,0,38,188,188,188,25,0,0,0,188,188,188,50,0,0,0,100,188,188,188,188,188,100,0,0;
    0,0,113,188,188,150,0,0,100,188,188,150,0,0,0,0,188,188,188,50,0,0,13,188,188,188,88,63,100,50,0,0;
    0,0,63,188,188,188,13,0,150,188,188,88,0,0,0,0,188,188,188,50,0,0,50,188,188,188,38,0,0,0,0,0;
    0,0,13,188,188,188,75,13,188,188,188,25,0,0,0,0,188,188,188,50,0,0,13,175,188,188,188,113,25,0,0,0;
    0,0,0,125,188,188,100,63,188,188,150,0,0,0,0,0,188,188,188,50,0,0,0,75,175,188,188,188,188,75,0,0;
    0,0,0,75,188,188,150,100,188,188,88,0,0,0,0,0,188,188,188,50,0,0,0,0,13,88,175,188,188,175,13,0;
    0,0,0,13,188,188,188,150,188,188,25,0,0,0,0,0,188,188,188,50,0,0,0,0,0,0,13,188,188,188,50,0;
    0,0,0,0,138,188,188,188,188,150,0,0,0,0,0,0,188,188,188,50,0,0,13,163,100,63,88,188,188,188,25,0;
    0,0,0,0,88,188,188,188,188,88,0,0,0,0,0,0,188,188,188,50,0,0,63,188,188,188,188,188,188,100,0,0;
    0,0,0,0,25,188,188,188,188,25,0,0,0,0,0,0,188,188,188,50,0,0,25,100,138,188,175,138,75,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
icon_image = [0,0,0,0,0,0,0,15,15,15,15,15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,15,31,31,47,47,47,47,31,15,15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,31,47,63,63,79,79,79,79,79,63,47,31,15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,31,63,79,95,95,111,111,111,127,111,111,95,79,63,47,31,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,31,62,79,94,110,126,142,143,158,157,157,158,159,143,127,95,79,47,31,0,0,0,0,0,0,0,0,0,0,0;0,15,46,78,110,126,141,157,157,174,205,205,220,220,207,191,175,158,126,110,78,46,31,0,0,0,0,0,0,0,0,0;15,31,62,94,126,143,158,173,189,221,239,253,252,252,255,255,237,220,189,174,143,111,78,46,15,0,0,0,0,0,0,0;31,47,78,110,127,159,174,189,220,255,255,254,252,252,255,255,252,252,253,239,207,174,126,94,63,15,0,0,0,0,0,0;31,62,94,110,143,159,174,205,237,255,255,253,252,253,255,254,252,252,255,255,255,236,189,158,111,63,15,0,0,0,0,0;31,62,94,127,143,159,173,205,254,255,255,252,252,254,255,253,252,252,255,255,254,252,252,206,159,111,63,15,0,0,0,0;31,62,94,111,143,158,173,205,254,255,255,252,252,254,255,253,252,236,238,255,253,252,252,254,223,159,110,46,0,0,0,0;31,46,78,111,143,158,173,205,239,255,254,252,252,254,255,254,236,205,189,189,237,252,252,254,255,223,141,94,31,0,0,0;15,46,78,95,127,142,173,189,223,255,254,252,252,255,255,255,253,205,173,142,159,205,252,254,255,255,190,126,62,15,0,0;0,31,63,95,111,143,159,175,207,239,255,255,255,255,255,255,255,255,191,143,127,143,207,255,255,255,239,159,95,47,0,0;0,0,31,63,95,111,143,159,175,207,223,255,255,255,255,255,255,216,245,175,127,127,159,223,255,255,255,191,127,63,15,0;0,0,6,20,73,95,111,127,143,175,110,88,101,137,255,255,157,108,118,212,159,127,143,207,255,255,255,223,143,79,31,0;0,0,0,6,36,63,95,111,127,137,67,74,81,151,239,255,108,108,108,186,223,175,175,207,255,255,255,239,175,111,47,15;0,0,0,0,9,31,63,79,95,85,54,61,67,161,207,239,186,108,147,245,255,239,223,223,255,255,255,255,191,127,63,31;0,0,0,0,0,15,31,47,63,52,40,47,68,143,175,207,239,255,255,255,255,255,255,239,255,255,255,255,191,127,79,47;0,0,0,0,0,0,0,15,31,22,27,33,69,111,143,175,88,108,108,216,255,255,255,255,206,147,108,118,119,110,79,47;0,0,0,0,0,0,0,0,0,6,20,27,73,95,111,143,74,94,108,216,255,255,255,177,108,108,108,108,88,99,95,63;0,0,0,0,0,0,0,0,0,0,6,17,47,79,95,127,67,81,101,216,255,255,245,108,108,108,186,206,133,121,95,63;0,0,0,0,0,0,0,0,0,0,0,11,31,47,79,95,54,74,88,216,255,255,216,108,108,108,226,239,191,127,79,63;0,0,0,0,0,0,0,0,0,0,0,0,0,31,47,79,47,61,74,175,255,255,245,118,108,108,108,136,147,111,79,47;0,0,0,0,0,0,0,0,0,0,0,0,0,15,31,63,33,47,61,148,223,239,255,196,118,108,94,81,61,73,63,31;0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,31,27,40,54,121,175,207,223,223,214,151,88,67,47,37,45,15;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,20,27,40,94,143,159,175,175,175,159,137,47,40,20,26,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,20,27,67,111,127,122,72,88,103,81,33,27,13,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,13,40,63,79,77,40,40,33,27,20,13,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,13,31,47,43,33,27,20,14,9,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,15,15,15,15,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
icon_image(:,:,2) = [143,159,191,207,223,239,255,255,255,255,255,255,255,255,239,223,207,191,175,159,127,111,95,79,47,31,15,0,0,0,0,0;159,175,207,223,255,255,255,255,255,255,255,255,255,255,255,255,239,223,191,175,159,143,111,95,79,47,31,0,0,0,0,0;191,207,239,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,239,223,191,175,159,127,111,79,63,27,10,0,0,0;207,239,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,239,223,191,175,143,111,95,50,20,10,0,0;239,231,165,165,243,195,165,171,219,243,195,165,171,231,255,255,255,255,231,195,177,165,155,153,154,159,82,61,41,20,10,0;255,213,165,165,171,165,165,165,171,183,165,165,165,177,255,255,255,201,165,165,165,165,165,165,202,178,103,82,61,41,20,13;255,201,165,165,201,255,189,165,165,195,249,195,155,155,239,255,207,165,165,213,249,213,165,165,249,223,173,103,82,82,63,47;255,189,165,171,249,255,201,165,165,249,223,173,124,133,191,182,134,134,188,255,255,195,165,177,255,255,200,124,103,118,95,63;255,165,165,201,255,255,189,165,189,239,207,134,103,118,143,119,82,82,136,175,223,177,165,195,255,255,201,155,124,159,111,95;243,165,165,219,255,255,165,165,213,223,191,110,93,119,111,59,41,30,47,95,130,124,155,213,255,255,189,165,165,191,143,111;225,165,165,231,255,237,165,165,225,223,175,103,93,109,79,35,0,0,0,0,45,72,113,192,255,255,165,165,177,188,171,143;213,165,165,255,255,225,165,165,243,223,173,103,93,103,79,13,0,0,0,0,0,20,51,144,223,255,171,165,165,155,187,175;195,165,183,255,255,213,165,165,255,239,178,124,103,127,79,31,0,0,0,0,0,0,10,72,159,239,219,177,165,171,223,191;255,255,255,255,255,255,255,255,255,255,239,223,191,159,111,79,31,0,0,0,0,0,0,15,95,191,255,255,255,255,255,223;239,255,255,255,255,255,255,255,255,255,255,255,239,207,175,127,95,59,0,0,0,0,0,0,47,143,239,255,255,255,255,239;183,186,198,198,247,255,255,255,255,255,213,198,198,209,223,191,135,86,50,0,0,0,0,0,15,111,207,255,255,255,255,255;183,161,186,198,232,255,255,255,255,251,198,198,198,228,255,239,148,123,86,42,0,0,0,0,0,79,175,255,255,255,255,255;159,151,161,186,217,255,255,255,255,232,198,198,198,247,255,255,214,148,133,94,31,0,0,0,0,63,159,255,255,255,255,255;143,138,148,161,189,255,255,255,255,221,198,198,209,255,255,255,255,223,191,143,95,47,0,0,15,79,159,239,255,255,255,255;111,118,123,148,161,225,255,255,255,202,198,198,228,255,255,255,198,198,161,165,127,79,31,15,29,79,136,189,213,232,255,255;95,111,101,123,136,183,223,239,244,198,198,198,247,255,255,255,198,198,186,180,159,111,63,41,49,99,148,198,198,225,255,255;63,79,96,99,123,143,191,207,211,198,198,209,255,255,255,255,198,198,198,210,191,143,109,74,86,123,185,236,225,240,255,255;47,63,73,86,99,111,172,191,170,173,186,228,255,255,255,255,198,198,198,240,223,191,150,111,123,136,213,255,255,255,255,255;15,31,62,61,86,99,130,157,136,161,173,232,255,255,255,255,198,198,198,240,255,223,188,151,148,161,198,221,247,255,255,255;0,15,31,54,61,74,112,132,123,136,157,223,239,255,255,255,198,198,198,240,255,255,239,217,189,198,198,198,198,232,255,255;0,0,15,28,49,61,78,98,111,123,156,191,223,239,255,255,198,198,198,240,255,255,255,255,251,228,202,198,198,202,251,255;0,0,0,15,24,49,61,78,86,99,154,175,191,223,239,255,198,198,198,240,255,255,255,255,255,255,251,198,198,198,240,255;0,0,0,0,13,24,36,61,74,91,127,159,175,191,223,239,198,198,198,240,255,255,251,206,225,236,228,198,198,198,247,239;0,0,0,0,0,12,24,36,61,85,111,127,159,175,191,223,186,198,198,240,255,255,236,198,198,198,198,198,198,225,223,207;0,0,0,0,0,0,12,24,36,61,95,111,127,143,175,191,161,186,198,240,255,255,247,225,213,198,202,213,232,223,207,191;0,0,0,0,0,0,0,15,31,47,63,95,111,127,143,159,191,207,223,239,255,255,255,255,255,255,255,239,223,207,175,159;0,0,0,0,0,0,0,0,15,31,47,79,95,111,127,143,159,191,207,223,223,239,239,255,255,239,239,223,207,175,159,143];
icon_image(:,:,3) = [255,255,255,255,255,255,255,239,239,239,239,239,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,239,223,207,191,175;255,255,255,255,255,239,223,223,207,207,207,207,223,239,239,255,255,255,255,255,255,255,255,255,255,255,255,255,239,181,167,191;255,255,255,255,223,207,191,191,175,175,175,175,175,191,207,223,239,255,255,255,255,255,255,255,255,255,255,159,47,14,195,207;255,255,255,223,191,175,159,159,143,143,143,127,143,143,159,175,191,207,223,255,255,255,255,255,255,255,255,111,15,15,239,223;255,191,13,11,153,59,8,15,69,97,35,6,12,71,95,111,127,159,131,77,41,15,15,79,175,255,15,15,15,15,15,111;255,134,12,10,17,7,7,6,12,20,3,3,2,6,47,63,79,41,7,8,10,12,13,15,191,207,15,15,15,15,15,143;239,97,11,9,55,111,29,5,4,12,14,0,0,0,0,0,7,2,4,44,104,80,10,12,224,255,191,15,15,159,255,255;223,64,10,17,119,95,34,4,2,0,0,0,0,0,0,0,0,0,0,15,47,29,7,29,191,239,143,15,15,207,255,255;223,11,9,62,111,95,24,3,5,0,0,0,0,0,0,0,0,0,0,0,0,3,4,35,143,191,104,15,15,255,255,255;195,11,9,79,111,95,5,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,26,95,143,59,14,79,255,255,255;153,11,9,107,111,77,5,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,31,95,8,12,47,191,239,255;125,12,10,143,111,65,5,3,13,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,31,13,9,13,15,239,255;89,12,43,159,127,62,5,4,31,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,39,23,11,29,255,255;255,223,191,159,143,111,95,79,47,15,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,15,95,159,207,255,255;255,255,223,191,159,143,111,95,79,47,31,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,63,127,191,239,255;219,188,176,153,169,159,143,127,111,79,51,35,11,0,0,0,0,0,0,0,0,0,0,0,0,0,0,31,111,175,223,255;242,188,188,176,185,191,159,143,127,109,70,58,46,41,15,0,0,0,0,0,0,0,0,0,0,0,0,15,79,143,207,239;255,192,188,188,197,223,191,175,159,128,94,82,70,76,47,15,0,0,0,0,0,0,0,0,0,0,0,0,63,127,191,223;255,215,188,188,192,239,223,207,191,148,117,105,100,111,79,47,15,0,0,0,0,0,0,0,0,0,0,0,63,127,175,207;255,233,188,188,188,237,255,239,223,156,141,129,140,143,111,79,35,0,0,0,0,0,0,0,0,0,0,0,38,99,175,207;255,255,192,188,188,219,255,255,242,176,153,141,169,159,143,111,58,23,0,0,0,0,0,0,0,0,0,0,35,95,159,191;255,255,215,188,188,201,255,255,219,188,176,176,207,175,159,127,70,46,11,0,0,0,0,0,0,0,0,0,54,103,159,191;255,255,233,188,188,188,251,255,201,188,188,210,223,207,175,159,94,58,35,0,0,0,0,0,0,0,0,15,63,127,175,191;255,255,251,188,188,188,228,251,188,188,188,246,255,223,207,175,105,82,58,44,0,0,0,0,0,0,0,40,92,143,175,207;255,255,255,210,188,188,219,233,188,188,201,255,255,239,223,191,129,105,82,73,31,15,0,0,0,0,23,46,82,142,191,223;239,255,255,228,188,188,201,219,188,188,224,255,255,255,239,223,141,117,94,103,79,47,31,31,31,41,47,70,105,132,204,239;223,239,255,251,188,188,188,201,188,188,246,255,255,255,255,239,153,141,117,133,111,95,79,79,79,95,109,105,117,153,207,255;191,207,223,255,206,188,188,188,188,201,255,255,255,255,255,255,176,153,141,163,143,127,125,86,109,116,126,129,141,164,246,255;175,191,207,223,224,188,188,188,188,224,255,255,255,255,255,255,188,176,164,192,191,175,145,117,117,129,141,153,164,219,255,255;159,175,191,207,215,188,188,188,188,246,255,255,255,255,255,255,188,188,188,222,223,207,200,178,167,153,168,193,228,255,255,255;143,159,175,191,207,223,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,239,239,239,239,239,255,255,255,255,255,255;143,143,159,175,191,207,239,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_image_24x24 = [19,10,10,10,13,23,29,29,28,13,10,10,10,10,10,10,10,10,10,10,10,10,10,19;8,0,0,15,39,48,62,63,63,52,32,16,2,0,0,0,0,0,0,0,0,0,0,8;10,0,19,59,82,92,105,105,114,108,95,73,49,30,3,0,0,0,0,0,0,0,0,10;10,15,61,89,110,131,144,157,168,171,170,154,128,95,61,28,5,0,0,0,0,0,0,10;18,38,83,124,143,164,177,212,234,243,238,228,209,175,144,105,65,25,3,0,0,0,0,10;39,60,103,135,165,185,224,253,255,254,254,255,255,248,233,196,143,92,44,2,0,0,0,10;46,80,114,145,167,201,246,255,255,254,254,255,254,254,255,255,225,172,110,44,2,0,0,10;45,81,117,146,166,202,255,255,255,254,255,255,253,248,255,255,254,240,183,111,37,0,0,10;42,67,106,145,165,201,246,255,254,254,255,255,236,202,202,242,254,254,249,176,91,18,0,10;26,62,96,128,159,186,229,255,254,254,255,255,247,193,141,157,222,254,255,233,143,58,6,10;11,30,74,108,141,166,198,234,255,255,255,255,255,243,167,118,152,231,255,255,193,101,28,10;10,3,34,75,107,132,155,189,185,201,219,255,243,234,231,145,126,195,255,255,228,132,50,12;10,0,7,36,73,104,124,143,136,154,198,252,212,209,234,221,180,203,255,255,247,171,80,28;10,0,0,6,33,61,82,96,102,118,168,220,220,212,243,255,247,234,254,255,255,183,98,51;10,0,0,0,3,16,38,54,67,86,131,174,223,255,255,255,255,254,255,255,255,195,108,60;10,0,0,0,0,0,0,21,44,69,100,134,148,195,234,255,255,252,231,212,209,166,111,74;10,0,0,0,0,0,0,1,14,40,76,107,122,167,232,255,255,228,209,212,219,157,104,73;10,0,0,0,0,0,0,0,1,10,43,76,97,134,195,255,255,209,209,228,225,158,93,62;10,0,0,0,0,0,0,0,0,0,18,50,68,103,152,219,245,228,206,193,160,118,76,42;10,0,0,0,0,0,0,0,0,0,2,19,48,79,116,162,193,202,184,151,116,78,46,22;10,0,0,0,0,0,0,0,0,0,0,2,20,48,76,118,133,142,137,124,82,50,13,10;10,0,0,0,0,0,0,0,0,0,0,0,1,15,34,62,81,73,70,58,35,15,0,10;8,0,0,0,0,0,0,0,0,0,0,0,0,0,3,17,30,26,25,18,5,0,0,8;19,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,10,10,10,10,19];
icon_image_24x24(:,:,2) = [153,179,209,234,250,255,255,255,255,255,250,234,213,188,164,135,111,89,57,33,10,10,10,19;178,206,240,255,255,255,255,255,255,255,255,255,249,229,200,172,146,114,88,54,20,4,0,8;213,246,255,255,255,255,255,255,255,255,255,255,255,255,252,230,202,169,131,98,59,22,4,10;246,236,235,246,235,233,244,247,236,233,247,255,255,254,244,236,228,206,172,137,96,55,23,13;255,232,232,233,236,232,232,235,236,228,230,255,254,235,232,236,236,232,223,177,135,98,57,40;252,232,233,252,249,232,232,241,215,183,184,206,204,210,241,255,241,232,244,233,178,141,105,73;247,232,240,255,246,232,236,221,169,141,137,129,105,106,158,219,232,232,249,240,224,188,141,103;243,232,244,255,241,232,239,201,152,130,110,63,22,14,36,105,167,220,254,235,232,230,179,136;238,232,249,255,238,232,242,203,152,125,92,27,0,0,0,10,60,145,239,233,232,235,207,173;236,232,254,255,235,232,251,228,183,145,105,42,0,0,0,0,0,48,166,241,235,235,244,206;251,255,255,255,255,255,255,253,243,212,162,108,55,11,0,0,0,1,83,214,255,255,255,240;219,233,237,244,255,255,255,255,241,236,222,188,139,74,8,0,0,0,32,167,253,255,255,253;190,208,233,240,255,255,255,253,237,237,246,247,185,135,61,0,0,0,8,122,244,255,255,255;159,182,204,232,251,255,255,248,237,237,251,255,232,180,128,54,8,0,13,111,233,255,255,255;124,146,175,204,243,255,255,244,237,240,255,255,255,230,176,117,46,8,35,131,229,255,255,255;96,114,139,166,201,231,255,238,237,246,255,255,237,234,201,155,91,51,75,152,228,240,248,255;64,83,110,135,159,198,223,227,237,249,255,255,237,237,238,201,150,114,128,177,237,240,249,255;31,50,81,106,129,161,179,202,225,255,255,255,237,237,247,246,203,170,175,216,254,255,254,255;10,19,51,74,99,130,152,168,203,242,255,255,237,237,247,255,249,233,229,237,238,247,254,255;10,2,18,49,74,96,120,144,176,211,241,255,237,237,247,255,255,255,248,241,237,237,244,255;10,0,1,18,42,72,93,114,155,180,212,242,237,237,247,255,255,253,255,255,246,237,237,245;10,0,0,1,17,40,69,94,123,154,178,213,222,237,247,255,255,238,238,242,237,237,226,210;8,0,0,0,2,18,37,64,99,120,146,172,189,217,239,255,254,244,240,237,236,229,201,178;19,10,10,10,10,11,29,48,82,107,126,145,173,203,224,235,243,253,253,243,232,206,173,154];
icon_image_24x24(:,:,3) = [255,255,255,255,252,241,235,235,236,251,255,255,255,255,255,255,255,255,255,255,246,224,201,184;255,255,255,240,215,206,192,191,191,202,221,239,253,255,255,255,255,255,255,255,243,215,208,200;255,255,235,194,172,162,149,149,139,146,159,181,205,225,251,255,255,255,255,255,202,193,237,225;255,192,153,149,114,94,97,89,68,63,77,100,126,157,171,182,189,202,218,202,193,193,193,234;247,163,129,100,88,67,58,32,16,8,12,27,44,62,82,120,152,173,211,205,193,202,214,246;217,147,117,115,83,51,22,2,0,0,0,0,0,5,18,58,95,123,186,223,193,222,255,255;200,132,117,109,78,39,6,0,0,0,0,0,0,0,0,0,23,61,135,176,191,234,255,255;191,131,121,108,75,39,0,0,0,0,0,0,0,0,0,0,0,11,70,113,164,243,255,255;182,142,139,108,72,39,7,0,0,0,0,0,0,0,0,0,0,0,5,59,123,188,243,255;192,145,156,126,74,51,24,0,0,0,0,0,0,0,0,0,0,0,0,18,86,155,244,255;254,225,179,146,113,89,56,20,0,0,0,0,0,0,0,0,0,0,0,0,61,153,226,255;252,229,195,162,147,122,99,65,36,8,0,0,0,0,0,0,0,0,0,0,27,122,204,252;255,238,225,197,181,150,130,107,81,61,33,3,0,0,0,0,0,0,0,0,8,83,174,236;255,245,234,228,217,193,172,146,119,101,79,34,4,0,0,0,0,0,0,0,0,71,156,213;255,251,234,234,241,238,216,185,158,141,123,80,31,0,0,0,0,0,0,0,0,59,146,204;255,255,237,234,241,255,255,212,184,170,154,120,67,15,0,0,0,0,0,0,0,54,131,190;255,255,243,234,236,255,251,232,217,206,179,147,96,46,1,0,0,0,0,0,3,64,140,191;255,255,248,234,234,251,245,234,235,244,211,178,125,83,39,0,0,0,0,0,29,96,159,202;254,255,254,234,234,245,241,234,243,255,236,204,157,117,85,35,9,3,4,17,56,120,176,221;236,254,255,241,234,243,237,234,247,255,253,235,180,145,123,92,61,52,54,73,103,146,193,243;201,222,250,245,234,237,234,234,254,255,255,253,210,179,164,136,121,108,117,130,157,177,218,255;177,195,220,248,234,234,234,241,255,255,255,255,233,217,208,192,173,154,157,176,194,217,243,255;155,175,196,219,232,234,234,247,255,255,255,255,234,234,241,237,223,214,208,213,231,245,254,255;149,157,177,196,227,254,255,255,255,255,255,255,255,255,255,255,254,254,254,255,255,255,255,255];

icon_zoom = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,254,254,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,254,254,254,255,255;255,255,254,254,255,254,254,254,255,194,172,156,143,140,144,167,221,255,255,254,254,254,254,254,255,255,254,254,254,254,254,254;255,255,254,254,254,254,254,254,180,136,124,113,104,100,101,110,130,197,254,254,254,254,254,254,255,255,254,254,254,254,254,255;255,254,254,254,254,255,254,153,145,161,188,205,212,196,173,127,103,106,159,255,255,255,254,254,255,255,255,254,254,255,255,255;255,254,254,254,255,255,177,139,172,248,255,255,254,251,227,177,131,96,107,211,255,254,254,254,255,255,255,254,254,255,255,255;255,254,254,254,255,179,141,166,254,255,255,244,231,222,215,221,205,111,96,121,221,254,254,254,255,255,254,254,254,255,255,255;255,254,254,255,255,156,135,198,255,255,255,234,223,215,206,198,182,137,101,95,168,254,254,255,255,255,254,254,254,255,255,255;255,254,254,255,255,137,127,231,255,245,230,222,214,205,197,183,166,158,104,88,145,254,254,255,255,255,254,254,254,255,255,255;255,254,254,255,255,126,121,240,252,232,222,214,205,197,190,181,168,162,100,87,136,254,254,255,255,255,254,254,254,254,255,255;254,254,254,255,255,119,113,231,241,224,215,208,198,192,186,194,202,168,98,86,135,254,254,255,255,255,255,254,254,254,255,255;255,255,255,255,255,119,109,200,223,214,203,198,190,190,214,237,224,160,96,88,137,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,134,109,150,192,209,190,183,188,206,239,245,197,107,89,92,156,255,255,255,255,255,255,255,255,255,255,255;228,205,205,205,248,163,111,109,163,215,148,119,140,201,218,181,112,41,43,117,199,255,252,255,255,255,255,255,255,255,255,255;245,205,205,205,235,254,142,100,107,129,102,109,126,168,161,113,42,37,37,72,159,255,255,255,255,255,255,255,255,255,255,255;255,208,205,205,222,255,221,119,98,77,62,69,68,105,107,94,64,37,50,84,99,181,244,255,255,255,255,255,255,255,255,255;255,225,205,205,208,255,255,211,144,72,41,38,47,87,88,92,103,94,88,87,87,97,204,255,255,255,255,255,255,255,255,255;255,238,205,205,205,242,255,255,214,107,74,61,81,105,111,136,115,70,43,74,87,87,99,179,238,218,205,208,218,235,255,255;255,255,208,205,205,228,255,255,244,205,205,205,248,255,254,248,180,119,66,83,88,87,87,74,153,205,205,205,205,228,255,255;255,255,225,205,205,215,255,255,228,198,195,199,236,236,235,232,178,177,115,99,94,88,84,37,48,182,231,238,228,242,255,255;255,255,238,205,205,205,252,255,215,205,202,222,241,241,239,237,187,202,173,155,119,94,74,37,37,49,171,255,255,255,255,255;255,255,252,205,205,205,235,252,205,205,205,248,255,255,254,254,204,197,202,226,171,117,93,41,37,37,50,181,248,255,255,255;255,255,255,222,205,205,228,238,205,205,215,255,255,255,255,255,205,205,204,230,232,167,112,74,41,37,37,55,205,235,255,255;255,255,255,235,205,205,215,228,205,205,231,255,255,255,255,255,205,205,205,242,255,228,169,120,88,64,40,55,205,208,251,255;255,255,255,252,205,205,205,215,205,205,248,255,255,255,255,255,205,205,205,242,255,255,241,171,113,98,103,108,191,205,235,255;255,255,255,255,218,205,205,205,205,215,255,255,255,255,255,255,205,205,205,242,255,255,244,181,154,137,153,164,186,205,246,255;255,255,255,255,231,205,205,205,205,231,255,255,255,255,255,255,205,205,205,242,255,255,230,187,150,122,148,178,184,227,255,255;255,255,255,255,248,205,205,205,205,248,255,255,255,255,255,255,205,205,205,242,255,255,248,228,218,205,208,218,235,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_zoom(:,:,2) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,249,251,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,243,228,224,251,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,236,224,224,255,255;255,247,224,224,251,234,224,226,243,190,152,126,115,132,144,167,221,255,247,234,228,224,224,232,245,255,224,224,224,224,224,236;255,241,224,224,226,224,224,224,152,112,94,83,74,75,101,110,130,179,224,224,224,224,224,224,247,249,224,224,224,224,224,241;255,236,224,224,236,255,232,123,115,141,186,185,180,166,173,127,87,76,129,241,253,241,224,224,253,255,247,224,224,243,255,255;255,232,224,226,253,255,158,109,142,246,255,241,224,225,227,173,101,66,89,211,255,234,224,228,255,255,241,224,224,249,255,255;255,224,224,236,255,179,119,136,232,255,255,224,201,200,214,205,175,82,92,121,221,228,224,234,255,255,236,224,224,255,255,255;251,224,224,243,255,156,105,168,241,255,255,208,193,201,206,176,152,107,101,95,161,224,224,241,255,255,232,224,232,255,255,255;245,224,224,247,255,131,97,201,245,245,230,192,184,193,197,160,136,128,96,86,121,224,224,243,255,255,224,224,228,247,253,255;241,224,224,255,255,116,91,210,248,232,215,184,175,192,190,167,138,132,70,57,114,224,224,247,255,255,226,224,224,224,253,255;234,224,230,255,255,105,83,199,241,224,203,178,168,192,186,194,184,140,70,70,133,224,224,247,255,255,243,228,224,226,255,255;255,255,255,255,255,119,110,200,221,214,204,198,190,190,214,237,224,160,96,88,137,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,134,110,150,190,208,190,183,188,206,239,245,197,115,91,92,156,255,255,255,255,255,255,255,255,255,255,255;245,236,236,236,252,163,111,109,162,215,169,149,171,225,218,181,133,72,72,123,199,255,252,255,255,255,255,255,255,255,255,255;251,236,236,236,247,254,142,100,107,131,133,139,157,183,161,113,73,68,68,86,159,255,255,255,255,255,255,255,255,255,255,255;255,237,236,236,242,255,221,119,98,90,93,100,99,109,107,94,78,68,73,86,99,181,244,255,255,255,255,255,255,255,255,255;255,243,236,236,237,255,255,211,144,91,72,69,71,87,88,92,103,94,88,87,87,98,204,255,255,255,255,255,255,255,255,255;255,249,236,236,236,250,255,255,214,135,104,91,96,105,111,135,144,100,73,83,87,87,99,179,249,241,236,237,241,247,255,255;255,255,237,236,236,245,255,255,250,236,236,236,252,255,254,248,211,148,97,91,88,87,87,90,184,236,236,236,236,245,255,255;255,255,243,236,236,239,255,255,245,229,226,223,236,235,232,230,209,205,145,107,94,88,86,68,79,212,246,249,245,250,255,255;255,255,249,236,236,236,254,255,239,236,232,238,241,239,237,237,218,232,203,162,119,94,83,68,68,80,177,255,255,255,255,255;255,255,254,236,236,236,247,254,236,236,236,252,255,255,254,254,232,228,232,234,173,117,95,70,68,68,81,199,252,255,255,255;255,255,255,242,236,236,245,249,236,236,239,255,255,255,255,255,236,236,234,239,235,166,111,86,70,68,68,85,236,247,255,255;255,255,255,247,236,236,239,245,236,236,246,255,255,255,255,255,236,236,236,250,255,227,169,119,90,78,69,86,236,237,253,255;255,255,255,254,236,236,236,239,236,236,252,255,255,255,255,255,236,236,236,250,255,255,241,169,112,98,105,139,222,236,247,255;255,255,255,255,241,236,236,236,236,239,255,255,255,255,255,255,236,236,236,250,255,255,246,207,170,147,167,195,216,234,251,255;255,255,255,255,246,236,236,236,236,246,255,255,255,255,255,255,236,236,236,250,255,255,241,218,181,153,178,209,215,241,255,255;255,255,255,255,252,236,236,236,236,252,255,255,255,255,255,255,236,236,236,250,255,255,252,245,241,236,237,241,247,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_zoom(:,:,3) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,238,244,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,222,183,172,244,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,205,172,172,255,255;255,233,172,172,244,200,172,178,222,183,117,74,67,118,144,167,221,255,233,200,183,172,172,194,227,255,172,172,172,172,172,205;255,216,172,172,178,172,172,172,103,71,42,31,23,30,101,110,130,148,172,172,172,172,172,172,233,238,172,172,172,172,172,216;255,205,172,172,205,255,194,71,63,107,182,150,128,114,173,127,59,24,77,216,250,216,172,172,250,255,233,172,172,222,255,255;255,194,172,178,250,255,127,57,90,243,255,216,172,180,227,166,50,15,58,211,255,200,172,183,255,255,216,172,172,238,255,255;255,172,172,205,255,179,81,84,194,255,255,190,149,162,214,181,123,30,85,121,221,183,172,200,255,255,205,172,172,255,255,255;244,172,172,222,255,156,53,116,216,255,255,163,141,177,206,138,100,55,101,95,147,172,172,216,255,255,194,172,194,255,255,255;227,172,172,233,255,120,46,149,227,245,230,140,133,172,197,122,84,76,82,83,79,172,172,222,255,255,172,172,183,233,250,255;216,172,172,255,255,98,39,158,241,232,201,133,123,181,190,143,86,81,19,5,76,172,172,233,255,255,178,172,172,172,250,255;200,172,189,255,255,81,31,147,241,224,182,126,116,192,186,194,156,91,21,43,130,172,172,233,255,255,222,183,172,178,255,255;255,255,255,255,255,119,110,200,221,214,204,198,190,190,214,237,224,160,96,88,137,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,134,110,150,190,208,190,183,188,206,239,245,197,114,91,92,156,255,255,255,255,255,255,255,255,255,255,255;243,232,232,232,252,163,111,109,162,215,166,146,168,223,218,181,130,69,68,122,199,255,252,255,255,255,255,255,255,255,255,255;251,232,232,232,246,254,142,100,107,131,129,136,154,181,161,113,69,64,64,85,159,255,255,255,255,255,255,255,255,255,255,255;255,234,232,232,240,255,221,119,98,88,89,96,96,108,107,94,77,64,70,86,99,181,244,255,255,255,255,255,255,255,255,255;255,241,232,232,234,255,255,211,144,89,68,65,69,87,88,92,103,94,88,87,87,98,204,255,255,255,255,255,255,255,255,255;255,247,232,232,232,249,255,255,214,132,100,87,94,105,111,135,141,96,69,82,87,87,99,179,247,238,232,234,238,246,255,255;255,255,234,232,232,243,255,255,249,232,232,232,252,255,254,248,207,145,93,90,88,87,87,88,180,232,232,232,232,243,255,255;255,255,241,232,232,237,255,255,243,226,222,221,236,235,232,230,205,201,142,106,94,88,86,64,75,209,244,247,243,249,255,255;255,255,247,232,232,232,254,255,237,232,229,237,241,239,237,237,214,229,200,161,119,94,82,64,64,77,176,255,255,255,255,255;255,255,254,232,232,232,246,254,232,232,232,252,255,255,254,254,229,224,229,233,173,117,95,66,64,64,77,197,252,255,255,255;255,255,255,240,232,232,243,247,232,232,237,255,255,255,255,255,232,232,231,238,235,166,111,85,66,64,64,82,232,246,255,255;255,255,255,246,232,232,237,243,232,232,244,255,255,255,255,255,232,232,232,249,255,227,169,119,90,77,66,83,232,234,253,255;255,255,255,254,232,232,232,237,232,232,252,255,255,255,255,255,232,232,232,249,255,255,241,169,112,98,105,135,218,232,246,255;255,255,255,255,238,232,232,232,232,237,255,255,255,255,255,255,232,232,232,249,255,255,246,204,168,146,166,191,213,231,251,255;255,255,255,255,244,232,232,232,232,244,255,255,255,255,255,255,232,232,232,249,255,255,240,214,177,150,175,205,212,239,255,255;255,255,255,255,252,232,232,232,232,252,255,255,255,255,255,255,232,232,232,249,255,255,252,243,238,232,234,238,246,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_zoom_24x24(:,:,1) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,254,252,252,254,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,229,167,145,129,128,150,213,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,217,151,143,153,152,137,113,112,180,249,255,255,255,255,255,255,255,255,255;255,255,255,255,225,141,176,240,250,253,236,179,116,91,198,255,255,255,255,255,255,255,255,255;255,255,255,240,152,162,253,255,249,232,220,221,182,98,111,219,255,255,255,255,255,255,255,255;255,255,255,231,127,197,255,247,231,217,206,193,172,125,87,175,255,255,255,255,255,255,255,255;255,255,255,225,109,215,254,228,216,206,194,181,172,133,87,162,255,255,255,255,255,255,255,255;255,255,255,223,99,203,242,217,207,195,189,202,206,132,87,159,255,255,255,255,255,255,255,255;255,255,255,225,104,159,212,205,190,187,213,249,204,113,87,168,255,255,255,255,255,255,255,255;253,239,239,227,129,107,177,203,160,187,224,203,136,81,108,207,255,254,255,255,255,255,255,255;255,243,239,242,208,101,103,136,136,156,160,115,72,71,85,169,252,255,255,255,255,255,255,255;255,248,239,239,250,180,107,87,83,81,88,89,79,73,83,96,184,255,255,255,255,255,255,255;255,252,239,239,248,254,193,121,94,91,106,129,147,100,87,87,96,183,255,255,255,255,255,255;255,255,242,239,244,255,255,236,231,237,245,244,207,135,89,87,87,93,189,240,239,242,249,255;255,255,246,239,240,255,252,235,227,233,237,233,223,208,139,99,87,78,86,182,242,242,250,255;255,255,250,239,239,252,248,239,239,252,252,250,235,237,225,153,98,71,71,80,197,255,254,255;255,255,254,239,239,248,244,239,246,255,255,255,239,239,247,229,147,88,71,71,88,244,254,255;255,255,255,244,239,246,242,239,249,255,255,255,239,239,248,255,231,155,88,76,108,231,245,255;255,255,255,248,239,242,239,239,254,255,255,255,239,239,248,255,255,221,156,144,188,225,236,255;255,255,255,253,239,239,239,244,255,255,255,255,239,239,248,255,252,230,195,190,216,228,244,255;255,255,255,255,243,239,239,249,255,255,255,255,239,239,248,255,254,246,242,239,243,248,254,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_zoom_24x24(:,:,2) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,253,251,253,255;255,255,255,255,255,255,255,255,254,252,252,254,255,255,255,255,255,255,255,255,248,247,254,255;255,249,248,252,248,247,226,165,138,122,125,150,213,255,251,249,247,248,250,248,247,247,247,252;255,247,247,247,249,210,143,136,147,144,131,113,112,173,242,249,249,247,250,249,247,248,250,254;254,247,247,254,223,134,169,239,249,245,229,179,111,84,193,255,250,247,251,251,247,251,255,255;252,247,250,240,148,154,247,255,246,225,214,219,175,90,110,219,248,247,253,250,247,252,255,255;251,247,251,231,122,190,251,247,227,210,201,191,165,118,87,173,247,247,255,248,247,253,255,255;249,247,253,225,104,208,252,228,211,198,192,180,164,126,81,155,247,247,255,247,247,248,253,255;249,247,255,223,93,194,240,217,202,187,188,202,202,126,81,157,248,248,255,252,248,248,255,255;255,255,255,225,104,160,211,205,190,187,213,249,204,113,87,168,255,255,255,255,255,255,255,255;254,249,249,233,129,107,176,202,167,196,231,203,139,86,109,207,255,254,255,255,255,255,255,255;255,250,249,250,208,101,103,137,145,165,165,115,81,81,90,169,252,255,255,255,255,255,255,255;255,252,249,249,252,180,107,90,92,90,90,89,86,82,85,96,184,255,255,255,255,255,255,255;255,254,249,249,252,254,193,125,103,99,107,129,146,99,87,87,96,183,255,255,255,255,255,255;255,255,250,249,251,255,255,245,240,242,245,244,216,142,93,87,87,93,194,249,249,250,253,255;255,255,251,249,249,255,254,245,238,236,235,232,232,217,143,99,87,83,95,191,249,250,253,255;255,255,253,249,249,254,252,249,248,252,251,250,244,246,229,154,98,81,81,85,197,255,255,255;255,255,255,249,249,252,251,249,251,255,255,255,249,249,251,230,146,93,81,81,97,248,255,255;255,255,255,251,249,251,250,249,253,255,255,255,249,249,252,255,230,154,92,84,118,241,251,255;255,255,255,252,249,250,249,249,255,255,255,255,249,249,252,255,255,222,156,144,193,234,248,255;255,255,255,254,249,249,249,251,255,255,255,255,249,249,252,255,252,239,204,197,226,237,250,255;255,255,255,255,250,249,249,253,255,255,255,255,249,249,252,255,255,251,250,249,250,252,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_zoom_24x24(:,:,3) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,251,244,251,255;255,255,255,255,255,255,255,255,254,252,252,254,255,255,255,255,255,255,255,255,237,234,252,255;255,238,237,246,237,235,219,160,128,109,121,150,213,254,245,238,234,237,242,237,234,234,234,247;255,234,234,235,238,196,130,125,136,130,120,113,111,162,229,238,238,234,241,238,234,237,241,252;252,234,235,252,220,121,155,237,247,231,218,179,103,71,185,255,242,234,245,245,234,244,255,255;248,234,241,240,143,141,237,255,242,211,204,219,161,77,109,219,237,234,250,241,234,248,255,255;244,234,245,231,114,177,245,247,220,196,194,189,151,105,87,170,234,234,254,237,234,251,255,255;239,234,250,225,94,194,247,228,201,185,187,179,151,113,71,142,234,235,255,235,234,237,251,255;238,234,254,223,81,180,238,217,191,174,186,202,198,114,70,154,237,237,255,246,237,237,254,255;255,255,255,225,104,160,211,205,190,187,213,249,204,113,87,168,255,255,255,255,255,255,255,255;254,248,248,232,129,107,176,202,166,195,230,203,138,85,109,207,255,254,255,255,255,255,255,255;255,249,248,249,208,101,103,137,144,164,164,115,81,80,89,169,252,255,255,255,255,255,255,255;255,252,248,248,252,180,107,90,91,89,90,89,85,81,85,96,184,255,255,255,255,255,255,255;255,254,248,248,252,254,193,125,102,98,107,129,146,99,87,87,96,183,255,255,255,255,255,255;255,255,249,248,250,255,255,244,239,241,245,244,215,141,93,87,87,93,193,249,248,249,252,255;255,255,251,248,249,255,254,244,237,235,235,232,231,216,142,99,87,83,94,191,249,249,253,255;255,255,253,248,248,254,252,248,247,252,251,250,243,245,229,154,98,80,80,85,197,255,255,255;255,255,255,248,248,252,250,248,251,255,255,255,248,248,250,230,146,92,80,80,96,248,255,255;255,255,255,250,248,251,249,248,252,255,255,255,248,248,252,255,230,154,92,83,117,240,250,255;255,255,255,252,248,249,248,248,255,255,255,255,248,248,252,255,255,222,156,144,192,233,246,255;255,255,255,254,248,248,248,250,255,255,255,255,248,248,252,255,252,238,203,196,225,236,250,255;255,255,255,255,249,248,248,252,255,255,255,255,248,248,252,255,255,251,249,248,249,252,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];

icon_plot = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,254,254,254,255,255;255,255,255,255,171,228,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,254,254,254,255,255;255,255,254,254,69,134,254,254,254,255,254,254,254,255,255,255,255,255,255,254,254,254,254,254,255,255,254,254,254,254,254,254;255,254,254,224,3,39,254,254,254,254,254,254,254,254,255,255,255,254,254,254,254,254,254,254,255,255,254,254,254,254,254,254;255,254,254,128,0,0,200,254,254,254,255,254,254,254,255,255,254,254,254,254,255,254,254,254,255,255,255,254,254,254,255,255;255,254,254,36,0,0,102,254,254,255,255,254,254,254,255,255,254,254,254,255,255,254,254,254,255,255,254,254,254,255,255,255;255,254,191,0,0,0,18,245,254,255,255,254,254,254,255,254,254,254,255,255,255,254,254,254,255,249,217,238,253,255,255,255;255,254,96,0,0,0,0,167,254,255,255,254,254,254,255,254,254,254,255,255,255,254,254,254,255,204,53,166,246,255,255,255;255,254,254,240,0,60,254,254,255,255,255,254,254,254,255,254,254,254,255,255,254,254,254,254,242,138,16,179,251,255,255,255;254,254,254,240,0,60,254,254,255,255,255,254,254,255,255,254,254,254,254,254,254,254,254,254,199,53,95,223,254,254,255,255;254,254,254,240,0,60,254,254,255,255,254,254,254,255,255,255,254,254,254,254,255,254,254,242,138,27,180,251,254,254,255,255;255,255,255,240,0,60,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,254,207,76,89,222,255,255,255,255,255;255,255,255,240,0,60,255,255,255,255,255,255,255,255,255,255,255,238,251,255,255,255,232,128,53,153,244,255,255,255,255,255;221,191,191,180,0,60,255,255,255,253,191,159,157,176,246,255,212,191,195,242,255,238,145,18,112,225,254,255,255,255,255,255;242,191,191,180,0,60,255,255,250,196,86,52,48,73,182,248,191,191,191,225,244,158,48,94,211,252,255,255,255,255,255,255;255,195,191,180,0,60,255,247,192,65,21,74,82,53,90,217,224,191,208,247,187,45,73,202,252,255,255,255,255,255,255,255;255,216,191,180,0,60,245,175,70,42,112,165,177,134,48,130,221,252,253,212,87,55,184,250,255,255,255,255,255,255,255,255;255,234,191,180,0,54,181,53,42,114,178,191,224,222,100,11,75,134,140,102,12,139,241,255,234,208,191,195,208,230,255,255;255,255,195,180,0,43,65,51,151,178,190,191,246,252,217,127,40,26,32,49,128,228,255,221,191,191,191,191,191,221,255,255;255,255,216,180,0,29,28,160,211,191,191,204,255,255,254,238,149,120,117,177,235,254,251,191,191,191,225,234,221,238,255,255;255,255,234,180,0,20,68,211,204,191,191,225,255,255,255,255,190,186,184,235,255,255,238,191,191,191,242,255,255,255,255,255;255,255,251,180,0,23,119,230,191,191,191,246,255,255,255,255,191,191,191,238,255,255,251,195,191,191,191,216,246,255,255,255;255,255,255,200,0,41,204,231,191,191,204,255,255,255,255,255,191,191,191,238,255,246,255,230,195,191,191,191,191,230,255,255;255,255,255,216,0,45,204,221,191,191,225,255,255,255,255,255,191,191,191,238,255,27,93,186,251,225,195,191,191,195,251,255;255,255,255,236,0,27,115,122,115,115,148,153,153,153,153,153,115,115,115,143,153,9,0,0,27,126,216,191,191,191,238,255;255,255,255,240,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,83,182,191,246,255;255,255,255,252,167,97,115,115,115,135,153,153,153,153,153,153,115,115,115,143,153,9,0,0,18,90,162,191,191,221,255,255;255,255,255,255,246,191,191,191,191,246,255,255,255,255,255,255,191,191,191,238,255,21,81,156,208,191,195,208,230,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,240,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_plot(:,:,2) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,248,250,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,240,222,217,250,255;255,255,255,255,171,228,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,232,217,217,255,255;255,245,217,217,68,122,217,219,240,250,230,217,219,245,255,255,255,255,245,230,222,217,217,227,242,255,217,217,217,217,217,232;255,237,217,191,3,33,217,217,219,224,217,217,217,222,255,255,255,232,217,217,217,217,217,217,245,248,217,217,217,217,217,237;255,232,217,110,0,0,179,217,217,230,253,230,217,217,255,255,235,217,217,237,253,237,217,217,253,255,245,217,217,240,255,255;255,227,217,31,0,0,93,217,217,253,255,237,217,222,255,250,217,217,232,255,255,230,217,222,255,255,237,217,217,248,255,255;255,217,163,0,0,0,16,209,227,255,255,230,217,227,255,237,217,217,250,255,255,222,217,230,255,249,198,203,216,255,255,255;250,217,82,0,0,0,0,143,237,255,255,222,217,237,255,227,217,217,255,255,245,217,217,237,255,204,47,142,220,255,255,255;242,217,217,231,0,58,217,217,242,255,255,217,217,240,255,227,217,217,245,253,224,217,217,240,242,138,14,153,219,245,253,255;237,217,217,240,0,57,217,217,250,255,245,217,217,248,255,237,217,217,217,217,227,217,217,244,199,53,82,191,217,217,253,255;230,217,224,240,0,56,217,217,255,255,240,217,217,255,255,255,235,219,219,235,253,217,217,233,138,27,170,219,217,219,255,255;255,255,255,240,0,60,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,254,207,76,89,222,255,255,255,255,255;255,255,255,240,0,60,255,255,255,255,255,255,255,255,255,255,255,248,253,255,255,255,232,128,53,153,244,255,255,255,255,255;242,230,230,216,0,60,255,255,255,253,217,191,189,203,246,255,238,230,232,250,255,238,145,18,112,225,254,255,255,255,255,255;250,230,230,216,0,60,255,255,250,197,104,62,58,79,182,248,230,230,230,243,244,158,48,94,211,252,255,255,255,255,255,255;255,232,230,216,0,60,255,247,192,69,25,89,99,54,90,217,242,230,237,249,187,45,73,202,252,255,255,255,255,255,255,255;255,240,230,216,0,60,245,175,70,46,134,198,204,134,48,130,221,252,253,212,87,55,184,250,255,255,255,255,255,255,255,255;255,247,230,216,0,56,181,53,42,136,214,230,242,222,100,11,90,161,169,106,12,139,241,255,247,237,230,232,237,245,255,255;255,255,232,216,0,47,65,51,156,215,229,230,252,252,217,127,49,32,39,51,128,228,255,242,230,230,230,230,230,242,255,255;255,255,240,216,0,34,28,160,232,230,230,235,255,255,254,238,179,144,141,185,235,254,253,230,230,230,243,247,242,248,255,255;255,255,247,216,0,24,68,211,235,230,230,243,255,255,255,255,229,224,222,245,255,255,248,230,230,230,250,255,255,255,255,255;255,255,253,216,0,28,127,232,230,230,230,252,255,255,255,255,230,230,230,248,255,255,253,232,230,230,230,240,252,255,255,255;255,255,255,224,0,50,223,244,230,230,235,255,255,255,255,255,230,230,230,248,255,246,255,245,232,230,230,230,230,245,255,255;255,255,255,231,0,54,235,242,230,230,243,255,255,255,255,255,230,230,230,248,255,27,93,186,253,243,232,230,230,232,253,255;255,255,255,238,0,32,138,141,138,138,151,153,153,153,153,153,138,138,138,149,153,9,0,0,27,126,217,230,230,230,248,255;255,255,255,240,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100,219,230,252,255;255,255,255,252,180,117,138,138,138,146,153,153,153,153,153,153,138,138,138,149,153,9,0,0,22,108,195,230,230,242,255,255;255,255,255,255,252,230,230,230,230,252,255,255,255,255,255,255,230,230,230,248,255,21,83,171,237,230,232,237,245,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,240,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_plot(:,:,3) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,235,242,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,214,166,153,242,255;255,255,255,255,171,228,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,194,153,153,255,255;255,227,153,153,65,99,153,160,214,242,187,153,160,227,255,255,255,255,227,187,166,153,153,181,221,255,153,153,153,153,153,194;255,207,153,135,2,23,153,153,160,173,153,153,153,166,255,255,255,194,153,153,153,153,153,153,227,235,153,153,153,153,153,207;255,194,153,77,0,0,143,153,153,187,248,187,153,153,255,255,201,153,153,207,248,207,153,153,248,255,227,153,153,214,255,255;255,181,153,23,0,0,78,153,153,248,255,207,153,166,255,242,153,153,194,255,255,187,153,166,255,255,207,153,153,235,255,255;255,153,115,0,0,0,13,148,181,255,255,187,153,181,255,207,153,153,242,255,255,166,153,187,255,249,166,143,152,255,255,255;242,153,58,0,0,0,0,101,207,255,255,166,153,207,255,181,153,153,255,255,227,153,153,207,255,204,38,100,175,255,255,255;221,153,153,214,0,55,153,153,221,255,255,153,153,214,255,181,153,153,227,248,173,153,153,214,242,138,10,108,164,227,248,255;207,153,153,240,0,52,153,153,242,255,227,153,153,235,255,207,153,153,153,153,181,153,153,226,199,53,60,134,153,153,248,255;187,153,173,240,0,49,153,153,255,255,214,153,153,255,255,255,201,160,160,201,248,153,153,215,138,27,152,164,153,160,255,255;255,255,255,240,0,60,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,254,207,76,89,222,255,255,255,255,255;255,255,255,240,0,60,255,255,255,255,255,255,255,255,255,255,255,247,253,255,255,255,232,128,53,153,244,255,255,255,255,255;239,226,226,213,0,60,255,255,255,253,215,188,185,200,246,255,236,226,228,249,255,238,145,18,112,225,254,255,255,255,255,255;249,226,226,213,0,60,255,255,250,197,102,61,57,78,182,248,226,226,226,241,244,158,48,94,211,252,255,255,255,255,255,255;255,228,226,213,0,60,255,247,192,69,25,88,97,54,90,217,240,226,234,249,187,45,73,202,252,255,255,255,255,255,255,255;255,237,226,213,0,60,245,175,70,46,132,195,201,134,48,130,221,252,253,212,87,55,184,250,255,255,255,255,255,255,255,255;255,245,226,213,0,56,181,53,42,133,210,226,240,222,100,11,89,159,166,106,12,139,241,255,245,234,226,228,234,243,255,255;255,255,228,213,0,46,65,51,155,211,225,226,251,252,217,127,48,31,38,50,128,228,255,239,226,226,226,226,226,239,255,255;255,255,237,213,0,34,28,160,229,226,226,232,255,255,254,238,176,142,138,184,235,254,253,226,226,226,241,245,239,247,255,255;255,255,245,213,0,23,68,211,232,226,226,241,255,255,255,255,225,220,218,244,255,255,247,226,226,226,249,255,255,255,255,255;255,255,253,213,0,28,126,232,226,226,226,251,255,255,255,255,226,226,226,247,255,255,253,228,226,226,226,237,251,255,255,255;255,255,255,222,0,49,220,242,226,226,232,255,255,255,255,255,226,226,226,247,255,246,255,243,228,226,226,226,226,243,255,255;255,255,255,229,0,53,232,239,226,226,241,255,255,255,255,255,226,226,226,247,255,27,93,186,253,241,228,226,226,228,253,255;255,255,255,238,0,32,136,139,136,136,151,153,153,153,153,153,136,136,136,148,153,9,0,0,27,126,217,226,226,226,247,255;255,255,255,240,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,98,215,226,251,255;255,255,255,252,179,115,136,136,136,145,153,153,153,153,153,153,136,136,136,148,153,9,0,0,21,106,191,226,226,239,255,255;255,255,255,255,251,226,226,226,226,251,255,255,255,255,255,255,226,226,226,247,255,21,83,169,234,226,228,234,243,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,240,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_plot_24x24(:,:,1) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,186,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,254,254,255,255;255,254,254,39,215,254,255,255,254,254,255,255,255,255,255,254,254,254,255,254,254,254,254,255;255,254,197,0,117,254,254,254,254,254,254,255,255,254,254,254,254,254,254,254,254,254,254,255;255,254,102,0,27,251,254,255,255,254,254,255,255,254,255,255,255,254,255,255,254,255,255,255;255,245,15,0,0,182,254,255,255,254,254,255,254,254,255,255,254,254,255,210,208,251,255,255;255,212,96,0,72,164,255,255,255,254,255,255,254,254,255,255,254,254,246,82,121,250,255,255;254,254,240,0,179,254,255,255,254,254,255,255,254,254,254,254,254,254,175,51,184,254,255,255;254,254,240,0,179,254,255,255,254,254,255,255,255,254,254,255,254,244,94,101,250,254,255,255;255,255,240,0,180,255,255,255,255,255,255,255,255,255,255,255,255,181,62,179,255,255,255,255;249,209,197,0,180,255,255,250,188,173,204,255,243,234,252,255,209,54,113,241,255,255,255,255;255,219,197,0,180,255,237,124,55,50,104,237,212,209,234,222,77,88,226,255,255,255,255,255;255,234,197,0,174,228,108,47,120,125,58,150,220,212,233,109,67,216,255,255,255,255,255,255;255,246,197,0,148,89,46,161,209,215,148,34,126,197,130,35,184,255,255,255,255,255,255,255;255,255,203,0,76,52,185,212,209,231,245,157,59,48,66,161,252,252,231,212,209,216,237,255;255,255,215,0,40,151,246,209,209,240,255,255,189,168,206,254,255,228,209,212,222,216,240,255;255,255,226,0,65,201,234,209,212,255,255,255,209,209,234,255,255,209,209,228,255,255,252,255;255,255,237,0,130,227,224,209,228,255,255,255,209,209,234,231,222,231,209,209,212,234,252,255;255,255,240,0,148,228,216,209,237,255,255,255,209,209,234,195,0,63,148,209,209,209,228,255;255,255,240,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,111,209,209,255;255,255,255,214,125,125,125,134,153,153,153,153,125,125,140,117,0,7,75,162,209,209,228,255;255,255,255,255,219,209,209,237,255,255,255,255,209,209,234,219,157,215,216,209,219,234,252,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_plot_24x24(:,:,2) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,250,243,250,255;255,255,255,186,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,235,232,252,255;255,236,235,38,199,233,244,247,236,233,247,255,255,254,244,236,232,235,241,235,232,232,232,247;255,232,180,0,108,232,232,235,236,232,235,255,254,235,232,236,236,232,240,236,232,235,240,252;252,232,93,0,26,229,232,252,252,232,235,255,241,232,241,255,241,232,244,244,232,243,255,255;247,224,14,0,0,166,238,255,247,232,238,255,232,232,254,255,235,232,249,199,190,243,255,255;243,194,92,0,68,150,244,255,243,232,243,250,232,233,255,249,232,232,245,76,110,245,255,255;238,232,234,0,168,232,247,255,238,232,247,254,232,232,238,233,232,233,175,47,168,235,250,255;236,232,239,0,166,232,252,255,236,232,252,255,247,235,236,249,235,226,94,97,231,235,254,255;255,255,240,0,180,255,255,255,255,255,255,255,255,255,255,255,255,181,62,179,255,255,255,255;253,237,223,0,180,255,255,250,207,196,225,255,250,247,254,255,209,54,113,241,255,255,255,255;255,241,223,0,180,255,237,126,62,57,111,237,238,237,247,222,77,88,226,255,255,255,255,255;255,247,223,0,177,228,108,50,136,142,59,150,238,238,239,109,67,216,255,255,255,255,255,255;255,251,223,0,156,89,46,172,237,239,148,34,126,197,130,35,184,255,255,255,255,255,255,255;255,255,226,0,83,52,185,238,237,246,245,157,67,55,70,161,252,254,246,238,237,240,248,255;255,255,230,0,45,151,251,237,237,249,255,255,214,191,217,254,255,244,237,238,242,240,249,255;255,255,234,0,73,205,247,237,238,255,255,255,237,237,247,255,255,237,237,244,255,255,254,255;255,255,239,0,147,239,243,237,244,255,255,255,237,237,247,231,222,246,237,237,238,247,254,255;255,255,240,0,167,244,240,237,248,255,255,255,237,237,247,195,0,63,155,230,237,237,244,255;255,255,240,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,119,237,237,255;255,255,255,217,142,142,142,146,153,153,153,153,142,142,148,117,0,8,84,177,237,237,244,255;255,255,255,255,241,237,237,248,255,255,255,255,237,237,247,219,158,230,240,237,241,247,254,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_plot_24x24(:,:,3) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,243,222,243,255;255,255,255,186,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,202,193,246,255;255,205,202,35,171,197,226,234,205,197,234,255,255,251,226,205,193,202,218,202,193,193,193,234;255,193,150,0,94,193,193,202,205,193,202,255,251,202,193,205,205,193,214,205,193,202,214,246;246,193,79,0,25,191,193,246,246,193,202,255,218,193,218,255,218,193,226,226,193,222,255,255;234,186,13,0,0,139,209,255,234,193,209,255,193,193,251,255,202,193,239,177,158,230,255,255;222,161,85,0,62,125,226,255,222,193,222,243,193,197,255,239,193,193,242,65,92,238,255,255;209,193,225,0,148,193,234,255,209,193,234,251,193,193,209,197,193,197,175,39,140,202,243,255;205,193,236,0,143,193,246,255,205,193,246,255,234,202,205,239,202,194,94,91,199,202,251,255;255,255,240,0,180,255,255,255,255,255,255,255,255,255,255,255,255,181,62,179,255,255,255,255;252,234,220,0,180,255,255,250,204,194,222,255,249,245,254,255,209,54,113,241,255,255,255,255;255,238,220,0,180,255,237,126,61,56,110,237,236,234,245,222,77,88,226,255,255,255,255,255;255,245,220,0,177,228,108,49,134,140,59,150,236,236,238,109,67,216,255,255,255,255,255,255;255,251,220,0,155,89,46,172,234,236,148,34,126,197,130,35,184,255,255,255,255,255,255,255;255,255,223,0,82,52,185,236,234,244,245,157,66,54,69,161,252,254,244,236,234,237,247,255;255,255,229,0,45,151,251,234,234,248,255,255,211,188,215,254,255,243,234,236,240,237,248,255;255,255,233,0,72,205,245,234,236,255,255,255,234,234,245,255,255,234,234,243,255,255,254,255;255,255,239,0,145,237,241,234,243,255,255,255,234,234,245,231,222,244,234,234,236,245,254,255;255,255,240,0,165,243,237,234,247,255,255,255,234,234,245,195,0,63,154,227,234,234,243,255;255,255,240,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,12,118,234,234,255;255,255,255,216,140,140,140,145,153,153,153,153,140,140,147,117,0,8,83,175,234,234,243,255;255,255,255,255,238,234,234,247,255,255,255,255,234,234,245,219,158,229,237,234,238,245,254,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];

icon_roi = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,204,134,75,29,3,29,75,134,204,255,255,255,255,255,255,255,255,255,255,255,255,255,254,254,255,255;255,255,255,255,245,138,51,121,180,226,252,226,180,121,51,138,245,255,255,255,255,255,255,255,255,255,255,255,254,254,255,255;255,255,254,254,122,117,219,254,255,255,254,254,254,255,220,117,122,255,255,254,254,254,254,254,255,255,254,254,254,254,254,255;255,255,254,163,91,254,254,254,254,254,254,254,254,254,255,255,91,164,254,254,254,254,254,254,255,255,254,254,254,254,254,255;255,255,254,85,170,255,254,254,254,254,255,254,254,254,255,255,170,85,254,255,255,255,254,254,255,255,255,254,254,255,255,255;255,254,254,25,230,255,255,254,254,255,255,255,254,254,255,255,229,25,255,255,255,254,254,254,255,255,255,254,254,255,255,255;255,254,254,25,230,255,254,254,254,255,255,254,254,254,255,255,229,25,255,255,255,254,254,254,255,255,255,254,254,255,255,255;255,254,254,85,170,255,254,254,255,255,255,254,254,255,255,254,169,85,255,255,255,121,33,3,3,33,121,254,254,255,255,255;255,254,254,164,91,255,254,254,255,255,255,254,254,255,255,254,91,163,255,255,84,24,104,202,202,104,24,84,254,255,255,255;255,254,254,255,122,117,219,254,255,255,255,254,254,255,220,117,122,254,254,115,28,214,254,255,255,255,214,28,115,254,255,255;254,254,254,255,245,138,51,121,180,226,252,225,179,121,51,138,245,254,225,21,128,254,254,255,255,255,255,127,21,225,255,255;255,255,255,255,255,255,204,134,75,29,3,29,75,134,204,255,255,255,115,28,255,255,255,255,255,255,255,255,28,115,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,243,47,81,255,255,255,255,255,255,255,255,81,48,255,255;231,209,209,209,249,255,255,255,255,255,222,209,209,219,255,255,224,209,12,143,255,255,255,255,255,255,255,255,148,15,255,255;246,209,209,209,237,255,255,255,255,252,209,209,209,234,255,255,209,209,2,194,255,255,255,255,255,255,255,255,211,2,255,255;255,212,209,209,224,255,255,255,255,237,209,209,209,249,255,255,234,209,0,248,255,255,255,255,255,255,255,255,251,0,255,255;255,228,209,209,212,255,255,255,255,228,209,209,219,255,255,255,255,255,2,211,255,255,255,255,255,255,255,255,211,2,255,255;255,240,209,209,209,243,255,255,255,212,209,209,234,255,255,255,209,209,12,141,255,255,255,255,240,222,209,212,129,14,255,255;255,255,235,235,235,235,235,235,235,235,235,235,235,235,235,235,209,209,39,77,255,255,255,231,209,209,209,209,66,43,255,255;255,255,235,209,209,219,255,255,231,209,209,219,255,255,255,235,209,209,94,27,255,255,252,209,209,209,234,240,25,110,255,255;255,255,235,209,209,209,252,255,219,209,209,234,255,255,255,235,209,209,185,20,128,255,243,209,209,209,246,128,21,226,255,255;255,255,235,209,209,209,237,252,209,209,209,249,255,255,255,235,209,209,209,110,28,215,252,212,209,209,176,25,112,255,255,255;255,255,235,224,209,209,231,240,209,209,219,255,255,255,255,235,209,209,209,243,84,24,104,188,168,85,20,69,209,237,255,255;255,255,235,237,209,209,219,231,209,209,234,255,255,255,255,235,209,209,209,243,255,121,33,3,3,30,101,209,209,212,252,255;255,255,235,252,209,209,209,219,209,209,249,255,255,255,255,235,209,209,209,243,255,255,255,255,255,255,252,209,209,209,243,255;255,255,235,235,235,235,235,235,235,235,235,235,235,235,235,235,209,209,209,243,255,255,252,216,231,240,234,209,209,209,249,255;255,255,255,255,234,209,209,209,209,234,255,255,255,255,255,255,209,209,209,243,255,255,240,209,209,209,209,209,209,231,255,255;255,255,255,255,249,209,209,209,209,249,255,255,255,255,255,255,209,209,209,243,255,255,249,231,222,209,212,222,237,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_roi(:,:,2) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,250,252,255;255,255,255,255,255,255,205,137,79,34,9,34,79,137,205,255,255,255,255,255,255,255,255,255,255,255,255,246,235,232,252,255;255,255,255,255,245,141,56,124,182,227,252,227,182,124,56,141,245,255,255,255,255,255,255,255,255,255,255,241,232,232,255,255;255,249,232,232,124,113,201,233,246,252,240,232,233,249,221,120,125,255,249,240,235,232,232,238,247,255,232,232,232,232,232,241;255,244,232,151,87,232,232,232,233,236,232,232,232,235,255,255,95,157,232,232,232,232,232,232,249,250,232,232,232,232,232,244;255,241,232,81,163,255,238,232,232,240,254,240,232,232,255,255,164,81,232,244,254,244,232,232,254,255,249,232,232,246,255,255;255,238,232,28,230,255,241,232,232,254,255,244,232,235,255,252,210,28,241,255,255,240,232,235,255,255,244,232,232,250,255,255;255,232,232,29,231,255,238,232,238,255,255,240,232,238,255,244,210,28,252,255,255,235,232,240,255,255,241,232,232,255,255,255;252,232,232,86,172,255,232,232,244,255,255,235,232,244,255,238,157,81,255,255,249,203,177,162,162,180,206,232,238,255,255,255;247,232,232,162,95,250,232,232,247,255,255,232,232,246,255,238,87,151,249,254,194,174,199,232,239,208,174,193,235,249,254,255;244,232,232,255,125,117,201,232,252,255,249,232,232,250,221,115,114,232,232,201,176,224,232,249,255,255,225,175,201,232,254,255;240,232,236,255,245,135,51,113,182,227,243,206,166,124,56,141,234,233,227,173,216,232,232,249,255,255,246,206,172,227,255,255;255,255,255,255,255,255,205,137,79,34,9,34,79,137,205,255,255,255,212,178,255,255,255,255,255,255,255,255,178,212,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,250,187,200,255,255,255,255,255,255,255,255,200,187,255,255;246,237,237,237,253,255,255,255,255,255,242,237,237,241,255,255,243,237,170,220,255,255,255,255,255,255,255,255,223,171,255,255;251,237,237,237,248,255,255,255,255,254,237,237,237,247,255,255,237,237,161,235,255,255,255,255,255,255,255,255,242,161,255,255;255,238,237,237,243,255,255,255,255,248,237,237,237,253,255,255,247,237,156,253,255,255,255,255,255,255,255,255,254,156,255,255;255,244,237,237,238,255,255,255,255,244,237,237,241,255,255,255,255,255,161,242,255,255,255,255,255,255,255,255,242,161,255,255;255,249,237,237,237,250,255,255,255,238,237,237,247,255,255,255,237,237,170,220,255,255,255,255,249,242,237,238,215,171,255,255;255,255,0,0,0,0,0,0,0,0,0,0,0,0,0,0,237,237,183,198,255,255,255,246,237,237,237,237,194,185,255,255;255,255,0,237,237,241,255,255,246,237,237,241,255,255,255,0,237,237,204,177,255,255,254,237,237,237,247,249,177,209,255,255;255,255,0,237,237,237,254,255,241,237,237,247,255,255,255,0,237,237,231,174,216,255,250,237,237,237,251,216,174,246,255,255;255,255,0,237,237,237,248,254,237,237,237,253,255,255,255,0,237,237,237,209,178,243,254,238,237,237,228,176,211,255,255,255;255,255,0,243,237,237,246,249,237,237,241,255,255,255,255,0,237,237,237,250,201,176,208,234,226,201,174,195,237,248,255,255;255,255,0,248,237,237,241,246,237,237,247,255,255,255,255,0,237,237,237,250,255,214,180,162,162,179,206,237,237,238,254,255;255,255,0,254,237,237,237,241,237,237,253,255,255,255,255,0,237,237,237,250,255,255,255,255,255,255,254,237,237,237,250,255;255,255,0,0,0,0,0,0,0,0,0,0,0,0,0,0,237,237,237,250,255,255,254,240,246,249,247,237,237,237,253,255;255,255,255,255,247,237,237,237,237,247,255,255,255,255,255,255,237,237,237,250,255,255,249,237,237,237,237,237,237,246,255,255;255,255,255,255,253,237,237,237,237,253,255,255,255,255,255,255,237,237,237,250,255,255,253,246,242,237,238,242,248,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_roi(:,:,3) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,243,246,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,230,202,193,246,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,218,193,193,255,255;255,239,193,193,251,236,202,197,230,246,214,193,197,239,255,255,255,255,239,214,202,193,193,209,234,255,193,193,193,193,193,218;255,226,193,215,234,193,193,193,197,205,193,193,193,202,255,255,255,231,193,193,193,193,193,193,239,243,193,193,193,193,193,226;255,218,193,234,230,255,209,193,193,214,251,214,193,193,255,255,233,234,193,226,251,226,193,193,251,255,239,193,193,230,255,255;255,209,193,249,251,255,218,193,193,251,255,226,193,202,255,246,199,249,218,255,255,214,193,202,255,255,226,193,193,243,255,255;255,193,193,251,255,255,209,193,209,255,255,214,193,209,255,226,199,249,246,255,255,202,193,214,255,255,218,193,193,255,255,255;246,193,193,247,255,255,193,193,226,255,255,202,193,226,255,209,214,234,255,255,239,112,56,36,37,64,119,193,209,255,255,255;234,193,193,245,255,243,193,193,234,255,255,193,193,230,255,209,233,215,239,251,92,51,101,190,210,126,51,88,202,239,251,255;226,193,193,255,255,245,202,193,246,255,239,193,193,243,255,242,225,193,193,107,55,169,193,239,255,255,172,53,107,193,251,255;214,193,205,255,255,239,243,226,255,255,230,200,211,255,255,255,223,197,179,51,145,193,193,239,255,255,230,120,49,179,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,135,60,255,255,255,255,255,255,255,255,60,135,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,249,77,106,255,255,255,255,255,255,255,255,106,77,255,255;244,234,234,234,252,255,255,255,255,255,240,234,234,238,255,255,241,234,47,162,255,255,255,255,255,255,255,255,164,48,255,255;251,234,234,234,247,255,255,255,255,254,234,234,234,245,255,255,234,234,36,209,255,255,255,255,255,255,255,255,218,36,255,255;255,236,234,234,241,255,255,255,255,247,234,234,234,252,255,255,245,234,33,251,255,255,255,255,255,255,255,255,252,33,255,255;255,243,234,234,236,255,255,255,255,243,234,234,238,255,255,255,255,255,36,218,255,255,255,255,255,255,255,255,218,36,255,255;255,248,234,234,234,249,255,255,255,236,234,234,245,255,255,255,234,234,47,160,255,255,255,255,248,240,234,236,155,47,255,255;255,255,0,0,0,0,0,0,0,0,0,0,0,0,0,0,234,234,73,104,255,255,255,244,234,234,234,234,100,75,255,255;255,255,0,234,234,238,255,255,244,234,234,238,255,255,255,0,234,234,126,59,255,255,254,234,234,234,245,248,59,133,255,255;255,255,0,234,234,234,254,255,238,234,234,245,255,255,255,0,234,234,212,54,147,255,249,234,234,234,251,147,54,230,255,255;255,255,0,234,234,234,247,254,234,234,234,252,255,255,255,0,234,234,234,133,60,221,254,236,234,234,203,59,134,255,255,255;255,255,0,241,234,234,244,248,234,234,238,255,255,255,255,0,234,234,234,249,109,57,126,204,195,117,55,102,234,247,255,255;255,255,0,247,234,234,238,244,234,234,245,255,255,255,255,0,234,234,234,249,255,141,64,37,37,63,132,234,234,236,254,255;255,255,0,254,234,234,234,238,234,234,252,255,255,255,255,0,234,234,234,249,255,255,255,255,255,255,254,234,234,234,249,255;255,255,0,0,0,0,0,0,0,0,0,0,0,0,0,0,234,234,234,249,255,255,254,237,244,248,245,234,234,234,252,255;255,255,255,255,245,234,234,234,234,245,255,255,255,255,255,255,234,234,234,249,255,255,248,234,234,234,234,234,234,244,255,255;255,255,255,255,252,234,234,234,234,252,255,255,255,255,255,255,234,234,234,249,255,255,252,244,240,234,236,240,247,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_roi_24x24(:,:,1) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,212,142,100,70,75,110,157,232,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,229,115,106,147,169,166,140,97,141,248,255,255,255,255,255,255,255,254,254,255,255;255,254,222,118,175,246,255,255,254,254,236,150,139,223,255,254,254,254,255,254,254,254,254,255;255,254,129,149,254,254,254,254,254,254,254,255,118,129,254,254,254,254,254,254,254,254,254,255;255,254,81,214,255,254,254,255,255,254,254,255,162,81,255,255,255,254,255,255,254,255,255,255;255,254,81,214,255,254,254,255,255,254,254,255,161,81,255,213,78,52,52,78,212,255,255,255;255,254,130,150,255,254,255,255,255,254,255,255,118,129,255,147,39,76,76,39,146,255,255,255;254,254,223,118,175,246,255,255,254,254,236,150,138,217,158,52,199,253,254,199,52,158,249,255;254,254,255,229,115,106,147,169,165,139,97,141,248,178,56,161,254,254,255,255,160,56,179,255;255,255,255,255,212,142,100,70,75,110,157,232,255,82,49,246,255,255,255,255,246,49,82,255;249,209,209,228,255,255,255,255,219,209,219,255,243,50,93,255,255,255,255,255,255,94,54,255;255,219,209,216,255,255,255,249,209,209,231,255,212,48,139,255,255,255,255,255,255,152,58,255;255,234,209,209,246,255,255,237,209,209,246,255,224,51,154,255,255,255,255,255,255,162,61,255;255,246,209,209,234,255,255,228,209,216,255,255,255,53,108,255,255,255,255,255,255,108,53,255;255,239,232,228,232,240,240,229,228,234,241,238,209,56,51,249,255,252,231,212,204,47,63,255;255,240,230,209,212,255,246,209,209,240,255,240,209,123,45,196,255,228,209,212,171,42,141,255;255,240,238,209,209,246,234,209,212,255,255,240,209,197,114,68,243,209,209,217,68,124,237,255;255,240,247,209,209,234,224,209,228,255,255,240,209,209,232,101,49,115,104,40,84,232,252,255;255,240,249,224,209,228,216,209,237,255,255,240,209,209,234,213,78,52,48,67,175,209,228,255;255,239,238,235,228,230,228,228,240,240,241,238,209,209,234,255,255,249,255,255,231,209,209,255;255,255,255,249,209,209,209,224,255,255,255,255,209,209,234,255,255,212,212,222,209,209,228,255;255,255,255,255,219,209,209,237,255,255,255,255,209,209,234,255,252,228,216,209,219,234,252,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_roi_24x24(:,:,2) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,213,145,104,74,79,113,159,233,255,255,255,255,255,255,255,255,250,243,250,255;255,255,255,230,118,110,150,171,168,143,101,144,248,255,255,255,255,255,255,255,235,232,252,255;255,236,206,117,164,226,244,247,236,233,229,152,142,223,244,236,232,235,241,235,232,232,232,247;255,232,121,140,236,232,232,235,236,232,235,255,121,123,232,236,236,232,240,236,232,235,240,252;252,232,78,212,249,232,232,252,252,232,235,255,155,78,241,255,241,232,244,244,232,243,255,255;247,232,80,215,246,232,238,255,247,232,238,255,150,78,254,242,192,178,181,193,223,247,255,255;243,232,127,152,241,232,244,255,243,232,243,250,111,122,255,217,182,189,196,182,207,250,255,255;238,232,219,121,166,225,247,255,238,232,229,152,129,203,215,184,220,233,255,221,183,213,248,255;236,232,254,230,109,100,148,171,156,130,100,144,240,218,184,222,235,235,255,246,214,184,231,255;255,255,255,255,213,145,104,74,79,113,159,233,255,199,188,252,255,255,255,255,252,188,199,255;253,237,237,244,255,255,255,255,241,237,241,255,250,185,204,255,255,255,255,255,255,204,186,255;255,241,237,240,255,255,255,253,237,237,246,255,238,179,218,255,255,255,255,255,255,223,183,255;255,247,237,237,251,255,255,248,237,237,251,255,243,179,223,255,255,255,255,255,255,226,183,255;255,251,237,237,247,255,255,244,237,240,255,255,255,185,209,255,255,255,255,255,255,209,185,255;255,48,34,68,67,70,70,65,65,68,72,33,237,189,189,253,255,254,246,238,236,188,192,255;255,64,166,237,238,255,251,237,237,249,255,70,237,212,185,237,255,244,237,238,227,183,219,255;255,64,170,237,237,251,247,237,238,255,255,70,237,234,210,194,251,237,237,241,194,214,250,255;255,64,173,237,237,247,243,237,244,255,255,70,237,237,246,205,188,211,206,184,198,246,254,255;255,64,174,243,237,244,240,237,248,255,255,70,237,237,247,242,198,183,181,194,227,237,244,255;255,48,36,71,65,66,65,65,70,70,72,33,237,237,247,255,255,253,255,255,246,237,237,255;255,255,255,253,237,237,237,243,255,255,255,255,237,237,247,255,255,238,238,242,237,237,244,255;255,255,255,255,241,237,237,248,255,255,255,255,237,237,247,255,254,244,240,237,241,247,254,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];
icon_roi_24x24(:,:,3) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,243,222,243,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,202,193,246,255;255,205,209,243,218,199,226,234,205,197,236,255,255,252,226,205,193,202,218,202,193,193,193,234;255,193,223,221,205,193,193,202,205,193,202,255,253,228,193,205,205,193,214,205,193,202,214,246;246,193,237,247,239,193,193,246,246,193,202,255,231,235,218,255,218,193,226,226,193,222,255,255;234,193,242,255,230,193,209,255,234,193,209,255,216,235,251,219,87,67,77,91,167,234,255,255;222,193,240,255,218,193,226,255,222,193,222,243,226,225,255,153,61,83,100,62,127,243,255,255;209,193,241,255,223,195,234,255,209,193,236,253,221,197,144,69,159,196,254,162,68,140,238,255;205,193,251,255,231,229,250,255,222,221,252,255,235,153,74,165,202,202,255,230,142,73,188,255;255,255,255,255,255,255,255,255,255,255,255,255,255,107,78,247,255,255,255,255,247,78,107,255;252,234,234,243,255,255,255,255,238,234,238,255,249,79,117,255,255,255,255,255,255,117,82,255;255,238,234,237,255,255,255,252,234,234,244,255,236,80,161,255,255,255,255,255,255,167,85,255;255,245,234,234,251,255,255,247,234,234,251,255,241,82,172,255,255,255,255,255,255,175,87,255;255,251,234,234,245,255,255,243,234,237,255,255,255,81,129,255,255,255,255,255,255,129,81,255;255,48,33,67,66,70,70,65,64,67,72,33,234,89,83,250,255,254,244,236,229,81,92,255;255,64,166,234,236,255,251,234,234,248,255,70,234,153,77,205,255,243,234,236,193,75,162,255;255,64,169,234,234,251,245,234,236,255,255,70,234,222,138,95,245,234,234,233,95,143,241,255;255,64,173,234,234,245,241,234,243,255,255,70,234,234,243,123,79,140,135,75,116,243,254,255;255,64,174,241,234,243,237,234,247,255,255,70,234,234,245,219,104,80,78,98,202,234,243,255;255,48,36,70,64,65,64,64,70,70,72,33,234,234,245,255,255,252,255,255,244,234,234,255;255,255,255,252,234,234,234,241,255,255,255,255,234,234,245,255,255,236,236,240,234,234,243,255;255,255,255,255,238,234,234,247,255,255,255,255,234,234,245,255,254,243,237,234,238,245,254,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255];

icon_profile = icon_roi;
icon_profile_24x24 = icon_roi_24x24;

icon_tooltips_24x24(:,:,1) = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,220,102,0;0,0,0,0,0,0,119,68,0,0,0,0,0,0,0,0,0,0,0,0,0,220,102,0;0,0,0,0,0,0,186,119,0,0,0,0,0,0,0,0,0,0,0,0,0,220,102,0;0,0,0,0,0,119,254,254,102,102,237,254,237,102,0,102,237,254,237,102,0,220,102,0;0,0,0,0,0,0,186,119,17,254,102,0,119,237,17,254,102,0,119,237,17,220,102,0;0,0,0,0,0,0,186,119,68,254,0,0,34,254,102,254,0,0,34,254,34,220,102,0;0,0,0,0,0,0,186,119,17,254,102,0,119,237,17,254,102,0,119,237,17,220,102,0;0,0,0,0,0,0,152,254,119,102,237,254,237,102,0,102,237,254,237,102,0,220,102,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;255,255,255,255,255,255,255,34,136,255,68,0,0,0,0,0,0,0,0,0,0,0,0,0;153,153,204,255,170,153,153,17,102,153,34,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,170,255,34,0,0,0,221,255,0,17,255,255,255,255,187,0,0,170,255,255,255,119;0,0,238,238,0,0,0,17,255,204,0,85,255,255,153,221,255,85,85,255,204,153,255,221;0,17,255,187,0,0,0,51,255,170,0,119,255,119,0,119,255,102,85,255,255,255,255,170;0,68,255,136,0,0,0,119,255,119,0,153,255,51,0,136,255,85,0,68,153,221,255,187;0,119,255,119,0,0,0,136,255,85,0,204,255,102,68,238,238,17,221,238,17,85,255,170;0,136,255,68,0,0,0,187,255,34,0,255,255,255,255,255,85,0,204,255,255,255,255,51;0,102,153,17,0,0,0,119,136,0,51,255,204,153,153,51,0,0,17,153,153,153,34,0;0,0,0,0,0,0,0,0,0,0,85,255,119,0,0,0,0,0,0,0,0,0,0,0];
icon_tooltips_24x24(:,:,2) = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,220,102,0;0,0,0,0,0,0,119,68,0,0,0,0,0,0,0,0,0,0,0,0,0,220,102,0;0,0,0,0,0,0,186,119,0,0,0,0,0,0,0,0,0,0,0,0,0,220,102,0;0,0,0,0,0,119,254,254,102,102,237,254,237,102,0,102,237,254,237,102,0,220,102,0;0,0,0,0,0,0,186,119,17,254,102,0,119,237,17,254,102,0,119,237,17,220,102,0;0,0,0,0,0,0,186,119,68,254,0,0,34,254,102,254,0,0,34,254,34,220,102,0;0,0,0,0,0,0,186,119,17,254,102,0,119,237,17,254,102,0,119,237,17,220,102,0;0,0,0,0,0,0,152,254,119,102,237,254,237,102,0,102,237,254,237,102,0,220,102,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;64,64,64,64,64,64,64,9,34,64,17,0,0,0,0,0,0,0,0,0,0,0,0,0;38,38,51,64,43,38,38,4,26,38,9,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,43,64,9,0,0,0,55,64,0,4,64,64,64,64,47,0,0,43,64,64,64,30;0,0,60,60,0,0,0,4,64,51,0,21,64,64,38,55,64,21,21,64,51,38,64,55;0,4,64,47,0,0,0,13,64,43,0,30,64,30,0,30,64,26,21,64,64,64,64,43;0,17,64,34,0,0,0,30,64,30,0,38,64,13,0,34,64,21,0,17,38,55,64,47;0,30,64,30,0,0,0,34,64,21,0,51,64,26,17,60,60,4,55,60,4,21,64,43;0,34,64,17,0,0,0,47,64,9,0,64,64,64,64,64,21,0,51,64,64,64,64,13;0,26,38,4,0,0,0,30,34,0,13,64,51,38,38,13,0,0,4,38,38,38,9,0;0,0,0,0,0,0,0,0,0,0,21,64,30,0,0,0,0,0,0,0,0,0,0,0];
icon_tooltips_24x24(:,:,3) = [255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,254,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,254,255,255;255,255,255,255,255,255,254,255,255,255,255,255,255,255,255,255,255,255,255,255,255,254,255,255;255,255,255,255,255,255,254,254,255,255,254,254,254,255,255,255,254,254,254,255,255,254,255,255;255,255,255,255,255,255,254,255,255,254,255,255,255,254,255,254,255,255,255,254,255,254,255,255;255,255,255,255,255,255,254,255,255,254,255,255,255,254,255,254,255,255,255,254,255,254,255,255;255,255,255,255,255,255,254,255,255,254,255,255,255,254,255,254,255,255,255,254,255,254,255,255;255,255,255,255,255,255,254,254,255,255,254,254,254,255,255,255,254,254,254,255,255,254,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255;64,64,64,64,64,64,64,230,153,64,204,255,255,255,255,255,255,255,255,255,255,255,255,255;140,140,102,64,128,140,140,242,179,140,230,255,255,255,255,255,255,255,255,255,255,255,255,255;255,255,128,64,230,255,255,255,89,64,255,242,64,64,64,64,115,255,255,128,64,64,64,166;255,255,77,77,255,255,255,242,64,102,255,191,64,64,140,89,64,191,191,64,102,140,64,89;255,242,64,115,255,255,255,217,64,128,255,166,64,166,255,166,64,179,191,64,64,64,64,128;255,204,64,153,255,255,255,166,64,166,255,140,64,217,255,153,64,191,255,204,140,89,64,115;255,166,64,166,255,255,255,153,64,191,255,102,64,179,204,77,77,242,89,77,242,191,64,128;255,153,64,204,255,255,255,115,64,230,255,64,64,64,64,64,191,255,102,64,64,64,64,217;255,179,140,242,255,255,255,166,153,255,217,64,102,140,140,217,255,255,242,140,140,140,230,255;255,255,255,255,255,255,255,255,255,255,191,64,166,255,255,255,255,255,255,255,255,255,255,255];

icon_update_24x24(:,:,1) = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,36,96,36,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,155,179,24,0;0,143,155,71,155,167,78,32,26,0,0,0,0,0,0,7,75,149,108,155,179,179,179,59;0,179,179,167,143,101,25,0,0,0,0,0,0,0,0,0,0,19,97,143,179,155,120,24;24,179,167,24,26,9,0,1,8,0,0,0,0,0,0,0,52,162,83,83,179,96,0,0;59,179,120,5,6,0,15,5,45,0,0,0,0,32,7,0,155,179,47,120,179,59,0,0;96,179,78,5,0,19,62,0,28,0,0,11,133,167,0,47,179,179,12,155,179,36,0,0;131,179,26,0,1,86,59,0,45,76,51,12,179,179,131,167,179,167,0,167,179,155,36,0;143,179,9,0,24,179,24,0,143,179,24,0,59,155,143,47,155,155,0,71,153,155,12,0;0,0,4,0,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,5,2,0,0;11,81,1,0,4,0,0,0,65,81,65,0,22,38,6,0,0,0,0,3,4,6,0,0;0,65,1,0,1,0,0,11,81,81,43,0,75,81,38,0,0,0,0,1,0,2,0,0;0,38,1,0,16,0,0,32,81,81,16,0,54,75,22,0,0,0,0,1,0,1,0,0;0,16,1,0,29,0,0,49,81,70,0,0,0,0,0,0,0,0,0,4,0,1,0,0;0,0,7,0,25,0,0,75,81,43,0,0,81,81,38,0,0,6,43,33,0,7,32,0;0,0,15,0,15,0,16,81,81,27,0,0,81,81,38,0,0,49,81,15,0,18,27,0;0,0,17,0,1,13,38,81,75,0,0,0,81,81,38,0,0,81,42,1,0,6,6,0;0,0,7,14,0,8,41,81,49,0,0,0,81,81,38,0,0,34,12,0,14,36,6,0;0,0,0,34,6,0,11,42,32,0,0,0,81,81,38,0,5,6,0,6,48,81,49,0;0,0,0,38,40,6,0,1,7,7,4,1,74,57,20,6,1,0,4,6,43,81,81,0;0,0,0,11,81,48,14,0,0,0,0,0,0,0,0,0,0,14,45,59,81,81,49,0;0,0,0,0,65,81,76,20,7,4,1,1,1,1,5,7,8,46,70,81,65,38,6,0;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
icon_update_24x24(:,:,2) = [54,54,54,54,54,54,54,54,162,135,66,54,54,54,54,54,54,54,54,54,54,54,54,54;54,54,54,54,54,54,54,54,178,230,222,158,81,54,54,54,54,54,54,54,70,96,70,54;54,54,54,54,54,54,54,54,95,230,230,230,230,178,100,54,54,54,54,54,122,132,64,54;54,117,122,85,122,127,101,162,211,230,230,230,230,230,230,225,189,126,101,122,132,132,132,80;54,132,132,127,117,176,217,230,230,230,230,230,230,230,230,230,230,220,129,117,132,122,106,64;64,132,127,64,162,226,230,230,211,230,230,230,230,230,230,217,168,141,90,90,132,96,54,54;80,132,106,134,225,230,219,148,115,230,230,230,230,213,125,60,122,132,75,106,132,80,54,54;96,132,101,209,230,221,131,54,191,230,230,181,157,127,54,75,132,132,59,122,132,70,54,54;111,132,160,230,230,184,80,54,189,188,101,59,132,132,111,127,132,127,54,127,132,122,70,54;117,132,199,230,215,132,64,54,117,132,64,54,80,122,117,75,122,122,54,85,124,122,59,54;54,54,221,230,166,54,54,54,54,54,54,54,54,54,54,54,54,54,54,69,126,79,54,54;68,162,229,230,115,54,54,54,140,162,140,54,83,104,61,54,54,54,54,91,222,199,54,54;54,140,229,230,72,54,54,68,162,162,112,54,155,162,104,54,54,54,54,71,230,228,54,54;54,104,229,230,91,54,54,97,162,162,75,54,126,155,83,54,54,54,54,72,230,228,54,54;54,75,229,230,147,54,54,119,162,148,54,54,54,54,54,54,54,54,54,115,230,228,54,54;54,54,225,230,191,54,54,155,162,112,54,54,162,162,104,54,54,61,112,201,230,225,97,54;54,54,209,230,219,54,75,162,162,90,54,54,162,162,104,54,54,119,162,219,230,214,90,54;54,54,166,230,230,158,104,162,155,54,54,54,162,162,104,54,54,162,197,230,230,150,61,54;54,54,74,221,230,219,156,162,119,54,54,54,162,162,104,54,54,146,223,230,220,113,61,54;54,54,54,173,227,230,222,197,97,54,54,54,162,162,104,54,148,215,230,227,192,162,119,54;54,54,54,104,199,227,230,230,208,166,115,72,169,185,183,207,230,230,224,153,112,162,162,54;54,54,54,68,162,192,221,230,230,230,230,230,230,230,230,230,230,220,188,133,162,162,119,54;54,54,54,54,140,162,167,169,199,221,228,228,229,229,223,199,154,127,148,162,140,104,61,54;54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54,54];
icon_update_24x24(:,:,3) = [207,207,207,207,207,207,207,207,80,112,192,207,207,207,207,207,207,207,207,207,207,207,207,207;207,207,207,207,207,207,207,207,61,0,10,85,175,207,207,207,207,207,207,207,179,133,179,207;207,207,207,207,207,207,207,207,158,0,0,0,0,61,153,207,207,207,207,207,87,69,189,207;207,96,87,152,87,78,133,79,22,0,0,0,0,0,0,6,32,84,124,87,69,69,69,161;207,69,69,78,96,44,16,0,0,0,0,0,0,0,0,0,0,7,93,96,69,87,114,189;189,69,78,189,81,7,0,1,26,0,0,0,0,0,0,16,60,63,143,143,69,133,207,207;161,69,114,118,9,0,15,102,123,0,0,0,0,12,121,200,87,69,170,114,69,161,207,207;133,69,133,31,0,13,105,207,39,0,0,55,51,78,207,170,69,69,198,87,69,179,207,207;105,69,83,0,1,39,161,207,36,29,138,198,69,69,105,78,69,78,207,78,69,87,179,207;96,69,42,0,19,69,189,207,96,69,189,207,161,87,96,170,87,87,207,152,86,87,198,207;207,207,14,0,81,207,207,207,207,207,207,207,207,207,207,207,207,207,207,190,127,179,207,207;205,193,3,0,140,207,207,207,196,193,196,207,203,200,206,207,207,207,207,166,13,43,207,207;207,196,3,0,187,207,207,205,193,193,199,207,194,193,200,207,207,207,207,188,0,5,207,207;207,200,3,0,184,207,207,201,193,193,204,207,198,194,203,207,207,207,207,187,0,3,207,207;207,204,3,0,135,207,207,198,193,195,207,207,207,207,207,207,207,207,207,140,0,3,207,207;207,207,13,0,78,207,207,194,193,199,207,207,193,193,200,207,207,206,199,77,0,13,201,207;207,207,42,0,31,207,204,193,193,202,207,207,193,193,200,207,207,198,193,31,0,41,202,207;207,207,97,0,1,100,200,193,194,207,207,207,193,193,200,207,207,193,95,1,0,99,206,207;207,207,191,29,0,22,143,193,198,207,207,207,193,193,200,207,207,144,22,0,29,186,206,207;207,207,207,113,10,0,22,95,201,207,207,207,193,193,200,207,102,23,0,10,110,193,198,207;207,207,207,200,90,10,0,1,32,81,140,187,174,131,79,32,1,0,11,96,199,193,193,207;207,207,207,205,193,110,29,0,0,0,0,0,0,0,0,0,0,29,111,197,193,193,198,207;207,207,207,207,196,193,179,97,43,14,3,3,3,3,14,43,99,184,195,193,196,200,206,207;207,207,207,207,207,207,207,207,207,207,207,207,207,207,207,207,207,207,207,207,207,207,207,207];

icon_manual = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0;0,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0;0,NaN,1,1,1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,NaN,0,0,NaN,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,0;0,NaN,NaN,NaN,1,1,1,1,1,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,0,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,1,1,NaN,NaN,NaN,NaN,0,NaN,0,NaN,NaN,NaN,0,0,0,NaN,0,NaN,NaN,0,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,1,1,NaN,NaN,NaN,0,NaN,NaN,NaN,0,0,NaN,0,NaN,0,NaN,NaN,0,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,1,NaN,NaN,0,NaN,NaN,0,0,NaN,NaN,0,NaN,0,NaN,NaN,0,NaN,NaN,0;0,NaN,0,0,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,0,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,0,NaN,NaN,0;0,NaN,0,0,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,0,NaN,NaN,0,NaN,0,0,NaN,NaN,NaN,0,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,NaN,0,NaN,0,NaN,0,0,NaN,0,NaN,NaN,0,0,0,NaN,NaN,0;0,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,0,NaN,0,NaN,0,0,NaN,NaN,NaN,0,0,NaN,NaN,NaN,NaN,0;0,NaN,0,0,NaN,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,0,NaN,0,0,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,0,NaN,NaN,NaN,NaN,NaN,0,NaN,0,NaN,0,0,0,0,NaN,NaN,NaN,1,1,NaN,NaN,0;0,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,NaN,0,NaN,0,0,NaN,NaN,NaN,NaN,NaN,1,1,NaN,NaN,NaN,NaN,0;0,NaN,0,0,NaN,NaN,0,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,0,NaN,NaN,NaN,NaN,1,1,1,1,NaN,NaN,NaN,1,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,1,1,1,NaN,NaN,NaN,NaN,1,1,NaN,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,0,NaN,1,1,NaN,NaN,NaN,1,1,1,NaN,NaN,NaN,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,0,NaN,NaN,NaN,NaN,1,1,1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0;0,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,1,1,1,1,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0;NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,0,0,0;NaN,NaN,NaN,0,0,0,0,0,0,NaN,NaN,NaN,0,0,NaN,NaN,0,NaN,NaN,0,0,NaN,NaN,0,0,0,0,0,0,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,0,0,NaN,0,NaN,0,0,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,0,0,0,0,0,0,NaN,NaN,0,0,NaN,NaN,0,0,0,0,0,0,NaN,NaN,0,0,0,0,0,0,0,0,0;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
icon_manual(:,:,2) = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0;0,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0;0,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,NaN,0,0,NaN,NaN,NaN,NaN,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,NaN,NaN,0;0,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,0.500000000000000,0.500000000000000,0,NaN,NaN,0.500000000000000,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,0,NaN,0,NaN,NaN,NaN,0.500000000000000,0.500000000000000,0.500000000000000,NaN,0,NaN,NaN,0.500000000000000,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,0,NaN,NaN,NaN,0.500000000000000,0.500000000000000,NaN,0,NaN,0,NaN,NaN,0.500000000000000,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,0,NaN,NaN,0.500000000000000,0.500000000000000,NaN,NaN,0,NaN,0,NaN,NaN,0.500000000000000,NaN,NaN,0;0,NaN,0,0,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,0.500000000000000,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,0.500000000000000,NaN,NaN,0;0,NaN,0,0,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,0.500000000000000,NaN,NaN,0,NaN,0,0,NaN,NaN,NaN,0.500000000000000,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,NaN,0,NaN,0.500000000000000,NaN,0,0,NaN,0,NaN,NaN,0.500000000000000,0.500000000000000,0.500000000000000,NaN,NaN,0;0,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,0,NaN,0.500000000000000,NaN,0,0,NaN,NaN,NaN,0.500000000000000,0.500000000000000,NaN,NaN,NaN,NaN,0;0,NaN,0,0,NaN,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,0.500000000000000,NaN,0,0,NaN,0.500000000000000,0.500000000000000,0.500000000000000,NaN,NaN,NaN,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,0,NaN,NaN,NaN,NaN,NaN,0,NaN,0.500000000000000,NaN,0.500000000000000,0.500000000000000,0.500000000000000,0.500000000000000,NaN,NaN,NaN,0,0,NaN,NaN,0;0,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,NaN,0,NaN,0.500000000000000,0.500000000000000,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,NaN,NaN,0;0,NaN,0,0,NaN,NaN,0,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,0,NaN,NaN,NaN,NaN,0,0,0,0,NaN,NaN,NaN,0,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,0,0,NaN,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,0,NaN,0,0,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,0,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0;0,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,0,0,0,0,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0;NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,0,0,0;NaN,NaN,NaN,0,0,0,0,0,0,NaN,NaN,NaN,0,0,NaN,NaN,0,NaN,NaN,0,0,NaN,NaN,0,0,0,0,0,0,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,0,0,NaN,0,NaN,0,0,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,0,0,0,0,0,0,NaN,NaN,0,0,NaN,NaN,0,0,0,0,0,0,NaN,NaN,0,0,0,0,0,0,0,0,0;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
icon_manual(:,:,3) = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0;0,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0;0,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,NaN,0,0,NaN,NaN,NaN,NaN,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,NaN,NaN,0;0,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,0.250000000000000,0.250000000000000,1,NaN,NaN,0.250000000000000,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,0,NaN,0,NaN,NaN,NaN,0.250000000000000,0.250000000000000,0.250000000000000,NaN,1,NaN,NaN,0.250000000000000,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,0,NaN,NaN,NaN,0.250000000000000,0.250000000000000,NaN,1,NaN,1,NaN,NaN,0.250000000000000,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,0,NaN,NaN,0.250000000000000,0.250000000000000,NaN,NaN,1,NaN,1,NaN,NaN,0.250000000000000,NaN,NaN,0;0,NaN,0,0,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,0.250000000000000,NaN,NaN,NaN,NaN,1,1,1,NaN,NaN,0.250000000000000,NaN,NaN,0;0,NaN,0,0,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,0.250000000000000,NaN,NaN,1,NaN,1,1,NaN,NaN,NaN,0.250000000000000,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,NaN,0,NaN,0.250000000000000,NaN,1,1,NaN,1,NaN,NaN,0.250000000000000,0.250000000000000,0.250000000000000,NaN,NaN,0;0,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,0,NaN,0.250000000000000,NaN,1,1,NaN,NaN,NaN,0.250000000000000,0.250000000000000,NaN,NaN,NaN,NaN,0;0,NaN,0,0,NaN,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,0.250000000000000,NaN,1,1,NaN,0.250000000000000,0.250000000000000,0.250000000000000,NaN,NaN,NaN,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,0,NaN,NaN,NaN,NaN,NaN,0,NaN,0.250000000000000,NaN,0.250000000000000,0.250000000000000,0.250000000000000,0.250000000000000,NaN,NaN,NaN,0,0,NaN,NaN,0;0,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,NaN,0,NaN,0.250000000000000,0.250000000000000,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,NaN,NaN,0;0,NaN,0,0,NaN,NaN,0,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,0,NaN,NaN,NaN,NaN,0,0,0,0,NaN,NaN,NaN,0,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,0,0,NaN,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,0,NaN,0,0,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,0;0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,0,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0;0,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,0,0,0,0,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0;NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,NaN,NaN,NaN,NaN;0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,0,0,0,NaN,NaN,NaN,NaN,NaN,0,0,0;NaN,NaN,NaN,0,0,0,0,0,0,NaN,NaN,NaN,0,0,NaN,NaN,0,NaN,NaN,0,0,NaN,NaN,0,0,0,0,0,0,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,NaN,NaN,0,0,NaN,0,NaN,0,0,NaN,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,0,0,0,0,0,0,NaN,NaN,0,0,NaN,NaN,0,0,0,0,0,0,NaN,NaN,0,0,0,0,0,0,0,0,0;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0,0,0,0,0,0,0,0,0,0,0,0,0,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
% arrow = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;...
%          NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN;...
%          NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN;...
%          NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN;...
%          NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN;...
%          NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN;...
%          NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN;...
%          NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN;...
%          NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN;...
%          NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN;...
%          NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;...
%          NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
% arrow = repmat(arrow,[1 1 3]);                      %Play arrow for play buttons
arrowAll = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
arrowAll = repmat(arrowAll,[1 1 3]);
arrowZoom = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    0.6,0.6,0.6,0  ,0  ,0  ,0  ,0  ,NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN,0  ,0  ,0  ,0  ,0.6,0.6,0.6,0.6;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
arrowZoom(:,:,2) = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    0.6,0.6,0.6,0  ,0  ,0  ,0  ,0  ,NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN,0  ,0  ,0  ,0  ,0.6,0.6,0.6,0.78;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
arrowZoom(:,:,3) = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    0.6,0.6,0.6,0  ,0  ,0  ,0  ,0  ,NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN,0  ,0  ,0  ,0  ,0.6,0.6,0.6,0.6;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
pausebt = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;...
    NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;...
    NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;...
    NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;...
    NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;...
    NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;...
    NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;...
    NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;...
    NaN,0  ,0  ,NaN,0  ,0  ,NaN,NaN;...
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
pausebt = repmat(pausebt,[1 1 3]);                  %Pause icon for play buttons
lockOpen = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,0  ,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,0  ,0  ,NaN,0  ,NaN;
    NaN,0  ,NaN,0  ,0  ,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
lockOpen = repmat(lockOpen,[1 1 3]);
lockClosed = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,0  ,0  ,NaN,0  ,NaN;
    NaN,0  ,NaN,0  ,0  ,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
lockClosed = repmat(lockClosed,[1 1 3]);
arrows =   [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ;
    NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
arrows = repmat(arrows,[1 1 3]);
arrow_lr(:,:) = [ ...
    NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN   2   2   2 NaN NaN   2   2   2 NaN NaN NaN NaN ; ...
    NaN NaN   2   2   2   1   2 NaN NaN   2   1   2   2   2 NaN NaN ; ...
    2   2   2   1   1   1   2   2   2   2   1   1   1   2   2   2 ; ...
    2   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2 ; ...
    2   2   2   1   1   1   2   2   2   2   1   1   1   2   2   2 ; ...
    NaN NaN   2   2   2   1   2 NaN NaN   2   1   2   2   2 NaN NaN ; ...
    NaN NaN NaN NaN   2   2   2 NaN NaN   2   2   2 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN ];
arrow_ud = arrow_lr';
marker =   [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
marker = repmat(marker, [1 1 3]);
menuBarIcon(:,:,1) = [0  ,0  ,0  ,0  ,0  ,NaN,0  ,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,NaN,NaN,0  ,0  ,0  ,NaN,NaN;
    0  ,0  ,0  ,0  ,0  ,NaN,0  ,NaN,0  ,NaN,0  ,NaN,NaN,NaN,0  ,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,NaN,0  ,0  ,0  ,0  ,0  ,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,NaN,0  ,NaN,NaN,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,NaN,NaN,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0.04,0.04,0.04,0.04,0.04,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.04;
    0.04,1  ,1  ,1  ,0.04,0.04,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.04,NaN;
    0.04,1  ,1  ,1  ,0.04,1  ,0.04,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN;
    0.04,1  ,1  ,1  ,0.04,0.04,0.04,0.04,NaN,0  ,1  ,NaN,1  ,0  ,0  ,0;
    0.04,1  ,1  ,1  ,1  ,1  ,1  ,0.04,NaN,0  ,1  ,NaN,1  ,NaN,1  ,NaN;
    0.04,1  ,1  ,1  ,1  ,1  ,1  ,0.04,NaN,0  ,NaN,1  ,NaN,0  ,0  ,0;
    0.04,1  ,1  ,1  ,1  ,1  ,1  ,0.04,NaN,0  ,1  ,NaN,0  ,0.501960784313726,0.501960784313726,0.501960784313726;
    0.04,1  ,1  ,1  ,1  ,1  ,1  ,0.04,NaN,0  ,NaN,0  ,0.501960784313726,0.501960784313726,0.501960784313726,0.501960784313726;
    0.04,1  ,1  ,1  ,1  ,1  ,1  ,0.04,NaN,0  ,0  ,0.501960784313726,0.501960784313726,0.501960784313726,0.501960784313726,0.501960784313726;
    0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0;];

menuBarIcon(:,:,2) = [0  ,0  ,0  ,0  ,0  ,NaN,0  ,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,NaN,NaN,0  ,0  ,0  ,NaN,NaN;
    0  ,0  ,0  ,0  ,0  ,NaN,0  ,NaN,0  ,NaN,0  ,NaN,NaN,NaN,0  ,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,NaN,0  ,0  ,0  ,0  ,0  ,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,NaN,0  ,NaN,NaN,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,NaN,NaN,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0.14,0.14,0.14,0.14,0.14,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.14;
    0.14,1  ,1  ,1  ,0.14,0.14,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.14,NaN;
    0.14,1  ,1  ,1  ,0.14,1  ,0.14,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN;
    0.14,1  ,1  ,1  ,0.14,0.14,0.14,0.14,NaN,0  ,1  ,NaN,1  ,0  ,0  ,0;
    0.14,1  ,1  ,1  ,1  ,1  ,1  ,0.14,NaN,0  ,1  ,NaN,1  ,NaN,1  ,NaN;
    0.14,1  ,1  ,1  ,1  ,1  ,1  ,0.14,NaN,0  ,NaN,1  ,NaN,0  ,0  ,0;
    0.14,1  ,1  ,1  ,1  ,1  ,1  ,0.14,NaN,0  ,1  ,NaN,0  ,0.501960784313726,0.501960784313726,0.501960784313726;
    0.14,1  ,1  ,1  ,1  ,1  ,1  ,0.14,NaN,0  ,NaN,0  ,0.501960784313726,0.501960784313726,0.501960784313726,0.501960784313726;
    0.14,1  ,1  ,1  ,1  ,1  ,1  ,0.14,NaN,0  ,0  ,0.501960784313726,0.501960784313726,0.501960784313726,0.501960784313726,0.501960784313726;
    0.14,0.14,0.14,0.14,0.14,0.14,0.14,0.14,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0;];
menuBarIcon(:,:,3) = [0  ,0  ,0  ,0  ,0  ,NaN,0  ,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,NaN,NaN,0  ,0  ,0  ,NaN,NaN;
    0  ,0  ,0  ,0  ,0  ,NaN,0  ,NaN,0  ,NaN,0  ,NaN,NaN,NaN,0  ,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,NaN,0  ,0  ,0  ,0  ,0  ,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,NaN,0  ,NaN,NaN,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,NaN,NaN,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0.42,0.42,0.42,0.42,0.42,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.42;
    0.42,1  ,1  ,1  ,0.42,0.42,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.42,NaN;
    0.42,1  ,1  ,1  ,0.42,1  ,0.42,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN;
    0.42,1  ,1  ,1  ,0.42,0.42,0.42,0.42,NaN,0  ,0  ,NaN,0  ,0  ,0  ,0;
    0.42,1  ,1  ,1  ,1  ,1  ,1  ,0.42,NaN,0  ,0  ,NaN,0  ,NaN,0  ,NaN;
    0.42,1  ,1  ,1  ,1  ,1  ,1  ,0.42,NaN,0  ,NaN,0  ,NaN,0  ,0  ,0;
    0.42,1  ,1  ,1  ,1  ,1  ,1  ,0.42,NaN,0  ,0  ,NaN,0  ,0  ,0  ,0;
    0.42,1  ,1  ,1  ,1  ,1  ,1  ,0.42,NaN,0  ,NaN,0  ,0  ,0  ,0  ,0;
    0.42,1  ,1  ,1  ,1  ,1  ,1  ,0.42,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0.42,0.42,0.42,0.42,0.42,0.42,0.42,0.42,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0;];
RGB1 =    [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    1  ,1  ,1  ,NaN,NaN,NaN,NaN,1  ,1  ,NaN,NaN,1  ,1  ,1  ,NaN,NaN;
    1  ,NaN,NaN,1  ,NaN,NaN,1  ,NaN,NaN,1  ,NaN,1  ,NaN,NaN,1  ,NaN;
    1  ,NaN,NaN,1  ,NaN,1  ,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,1;
    1  ,NaN,NaN,1  ,NaN,1  ,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,NaN;
    1  ,1  ,1  ,NaN,NaN,1  ,NaN,1  ,1  ,1  ,NaN,1  ,NaN,NaN,1  ,NaN;
    1  ,NaN,1  ,NaN,NaN,1  ,NaN,NaN,NaN,1  ,NaN,1  ,NaN,NaN,NaN,1;
    1  ,NaN,NaN,1  ,NaN,NaN,1  ,NaN,NaN,1  ,NaN,1  ,NaN,NaN,NaN,1;
    1  ,NaN,NaN,NaN,1  ,NaN,NaN,1  ,1  ,1  ,NaN,1  ,1  ,1  ,1  ,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
RGB = NaN(10,16,3);
RGB(:,1:5,1) = RGB1(:,1:5);
RGB(:,6:11  ,1) = RGB1(:,6:11)*0;
RGB(:,12:16,1) = RGB1(:,12:16)*0;
RGB(:,6:11  ,2) = 0.8 * RGB1(:,6:11);
RGB(:,12:16,3) = RGB1(:,12:16);
RGB2(:,:,1) = [1  ,1  ,1  ,1  ,1  ,0,0,0,0,0,0,0,0,0,0,0;
    1  ,1  ,1  ,1  ,1  ,0,0,0,0,0,0,0,0,0,0,0;
    1  ,1  ,1  ,1  ,1  ,0,0,0,0,0,0,0,0,0,0,0;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    1  ,1  ,1  ,1  ,1  ,0,0,0,0,0,0,0,0,0,0,0;
    1  ,1  ,1  ,1  ,1  ,0,0,0,0,0,0,0,0,0,0,0;
    1  ,1  ,1  ,1  ,1  ,0,0,0,0,0,0,0,0,0,0,0;];
RGB2(:,:,2) =  [0,0,0,0,0,1  ,1  ,1  ,1  ,1  ,1  ,0,0,0,0,0;
    0,0,0,0,0,1  ,1  ,1  ,1  ,1  ,1  ,0,0,0,0,0;
    0,0,0,0,0,1  ,1  ,1  ,1  ,1  ,1  ,0,0,0,0,0;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0,0,0,0,0,1  ,1  ,1  ,1  ,1  ,1  ,0,0,0,0,0;
    0,0,0,0,0,1  ,1  ,1  ,1  ,1  ,1  ,0,0,0,0,0;
    0,0,0,0,0,1  ,1  ,1  ,1  ,1  ,1  ,0,0,0,0,0;];
RGB2(:,:,3) =   [0,0,0,0,0,0,0,0,0,0,0,1  ,1  ,1  ,1  ,1;
    0,0,0,0,0,0,0,0,0,0,0,1  ,1  ,1  ,1  ,1;
    0,0,0,0,0,0,0,0,0,0,0,1  ,1  ,1  ,1  ,1;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0,0,0,0,0,0,0,0,0,0,0,1  ,1  ,1  ,1  ,1;
    0,0,0,0,0,0,0,0,0,0,0,1  ,1  ,1  ,1  ,1;
    0,0,0,0,0,0,0,0,0,0,0,1  ,1  ,1  ,1  ,1;];
objBt(:,:,1) =     [NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN;
    NaN,NaN,1  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN;
    NaN,NaN,1  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN;
    NaN,NaN,1  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN;
    NaN,NaN,1  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN;
    NaN,NaN,1  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN;
    NaN,NaN,1  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN;
    NaN,NaN,1  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN;
    NaN,NaN,1  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN;
    0  ,0  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,0  ,0  ,0  ,0;
    NaN,NaN,1  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN;
    NaN,NaN,1  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN;
    NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
objBt(:,:,2) =     [NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
objBt(:,:,3) = objBt(:,:,2);

tifParBt =         [NaN,NaN,0  ,0  ,0  ,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,0  ,0;
    NaN,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0  ,NaN,0  ,0  ,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,0  ,NaN,NaN,0  ,0  ,NaN,NaN,NaN,0  ,NaN,0  ,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,0  ,NaN,0  ,NaN,NaN,0  ,NaN,NaN,0  ,0  ,NaN,NaN,NaN;
    NaN,0  ,0  ,0  ,0  ,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,0  ,NaN,0  ,NaN,NaN,NaN,NaN;];
tifParBt = repmat(tifParBt, [1 1 3]);
meanBt{1}(:,:,1) =  [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,0  ,1  ,1  ,1  ,0  ,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    1  ,0.5,0.5,0.5,1  ,0.5,0.5,0.5,0  ,1  ,1  ,1  ,0  ,0.5,0.5,0.5,1  ,0.5,0.5,0.5,1;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,0  ,1  ,1  ,1  ,0  ,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
meanBt{1}(:,:,2) =  [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,0  ,0  ,0  ,0  ,0  ,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    0  ,0.5,0.5,0.5,0  ,0.5,0.5,0.5,0  ,0  ,0  ,0  ,0  ,0.5,0.5,0.5,0  ,0.5,0.5,0.5,0;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,0  ,0  ,0  ,0  ,0  ,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
meanBt{2}(:,:,1) =  [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,0.5,0.5,0.5,NaN;
    1  ,0.5,0.5,0.5,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,0.5,0.5,0.5,0;
    NaN,0.5,0.5,0.5,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
meanBt{2}(:,:,2) =  [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.5,0.5,0.5,NaN;
    0  ,0.5,0.5,0.5,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.5,0.5,0.5,0;
    NaN,0.5,0.5,0.5,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN,0.5,0.5,0.5,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
meanBt{3}(:,:,1) =  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;];
meanBt{3}(:,:,2) =  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;];
meanBt{1}(:,:,3) = meanBt{1}(:,:,2);
meanBt{2}(:,:,3) = meanBt{2}(:,:,2);
meanBt{3}(:,:,3) = meanBt{3}(:,:,2);
meanBt{4}(:,:,1) =  [0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0  ,0.3,0.3,0.3,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.3,0.3,0.3,0;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0  ,0.3,0.3,0.3,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,0.3,0.3,0.3,0;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0  ,0.3,0.3,0.3,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.3,0.3,0.3,0;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;];
meanBt{4}(:,:,2) =  [0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0  ,0.3,0.3,0.3,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,0.3,0.3,0.3,0;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0  ,0.3,0.3,0.3,0  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,0  ,0.3,0.3,0.3,0;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0  ,0.3,0.3,0.3,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,0  ,0.3,0.3,0.3,0;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0;
    0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0  ,0.3,0.3,0.3,0;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;];
meanBt{4}(:,:,3) =  [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,1  ,1  ,1  ,0,0,0,0,0,0,0,0,0,0,0,0,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,0,0,0,0,0,0,0,0,0,0,0,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,0,0,0,0,0,0,0,0,0,0,0,0,1  ,1  ,1  ,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,1  ,1  ,1  ,0,0,0,0,0,0,0,0,0,0,0,0,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,0,0,0,0,0,0,0,0,0,0,0,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,0,0,0,0,0,0,0,0,0,0,0,0,1  ,1  ,1  ,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,1  ,1  ,1  ,0,0,0,0,0,0,0,0,0,0,0,0,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,0,0,0,0,0,0,0,0,0,0,0,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,0,0,0,0,0,0,0,0,0,0,0,0,1  ,1  ,1  ,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0,1  ,1  ,1  ,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;];
meanBt{5}(:,:,1) =   [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN;
    NaN,NaN,NaN,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,NaN,NaN,NaN;
    NaN,NaN,NaN,1  ,1  ,0  ,1  ,0  ,1  ,0  ,0  ,0  ,1  ,0  ,1  ,0  ,1  ,1  ,NaN,NaN,NaN;
    NaN,NaN,NaN,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,NaN,NaN,NaN;
    NaN,NaN,NaN,1  ,1  ,0  ,1  ,0  ,1  ,0  ,0  ,0  ,1  ,0  ,1  ,0  ,1  ,1  ,NaN,NaN,NaN;
    NaN,NaN,NaN,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,NaN,NaN,NaN;
    NaN,NaN,NaN,1  ,1  ,0  ,1  ,0  ,1  ,0  ,0  ,0  ,1  ,0  ,1  ,0  ,1  ,1  ,NaN,NaN,NaN;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    NaN,NaN,NaN,1  ,1  ,0  ,1  ,0  ,1  ,0  ,0  ,0  ,1  ,0  ,1  ,0  ,1  ,1  ,NaN,NaN,NaN;
    NaN,NaN,NaN,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,NaN,NaN,NaN;
    NaN,NaN,NaN,1  ,1  ,0  ,1  ,0  ,1  ,0  ,0  ,0  ,1  ,0  ,1  ,0  ,1  ,1  ,NaN,NaN,NaN;
    NaN,NaN,NaN,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,NaN,NaN,NaN;
    NaN,NaN,NaN,1  ,1  ,0  ,1  ,0  ,1  ,0  ,0  ,0  ,1  ,0  ,1  ,0  ,1  ,1  ,NaN,NaN,NaN;
    NaN,NaN,NaN,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,0  ,1  ,NaN,NaN,NaN;
    NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
meanBt{5}(:,:,2) =   [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
meanBt{5}(:,:,3) = meanBt{5}(:,:,2);
meanBt{6}        = RGB;
histIcon(:,:,1) =  [NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,1  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,1  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,1  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,1  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,1  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,1  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN;
    NaN,0  ,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,1  ,NaN;
    NaN,0  ,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,0  ,0  ,0  ,NaN,NaN,1  ,NaN;
    NaN,0  ,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,0  ,0  ,0  ,NaN,NaN,1  ,NaN;
    NaN,0  ,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,1  ,NaN;
    NaN,0  ,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,0;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,0;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,0;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
histIcon(:,:,2) =  [NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,1  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,1  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,0  ,0  ,0  ,0  ,1  ,0  ,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,0  ,0  ,0  ,1  ,0  ,0  ,0  ,NaN,NaN,0  ,0  ,0  ,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,0  ,0  ,1  ,0  ,0  ,0  ,0  ,NaN,NaN,0  ,0  ,0  ,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,0  ,0  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,0  ,NaN;
    NaN,0  ,NaN,0  ,0  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
    NaN,0  ,0  ,0  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    NaN,0  ,0  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0  ,0  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
histIcon(:,:,3) =  [NaN,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN;
    NaN,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,NaN;
    NaN,1  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0  ,NaN;
    NaN,1  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,NaN;
    NaN,1  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,1  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,1  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,1  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,1  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,1  ,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,NaN;
    NaN,1  ,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,0  ,0  ,0  ,NaN,NaN,0  ,NaN;
    NaN,1  ,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,0  ,0  ,0  ,NaN,NaN,0  ,NaN;
    NaN,1  ,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,0  ,NaN;
    NaN,1  ,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
    NaN,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    NaN,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
histIcon(:,1  ,:) = [];
histIcon(:,2,:) = [];
exportBt(:,:,1) =  [NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,1  ,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,1  ,NaN,1  ,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,0.3,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,1  ,NaN,NaN,NaN,1  ,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN,0.3,0.3,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,1  ,NaN,NaN,NaN,1  ,NaN,NaN,1  ,NaN,NaN,NaN,0.3,0.3,0.3,0.5,0.3,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,1  ,1  ,NaN,NaN,NaN,NaN,0.3,0.5,0.5,0.5,0.5,0.3,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.3,0.5,0.5,0.5,0.5,0.5,0.3,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.3,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.3;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.3,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.3,0.5,0.5,0.5,0.5,0.5,0.3,NaN,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,0.3,0.5,0.5,0.5,0.5,0.3,NaN,NaN,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,0.3,0.3,0.3,0.5,0.3,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,0  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,0.3,0.3,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,1  ,1  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,0.3,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,1  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,0  ,0  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
exportBt(:,:,2) =  [NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,0  ,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0.3,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,0.3,0.3,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,NaN,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,0.3,0.3,0.3,0.5,0.3,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,0.3,0.5,0.5,0.5,0.5,0.3,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.3,0.5,0.5,0.5,0.5,0.5,0.3,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.3,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.3;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.3,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.3,0.5,0.5,0.5,0.5,0.5,0.3,NaN,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,0.3,0.5,0.5,0.5,0.5,0.3,NaN,NaN,NaN;
    NaN,0  ,0  ,1  ,1  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,NaN,0.3,0.3,0.3,0.5,0.3,NaN,NaN,NaN,NaN;
    NaN,0  ,1  ,1  ,0  ,0  ,0  ,1  ,1  ,1  ,0  ,0  ,NaN,NaN,NaN,0.3,0.3,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,1  ,0  ,1  ,1  ,0  ,0  ,1  ,1  ,1  ,0  ,NaN,NaN,NaN,0.3,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,1  ,0  ,0  ,0  ,0  ,1  ,1  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,1  ,1  ,0  ,1  ,1  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,0  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
exportBt(:,:,3) =  [NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,0  ,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0.3,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,0.3,0.3,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,NaN,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,0.3,0.3,0.3,0.5,0.3,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,0.3,0.5,0.5,0.5,0.5,0.3,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.3,0.5,0.5,0.5,0.5,0.5,0.3,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.3,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.3;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0.3,0.5,0.5,0.5,0.5,0.5,0.5,0.3,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.3,0.5,0.5,0.5,0.5,0.5,0.3,NaN,NaN;
    NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,0.3,0.5,0.5,0.5,0.5,0.3,NaN,NaN,NaN;
    NaN,1  ,1  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,1  ,NaN,0.3,0.3,0.3,0.5,0.3,NaN,NaN,NaN,NaN;
    NaN,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,NaN,NaN,NaN,0.3,0.3,NaN,NaN,NaN,NaN,NaN;
    NaN,1  ,0  ,0  ,0  ,1  ,0  ,0  ,0  ,0  ,0  ,1  ,NaN,NaN,NaN,0.3,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,1  ,1  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,1  ,1  ,1  ,0  ,0  ,0  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
exportBt(:,1  ,:) = [];
exportBt(:,7,:) = [];
exportBt(:,end,:) = [];
roiRectBt =          [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
roiRectBt = repmat(roiRectBt,[1 1 3]);
roiEllipseBt =       [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN;
    NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
roiEllipseBt = repmat(roiEllipseBt,[1 1 3]);
roiPolygonBt =     [NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0  ,0  ,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN;
    NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0  ,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
roiPolygonBt = repmat(roiPolygonBt,[1 1 3]);
roiExportBt(:,:,1) =[0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN,NaN;
    0.2,0.2,0.2,1  ,1  ,0.2,0.2,0.2,0.2,0.2,0.2,0.2,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,NaN,NaN,NaN,NaN;
    0.2,0.2,1  ,0.2,0.2,1  ,0.2,0.2,0.2,0.2,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN;
    0.2,1  ,0.2,0.2,0.2,0.2,1  ,0.2,0.2,0.2,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN;
    0.2,1  ,0.2,0.2,0.2,0.2,0.2,1  ,0.2,0.2,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN;
    0.2,1  ,0.2,0.2,0.2,0.2,0.2,0.2,1  ,0.2,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1;
    0.2,1  ,0.2,0.2,0.2,0.2,0.2,0.2,1  ,0.2,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN;
    0.2,0.2,1  ,0.2,0.2,0.2,0.2,0.2,1  ,0.2,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN;
    0.2,0.2,0.2,1  ,0.2,0.2,0.2,0.2,1  ,0.2,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN;
    0.2,0.2,0.2,0.2,1  ,1  ,1  ,1  ,0.2,0.2,0.2,0.2,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,NaN,NaN,NaN,NaN;
    0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN,NaN;];
roiExportBt(:,:,2) = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN;
    0.2,0.2,0.2,1  ,1  ,0.2,0.2,0.2,0.2,0.2,0.2,0.2,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN;
    0.2,0.2,1  ,0.2,0.2,1  ,0.2,0.2,0.2,0.2,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN;
    0.2,1  ,0.2,0.2,0.2,0.2,1  ,0.2,0.2,0.2,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    0.2,1  ,0.2,0.2,0.2,0.2,0.2,1  ,0.2,0.2,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
    0.2,1  ,0.2,0.2,0.2,0.2,0.2,0.2,1  ,0.2,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0.2,1  ,0.2,0.2,0.2,0.2,0.2,0.2,1  ,0.2,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
    0.2,0.2,1  ,0.2,0.2,0.2,0.2,0.2,1  ,0.2,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    0.2,0.2,0.2,1  ,0.2,0.2,0.2,0.2,1  ,0.2,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN;
    0.2,0.2,0.2,0.2,1  ,1  ,1  ,1  ,0.2,0.2,0.2,0.2,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN;
    0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN;];
roiExportBt(:,:,3) = [1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN;
    1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN;
    1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN;
    1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
    1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
    1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN;
    1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN;
    1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN;];
roiImportBt(:,:,1) = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.2,0.2,0.2,1  ,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2;
    NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.2,0.2,0.2,1  ,1  ,0.2,0.2,0.2,0.2,0.2,0.2,0.2;
    NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0.2,0.2,0.2,0.2,0.2,0.2;
    NaN,0  ,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0.2,0.2,0.2,0.2,0.2;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0.2,0.2,0.2,0.2;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0.2,0.2,0.2;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0.2,0.2,0.2,0.2;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0.2,0.2,0.2,0.2,0.2;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,0  ,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,0.2,0.2,0.2,0.2,0.2,0.2;
    NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0.2,0.2,0.2,1  ,1  ,0.2,0.2,0.2,0.2,0.2,0.2,0.2;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.2,0.2,0.2,1  ,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2;];
roiImportBt(:,:,2) = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.2,0.2,0.2,0  ,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2;
    NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.2,0.2,0.2,0  ,0  ,0.2,0.2,0.2,0.2,0.2,0.2,0.2;
    NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.2,0.2,0.2,0.2,0.2,0.2;
    NaN,0  ,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.2,0.2,0.2,0.2,0.2;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.2,0.2,0.2,0.2;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.2,0.2,0.2;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.2,0.2,0.2,0.2;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.2,0.2,0.2,0.2,0.2;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,0  ,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0.2,0.2,0.2,0.2,0.2,0.2;
    NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0.2,0.2,0.2,0  ,0  ,0.2,0.2,0.2,0.2,0.2,0.2,0.2;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.2,0.2,0.2,0  ,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2;];
roiImportBt(:,:,3) = [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,0  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1;
    NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,0  ,0  ,1  ,1  ,1  ,1  ,1  ,1  ,1;
    NaN,NaN,0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,1  ,1  ,1;
    NaN,0  ,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,1  ,1;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,1;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,1;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,1  ,1;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,0  ,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,1  ,1  ,1  ,1  ,1  ,1;
    NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,0  ,0  ,1  ,1  ,1  ,1  ,1  ,1  ,1;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,0  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1;];
roiExportDataBt(:,:,1) =   [0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN,NaN;
    0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,NaN,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN;
    0  ,1  ,1  ,1  ,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN;
    0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1;
    0  ,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN;
    0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN;
    0  ,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,1  ,1  ,NaN,NaN,NaN,NaN;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,1  ,NaN,NaN,NaN,NaN,NaN;];
roiExportDataBt(:,:,2) =   [0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN;
    0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN;
    0  ,0  ,0  ,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
    0  ,0.5,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0  ,NaN,0.5,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
    0  ,NaN,NaN,0.5,NaN,NaN,NaN,NaN,NaN,0.5,0.5,0.5,0.5,0.5,0.5,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    0  ,NaN,NaN,NaN,0.5,0.5,0.5,0.5,0.5,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN;];
roiExportDataBt(:,:,3) =   [0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN;
    0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN;
    0  ,0  ,0  ,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
    0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0;
    0  ,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
    0  ,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN;
    0  ,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN;
    0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN;
    0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN;];
roiDeleteBt(:,:,1) =       [NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,NaN,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,NaN,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1  ,1  ,1  ,NaN,NaN,NaN;];
roiDeleteBt(:,:,2) = NaN;
roiDeleteBt(:,:,3) = NaN;
roiShiftBt =               [NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN;
    NaN,NaN,0  ,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,0  ,NaN,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN;
    NaN,NaN,0  ,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,0  ,NaN,NaN;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN,0  ,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
roiShiftBt = repmat(roiShiftBt, [1 1 3]);
roiRotateBt =              [NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
roiRotateBt = repmat(roiRotateBt, [1 1 3]);
roiScaleBt =               [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,0  ,NaN,0  ,NaN,0  ,NaN,0  ,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    NaN,0  ,NaN,0  ,NaN,0  ,NaN,0  ,NaN,0  ,NaN,0  ,NaN,0  ,NaN,0  ,NaN,0  ,NaN,0  ,NaN,0  ,NaN;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
roiScaleBt = repmat(roiScaleBt, [1 1 3]);
invertBt =                 [NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;
    1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0;
    1  ,1  ,1  ,1  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0  ,0  ,0  ,0;
    0.8,0.8,0.8,0.8,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.2,0.2,0.2,0.2;
    0.8,0.8,0.8,0.8,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.2,0.2,0.2,0.2;
    0.7,0.7,0.7,0.7,NaN,NaN,NaN,NaN,0  ,0  ,NaN,NaN,0.4,0.4,0.4,0.4;
    0.7,0.7,0.7,0.7,NaN,NaN,0  ,0  ,0  ,0  ,0  ,NaN,0.4,0.4,0.4,0.4;
    0.5,0.5,0.5,0.5,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.5,0.5,0.5,0.5;
    0.5,0.5,0.5,0.5,NaN,0  ,0  ,0  ,0  ,0  ,NaN,NaN,0.5,0.5,0.5,0.5;
    0.4,0.4,0.4,0.4,NaN,NaN,0  ,0  ,NaN,NaN,NaN,NaN,0.7,0.7,0.7,0.7;
    0.4,0.4,0.4,0.4,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.7,0.7,0.7,0.7;
    0.2,0.2,0.2,0.2,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.8,0.8,0.8,0.8;
    0.2,0.2,0.2,0.2,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,0.8,0.8,0.8,0.8;
    0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1;
    0  ,0  ,0  ,0  ,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1  ,1  ,1  ,1;
    NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN;];
invertBt = repmat(invertBt, [1 1 3]);
invertBt(2:end-1,1:4,:) = permute(repmat(hot(14),[1 1 4]),[1 3 2]);
invertBt(2:end-1,end-3:end,:) = permute(repmat(1-hot(14),[1 1 4]),[1 3 2]);
flipBt = invertBt;
flipBt(2:end-1,end-3:end,:) = flipdim(permute(repmat(hot(14),[1 1 4]),[1 3 2]),1);
hand = load([matlabroot,'/toolbox/matlab/icons/pan.mat']);   %#ok Icon for panning hand
hand = struct2cell(hand);
hand = hand{1}(:,:,1);
hand(hand == 2) = NaN;
hand(hand == 1) = 2;
hand(hand == 0) = 1;
hand2 = ...
    [NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
    NaN   NaN   NaN   NaN     0     0   NaN     0     0   NaN     0     0   NaN   NaN   NaN   NaN;
    NaN   NaN   NaN     0     1     0     0     1     1     0     1     1     0   NaN   NaN   NaN;
    NaN   NaN   NaN     0     1     1     0     1     1     0     1     1     0     0   NaN   NaN;
    NaN   NaN   NaN     0     1     1     0     1     1     0     1     1     0     1     0   NaN;
    NaN   NaN     0   NaN     0     1     1     1     1     1     1     1     1     1     0   NaN;
    NaN     0     1     0     0     1     1     1     1     1     1     1     1     1     0   NaN;
    NaN     0     1     1     0     1     1     1     1     1     1     1     1     1     0   NaN;
    NaN   NaN     1     1     1     1     1     1     1     1     1     1     1     1     0   NaN;
    NaN   NaN     0     1     1     1     1     1     1     1     1     1     1     0   NaN   NaN;
    NaN   NaN     0     1     1     1     1     1     1     1     1     1     1     0   NaN   NaN;
    NaN   NaN   NaN     0     1     1     1     1     1     1     1     1     1     0   NaN   NaN;
    NaN   NaN   NaN   NaN     0     1     1     1     1     1     1     1     0   NaN   NaN   NaN;
    NaN   NaN   NaN   NaN   NaN     0     1     1     1     1     1     1     0   NaN   NaN   NaN;
    NaN   NaN   NaN   NaN   NaN     0     1     1     1     1     1     1     0   NaN   NaN   NaN;];

hand2(hand2 == 1) = 2;
hand2(hand2 == 0) = 1;

icon_linkWins(:,:,1) = [ ...
    0  11  11  11  11  11  11 255 255 NaN NaN NaN NaN NaN NaN NaN ; ...
    80  11  11  11  11  11  11 255 255 NaN NaN NaN NaN NaN NaN NaN ; ...
    80 240 240 240 240 240 240 240  80 NaN NaN NaN NaN NaN NaN NaN ; ...
    80 240 240 240 240 240 255 255 255 255 NaN NaN NaN NaN NaN NaN ; ...
    80 240 240 240 240 255 255 255 255 255 255 NaN NaN NaN NaN NaN ; ...
    80  80  80  80  80 255 255 255 255 255 255 NaN NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN 255 255 NaN NaN NaN NaN 255 255 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN 255 255 NaN NaN NaN NaN 255 255 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN 255 255 255 255 255 255 255 255 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN 255 255 255 255 255 255 255 255 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN 255 255 255 255 255 255 255 255  11  11 255 255 ; ...
    NaN NaN NaN NaN 255 255 255 255 255 255 255 255  11  11 255 255 ; ...
    NaN NaN NaN NaN 255 255 255 255 255 255 255 255 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80 240 240 240 240 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80 240 240 240 240 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80  80  80  80  80  80  80  80  80 ; ...
    ]/255;
icon_linkWins(:,:,2) = [ ...
    0 132 132 132 132 132 132   0   0 NaN NaN NaN NaN NaN NaN NaN ; ...
    80 132 132 132 132 132 132   0   0 NaN NaN NaN NaN NaN NaN NaN ; ...
    80 240 240 240 240 240 240 240  80 NaN NaN NaN NaN NaN NaN NaN ; ...
    80 240 240 240 240 240   0   0   0   0 NaN NaN NaN NaN NaN NaN ; ...
    80 240 240 240 240   0   0   0   0   0   0 NaN NaN NaN NaN NaN ; ...
    80  80  80  80  80   0   0   0   0   0   0 NaN NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0 NaN NaN NaN NaN   0   0 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0 NaN NaN NaN NaN   0   0 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 132 132   0   0 ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 132 132   0   0 ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80 240 240 240 240 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80 240 240 240 240 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80  80  80  80  80  80  80  80  80 ; ...
    ]/255;
icon_linkWins(:,:,3) = [ ...
    0 199 199 199 199 199 199   0   0 NaN NaN NaN NaN NaN NaN NaN ; ...
    80 199 199 199 199 199 199   0   0 NaN NaN NaN NaN NaN NaN NaN ; ...
    80 240 240 240 240 240 240 240  80 NaN NaN NaN NaN NaN NaN NaN ; ...
    80 240 240 240 240 240   0   0   0   0 NaN NaN NaN NaN NaN NaN ; ...
    80 240 240 240 240   0   0   0   0   0   0 NaN NaN NaN NaN NaN ; ...
    80  80  80  80  80   0   0   0   0   0   0 NaN NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0 NaN NaN NaN NaN   0   0 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0 NaN NaN NaN NaN   0   0 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 199 199   0   0 ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 199 199   0   0 ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80 240 240 240 240 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80 240 240 240 240 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80  80  80  80  80  80  80  80  80 ; ...
    ]/255;
icon_unlinkWins(:,:,1) = [ ...
    0  11  11  11  11  11  11  11 255 255 255 NaN NaN NaN NaN NaN ; ...
    80  11  11  11  11  11  11 255 255 255 255 255 255 NaN NaN NaN ; ...
    80 240 240 240 240 240 240 255 255 NaN NaN 255 255 NaN NaN NaN ; ...
    80 240 240 240 240 240 255 255  80 NaN NaN NaN 255 255 NaN NaN ; ...
    80 240 240 240 240 240 240 240  80 NaN NaN NaN 255 255 NaN NaN ; ...
    80  80  80  80  80  80  80  80  80 NaN NaN NaN 255 255 NaN NaN ; ...
    NaN NaN NaN NaN 255 NaN NaN NaN NaN NaN NaN 255 255 255 NaN NaN ; ...
    NaN NaN NaN NaN 255 255 NaN NaN NaN NaN 255 255 255 NaN NaN NaN ; ...
    NaN NaN NaN NaN 255 255 255 255 255 255 255 255 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN 255 255 255 255 255 255 255 255 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN 255 255 255 255 255 255 255 255  11  11 255 255 ; ...
    NaN NaN NaN NaN 255 255 255 255 255 255 255 255  11  11 255 255 ; ...
    NaN NaN NaN NaN 255 255 255 255 255 255 255 255 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80 240 240 240 240 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80 240 240 240 240 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80  80  80  80  80  80  80  80  80 ; ...
    ]/255;
icon_unlinkWins(:,:,2) = [ ...
    0 132 132 132 132 132 132 132   0   0   0 NaN NaN NaN NaN NaN ; ...
    80 132 132 132 132 132 132   0   0   0   0   0   0 NaN NaN NaN ; ...
    80 240 240 240 240 240 240   0   0 NaN NaN   0   0 NaN NaN NaN ; ...
    80 240 240 240 240 240 255 255  80 NaN NaN NaN   0   0 NaN NaN ; ...
    80 240 240 240 240 240 240 240  80 NaN NaN NaN   0   0 NaN NaN ; ...
    80  80  80  80  80  80  80  80  80 NaN NaN NaN   0   0 NaN NaN ; ...
    NaN NaN NaN NaN   0 NaN NaN NaN NaN NaN NaN   0   0   0 NaN NaN ; ...
    NaN NaN NaN NaN   0   0 NaN NaN NaN NaN   0   0   0 NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 132 132   0   0 ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 132 132   0   0 ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80 240 240 240 240 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80 240 240 240 240 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80  80  80  80  80  80  80  80  80 ; ...
    ]/255;
icon_unlinkWins(:,:,3) = [ ...
    0 199 199 199 199 199 199 199   0   0   0 NaN NaN NaN NaN NaN ; ...
    80 199 199 199 199 199 199   0   0   0   0   0   0 NaN NaN NaN ; ...
    80 240 240 240 240 240 240   0   0 NaN NaN   0   0 NaN NaN NaN ; ...
    80 240 240 240 240 240 255 255  80 NaN NaN NaN   0   0 NaN NaN ; ...
    80 240 240 240 240 240 240 240  80 NaN NaN NaN   0   0 NaN NaN ; ...
    80  80  80  80  80  80  80  80  80 NaN NaN NaN   0   0 NaN NaN ; ...
    NaN NaN NaN NaN   0 NaN NaN NaN NaN NaN NaN   0   0   0 NaN NaN ; ...
    NaN NaN NaN NaN   0   0 NaN NaN NaN NaN   0   0   0 NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 NaN NaN NaN NaN ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 199 199   0   0 ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 199 199   0   0 ; ...
    NaN NaN NaN NaN   0   0   0   0   0   0   0   0 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80 240 240 240 240 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80 240 240 240 240 240 240 240  80 ; ...
    NaN NaN NaN NaN NaN NaN NaN  80  80  80  80  80  80  80  80  80 ; ...
    ]/255;
icon_Record(:,:,1)=[
    NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN ,0.8500 ,0.8500 ,0.8500 ,0.8500    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN    ,NaN ,0.8500 ,0.8500 ,0.8500 ,0.8500 ,0.8500 ,0.8500 ,0.8500 ,0.8500    ,NaN    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN ,0.8500 ,0.8500 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,0.8500 ,0.8500 ,0.8500    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN ,0.8500 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,0.8500 ,0.8500    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN ,0.8500 ,0.8500 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,0.8500 ,0.8500    ,NaN    ,NaN;
    NaN    ,NaN ,0.8500 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,0.8500 ,0.8500    ,NaN    ,NaN;
    NaN    ,NaN ,0.8500 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,0.8500 ,0.8500    ,NaN    ,NaN;
    NaN    ,NaN ,0.8500 ,0.8500 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,0.8500 ,0.8500    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN ,0.8500 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,0.8500 ,0.8500    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN ,0.8500 ,0.8500 ,0.8500 ,1.0000 ,1.0000 ,1.0000 ,1.0000 ,0.8500 ,0.8500 ,0.8500    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN    ,NaN ,0.8500 ,0.8500 ,0.8500 ,0.8500 ,0.8500 ,0.8500 ,0.8500 ,0.8500    ,NaN    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN ,0.8500 ,0.8500 ,0.8500 ,0.8500    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN
    ];
icon_Record(:,:,2)=[
    NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN ,0.1600 ,0.1600 ,0.1600 ,0.1600    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN    ,NaN ,0.1600 ,0.1600 ,0.1600 ,0.1600 ,0.1600 ,0.1600 ,0.1600 ,0.1600    ,NaN    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN ,0.1600 ,0.1600      ,0      ,0      ,0      ,0      ,0 ,0.1600 ,0.1600 ,0.1600    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN ,0.1600      ,0      ,0      ,0      ,0      ,0      ,0      ,0 ,0.1600 ,0.1600    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN ,0.1600 ,0.1600      ,0      ,0      ,0      ,0      ,0      ,0      ,0      ,0 ,0.1600 ,0.1600    ,NaN    ,NaN;
    NaN    ,NaN ,0.1600      ,0      ,0      ,0      ,0      ,0      ,0      ,0      ,0      ,0 ,0.1600 ,0.1600    ,NaN    ,NaN;
    NaN    ,NaN ,0.1600      ,0      ,0      ,0      ,0      ,0      ,0      ,0      ,0      ,0 ,0.1600 ,0.1600    ,NaN    ,NaN;
    NaN    ,NaN ,0.1600 ,0.1600      ,0      ,0      ,0      ,0      ,0      ,0      ,0      ,0 ,0.1600 ,0.1600    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN ,0.1600      ,0      ,0      ,0      ,0      ,0      ,0      ,0 ,0.1600 ,0.1600    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN ,0.1600 ,0.1600 ,0.1600      ,0      ,0      ,0      ,0 ,0.1600 ,0.1600 ,0.1600    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN    ,NaN ,0.1600 ,0.1600 ,0.1600 ,0.1600 ,0.1600 ,0.1600 ,0.1600 ,0.1600    ,NaN    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN ,0.1600 ,0.1600 ,0.1600 ,0.1600    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN;
    NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN    ,NaN
    ];
icon_Record(:,:,3)=[
    NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN;
    NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN;
    NaN ,NaN ,NaN ,NaN ,NaN ,NaN   ,0   ,0   ,0   ,0 ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN;
    NaN ,NaN ,NaN ,NaN   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0 ,NaN ,NaN ,NaN ,NaN;
    NaN ,NaN ,NaN   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0 ,NaN ,NaN ,NaN;
    NaN ,NaN ,NaN   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0 ,NaN ,NaN ,NaN;
    NaN ,NaN   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0 ,NaN ,NaN;
    NaN ,NaN   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0 ,NaN ,NaN;
    NaN ,NaN   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0 ,NaN ,NaN;
    NaN ,NaN   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0 ,NaN ,NaN;
    NaN ,NaN ,NaN   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0 ,NaN ,NaN ,NaN;
    NaN ,NaN ,NaN   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0 ,NaN ,NaN ,NaN;
    NaN ,NaN ,NaN ,NaN   ,0   ,0   ,0   ,0   ,0   ,0   ,0   ,0 ,NaN ,NaN ,NaN ,NaN;
    NaN ,NaN ,NaN ,NaN ,NaN ,NaN   ,0   ,0   ,0   ,0 ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN;
    NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN;
    NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN
    ];


icon_RecordAll(:,:,1) = [
 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   1 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN NaN   1   1   1   1   1 NaN NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN   1   1   1   1   1   1   1 NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN   1   1   1   1   1   1   1 NaN NaN NaN NaN NaN NaN NaN NaN;
   1   1   1   1   1   1 NaN   1   1   1   1   1   1   1   1   1 NaN   1   1   1   1   1   1;
 NaN NaN NaN NaN NaN NaN NaN NaN   1   1   1   1   1   1   1 NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN   1   1   1   1   1   1   1 NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN NaN   1   1   1   1   1 NaN NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   1 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];

icon_RecordAll(:,:,2) = [
 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN NaN   0   0   0   0   0 NaN NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN   0   0   0   0   0   0   0 NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN   0   0   0   0   0   0   0 NaN NaN NaN NaN NaN NaN NaN NaN;
   0   0   0   0   0   0 NaN   0   0   0   0   0   0   0   0   0 NaN   0   0   0   0   0   0;
 NaN NaN NaN NaN NaN NaN NaN NaN   0   0   0   0   0   0   0 NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN   0   0   0   0   0   0   0 NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN NaN   0   0   0   0   0 NaN NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];

icon_RecordAll(:,:,3) = icon_RecordAll(:,:,2);

icon_RecordZoom(:,:,1) = [
    NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   1 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
NaN NaN   1 NaN NaN NaN NaN NaN NaN   1   1   1   1   1 NaN NaN NaN NaN NaN NaN   1 NaN NaN;
NaN NaN   1 NaN NaN NaN NaN NaN   1   1   1   1   1   1   1 NaN NaN NaN NaN NaN   1 NaN NaN;
NaN NaN   1 NaN NaN NaN NaN NaN   1   1   1   1   1   1   1 NaN NaN NaN NaN NaN   1 NaN NaN;
  1   1   1   1   1   1 NaN   1   1   1   1   1   1   1   1   1 NaN   1   1   1   1   1   1;
NaN NaN   1 NaN NaN NaN NaN NaN   1   1   1   1   1   1   1 NaN NaN NaN NaN NaN   1 NaN NaN;
NaN NaN   1 NaN NaN NaN NaN NaN   1   1   1   1   1   1   1 NaN NaN NaN NaN NaN   1 NaN NaN;
NaN NaN   1 NaN NaN NaN NaN NaN NaN   1   1   1   1   1 NaN NaN NaN NaN NaN NaN   1 NaN NaN;
NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   1 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];

icon_RecordZoom(:,:,2) = [
    NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
NaN NaN   0 NaN NaN NaN NaN NaN NaN   0   0   0   0   0 NaN NaN NaN NaN NaN NaN   0 NaN NaN;
NaN NaN   0 NaN NaN NaN NaN NaN   0   0   0   0   0   0   0 NaN NaN NaN NaN NaN   0 NaN NaN;
NaN NaN   0 NaN NaN NaN NaN NaN   0   0   0   0   0   0   0 NaN NaN NaN NaN NaN   0 NaN NaN;
  0   0   0   0   0   0 NaN   0   0   0   0   0   0   0   0   0 NaN   0   0   0   0   0   0;
NaN NaN   0 NaN NaN NaN NaN NaN   0   0   0   0   0   0   0 NaN NaN NaN NaN NaN   0 NaN NaN;
NaN NaN   0 NaN NaN NaN NaN NaN   0   0   0   0   0   0   0 NaN NaN NaN NaN NaN   0 NaN NaN;
NaN NaN   0 NaN NaN NaN NaN NaN NaN   0   0   0   0   0 NaN NaN NaN NaN NaN NaN   0 NaN NaN;
NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN   0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];

icon_RecordZoom(:,:,3) = icon_RecordZoom(:,:,2);

icon_Stop = [
 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN   0   0   0   0   0   0   0 NaN NaN;
 NaN NaN   0   0   0   0   0   0   0 NaN NaN;
 NaN NaN   0   0   0   0   0   0   0 NaN NaN;
 NaN NaN   0   0   0   0   0   0   0 NaN NaN;
 NaN NaN   0   0   0   0   0   0   0 NaN NaN;
 NaN NaN   0   0   0   0   0   0   0 NaN NaN;
 NaN NaN   0   0   0   0   0   0   0 NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN;
 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN];
 icon_Stop = repmat(icon_Stop,[1 1 3]);
copyRoiBt = [
    NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN   NaN     0     0     0     0     0     0     0     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN     1     1     1     1     1     1     1     0     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN   NaN     0     0     0     0     0     0     0     0     0     0   NaN   NaN     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1     1     1     1     1     1     1     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];

copyRoiBt(:,:,2) = [
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN   NaN     0     0     0     0     0     0     0     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN     0     0     0     0     0     0     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN   NaN     0     0     0     0     0     0     0     0     0     0   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0     0     0     0     0     0     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];

copyRoiBt(:,:,3) = [
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN   NaN     0     0     0     0     0     0     0     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN     0     0     0     0     0     0     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN   NaN     0     0     0     0     0     0     0     0     0     0   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     0     0     0     0     0     0     0     0     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN];




% hand3 = ...
%  [NaN   NaN   NaN   NaN   NaN   NaN   NaN    0     0   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
%  NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
%  NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
%  NaN   NaN   NaN   NaN   NaN   NaN   NaN     1     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN;
%  NaN   NaN   NaN   NaN     0     0     0     1     1   NaN     0     0   NaN   NaN   NaN   NaN;
%  NaN   NaN   NaN   NaN     1     1     0     1     1     0     1     1     0     0   NaN   NaN;
%  NaN   NaN     0   NaN     1     1     1     1     1     1     1     1     1     1     0   NaN;
%  NaN     0     1     0     1     1     1     1     1     1     1     1     1     1     0   NaN;
%  NaN     0     1     1     1     1     1     1     1     1     1     1     1     1     0   NaN;
%  NaN   NaN     1     1     1     1     1     1     1     1     1     1     1     1     0   NaN;
%  NaN   NaN     0     1     1     1     1     1     1     1     1     1     1     0   NaN   NaN;
%  NaN   NaN     0     1     1     1     1     1     1     1     1     1     1     0   NaN   NaN;
%  NaN   NaN   NaN     0     1     1     1     1     1     1     1     1     1     0   NaN   NaN;
%  NaN   NaN   NaN   NaN     0     1     1     1     1     1     1     1     0   NaN   NaN   NaN;
%  NaN   NaN   NaN   NaN   NaN     0     1     1     1     1     1     1     0   NaN   NaN   NaN;
%  NaN   NaN   NaN   NaN   NaN     0     1     1     1     1     1     1     0   NaN   NaN   NaN;];
%
% hand3(hand3 == 1) = 2;
% hand3(hand3 == 0) = 1;

% icon_book(:,:,1) = [ ...
% 	NaN NaN NaN 128 NaN NaN NaN NaN NaN NaN NaN 128 NaN NaN NaN ; ...
% 	NaN 128 128 NaN 128 NaN NaN NaN NaN NaN 128 NaN 128 128 NaN ; ...
% 	  0 128 NaN NaN NaN 128 NaN NaN NaN 128 NaN NaN NaN 128   0 ; ...
% 	  0 128 NaN   0 NaN NaN 128 NaN 128 NaN NaN   0 NaN 128   0 ; ...
% 	  0 128 NaN NaN   0 NaN NaN 128 NaN NaN   0 NaN NaN 128   0 ; ...
% 	  0 128 NaN   0 NaN NaN NaN 128 NaN   0 NaN NaN NaN 128   0 ; ...
% 	  0 128 NaN NaN   0 NaN NaN 128 NaN NaN   0 NaN NaN 128   0 ; ...
% 	  0 128 NaN NaN NaN   0 NaN 128 NaN   0 NaN NaN NaN 128   0 ; ...
% 	  0 128 128 128 NaN NaN NaN 128 NaN NaN NaN 128 128 128   0 ; ...
% 	  0 128 128 NaN 128 NaN NaN 128 NaN NaN 128 NaN 128 128   0 ; ...
% 	  0   0 128 128 NaN 128 NaN 128 NaN 128 NaN 128 128   0   0 ; ...
% 	NaN NaN   0   0 128 128 128 128 128 128 128   0   0 NaN NaN ; ...
% 	NaN NaN NaN NaN   0   0 128 128 128   0   0 NaN NaN NaN NaN ; ...
% 	NaN NaN NaN NaN NaN NaN   0   0   0 NaN NaN NaN NaN NaN NaN ; ...
% 	]/255;
% icon_book(:,:,2) = [ ...
% 	NaN NaN NaN   0 NaN NaN NaN NaN NaN NaN NaN   0 NaN NaN NaN ; ...
% 	NaN   0   0 NaN   0 NaN NaN NaN NaN NaN   0 NaN   0   0 NaN ; ...
% 	  0   0 NaN NaN NaN   0 NaN NaN NaN   0 NaN NaN NaN   0   0 ; ...
% 	  0   0 NaN   0 NaN NaN   0 NaN   0 NaN NaN   0 NaN   0   0 ; ...
% 	  0   0 NaN NaN   0 NaN NaN   0 NaN NaN   0 NaN NaN   0   0 ; ...
% 	  0   0 NaN   0 NaN NaN NaN   0 NaN   0 NaN NaN NaN   0   0 ; ...
% 	  0   0 NaN NaN   0 NaN NaN   0 NaN NaN   0 NaN NaN   0   0 ; ...
% 	  0   0 NaN NaN NaN   0 NaN   0 NaN   0 NaN NaN NaN   0   0 ; ...
% 	  0   0   0   0 NaN NaN NaN   0 NaN NaN NaN   0   0   0   0 ; ...
% 	  0   0   0 NaN   0 NaN NaN   0 NaN NaN   0 NaN   0   0   0 ; ...
% 	  0   0   0   0 NaN   0 NaN   0 NaN   0 NaN   0   0   0   0 ; ...
% 	NaN NaN   0   0   0   0   0   0   0   0   0   0   0 NaN NaN ; ...
% 	NaN NaN NaN NaN   0   0   0   0   0   0   0 NaN NaN NaN NaN ; ...
% 	NaN NaN NaN NaN NaN NaN   0   0   0 NaN NaN NaN NaN NaN NaN ; ...
% 	]/255;
% icon_book(:,:,3) = [ ...
% 	NaN NaN NaN   0 NaN NaN NaN NaN NaN NaN NaN   0 NaN NaN NaN ; ...
% 	NaN   0   0 NaN   0 NaN NaN NaN NaN NaN   0 NaN   0   0 NaN ; ...
% 	  0   0 NaN NaN NaN   0 NaN NaN NaN   0 NaN NaN NaN   0   0 ; ...
% 	  0   0 NaN   0 NaN NaN   0 NaN   0 NaN NaN   0 NaN   0   0 ; ...
% 	  0   0 NaN NaN   0 NaN NaN   0 NaN NaN   0 NaN NaN   0   0 ; ...
% 	  0   0 NaN   0 NaN NaN NaN   0 NaN   0 NaN NaN NaN   0   0 ; ...
% 	  0   0 NaN NaN   0 NaN NaN   0 NaN NaN   0 NaN NaN   0   0 ; ...
% 	  0   0 NaN NaN NaN   0 NaN   0 NaN   0 NaN NaN NaN   0   0 ; ...
% 	  0   0   0   0 NaN NaN NaN   0 NaN NaN NaN   0   0   0   0 ; ...
% 	  0   0   0 NaN   0 NaN NaN   0 NaN NaN   0 NaN   0   0   0 ; ...
% 	  0   0   0   0 NaN   0 NaN   0 NaN   0 NaN   0   0   0   0 ; ...
% 	NaN NaN   0   0   0   0   0   0   0   0   0   0   0 NaN NaN ; ...
% 	NaN NaN NaN NaN   0   0   0   0   0   0   0 NaN NaN NaN NaN ; ...
% 	NaN NaN NaN NaN NaN NaN   0   0   0 NaN NaN NaN NaN NaN NaN ; ...
% 	]/255;


colorbarIcon = load([matlabroot,'/toolbox/matlab/icons/colorbar.mat']);     %#ok Colorbar Icon
colorbarIcon = struct2cell(colorbarIcon);
colorbarIcon = permute(colorbarIcon{1},[2,1  ,3]);
%% Gui Properties
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
%Gui Figure
set(gui, 'Units', 'Pixel', 'Position', customConfig.winPos.gui, 'Name', ['matVis: ', allNames],...
    'MenuBar', 'none', 'Resize', 'off', 'NumberTitle', 'off',...
    'KeyPressFcn', @keyPress, 'CloseRequestFcn', @closeGui, 'Visible','off',...
    'WindowButtonDownFcn', @buttonDownGui,'WindowButtonMotionFcn', '','NextPlot', 'new', 'HandleVisibility', 'off', 'WindowStyle','normal');
% In order to avoid overflowing of axis in uipanels (a knwon Matlab bug,
% the renderer of the figure has to be changed (opengl can give problems
% with axis labels being flipped upside down with some graphic cards).
set(gui, 'Renderer',' zbuffer');
% Set 'WindowButtonMotionFcn' later to avoid errors due to non-existing
% image and zoom windows
try
    gui_jf = get(gui,'JavaFrame');
    gui_jf.setFigureIcon(javax.swing.ImageIcon(im2java(uint8(icon_matVis))));
catch    %#ok
end

panel_positionControls = uipanel(gui, 'units','pixel','Position', [-1 customConfig.winPos.gui(4)-150-(nDim-1)*40 customConfig.winPos.gui(4)+2 130+(nDim-1)*40],...
    'BackgroundColor',get(gui, 'Color'),'BorderType','none');

% Position and plot controls
txt_titles(1) = uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [5 customConfig.winPos.gui(4)-25 customConfig.winPos.gui(3)-10 20], 'String', 'Position and plot controls', ...
    'HorizontalAlignment', 'center','FontWeight','bold','FontSize',10);
txt_plots = uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [290 customConfig.winPos.gui(4)-30 25 15], 'String', 'Plots', ...
    'HorizontalAlignment', 'left');
%Dimension Controls
for i = 1:nDim
    %Text for dimension names
    txt_dimName(i) = uicontrol('Parent',panel_positionControls , 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
        'Position', [5 (nDim-i+2)*40+18 40 17], 'String', dimNames(i), 'FontWeight','bold',...
        'HorizontalAlignment', 'left','ForegroundColor', [(i==1) 0.8*(i==2) 0] );  %#ok
    %Slider for current values
    sld(i) = uicontrol('Parent',panel_positionControls , 'Style', 'Slider', 'Callback', @sliderCallback, 'Units', 'Pixel', ...
        'Position', [40 (nDim-i+2)*40+22 105 16], 'Min', 1  , 'Max', dim(i), ...
        'SliderStep', [1/(dim(i)-1) 10/(dim(i)-1)],'Value',currPos(i),'Userdata', i,'Tooltip',['Set position for dimension ''',dimNames{i},'''. You can use the mouse scroll wheel to move the slider position.'],'Tag',['Set position for dimension ''',dimNames{i},'''. You can use the mouse scroll wheel to move the slider position.']);  %#ok
    sldPos(:,i) = get(sld(i),'Position')+get(panel_positionControls,'Position').*[1 1 0 0];
    if usejava('awt')  % java enabled -> use it to update while dragging
      if oldMATLAB
        hListeners(i) = handle.listener(sld(i), 'ActionEvent', @sliderCallback);  %#ok
      else
        hListeners(i) = addlistener(sld(i), 'ContinuousValueChange', @sliderCallback);  %#ok
      end
        if i == nDim
            set(sld, 'Callback', '');
            setappdata(gui,'sliderListeners', hListeners); % -> hListeners can be list of handles
        end
    end
    %etxt windows for current values
    etxt(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Edit', 'Callback', {@etxtCallback,i}, 'Units', 'Pixel', ...
        'Position', [200 (nDim-i+2)*40+22 50 16], 'String', num2str(currPos(i)),...
        'Tag',['Set position for dimension ''',dimNames{i},'''']);  %#ok
    %Text for lengths of dimensions
    dimSize(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
        'Position', [250 (nDim-i+2)*40+20 75 17], 'String', [' / ', num2str(dim(i))], ...
        'HorizontalAlignment', 'left','Tooltip',['Size of data set along dimension ''',dimNames{i},''''],'Tag',['Size of data set along dimension ''',dimNames{i},'''']);  %#ok
    % Slider for Zoom Values
    sld_down(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Slider', 'Units', 'Pixel', ...
        'Position', [40 (nDim-i+2)*40+15 105 8], 'Min', 1  , 'Max', dim(i), ...
        'SliderStep', [1/(dim(i)-1) 10/(dim(i)-1)],'Value',1,'Userdata', i,...
        'Callback',{@updateZoom,i,'sld'},'Tag',['Minimum of zoom range for dimension ''',dimNames{i},'''']);  %#ok
    sld_up(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Slider', 'Units', 'Pixel', ...
        'Position', [40 (nDim-i+2)*40+8 105 8], 'Min', 1  , 'Max', dim(i), ...
        'SliderStep', [1/(dim(i)-1) 10/(dim(i)-1)],'Value',dim(i),'Userdata', i,...
        'Callback',{@updateZoom,i,'sld'},'Tag',['Maximum of zoom range for dimension ''',dimNames{i},'''']);  %#ok
    % Etxt for zoom values
    etxt_down(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Edit', 'Callback', {@etxtCallback,i}, 'Units', 'Pixel', ...
        'Position', [200 (nDim-i+2)*40+8 26 15], 'String', zoomVal(i,1),'FontSize',7,...
        'Callback',{@updateZoom,i,'etxt'},'Tooltip',['Minimum of zoom range for dimension ''',dimNames{i},''''],'Tag',['Minimum of zoom range for dimension ''',dimNames{i},'''']);  %#ok
    etxt_up(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Edit', 'Callback', {@etxtCallback,i}, 'Units', 'Pixel', ...
        'Position', [224 (nDim-i+2)*40+8 26 15], 'String', zoomVal(i,2),'FontSize',7,...
        'Callback',{@updateZoom,i,'etxt'},'Tooltip',['Maximum of zoom range for dimension ''',dimNames{i},''''],'Tag',['Maximum of zoom range for dimension ''',dimNames{i},'''']);  %#ok
    %Play button
    btPlayAll(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Togglebutton', 'Callback', {@playCallback,i,'all',1}, 'Units', 'Pixel', ...
        'Position', [150 (nDim-i+2)*40+23 30 15], 'Value',0  , 'CData', arrowAll,...
        'BackgroundColor', get(gui, 'Color'),'Tooltip',['Play along dimension ''',dimNames{i},'''. Right click for reverse play.'],'Tag',['Play along dimension ''',dimNames{i},'''. Right click for reverse play.']);  %#ok
    btPlayZoom(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Togglebutton', 'Callback', {@playCallback,i,'zoom',1}, 'Units', 'Pixel', ...
        'Position', [150 (nDim-i+2)*40+8 30 15], 'Value',0  , 'CData', arrowZoom,...
        'BackgroundColor', get(gui, 'Color'),'Tooltip',['Play across zoom range of dimension ''',dimNames{i},'''. Right click for reverse play.'],'Tag',['Play across zoom range of dimension ''',dimNames{i},'''. Right click for reverse play.']);  %#ok
    % Zoom Width
    txt_zoomWidth(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
        'Position', [250 (nDim-i+2)*40+5 75 17], 'String', [' = ', num2str(zoomVal(i,2))], ...
        'HorizontalAlignment', 'left','FontSize',7,'Tooltip',['Width of zoom range of  ''',dimNames{i},''''],'Tag',['Width of zoom range of  ''',dimNames{i},'''']);  %#ok
    if customDimScale
        txt_zoomDownScaled(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
            'Position', [200 (nDim-i+2)*40+8-10 26 10], 'String', num2str(pxWidth(i)*round(dimScale(i,1)/pxWidth(i)),'%6.1f'), ...
            'HorizontalAlignment', 'center','FontSize',7,'Tooltip',['Minimum of zoom range for dimension ''',dimNames{i},''' in dimension units'],'Tag',['Minimum of zoom range for dimension ''',dimNames{i},''' in dimension units']);  %#ok
        txt_zoomUpScaled(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
            'Position', [224 (nDim-i+2)*40+8-10 26 10], 'String', num2str(pxWidth(i)*round(dimScale(i,2)/pxWidth(i)),'%6.1f'), ...
            'HorizontalAlignment', 'center','FontSize',7,'Tooltip',['Maximum of zoom range for dimension ''',dimNames{i},''' in dimension units'],'Tag',['Maximum of zoom range for dimension ''',dimNames{i},''' in dimension units']);  %#ok
        txt_zoomWidthScale(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
            'Position', [250 (nDim-i+2)*40+8-10 75 10], 'String', [' = ', num2str(pxWidth(i)*zoomVal(i,2),'%6.2f') ' ' dimUnits{i}], ...
            'HorizontalAlignment', 'left','FontSize',7,'Tooltip',['Width of zoom range of  ''',dimNames{i},''' in dimension units.'],'Tag',['Width of zoom range of  ''',dimNames{i},''' in dimension units.']);  %#ok
        txt_pxScale(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
            'Position', [40 (nDim-i+2)*40+8-10 105 10], 'String', [num2str(pxWidth(i),'%6.2f') ' '  dimUnits{i} '/px'], ...
            'HorizontalAlignment', 'center','FontSize',7,'Tooltip',['Pixel size of ''',dimNames{i},''''],'Tag',['Pixel size of ''',dimNames{i},'''']);  %#ok
    end
    % Checkboxes for displaying plots
    cbPlots(i) = uicontrol('Parent',panel_positionControls, 'Style', 'checkbox', 'Callback', {@cbCallback,i}, 'Units', 'Pixel', ...
        'Position', [295 (nDim-i+2)*40+23 20 17], 'Value', any(i == customConfig.plotDim), 'String', '',...
        'BackgroundColor', get(gui, 'Color'),'Tooltipstring', 'Show/hide plot along respective dimension.',...
        'Tag',['Show/hide plot along dimension''',dimNames{i},'''']);  %#ok
    % Toggle button to lock plot x-axes to zoom values
    tb_lockPlots2Zoom(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Togglebutton', 'Callback', {@updateZoom,i,'lockBt'}, 'Units', 'Pixel', ...
        'Position', [295 (nDim-i+2)*40+9 15 15], 'Value',0  , 'CData', lockOpen,...
        'BackgroundColor', get(gui, 'Color'),'Tooltipstring', sprintf('Locks x-lim of plot axes to zoom interval. Requires pressed ''x-lim button'' (under plot controls)\nIn RGB-Stretch mode it limits the stretch range.'),...
        'Tag', ['Lock x-lim of plot axes to zoom interval.' char(10) 'Requires pressed ''x-lim button'' (under plot controls).' char(10) 'In RGB-Stretch mode it limits the stretch range.']);  %#ok
end
if customDimScale
    set(panel_positionControls, 'Children',[tb_lockPlots2Zoom cbPlots txt_zoomWidthScale txt_zoomDownScaled txt_zoomUpScaled txt_pxScale txt_zoomWidth btPlayZoom btPlayAll sld_up sld_down sld txt_dimName dimSize flipdim(reshape((cat(1,etxt_down,etxt_up)),[1 2*nDim]),2) etxt(nDim:-1:1)]);
else
    set(panel_positionControls, 'Children',[tb_lockPlots2Zoom cbPlots txt_zoomWidth btPlayZoom btPlayAll sld_up sld_down sld txt_dimName dimSize flipdim(reshape((cat(1,etxt_down,etxt_up)),[1 2*nDim]),2) etxt(nDim:-1:1)]);
end
sldPlaySpeed = uicontrol('Parent',panel_positionControls, 'Style', 'Slider', 'Units', 'Pixel', ...
    'Position', [185 88 10 nDim*dimHeight-10], 'Min', -10, 'Max', 10, ...
    'SliderStep', [0.05 0.2],'Value',0, ...
    'Callback', @playSpeedCallback, 'TooltipString','Playspeed: 0',...
    'Tag',['Speed of playback. Current value: 0.' char(10) 'Negative values: faster speed by skipping over images.' char(10) 'Positive values: slower speed by introducing pauses.'] );

% Find Min/Max Buttons
uicontrol('Parent',panel_positionControls, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [1 45 38 25], 'String', {'Set/Get';'Position'}, ...
    'HorizontalAlignment', 'center','FontSize',7);  %#ok
for i=1:4
    bt_setJumpPos(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Pushbutton', 'Callback', {@setJumpPos,i}, 'Units', 'Pixel', 'FontSize',7,...
        'Position', [40+(i-1)*20 60 20 15], 'ForeGroundColor','r', 'String', num2str(i),...
        'BackgroundColor', get(gui, 'Color'),'Tooltipstring', 'Click to save current position.',...
        'Tag', 'Click to save current position.');  %#ok
    bt_getJumpPos(i) = uicontrol('Parent',panel_positionControls, 'Style', 'Pushbutton', 'Callback', {@jumpPos,i}, 'Units', 'Pixel', 'FontSize',7,...
        'Position', [40+(i-1)*20 45 20 15], 'ForeGroundColor',[0 .7 0], 'String', num2str(i),...
        'BackgroundColor', get(gui, 'Color'),'Tooltipstring', 'Click to go saved position','Enable','off',...
        'Tag','Click to go to saved position');  %#ok
end
% Lock to zoom
bt_savedZoom = uicontrol('Parent',panel_positionControls, 'Style', 'Togglebutton', 'Units', 'Pixel', 'FontSize',7,...
    'Position', [120 45 38 30], 'String', {'+Zoom'}, 'Value',0,...
    'Tooltipstring', ['Pressed: Apply saved zoom setting together with position.' char(10) 'Unpressed: Only jump to position, leave zoom unchanged.'],...
    'Tag', ['Pressed: Apply saved zoom setting together with position.' char(10) 'Unpressed: Only jump to position, leave zoom unchanged.']);
% Go to global min/max
bt_jumpMin = uicontrol('Parent',panel_positionControls, 'Style', 'Pushbutton', 'Callback', {@jumpPos,'min'}, 'Units', 'Pixel', ...
    'Position', [160 60 50 15],  'String', 'Min-global', 'FontSize',7,...
    'BackgroundColor', get(gui, 'Color'),'Tooltipstring', ['Jumps to (first) position with minimum value: [' num2str(minPos) ']'],...
    'Tag',['Jump to (first) position with the minimum value of the data set: [' num2str(minPos) ']']);  %#ok
bt_jumpMax = uicontrol('Parent',panel_positionControls, 'Style', 'Pushbutton', 'Callback', {@jumpPos,'max'}, 'Units', 'Pixel', 'FontSize',7, ...
    'Position', [160 45 50 15],  'String', 'Max-global', ...
    'BackgroundColor', get(gui, 'Color'),'Tooltipstring', ['Jumps to (first) position with maximum value: [' num2str(maxPos) ']'],...
    'Tag',['Jump to (first) position with the maximum value of the data set: [' num2str(minPos) ']']);  %#ok
% Push button: find min/max in zoom/image
tb_findLocalMinImg = uicontrol('Parent',panel_positionControls, 'Style', 'pushbutton', 'Callback', {@findExtremum,'min','image'}, 'Units', 'Pixel', ...
    'Position', [213 60 50 15],  'String', 'Min-Img','TooltipString', 'Find minimum in current image',...
    'BackgroundColor', get(gui, 'Color'), 'Value', customConfig.plotZoom,'Tag', 'Find minimum in current image');
tb_findLocalMaxImg = uicontrol('Parent',panel_positionControls, 'Style', 'pushbutton', 'Callback', {@findExtremum,'max','image'}, 'Units', 'Pixel', ...
    'Position', [213 45 50 15],  'String', 'Max-Img','TooltipString', 'Find maximum in current image',...
    'BackgroundColor', get(gui, 'Color'), 'Value', customConfig.plotZoom,'Tag', 'Find maximum in current image');
tb_findLocalMinZm = uicontrol('Parent',panel_positionControls, 'Style', 'pushbutton', 'Callback', {@findExtremum,'min','zoom'}, 'Units', 'Pixel', ...
    'Position', [266 60 50 15],  'String', 'Min-Zoom','TooltipString', 'Find minimum in current zoom',...
    'BackgroundColor', get(gui, 'Color'), 'Value', customConfig.plotZoom,'Tag', 'Find minimum in current zoom');
tb_findLocalMaxZm = uicontrol('Parent',panel_positionControls, 'Style', 'pushbutton', 'Callback', {@findExtremum,'max','zoom'}, 'Units', 'Pixel', ...
    'Position', [266 45 50 15],  'String', 'Max-Zoom','TooltipString', 'Find maximum in current zoom',...
    'BackgroundColor', get(gui, 'Color'), 'Value', customConfig.plotZoom,'Tag', 'Find maximum in current zoom');

%Toggle button for applying zoom plots
tbPlotsXLim = uicontrol('Parent',panel_positionControls, 'Style', 'togglebutton', 'Callback', @updatePlots, 'Units', 'Pixel', ...
    'Position', [40 10 24 24],  'String', '','TooltipString', ['Set x-limits of plots.' char(10) '- Left click: Toggle auto / fixed.' char(10) '- Right click: Set fixed limits.'],...
    'BackgroundColor', get(gui, 'Color'), 'Value', customConfig.plotZoom, 'CData', permute(arrows,[2,1  ,3]),'Tag', ['Set x-limits of plots.' char(10) 'Left click: Toggle auto / fixed.' char(10) 'Right click: Set fixed limits.']);
%Toggle button to (y-)scale plots (otherwise scaled to min-max of complete
%data set
tbPlotsYLim = uicontrol('Parent',panel_positionControls, 'Style', 'togglebutton', 'Callback', @updatePlots, 'Units', 'Pixel', ...
    'Position', [70 10 24 24],  'String', '','TooltipString', ['Set y-limits of plots.' char(10) '- Left click: Toggle auto / fixed.' char(10) '- Right click: Set fixed limits.'],...
    'BackgroundColor', get(gui, 'Color'), 'Value', customConfig.plotScale, 'CData', arrows,...
    'Tag', ['Set y-limits of plots.' char(10) 'Left click: Toggle auto / fixed.' char(10) 'Right click: Set fixed limits.']);
%Toggle button to display marker in plots
tbMarker = uicontrol('Parent',panel_positionControls, 'Style', 'togglebutton', 'Callback', @updatePlots, 'Units', 'Pixel', ...
    'Position', [100 10 24 24],  'String', '', 'TooltipString', 'Show / hide marker in plots', ...
    'BackgroundColor', get(gui, 'Color'), 'Value', customConfig.marker, 'CData', marker, 'Tag', 'Show / hide marker in plots');
%Mean Plots
btMean = uicontrol('Parent',panel_positionControls, 'Style', 'Pushbutton', 'Units', 'Pixel', 'UserData', customConfig.plotMean, ...
    'Position', [130 10 24 24], 'CData', meanBt{customConfig.plotMean+1}, 'Callback', @meanPlots,...
    'ToolTipString', ['Average plot values over 1, 9, 25 pixels or the complete zoom area.' char(10) 'In stretch RGB mode all plots along the RGB dimension can be displayed alternatively.'],...
    'Tag', ['Average plot values over 1, 9, 25 pixels or the complete zoom area. In stretch RGB mode all plots along the' char(10) 'RGB dimension can be displayed alternatively.']);

% Features
% uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
%     'Position', [5 110 customConfig.winPos.gui(3)-10 20], 'String', 'Features', ...
%     'HorizontalAlignment', 'center','FontWeight','bold','FontSize',10);

%Button for Data Export
btExport = uicontrol('Parent',panel_positionControls, 'Style', 'Pushbutton', 'Units', 'Pixel', 'CData', exportBt,  ...
    'Position', [200 10 24 24], 'Callback', @exportData,  'Value', 0,...
    'ToolTipString', 'Export subset of data','Tag','Export subset of data either as a new variable in the workspace or into a new matVis.');  %[115 85 24 24]

%Toggle button for display of menu bars in all windows (excpett gui)
tbMenuBars = uicontrol('Parent', panel_positionControls, 'Style', 'Togglebutton', 'Units', 'Pixel', 'Value', customConfig.menuBarVis, ...
    'Position', [230 10 24 24], 'CData', menuBarIcon, 'Callback', @toggleMenuBars, ...
    'TooltipString', 'Toggle display of menu bars in Image/Zoom/Plot windows', ...
    'Tag', 'Toggle display of menu bars in Image/Zoom/Plot windows');


%Button to link/unlink window position/size
tbLinkWin = uicontrol('Parent',panel_positionControls, 'Style', 'Togglebutton', 'Units', 'Pixel', 'CData', icon_linkWins,  ...
    'Position', [260 10 24 24], 'Callback', @linkWins,  'Value', 1, 'Tag', 'Link window positions and sizes'); %[145 85 24 24]
if numel(data) == 1
    set(tbLinkWin, 'Visible','off', 'Value', 0);
end

%Button for display of CustomTif Parameter
tbTifPar = uicontrol('Parent',panel_positionControls, 'Style', 'Togglebutton', 'Units', 'Pixel', 'CData', tifParBt,  ...
    'Position', [290 10 24 24], 'Callback', @showTifPar,  'Value', 0, 'Tag', 'Display parameters of custom tif files used in the Schild lab (University of Göttingen).'); %[145 85 24 24]
if ~isCustomTif
    set(tbTifPar, 'Visible', 'off');
end

%Separator line
sepLine(1) = uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [0 350 customConfig.winPos.gui(3) 9], 'String','_______________________________________________________________________________________________________________', ...
    'HorizontalAlignment', 'left', 'FontWeight', 'Bold','Fontsize',6);

% Image controls
txt_titles(end+1) = uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [5 327 customConfig.winPos.gui(3)-10 20], 'String', 'Image controls', ...
    'HorizontalAlignment', 'center','FontWeight','bold','FontSize',10);
panel_imageControls = uipanel(gui, 'units','pixel','Position', [-1 135 customConfig.winPos.gui(4)+2 190],...
    'BackgroundColor',get(gui, 'Color'),'BorderType','none');  %
% xText
uicontrol('Parent', panel_imageControls, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [20 175 10 20], 'String', 'x', ...
    'HorizontalAlignment', 'left', 'Tag', 'Select dimension of data set to be used as the ''x-dimension'' of the displayed images');
% XDir control
bt_XDir = uicontrol('Parent', panel_imageControls, 'Style', 'pushbutton', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [42 181 13 13], 'String', '>', 'Callback',{@setAxisDir,'x'}, ...
    'HorizontalAlignment', 'left', 'Tag', 'x-axis direction of images', 'Tooltip', 'x-axis direction of images');
% yText
uicontrol('Parent', panel_imageControls, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [70 175 10 20], 'String', 'y', ...
    'HorizontalAlignment', 'left', 'Tag', 'Select dimension of data set to be used as the ''y-dimension'' of the displayed images');
% YDir control
bt_YDir = uicontrol('Parent', panel_imageControls, 'Style', 'pushbutton', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [92 181 13 13], 'String', 'v', 'Callback',{@setAxisDir,'y'}, ...
    'HorizontalAlignment', 'left', 'Tag', 'y-axis direction of images', 'Tooltip', 'y-axis direction of images');
% Switch position of pop1 and pop2 to place pop1 under the "y" direction
% and pop2 under the "x" direction. (Changed in v1.021
pop(1)  = uicontrol('Parent', panel_imageControls, 'Style', 'popupmenu', 'Callback', @popCallback, 'Units', 'Pixel', ...
    'Position', [56 160 50 20], 'String', dimNames, 'Value',xySel(1),'FontSize',7, 'Tag', 'Select dimension of data set to be used as the ''x-dimension'' of the displayed images');
pop(2)  = uicontrol('Parent', panel_imageControls, 'Style', 'popupmenu', 'Callback', @popCallback, 'Units', 'Pixel', ...
    'Position', [5 160 50 20], 'String', dimNames, 'Value',xySel(2),'FontSize',7, 'Tag', 'Select dimension of data set to be used as the ''y-dimension'' of the displayed images');
if withDipimage || withImageProcessingTB
    withFilter = 1;
    % Filter types
    if withDipimage
        if withImageProcessingTB
            s = {'none','gauss','median','mean','unsharp-ml','unsharp-dip'};
        else
            s = {'none','gauss','median','mean','unsharp-dip'};
            
        end
    else
        s = {'none','gauss','median','mean','unsharp-ml'};
    end
    % Filter
    uicontrol('Parent', panel_imageControls, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
        'Position', [125 173 30 20], 'String', 'Filter', ...
        'HorizontalAlignment', 'left', 'Tag', 'Select filter function.');
    popFilter = uicontrol('Parent', panel_imageControls, 'Style', 'popupmenu', 'Callback', @updateFilter, 'Units', 'Pixel', ...
        'Position', [108 160 60 20], 'String',s, 'Value',xySel(1),'FontSize',7, 'Tag', 'Select filter function','Tooltipstring','Select filter function');
    % Filter size
    txtFilterPar = uicontrol('Parent', panel_imageControls, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
        'Position', [167 173 32 20], 'String', 'Par', ...
        'HorizontalAlignment', 'center', 'Tag', 'Enter filter size (0 = no filter applied)');
    if withDipimage
        set(popFilter, 'Tag', 'Select filter function. Functions from DIPImage will be used.','Tooltip', 'Select filter function. Functions from DIPImage will be used.')
    else
        set(popFilter, 'Tag', 'Select filter function. Functions from Matlab''s Image Processing Toolbox will be used.', 'Tooltip', 'Select filter function. Functions from Matlab''s Image Processing Toolbox will be used.')
    end
    
    editFilterPar = uicontrol('Parent', panel_imageControls, 'Style', 'edit', 'Callback', @updateFilter, 'Units', 'Pixel', ...
        'Position', [168 160 30 20], 'String', '1','FontSize',7, 'Tag', 'Enter filter size in pixels (0 = no filter applied). Unsharp filter: alpha = size /10 for size <= 10, alpha = 1 for size > 10.','Tooltipstring','Enter filter size in pixels (0 = no filter applied)');
else
    withFilter = 0;
end

%Projection popup and button
projText = uicontrol('Parent', panel_imageControls, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [200 173 50 20], 'String', 'Projection', ...
    'HorizontalAlignment', 'left', 'Tag', 'Display projection of data along a certain dimension.');
projDimText = uicontrol('Parent', panel_imageControls, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [255 173 50 20], 'String', 'Proj. dim', ...
    'HorizontalAlignment', 'left', 'Tag', 'Display projection of data along a certain dimension.');
projMethodPop  = uicontrol('Parent', panel_imageControls, 'Style', 'popupmenu', 'Callback', {@projCallback,'method'}, 'Units', 'Pixel','FontSize',7, ...
    'Position', [200 160 50 20], 'String', {'none';'max';'min';'mean';'std';'var'; 'tile'}, 'Value',projMethod+1,...
    'Tag', 'Select type of projection.','Tooltipstring', 'Select type of projection.');
projDimPop  = uicontrol('Parent', panel_imageControls,  'Style', 'popupmenu','Callback', {@projCallback,'dim'}, 'Units', 'Pixel', ...
    'Position', [250 160 50 20],  'Value',1 ,'Tooltipstring', 'Select dimension along which data will be projected' ,...
    'Enable','on', 'String','Dim 3','FontSize',7,'Tag', 'Select dimension along which to project.');
if nDim < 3
    set([projMethodPop projText projDimText projDimPop], 'Enable', 'off');
else
    set(projDimPop, 'String', dimNames(3:end));  % ToDo: Adjust so it can be used with start parameter
end
bt_zoomProj = uicontrol('Parent', panel_imageControls, 'Style', 'Togglebutton', 'Units', 'Pixel', 'Callback', @zoomProj,...
    'Position', [303 163 15 15], 'Value',0  , 'CData', lockOpen,...
    'BackgroundColor', get(gui, 'Color'),'Tooltipstring', 'If pressed only data inside zoom interval of the respective dimension are included in the projection.',...
    'Tag', 'If pressed only data inside zoom interval of the respective dimension are included in the projection.');

%Button Group Colormap
bg_colormap = uibuttongroup('Parent', panel_imageControls, 'Title', '', 'Units', 'Pixel', ...
    'Position', [0 140 customConfig.winPos.gui(3) 24], 'BackgroundColor', get(gui, 'Color'),...
    'TitlePosition', 'centertop', 'BorderType', 'none', 'HighlightColor', 'w',...
    'SelectionChangeFcn',@updateColormap,'Tag', 'Select contrast adjustment. Global: min/max of complete data set. Image: min/max of currently displayed image. Zoom: min /max of current zoom in displayed image. Manual: Select min/max using slider or edit control below.');
% Contrast Text
% uicontrol('Style', 'Text', 'Parent', bg_colormap, 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
%     'Position', [10 3 50 15], 'String', 'Contrast', ...
%     'HorizontalAlignment', 'left');
cmGlobal = uicontrol('Style', 'Radio', 'String', 'Global', 'Parent', bg_colormap, ...
    'Units', 'pixel', 'Position', [10 3 60 20], 'Value', strcmp(customConfig.colormapMode,'Global'),... %'Tag', 'yScaleOff',...
    'BackgroundColor', get(gui, 'Color'), 'HorizontalAlignment', 'left','Tag', 'Select contrast adjustment. Global: min/max of complete data set. Image: min/max of currently displayed image. Zoom: min /max of current zoom in displayed image. Manual: Select min/max using slider or edit control below.');
cmImage = uicontrol('Style', 'Radio', 'String', 'Image', 'Parent', bg_colormap, ...
    'Units', 'pixel', 'Position', [65 3 60 20], ...  %'Tag', 'yScaleLocal'
    'BackgroundColor', get(gui, 'Color'), 'Value', strcmp(customConfig.colormapMode,'Local'), 'HorizontalAlignment', 'left');
cmZoom = uicontrol('Style', 'Radio', 'String', 'Zoom', 'Parent', bg_colormap, ...
    'Units', 'pixel', 'Position', [120 3 60 20], ...  %'Tag', 'yScaleLocal'
    'BackgroundColor', get(gui, 'Color'), 'Value', strcmp(customConfig.colormapMode,'Local'), 'HorizontalAlignment', 'left');
cmManual = uicontrol('Style', 'Radio', 'String', 'Manual', 'Parent', bg_colormap, ...
    'Units', 'pixel', 'Position', [175 3 60 20], ... %'Tag', 'yScaleGlobal',...
    'BackgroundColor', get(gui, 'Color'), 'HorizontalAlignment', 'left', 'Value', strcmp(customConfig.colormapMode,'Manual'));
cmThresh = uicontrol('Style', 'Radio', 'String', '', 'Parent', bg_colormap, ...
    'Units', 'pixel', 'Position', [233 3 60 20], ... %'Tag', 'yScaleGlobal',...
    'BackgroundColor', get(gui, 'Color'), 'HorizontalAlignment', 'left', 'Value', strcmp(customConfig.colormapMode,'Thresh'));
popThresh = uicontrol('Parent', panel_imageControls, 'Style', 'popupmenu', 'Callback', @updateColormap, 'Units', 'Pixel', ...
    'Position', [250 137 60 24], 'String', {'triangle';'otsu';'isodata';'max entropy'},...
    'Value',1, 'TooltipString', 'Choose threshold method','FontSize',7, 'Tag','Choose threshold method.');
% Add additional thresholding algorithms if dipimages' threshold function
% is available
if withDipimage
    set(popThresh,'String',{'triangle';'otsu';'isodata';'min error';'background';'max entropy'});
end
cmStretchRGBMean = uicontrol('Style', 'Radio', 'String', 'MeanStretch', 'Parent', bg_colormap,  ...
    'Units', 'pixel', 'Position', [130 3 80 20], 'Tooltipstring', 'Stretch RGB map along dimension',...
    'BackgroundColor', get(gui, 'Color'), 'HorizontalAlignment', 'left', 'Value', 0,...
    'Visible', 'off');
cmStretchRGBMax = uicontrol('Style', 'Radio', 'String', 'MaxStretch', 'Parent', bg_colormap,  ...
    'Units', 'pixel', 'Position', [215 3 75 20], 'Tooltipstring', 'Stretch RGB map along dimension',...
    'BackgroundColor', get(gui, 'Color'), 'HorizontalAlignment', 'left', 'Value', 0,...
    'Visible', 'off');
% Popup for histogram / contrasts par selection
clear s
if nMat > 1 %|| withAlpha
    for i=1:nMat
        if isempty(varName{i})
            s{i} = ['Data' num2str(i)];
        else
            s{i} = varName{i};
        end
    end
    %     if withAlpha
    %         for i=1:nMat
    %             s{i+nMat} = [s{i} '-Alpha'];
    %         end
    %     end
    %Popup for selection of data set
    popContrastSel = uicontrol('Style','popup','FontSize', 7,'Parent',panel_imageControls,'String',s,'value',currContrastSel,'Position',[6 113 47 30],'Value',1,'Callback',@updateContrastSel);
    % Lock contrast adjustments across data sets (manual mode: linear scaling)
    tb_linkContrastAdjustments = uicontrol('Parent',panel_imageControls, 'Style', 'Togglebutton','Position',[54 125 15 15],'Value',linkContrastSettings,  ...
        'CData', lockOpen, 'Callback', @updateLinkContrast,...
        'BackgroundColor', get(gui, 'Color'),'Tooltipstring', 'Link contrast adjustments across data sets in manual mode. ',...
        'Tag', ['Link contrast adjustments across data sets in manual mode. (Linear scaling) ']);  %#ok
    if linkContrastSettings
        set(tb_linkContrastAdjustments, 'CData', lockClosed);
    end
end
%Values of slider limits
uicontrol('Parent', panel_imageControls, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [6 80 47 40], 'String', {'Slider';'min / max'}, ...
    'HorizontalAlignment', 'center','FontSize',7, 'Tag', 'Set minimum and maximum of contrast sliders.  The range of the histogram and the display of the colormap in between the sliders will be updated accordingly.');
sldLimMin = uicontrol('Parent', panel_imageControls, 'Style', 'Edit', 'Units', 'Pixel', ...
    'Position', [7 80 45 12], 'String', num2str(cmMinMax(1,1)),...
    'HorizontalAlignment', 'right','Callback', @updateSldLim, 'Tag', 'Set minimum of contrast sliders. The range of the histogram and the display of the colormap in between the sliders will be updated accordingly.');
sldLimMax = uicontrol('Parent', panel_imageControls, 'Style', 'Edit', 'Units', 'Pixel', ...
    'Position', [7 69 45 12], 'String', num2str(cmMinMax(1,2)), ...
    'HorizontalAlignment', 'right','Callback', @updateSldLim, 'Tag', 'Set maximum of contrast sliders. The range of the histogram and the display of the colormap in between the sliders will be updated accordingly.');
if ~isinteger(data{1})
    set(sldLimMin, 'String', num2str(cmMinMax(1  ,1),'%6.3f'));
    set(sldLimMax, 'String', num2str(cmMinMax(1  ,2),'%6.3f'));
end
%Slider Colormap
sldMin = uicontrol('Parent', panel_imageControls, 'Style', 'Slider', 'Callback', @updateColormap, 'Units', 'Pixel', ...
    'Position', [56 82 208 10], 'Min', minVal(1), 'Max', maxVal(1), 'Enable', 'off', ...
    'Value',cmMinMax(1,1),'Tooltipstring','Set minimum of colormap ', 'Tag', 'Set minimum (''black point'') of colormap. This value is linked to the edit box to the right.');
sldMax = uicontrol('Parent', panel_imageControls, 'Style', 'Slider', 'Callback', @updateColormap, 'Units', 'Pixel', ...
    'Position', [56 69 208 10], 'Min', minVal(1), 'Max', maxVal(1), 'Enable', 'off', ...
    'Value',cmMinMax(1,2),'Tooltipstring','Set maximum of colormap ', 'Tag', 'Set maximum (''white point'') of colormap. This value is linked to the edit box to the right.');
if isinteger(data)
    set(sldMax,'SliderStep', [1/double(maxVal-minVal-1) 10/double(maxVal-minVal-10)]);
    set(sldMin,'SliderStep', [1/double(maxVal-minVal-1) 10/double(maxVal-minVal-10)]);
else
    set(sldMax, 'SliderStep', [0.005 0.05]);
    set(sldMin, 'SliderStep', [0.005 0.05]);
end
%Values of colormap slider
uicontrol('Parent', panel_imageControls, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [267 80 47 40], 'String', {'Contrast';'min / max'}, ...
    'HorizontalAlignment', 'center','FontSize',7, 'Tag', 'Set minimum (''black point'') and maximum (''white point'') of colormap. These values are linked to the sliders to the left.');
valSldMin = uicontrol('Parent', panel_imageControls, 'Style', 'Edit', 'Units', 'Pixel', ...
    'Position', [268 80 45 12], 'String', num2str(cmMinMax(1,1)), 'Callback', @setColormapFromEdit,...
    'HorizontalAlignment', 'right', 'Tag', 'Set minimum (''black point'') of colormap. This value is linked to the slider to the left.');
valSldMax = uicontrol('Parent', panel_imageControls, 'Style', 'Edit', 'Units', 'Pixel', ...
    'Position', [268 69 45 12], 'String', num2str(cmMinMax(1,2)), 'Callback', @setColormapFromEdit,...
    'HorizontalAlignment', 'right', 'Tag', 'Set maximum (''white point'') of colormap. This value is linked to the slider to the left.');

% Small histogram on top of sliders
histAxGuiPos = [72 93 176 47]; %82
% histAxGui = axes('Parent', gui,'Units','pixel','Position', histAxGuiPos,'Tag','Histogram of current image. Middle click to change between lin / log y-scale.');
histAxGui = axes('Parent', panel_imageControls,'Units','pixel','Position', histAxGuiPos,'Tag','Histogram of current image. Middle click to change between lin / log y-scale. Light rectangle in the background indicates current range of colormap.');

if ~oldMATLAB
    histAxBg = get(gui, 'Color')*1.1;histAxBg(histAxBg>1)=1; % Avoid >1
    histAxBg = rectangle('Parent',histAxGui, 'Position',[0 0 176 44], ...
        'FaceColor',histAxBg, 'EdgeColor','none'); % 'FaceColor',get(gui, 'Color')*1.1, 'EdgeColor','none');
else
    histAxBg = rectangle('Parent',histAxGui, 'Position',[0 0 176 44], ...
    'FaceColor',get(gui, 'Color')*1.1, 'EdgeColor','none');
end
cb_showHist = uicontrol('Parent', panel_imageControls, 'Style', 'checkbox','Units','pixel',...
    'Position',[histAxGuiPos(1)+histAxGuiPos(3)-15 histAxGuiPos(2)+histAxGuiPos(4)-17 15 17],...
    'BackgroundColor',get(gui,'Color'),'Value', 0, 'String', '','Visible','off',...
    'Callback',@updateGuiHistVal,...
    'Tooltipstring', 'Show/hide histogram for RGB stretch mode.','Tag','Show/hide histogram for RGB stretch mode.');

hold(histAxGui,'on');
histXData = linspace(minVal(1),maxVal(1),256);
histXDataBin = histXData(1:end-1)+.5*diff(histXData(1:2));
histAxPlot = bar(histXDataBin, zeros(1,numel(histXDataBin)),'Parent',histAxGui, 'BaseValue',0,'BarWidth',1,'EdgeColor','k');
set(histAxGui,'XLim',[minVal(1) maxVal(1)],'Color','none','XTick',[],'YTick',[],'YColor',get(gui,'Color'));
set([histAxGui get(histAxGui,'Children')'], 'ButtonDownFcn', @buttonDownGuiHist);
txt_linlog = text(1,1, 'lin','Parent',histAxGui, 'FontSize',8, 'units','pixel','Position',[-13 15],'Rotation',90);
tb_playHist = uicontrol('Parent', panel_imageControls, 'Style', 'Togglebutton','Units','pixel',...
    'Position',[histAxGuiPos(1)+histAxGuiPos(3)+1 histAxGuiPos(2)-1 15 24],'CData',arrowAll(:,9:16,:),'BackgroundColor',get(gui,'Color'),...
    'Tooltipstring','Toggle update of histogram during playback','Tag','Toggle update of histogram during playback. As updating the histogram takes considerable time, this can be useful for speeding up the playback.','Value',defaultConfig.playHist);
tb_moveHist = uicontrol('Parent', panel_imageControls, 'Style', 'Togglebutton','Units','pixel','BackgroundColor',get(gui,'Color'),...
    'Position',[histAxGuiPos(1)+histAxGuiPos(3)+1 histAxGuiPos(2)+23 15 24],'CData',repmat(floor(arrow_lr(:,[3:7 10:14])/2),[1 1 3]),...
    'Tooltipstring','Toggle update of histogram when dragging the current position indicators.','Tag','Toggle update of histogram during change of current position (either using the sliders or the position lines in the plots). As updating the histogram takes considerable time, this can be useful for speeding up the image refresh.','Value',defaultConfig.moveHist,...
    'Callback','');  %@toggleMoveHist
% axis(histAxGui,'off');
set(txt_linlog, 'Visible','on','Position',get(txt_linlog,'Position')+[5 3 0]);
% Display image for contrast slider
contrastAx = axes('Parent', panel_imageControls,'Units','pixel','Position', [72 79 176 3]);
contrastSldIm = imagesc(histXDataBin,1,linspace(minVal(1),maxVal(1),numel(histXDataBin)) ,'Parent',contrastAx);  %#ok
set(contrastAx, 'Tag','Current colormap.','Visible','off');

%Slider RGB Stretch
sldMin_RGB = uicontrol('Parent', panel_imageControls, 'Style', 'Slider', 'Callback', @updateRgbStretchPar, 'Units', 'Pixel', ...
    'Position', [56 55 208 10], 'Min', -.5, 'Max', 1.5, 'Visible', 'off', ...
    'Value',rgbStretchSldVal(1),'SliderStep', [1/200, 1/20], 'Tag', 'Set lower end of RGB spectrum used for RGB stretch. Values <= zero correspond to blue, values >= 1 correspond to red. This value is linked to the edit box to the right.');
sldMax_RGB = uicontrol('Parent', panel_imageControls, 'Style', 'Slider', 'Callback', @updateRgbStretchPar, 'Units', 'Pixel', ...
    'Position', [56 42 208 10], 'Min', -.5, 'Max', 1.5, 'Visible', 'off', ...
    'Value',rgbStretchSldVal(2),'SliderStep', [1/200, 1/20], 'Tag', 'Set upper end of RGB spectrum used for RGB stretch. Values <= zero correspond to blue, values >= 1 correspond to red. This value is linked to the edit box to the right.');
%Values of colormap slider
valSldMin_RGB = uicontrol('Parent', panel_imageControls, 'Style', 'Edit', 'Units', 'Pixel', ...
    'Position', [270 53 40 12], 'String', num2str(rgbStretchSldVal(1)), 'Callback', @updateRgbStretchPar,...
    'HorizontalAlignment', 'right', 'Visible', 'off', 'Tag','Set lower end of RGB spectrum used for RGB stretch. Values <= zero correspond to blue, values >= 1 correspond to red. This value is linked to the slider to the left.');
valSldMax_RGB = uicontrol('Parent', panel_imageControls, 'Style', 'Edit', 'Units', 'Pixel', ...
    'Position', [270 42 40 12], 'String', num2str(rgbStretchSldVal(2)), 'Callback', @updateRgbStretchPar,...
    'HorizontalAlignment', 'right', 'Visible', 'off', 'Tag','Set upper end of RGB spectrum used for RGB stretch. Values <= zero correspond to blue, values >= 1 correspond to red. This value is linked to the slider to the left.');
% Display image for RGB Stretch slider
rgbTmp = zeros(101,3);
for i=1:101
    rgbTmp(i,:) = rgbValue(-.5+(i-1)/100*2,1);
end
rgbTmp = permute(repmat(rgbTmp,[1 1 30]),[3 1 2]);
rgbAx = axes('Parent', panel_imageControls,'Units','pixel','Position', [72 52 176 3],'Visible','off');
rgbSldIm = imagesc(rgbTmp,'Parent',rgbAx,'Visible','off');
axis(rgbAx,'off');

%Gamma Value Slider and etxt
sldGamma = uicontrol('Parent', panel_imageControls, 'Style', 'Slider', 'Callback', {@updateGamma,'sld','data'}, 'Units', 'Pixel', ...
    'Position', [56 55 208 10], 'Min', 0, 'Max', 5, 'SliderStep', [0.01 .05], ...
    'Value',customConfig.gamma(1), 'ToolTipString', 'Set Gamma Value','Tag', 'Gamma value to provide non-linear colormaps.');
valSldGamma = uicontrol('Parent', panel_imageControls, 'Style', 'Edit', 'Units', 'Pixel', ...
    'Position', [268 53 45 12], 'String', num2str(customConfig.gamma(1),'%6.3f'), 'Callback', {@updateGamma,'etxt','data'},...
    'HorizontalAlignment', 'right', 'ToolTipString', 'Set Gamma Value','Tag', 'Gamma value to provide non-linear colormaps. This value is linked to the slider to the left.');
strGamma = text('Units', 'Pixel', 'Parent',rgbAx,...
    'Position', [-32 8 1], 'String', '$\gamma$','Interpreter','Latex','Tag', 'Gamma value to provide non-linear colormaps. This value is linked to the edit box to the right.');

if withAlpha
    set(sldGamma, 'Position',[86 55 178 10]); % Make space for gamma contrast edit fields
    set(strGamma, 'Position',[-2 8 1]);
    text('Units', 'Pixel', 'Parent',rgbAx,...
    'Position', [-60 8 1], 'String', '$\alpha_{\mathrm{min}}$','Interpreter','Latex','Tag', 'Minimum of alpha contrast');
    edt_alphaMin = uicontrol('Parent', panel_imageControls, 'Style', 'Edit', 'Units', 'Pixel', ...
        'Position', [7 42 31 12], 'String', num2str(cmMinMax(nMat+1,1),'%6.2f'), 'Callback', {@updateAlphaContrast,'etxt'},...
        'HorizontalAlignment', 'right', 'ToolTipString', 'Set minimum of alpha contrast','Tag', 'Set minimum of alpha contrast');
    text('Units', 'Pixel', 'Parent',rgbAx,...
    'Position', [-30 8 1], 'String', '$\alpha_{\mathrm{max}}$','Interpreter','Latex','Tag', 'Maximum of alpha contrast');
    edt_alphaMax = uicontrol('Parent', panel_imageControls, 'Style', 'Edit', 'Units', 'Pixel', ...
        'Position', [37 42 31 12], 'String', num2str(cmMinMax(nMat+1,2),'%6.2f'), 'Callback', {@updateAlphaContrast,'etxt'},...
        'HorizontalAlignment', 'right', 'ToolTipString', 'Set maximum of alpha contrast','Tag', 'Set maximum of alpha contrast');
    % Gamma for alpha data: Gamma Value Slider and etxt
    sldGamma(2) = uicontrol('Parent', panel_imageControls, 'Style', 'Slider', 'Callback', {@updateGamma,'sld','alpha'}, 'Units', 'Pixel', ...
        'Position', [86 42 178 10], 'Min', 0, 'Max', 5, 'SliderStep', [0.01 .05], ...
        'Value',customConfig.gamma(min(nMat+1,numel(customConfig.gamma))), 'ToolTipString', 'Set Gamma Value','Tag', 'Gamma value to provide non-linear colormaps.');
    valSldGamma(2) = uicontrol('Parent', panel_imageControls, 'Style', 'Edit', 'Units', 'Pixel', ...
        'Position', [268 42 45 12], 'String', num2str(customConfig.gamma(min(nMat+1,numel(customConfig.gamma))),'%6.3f'), 'Callback', {@updateGamma,'etxt','alpha'},...
        'HorizontalAlignment', 'right', 'ToolTipString', 'Set Gamma Value','Tag', 'Gamma value to provide non-linear colormaps. This value is linked to the slider to the left.');
    strGamma(2) = text('Units', 'Pixel', 'Parent',rgbAx,...
        'Position', [-2 -5 1], 'String', '$\gamma_{\alpha}$','Interpreter','Latex','Tag', 'Gamma value to provide non-linear colormaps. This value is linked to the edit box to the right.');
end

%popup for colormap (look up table)
popLut = uicontrol('Parent', panel_imageControls, 'Style', 'popupmenu', 'Callback', @updateColormap, 'Units', 'Pixel', ...
    'Position', [7 14 80 15], 'String', {'Gray';'Gray (Range)';'4 x Gray';'Parula'; 'Jet'; 'HSV'; 'Hot'; 'Cool';...
    'Red 1';'Red 2';'Green 1';'Green 2';'Blue 1';'Blue 2'; ...
    'Rainbow1';'Rainbow2';'Rainbow3';'Rainbow4';...
    'Blue-Gray-Yellow (0 centered)';'Blue-Gray-Red (0 centered)';'Magenta-Gray-Green (0 centered)'},...
    'Value',customConfig.colormap, 'TooltipString', 'Choose colormap','FontSize',7, 'Tag','Select colormap.');

if ~isempty(defaultColormap{1})
    set(popLut, 'String', [get(popLut, 'String');'Default']);
end

btGap = 4;
btPos = 90;
btWidth = 24;

%Invert Colormap
tbInvert = uicontrol('Parent', panel_imageControls, 'Style', 'Togglebutton', 'Units', 'Pixel', 'Value', 0,  ...
    'Position', [btPos 8 24 24], 'Enable', 'on', 'CData', invertBt, 'Callback', @updateInvertMode,...
    'ToolTipString', 'Invert Colormap','Tag','Invert colormap');
btPos = btPos + btGap + btWidth;
%Flip Colormap
tbFlip = uicontrol('Parent', panel_imageControls, 'Style', 'Togglebutton', 'Units', 'Pixel', 'Value', 0,  ...
    'Position', [btPos 8 24 24], 'Enable', 'on', 'CData', flipBt, 'Callback', @updateFlipMode,...
    'ToolTipString', 'Flip colormap','Tag','Flip colormap');
btPos = btPos + btGap + btWidth;
%Toggle button for RGB display (disabled for 2D data)
tbSwitchRGB = uicontrol('Parent', panel_imageControls, 'Style', 'Togglebutton', 'Units', 'Pixel',  ...
    'Position', [btPos 8 24 24],  'Callback', @switchRGB, 'CData', RGB,...
    'TooltipString', 'Display three consecutive slices as RGB image','Tag', 'Display three consecutive slices as RGB image or use ''stretch modes'' to create color-coded projections along the RGB dimension.');
if nDim < 3 || withAlpha
    set(tbSwitchRGB, 'Enable', 'off');
end
btPos = btPos + btGap + btWidth;

%toggle button for display of colorbar
tbColorbar = uicontrol('Parent', panel_imageControls, 'Style', 'Togglebutton', 'Units', 'Pixel', 'CData', colorbarIcon, ...
    'Position', [btPos 8 24 24], 'String', '', 'Value', 0, 'Callback', @showColorbar, 'TooltipString', ' Show / hide colorbar ', 'Tag','Show / hide colorbar in the image window.');
btPos = btPos + btGap + btWidth;
%Toggle button for equal axes property
tbAspRatio = uicontrol('Parent', panel_imageControls, 'Style', 'Togglebutton', 'Units', 'Pixel',  ...
    'Position', [btPos 8 24 24], 'String', '1:1', 'Value', 1  , 'FontSize', 7, ...
    'Callback', @updateAspRatio,'Value', customConfig.aspectRatio, 'TooltipString', ['Switch axes between filling the figure or having a fixed aspect ratio.' char(10) '- Left click: Toggle fill / fixed.' char(10) '- Right click: Set aspect ratio.'],...
    'Tag', ['Switch axes between filling the figure or having a fixed aspect ratio.' char(10) '- Left click: Toggle fill / fixed.' char(10) '- Right click: Set aspect ratio.']);
btPos = btPos + btGap + btWidth;
%Toggle button for display of objects
tbShowObjects = uicontrol('Parent', panel_imageControls, 'Style', 'Togglebutton', 'Units', 'Pixel',  ...
    'Position', [btPos 8 24 24], 'CData', objBt, 'Callback', @toggleShowObjects, ...
    'TooltipString', 'Toggle display of objects in Image/Zoom windows', 'Value', customConfig.lineVis, ...
    'Tag', 'Toggle display of objects (position lines and zoom indicator) in Image/Zoom windows');
btPos = btPos + btGap + btWidth;
%100 % button
bt_100pct = uicontrol('Parent', panel_imageControls, 'Style', 'pushbutton', 'Units', 'Pixel',  ...
    'Position', [btPos 8 24 24], 'String', '100%','FontSize',6, 'Callback', @set100Pct, ...
    'TooltipString', 'Set display size to 100 % (1 image pixel equals 1 screen pixel). Right click to specify another pixel scale.', 'Value', 0, ...
    'Tag', 'Set display size of zoom window to 100 % (1 image pixel equals 1 screen pixel). Right click to specify another pixel scale.'); %#ok

%Separator line
sepLine(2) = uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [0 133 customConfig.winPos.gui(3) 10], 'String','__________________________________________________________________________________', ...
    'HorizontalAlignment', 'left', 'FontWeight', 'Bold','Fontsize',6);

% Windows, config & help
txt_titles(end+1) = uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [5 105 customConfig.winPos.gui(3)-10 20], 'String', 'Windows, configuration & help', ...
    'HorizontalAlignment', 'center','FontWeight','bold','FontSize',10);

%Toggle Button for Window Visibility & Window Controls
panel_windowButtons = uipanel(gui, 'units','pixel','Position', [-1 75 customConfig.winPos.gui(4)+2 30],...
    'BackgroundColor',get(gui, 'Color'),'BorderType','none'); %
%For image windows
tbWin(1) = uicontrol('Parent',panel_windowButtons, 'Style', 'Togglebutton', 'Units', 'Pixel', 'TooltipString', 'Show / hide Image window(s)', ...
    'Position', [7 1 28 28], 'CData', double(icon_image_24x24)/255, 'Value', customConfig.winVis.imageWin, 'Callback',{@windowVisibility,1},...
    'Tag', 'Show / hide image window(s)','BackgroundColor',get(gui,'Color'));
%For zoom windows
tbWin(2) = uicontrol('Parent',panel_windowButtons, 'Style', 'Togglebutton', 'Units', 'Pixel', 'TooltipString', 'Show / hide Zoom window(s)',...
    'Position', [36 1 28 28], 'CData', double(icon_zoom_24x24)/255, 'Value', customConfig.winVis.zoomWin, 'Callback',{@windowVisibility,2},...
    'Tag', 'Show / hide zoom window(s)','BackgroundColor',get(gui,'Color'));
%For plot windows
tbWin(3) = uicontrol('Parent',panel_windowButtons, 'Style', 'Togglebutton', 'Units', 'Pixel', 'TooltipString', 'Show / hide Plots window(s)',...
    'Position', [65 1 28 28], 'CData', double(icon_plot_24x24)/255, 'Value', customConfig.winVis.plotWin, 'Callback',{@windowVisibility,3},...
    'Tag', 'Show / hide plot window(s)','BackgroundColor',get(gui,'Color'));
%Histogram
tbHist = uicontrol('Parent',panel_windowButtons, 'Style', 'Togglebutton', 'Units', 'Pixel', 'Value', 0,  ...
    'Position', [94 1 28 28], 'Enable', 'on', 'CData', histIcon, 'Callback',@showHist ,...
    'ToolTipString', 'Show histogram of complete data set. Histogram not yet calculated. Press button to calculate and display.', 'FontSize',7,...
    'Tag', 'Show / hide histogram window','BackgroundColor',get(gui,'Color'));  %
if withAlpha
    set(tbHist, 'Callback', @update2DHist);
end
% ROI Manager
tbRoi = uicontrol('Parent',panel_windowButtons, 'Style', 'Togglebutton', 'Units', 'Pixel', 'Callback', @roiGui,  ...
    'Position', [123 1 28 14],  'Value', 0, 'Enable','on','CData',double(icon_roi_24x24(:,13:24,:))/255,...
    'ToolTipString', 'Show / hide ROI Manager',...
    'Tag', 'Show / hide ROI Manager','BackgroundColor',get(gui,'Color'));  %[175 85 24 24]
tbProfile = uicontrol('Parent',panel_windowButtons, 'Style', 'Togglebutton', 'Units', 'Pixel', 'Callback', @profileGui,  ...
    'Position', [123 15 28 14],  'Value', 0, 'Enable','on','CData',double(icon_roi_24x24(:,1:12,:))/255,...
    'ToolTipString', 'Show / hide profile Manager',...
    'Tag', 'Show / hide profile Manager','BackgroundColor',get(gui,'Color'));  %[175 85 24 24]
% Record button
tbRecord = uicontrol('Parent',panel_windowButtons, 'Style', 'Togglebutton', 'Units', 'Pixel', 'Callback', @vidGenerator,  ...
    'Position', [152 1 28 28],  'Value', 0, 'Enable','on','CData',double(icon_Record),...
    'ToolTipString', 'Show /hide video recording tool','String','REC','FontSize',8,'Foregroundcolor','w','FontWeight','bold',...
    'Tag', 'Show / hide video recording tool','BackgroundColor',0.2*[1 1 1],'tag','Open video generator tool to record videos');  %[175 85 24 24]
% tooltips button
tbTooltips = uicontrol('Parent',panel_windowButtons, 'Style', 'Togglebutton', 'Units', 'Pixel',...
    'Position', [181 1 28 28],  'CData', double(icon_tooltips_24x24)/255,...
    'Tooltipstring', 'Show / hide tooltips','ForegroundColor',get(gui,'Color'),'Value',1, 'Callback', @toggleTooltipDisplay,...
    'Tag', 'Show / hide tooltips','BackgroundColor',get(gui,'Color'));
% matVisGuide button
bt_matVisGuide = uicontrol('Parent',panel_windowButtons, 'Style', 'Pushbutton', 'Units', 'Pixel', 'Callback', @openMatVisGuide,...
    'Position', [210 1 28 28],  'CData', icon_manual(5:28,5:28,:),...
    'Tooltipstring', 'Open matVis Guide','ForegroundColor',get(gui,'Color'),...
    'Tag','Open matVis Guide. Opens a pdf-file in a webbrowser. Only available with internet connection.');  %#ok
% update matVis
bt_updateMatVis = uicontrol('Parent',panel_windowButtons, 'Style', 'Pushbutton', 'Units', 'Pixel', 'Callback', @updateMatVis,...
    'Position', [239 1 28 28],  'CData', double(icon_update_24x24)/255,...
    'Tooltipstring', 'Check online for matVis update.','ForegroundColor',get(gui,'Color'),...
    'Tag','Check online whether a new version of matVis is available. Older versions can either be overwritten or backed up with the respective version number.'); %#ok
% % Check for availability of website hosting the latest matVis version and
% % the matVisGuide
% [w flag] = urlread('http://www.colors-and-contrasts.com/');
% if ~flag
%     set([bt_matVisGuide bt_updateMatVis], 'Enable','off');
% end
%Toggle default and custom window position / visibility
btConfig = uicontrol('Parent',panel_windowButtons, 'Style', 'Pushbutton', 'Units', 'Pixel', 'TooltipString', 'Switch configuration: Custom - Default',...
    'Position', [272 15 45 15], 'String', 'CustConf','Callback',@switchConfig, 'FontSize',7,...
    'UserData', 0,'Tag','Switches matVis configuration between default and custom configuration. The matVis configuration includes position and visibility of windows, colormap, contrast mode, aspect ratio and many other parameters.');
% btSaveConfig
uicontrol('Parent',panel_windowButtons, 'Style', 'Pushbutton', 'Units', 'Pixel', 'TooltipString', 'Save Configuration', 'FontSize',7,...
    'Position', [272 1 45 15], 'String', 'SaveConf', 'Callback',@saveConfig,'Tag','Saves the current matVis configuration as the custom configuration. The config file (''matVisConfig.mat'') is saved in the same directory as the matVis.m file.');

%Separator line
sepLine(3) = uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [0 65 customConfig.winPos.gui(3) 10], 'String','__________________________________________________________________________________', ...
    'HorizontalAlignment', 'left', 'FontWeight', 'Bold','Fontsize',6);

%Status line
txt_tooltips = uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [5 0 customConfig.winPos.gui(3)-10 60], 'String','Enable to show tooltips.', ...
    'HorizontalAlignment', 'left', 'Fontsize',8);
%% Image, Zoom and Plot Windows
%Default window positions
%Image, zoom and plot windows
for i = 1:nMat
    imageWin(i) = figure('Units', 'Pixel', 'Position', customConfig.winPos.imageWin(i,:),'Name', ['Image (',varName{i},')'], ...
        'MenuBar', 'none', 'NumberTitle', 'off', 'WindowButtonMotionFcn', @mouseMotion,...
        'WindowButtonDownFcn', @buttonDownCallback, 'CloseRequestFcn', {@windowVisibility, 0},...
        'Visible', 'off', 'HandleVisibility', 'off', 'WindowStyle','normal', 'PaperPositionMode','auto');   %#ok
    try
        im_jf = get(imageWin(i),'JavaFrame');
        im_jf.setFigureIcon(javax.swing.ImageIcon(im2java(uint8(icon_image))));
    catch    %#ok
    end
    %winVis{imageWin(i)} =   'on';
    zoomWin(i) = figure('Units', 'Pixel', 'Position', customConfig.winPos.zoomWin(i,:),'Name', ['Zoom (', varName{i},')'], ...
        'MenuBar', 'none', 'NumberTitle', 'off', 'WindowButtonMotionFcn', @mouseMotion,...
        'WindowButtonDownFcn', @buttonDownCallback, 'CloseRequestFcn', {@windowVisibility, 0},...
        'Visible', 'off', 'HandleVisibility', 'off','KeyPressFcn',@zoomKeyPressFcn,...
        'KeyReleaseFcn',@zoomKeyReleaseFcn, 'WindowStyle','normal', 'PaperPositionMode','auto');    %#ok
    try
        zoom_jf = get(zoomWin(i),'JavaFrame');
        zoom_jf.setFigureIcon(javax.swing.ImageIcon(im2java(uint8(icon_zoom))));
    catch    %#ok
    end
    %winVis{zoomWin(i)} =   'on';
    plotWin(i) = figure('Units', 'Pixel', 'Position', customConfig.winPos.plotWin(i,:),'Name', ['Plots (', varName{i}, ')'], ...
        'MenuBar', 'none', 'NumberTitle', 'off', 'WindowButtonMotionFcn', @mouseMotion,...
        'WindowButtonDownFcn', @buttonDownCallback, 'CloseRequestFcn', {@windowVisibility, 0},...
        'Visible', 'off', 'HandleVisibility', 'off', 'WindowStyle','normal', 'PaperPositionMode','auto');    %#ok
    try
        plot_jf = get(plotWin(i),'JavaFrame');
        xx = javax.swing.ImageIcon(im2java(uint8(icon_plot)));
        plot_jf.setFigureIcon(xx);
    catch    %#ok
    end
    %winVis{plotWin(i)} =   'on';
end
try
    % Avoid that selfmade icons appear in figure container (dirty solution ...)                                                                                          -
    f = figure;
    fake_jf = get(f,'JavaFrame');
    xxx = fake_jf.getFigureIcon;
    fake_jf.setFigureIcon(xx);
    fake_jf.setFigureIcon(xxx);
    delete(f);
catch    %#ok
end
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
%Set Scrollwheel Callback if possible (starting from Matlab 2007a)
try
    set([imageWin,zoomWin], 'WindowScrollWheelFcn', @scrollWheelCallback);
catch    %#ok
end
try
    set(gui, 'WindowScrollWheelFcn', @scrollWheelCallbackGui);
catch    %#ok
end
%Set background to black if alpha map is used
if withAlpha
    set(imageWin, 'Color', 'k');
    set(zoomWin, 'Color', 'k');
end
%% Histogram Window
histWin = figure('Units', 'Pixel', 'Name', ['Histogram (',varName{1},') - Click: lin / log'], ...
    'MenuBar', 'none', 'NumberTitle', 'off', 'Visible', 'off', 'HandleVisibility', 'off','WindowStyle','normal', ...
    'WindowButtonDownFcn', {@updateHist,1}, 'CloseRequestFcn', {@showHist,1}, 'UserData', 0);
histAx(1) = axes('Parent',histWin, 'Box', 'on');
histAx(2) = axes('Parent',histWin ,'Position',get(histAx,'Position'),'Visible','off'); %Axis for contrast lines  ,'Position',get(histAx,'Position')
hold(histAx(1), 'on');
hold(histAx(2), 'on');
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
hist_jf = get(histWin,'JavaFrame');
histIcon(isnan(histIcon)) = 1;
hist_jf.setFigureIcon(javax.swing.ImageIcon(im2java(uint8(255*histIcon))));
warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
%% Start
s = get(tbWin, 'Value');
set(tbWin, 'Value', 1); %Draw all windows to initialize all objects and axes
drawImages;
if updateProj   %Update value of projection dimenions popmenu if projetion dimension is given as start par
    s = get(projDimPop, 'String');
    set(projDimPop, 'Value',find(strcmp(s, dimNames{projDim})));
    if projMethod == 6
        projCallback;
    end
end
if updateRGB   %Update value of RGB dimenions popmenu if RGB dimension is given as start par
    s = get(projDimPop, 'String');
    rgbCount = mod(rgbDim , nDim-1);
    switchRGB;
end
updateColormap;
drawObjects;
drawPlots;
if ~get(tb_moveHist, 'Value')
    updateGuiHistVal;
end
adjustGuiSize(gui);
setConfig(customConfig);
% Set window resize callback
set(imageWin, 'ResizeFcn',{@resizeWin, 'image'});
set(zoomWin, 'ResizeFcn',{@resizeWin, 'zoom'});
set(plotWin, 'ResizeFcn',{@resizeWin, 'plot'});
warning off MATLAB:gui:latexsup:UnsupportedFont  % Supress warning due to supposedly missing latex font cmss10: bug since Matlab 2011b
figure(gui);
% warning on MATLAB:gui:latexsup:UnsupportedFont
currWin = gui;
set(gui, 'WindowButtonMotionFcn', @mouseMotion,'Visible','on'); % Setting the MotionFcn earlier would cause errors due to non-existing image and zoom windows
clear varargin;
if nargout == 1
    out = [];
    updateOutputStrct;
    varargout{1} = out;
end
% timerHist = timer('TimerFcn',@updateGlobalHist,'StartDelay',.1);
% start(timerHist);
%% Callback Functions
%Callback for slider for current position
    function sliderCallback(hObject, event)    %#ok
        if debugMatVis, debugMatVisFcn(1); end
        dimNum = get(hObject, 'Userdata');
        currPos(dimNum) = round(get(sld(dimNum), 'Value'));
        updateSelection(dimNum);
        if debugMatVis, debugMatVisFcn(2); end
    end

%Callback for textboxes for current position
    function etxtCallback(hObject, event, dimNum)   %#ok
        if debugMatVis, debugMatVisFcn(1); end
        currPos(dimNum) = str2double(get(etxt(dimNum), 'String'));
        currPos(dimNum) = min(max(currPos(dimNum),1),dim(dimNum));
        updateSelection(dimNum);
        if debugMatVis, debugMatVisFcn(2); end
    end

% Set position from small position buttons to the right
    function jumpPos(hO,ev,in)  %#ok
        if debugMatVis, debugMatVisFcn(1); end
        switch ischar(in(1))
            case 1
                switch in
                    case 'min'
                        currPos = minPos;
                    case 'max'
                        currPos = maxPos;
                end
            case 0
                currPos = savedPos(:,in)';
                if get(bt_savedZoom, 'Value')
                    zoomVal = savedZoom(:,:,in);
                    zoomValXY([1,3]) = zoomVal(xySel(2),:);
                    zoomValXY([2,4]) = zoomVal(xySel(1),:);
                    updateZoom;
                end
        end
        updateSelection(1:numel(dim));
        if debugMatVis, debugMatVisFcn(2); end
    end
    function findExtremum(hO,ev,mnmx, zmim)
        % Find and go to local maximum or minimum (either within currently displayed
        % image or currently selected zoom region inside current image,
        % depending on user input)
        if withAlpha
            selectedRegion = currAlphaMap{1};
        else
            selectedRegion = currIm{1};
        end
        switch zmim
            case 'zoom'
                selectedRegion = selectedRegion(zoomValXY(2):zoomValXY(2)+zoomValXY(4)-1,zoomValXY(1):zoomValXY(1)+zoomValXY(3)-1);
                regionOffset = zoomValXY([2 1])-1;
            case 'image'
                regionOffset = [0 0];
        end
        selectedSize = size(selectedRegion);
        switch mnmx
            case 'min'
                [extrVal, extrPos] = min(selectedRegion(:));
            case 'max'
                [extrVal, extrPos] = max(selectedRegion(:));
        end
        [extrPosX, extrPosY] = ind2sub(selectedSize,extrPos);
        currPos(xySel(1)) =  extrPosX + regionOffset(1);
        currPos(xySel(2)) =  extrPosY + regionOffset(2);
        updateSelection(1:nDim);
    end
    function setJumpPos(hO, ev, num)  %#ok
        savedPos(:,num) = currPos;
        savedZoom(:,:,num) = zoomVal;
        set(bt_getJumpPos(num), 'Enable','on', 'ToolTipString', ['Click to go to saved position: [' num2str(currPos) ']'],...
            'Tag',['Click to go to saved position: [' num2str(currPos) ']']);
    end

%Callback for keyboard: number keys change selected position of
%respective dimension, without modifier advance, with control go back.
%With shift, go to position according to saved positions
    function keyPress(src,evnt)   %#ok
        if debugMatVis, debugMatVisFcn(1); end
        if length(evnt.Modifier) == 1
            if strcmp(evnt.Modifier{:}, 'alt')  % Jump to saved position
                if size(evnt.Key,2) > 6  % Numpad Keys ('numpad1', 'numpad2',...)
                    evnt.Key
                    posNum = str2double(evnt.Key(7));
                else
                    posNum = str2double(evnt.Key);
                end
                if size(savedPos,2) >= posNum
                    jumpPos(1,1,posNum);
                end
            elseif strcmp(evnt.Modifier{:}, 'control')   % Go one step back
                if size(evnt.Key,2) > 6
                    dimNum = str2double(evnt.Key(7));
                else
                    dimNum = str2double(evnt.Key);
                end
                if sum(dimNum == 1:nDim)> 0
                    currPos(dimNum) = max([1  ,currPos(dimNum)-1]);
                    updateSelection(dimNum);
                end
            end
        elseif sum(str2double(evnt.Character) == 1:nDim) > 0  % Advance current position
            dimNum = str2double(evnt.Character);
            currPos(dimNum) = min([currPos(dimNum)+1  ,dim(dimNum)]);
            updateSelection(dimNum);
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Callback for change of xy dimension popups
    function popCallback(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        xNew = get(pop(1), 'Value');
        yNew = get(pop(2), 'Value');
        if xNew == xySel(1) && yNew == xySel(2)
          if debugMatVis, debugMatVisFcn(2); end
          return
        end
        if xNew ~= yNew
            %Selected Rois will be lost after changing xySel
            if ~isempty(roiList)
                q = questdlg('All Roi information will be lost after changing x/y dimension. Continue anyway?','','Continue','Cancel','Cancel');
                if strcmp(q,'Cancel')
                  if debugMatVis, debugMatVisFcn(2); end
                  return
                end
            end
            %Update string of projection popup and projection display
            n = dimNames;
            n([xNew yNew]) = [];
            set(projDimPop, 'String',  n, 'Value',1);
            %Exchange projDim with xyDim
            if projDim == xNew
                projDim = xySel(1);
            else
                projDim = xySel(2);
            end
            % Update zoom values
            zoomValXY([1,3]) = zoomVal(yNew,:);
            zoomValXY([2,4]) = zoomVal(xNew,:);
            xySel = [xNew, yNew];
            set(txt_dimName, 'ForegroundColor', 'k');
            if ~isempty(roiWin)
                resetRois;
            end
            if ~isempty(roiPopDataExport)
                n = dimNames;
                n(xySel) = [];
                set(roiPopDataExport, 'String', n);
            end
        else
            %Exchange x and y
            for ii=1:nMat
                currIm{ii} = permute(currIm{ii},[2 1 3]);
            end
            set(pop(1), 'Value', xySel(2));
            set(pop(2), 'Value', xySel(1));
            xySel(1) = get(pop(1), 'Value');
            xySel(2) = get(pop(2), 'Value');
            zoomValXY = [zoomValXY(2) zoomValXY(1) zoomValXY(4) zoomValXY(3)];
            if customDimScale
                set(roiImage,'XData',dimScale(xySel(2),:),'YData',dimScale(xySel(1),:),'CData',currIm{1},'AlphaData',1);
            else
                set(roiImage,'CData',currIm{1},'AlphaData',1);
            end
            if nRois > 0
                for ii=1:nRois
                    x =  roiList(ii).index.x;
                    roiList(ii).index.x = roiList(ii).index.y;
                    roiList(ii).index.y = x;
                    roiList(ii).rectangle = roiList(ii).rectangle([2,1  ,4,3]);
                    roiList(ii).corners = circshift(roiList(ii).corners,[1 0]);
                    roiList(ii).mask = roiList(ii).mask';
                end
                drawRois;
            end
        end
        set(txt_dimName(xySel(1)), 'ForegroundColor', 'r');
        set(txt_dimName(xySel(2)), 'ForegroundColor', [0 0.8 0]);
        if rgbCount
            rgbCount = nDim - 2;
        end
        if get(tbSwitchRGB, 'Value')
            switchRGB;
        end
        if projMethod
            projCallback;
        else
            updateImages;
            drawPlots;
        end
        if customDimScale
            set([imHandle zoomHandle], 'XData', dimScale(xySel(2),:), 'YData', dimScale(xySel(1),:));
            set([imAx zoomAx], 'XLim', dimScale(xySel(2),:), 'YLim', dimScale(xySel(1),:));
        else
            set(imAx, 'XLim', .5+[0 size(currIm{1},2)], 'YLim', .5+[0 size(currIm{1},1)]);
        end
        for ii=1:nMat
            if ~withDimUnits
                xlabel(imAx(ii), dimNames{xySel(2)});
                ylabel(imAx(ii), dimNames{xySel(1)});
            else
                xlabel(imAx(ii),[dimNames{xySel(2)} ' [' dimUnits{xySel(2)} ']']);
                ylabel(imAx(ii),[dimNames{xySel(1)} ' [' dimUnits{xySel(1)} ']']);
            end
        end
        updateZoom;
        updateObjects;
        if debugMatVis, debugMatVisFcn(2); end
    end

%Callback for checkboxes for dimensions to be plotted
    function cbCallback(hObject, event, dimNum)  %#ok
        if debugMatVis, debugMatVisFcn(1); end
        if nargin == 0
            for ii = 1:nDim
                plotSel(ii) = get(cbPlots(ii), 'Value');
            end
        else
            plotSel(dimNum) = get(cbPlots(dimNum), 'Value');
        end
        drawPlots;
        if debugMatVis, debugMatVisFcn(2); end
    end
    function projCallback(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        %Update projection method
        prevProjMethod = projMethod;
        projMethod = get(projMethodPop, 'Value') - 1;
        if rgbCount == get(projDimPop,'Value')
          projMethod =0;
          set(projMethodPop, 'Enable','off')
          %set(projMethodPop, 'Value', 1)
        else
          set(projMethodPop, 'Enable','on')
        end
        set([sld etxt btPlayAll btPlayZoom cbPlots tb_lockPlots2Zoom], 'Enable', 'on');
        if projMethod
            if isPlaying
                isPlaying = 0;
                set(btPlayAll, 'CData', arrow, 'Value', 0);
            end
            %set(tbSwitchRGB, 'Enable', 'off');
            set(btMean, 'Enable', 'off', 'CData', squeeze(meanBt{1}));
            %Update projection dimension
            notXY = setdiff(1:nDim, xySel);
            projDim = notXY(get(projDimPop, 'Value'));
            set([sld(projDim) etxt(projDim) btPlayAll(projDim) btPlayZoom(projDim) tb_lockPlots2Zoom(projDim)], 'Enable', 'off');
            %Allow no plots other than xy dimension (plot dimensions)
            set(cbPlots(notXY),  'Value', 0 ,'Enable', 'off'); %necessary when call from WindowCloseRequestFcn
            plotSel(notXY) = 0;
        else
            %set(tbSwitchRGB, 'Enable', 'on');
            set(btMean, 'Enable', 'on', 'CData', squeeze(meanBt{get(btMean, 'UserData')+1}));
            set([tbWin(3) cbPlots tbShowObjects], 'Enable', 'on');
        end
        updateImages;  %drawImages;
        %Disable position lines and plots in tile mode
        if projMethod == 6
            set([lineHorIm lineVertIm lineHorZoom lineVertZoom], 'Visible', 'off');
            set(plotWin, 'Visible', 'off');
            set([tbWin(3) cbPlots tbShowObjects], 'Enable', 'off');
            set(cbPlots, 'Value', 0);
            set(zoomReg, 'Visible', 'on');
            set([sld(xySel) etxt(xySel) btPlayAll(xySel) btPlayZoom(xySel)],'Enable','off');
            set([sld_up(xySel(1)) sld_down(xySel(1))], 'Max', size(currIm{1},1));
            set([sld_up(xySel(2)) sld_down(xySel(2))], 'Max', size(currIm{1},2));
            if customDimScale
                set(imAx, 'XLim', dimScale(1,:), 'YLim',dimScale(1,:));
            else
                set(imAx, 'XLim', .5+[0 size(currIm{1},2)], 'YLim', .5+[0 size(currIm{1},1)]);
            end
        else
            if prevProjMethod == 6
                set([sld(xySel) etxt(xySel) btPlayAll(xySel) btPlayZoom(xySel)],'Enable','on');
                set([sld_down(xySel(1)) sld_down(xySel(2))], 'Value', 1);
                set(sld_up(xySel(1)), 'Value',dim(xySel(1)));
                set(sld_up(xySel(2)), 'Value',dim(xySel(2)));
                set([sld_up(xySel(1)) sld_down(xySel(1))], 'Max', dim(xySel(1)));
                set([sld_up(xySel(2)) sld_down(xySel(2))], 'Max', dim(xySel(2)));
                set([tbWin(3) tbShowObjects], 'Enable', 'on');
                updateZoom(1,1,xySel(1),'sld');
                updateZoom(1,1,xySel(2),'sld');
                if customDimScale
                    set(imAx, 'XLim', dimScale(1,:), 'YLim',dimScale(1,:));
                else
                    set(imAx, 'XLim', .5+[0 size(currIm{1},2)], 'YLim', .5+[0 size(currIm{1},1)]);
                end
            end
            updateObjects;
            drawPlots;
        end
        if any(get(bg_colormap, 'SelectedObject') == [cmImage cmZoom])
            updateColormap
        end
        if ~rgbCount
          updateGuiHistVal;
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Callback for play buttons
    function playCallback(hObject, event, pD, allZoom, increase)  %#ok
        if debugMatVis, debugMatVisFcn(1); end
        playDim = pD;
        set([btPlayAll btPlayZoom], 'Enable', 'off');
        switch allZoom
            case 'zoom'
                startZoom = zoomVal(playDim,1);
                stopZoom = zoomVal(playDim,1) + zoomVal(playDim,2) - 1;
                set(btPlayZoom(playDim), 'Enable', 'on', 'CData', pausebt, 'Callback', @pauseCallback);
                if currPos(playDim) < startZoom
                    currPos(playDim) = startZoom;
                end
            case 'all'
                startZoom = 1;
                stopZoom = dim(playDim);
                set(btPlayAll(playDim), 'Enable', 'on', 'CData', pausebt, 'Callback', @pauseCallback);
        end
        isPlaying = 1;
        if movdata.rec
            if strcmp(allZoom, 'zoom')
                currPos(playDim) = zoomVal(playDim,1);
                endPos = zoomVal(playDim,1)+zoomVal(playDim,2)-1;
            else
                currPos(playDim) = 1;
                endPos = dim(playDim);
            end
            updateSelection(playDim);
        end
        tic;
        playTime = nan(1,10);
        playCt = 0;
        while isPlaying == 1 
            v = get(sldPlaySpeed, 'Value');
            if round(v) < 0 % faster playback: frame skipping
                inc = ((-1*round(v))+1)*(-1)^(increase+1);
                ps = 0;
            elseif round(v) > 0 % slower playback: pause
                inc = (-1)^(increase+1);
                ps = (v/10)^2;
            else
                inc = (-1)^(increase+1);
                ps = 0;
            end
            playCt = playCt + 1;
            playTime(mod(playCt-1,10)+1) = toc;
            playTimeDiff = diff(sort(playTime));
            if any(~isnan(playTimeDiff))
                if round(v) < 0  % frame skipping
                    set(gui,'Name',['fps: ' num2str(1/mean(playTimeDiff,'omitnan'),'%6.2f') ' / ' num2str(inc/mean(playTimeDiff,'omitnan'),'%6.2f') ' (every ' num2str(inc) '. frame)']);  % Used to be nanmean
                else
                    set(gui,'Name',['fps: ' num2str(1/mean(playTimeDiff,'omitnan'),'%6.2f') ' (' num2str(ps) ' s pause)']);  % Used to be nanmesn
                end
            end
            currPos(playDim) = currPos(playDim)+inc;
            if currPos(playDim) > stopZoom
                currPos(playDim) = startZoom;
            end
            if currPos(playDim) < startZoom
                currPos(playDim) = stopZoom;
            end
            updateSelection(playDim);
            if ~any(xySel == playDim) && any([cmImage cmZoom] == get(bg_colormap, 'SelectedObject'))
                updateColormap;
            end
            drawnow;
            pause(ps);
            if movdata.rec && currPos(playDim) == endPos
                set(hObject,'Value',0)
                pauseCallback;
            end
        end
        if debugMatVis, debugMatVisFcn(1); end
    end

%Callback for pause button
    function pauseCallback(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        isPlaying = 0;
        set(gui, 'Name', ['matVis: ', allNames]);
        for ii = 1:nDim
            set(btPlayAll(ii), 'Enable', 'on', 'Callback', {@playCallback,ii,'all',1}, 'CData', arrowAll);
            set(btPlayZoom(ii), 'Enable', 'on', 'Callback', {@playCallback,ii,'zoom',1}, 'CData', arrowZoom);
        end
        if movdata.rec
            set(btPlayAll,'CData',icon_RecordAll);
            set(btPlayZoom,'CData',icon_RecordZoom);
        end
        % Update GUI hist
        if ~get(tb_playHist,'Value')
            updateGuiHistVal([],[],'forceUpdateGuiHist')
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Callbakck for colorbar-display toggle button
    function showColorbar(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if get(tbColorbar, 'Value') == 1
            for ii=1:nMat
                cb_axes(ii) = colorbar('peer', imAx(ii));      %#ok
            end
            if withAlpha
                set(cb_axes, 'XColor', 'w', 'YColor', 'w');
            end
            if any(currGamma ~=1)
                updateImages;
            end
        else
            for ii=1:nMat
                colorbar('peer', imAx(ii),'off');
            end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Callback for play speed slider
    function playSpeedCallback(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        set(sldPlaySpeed, 'Value', round(get(sldPlaySpeed, 'Value')));
        set(sldPlaySpeed,'TooltipString',['Playspeed: ',num2str(get(sldPlaySpeed,'Value'))],...
            'Tag', ['Speed of playback. Current value: ' num2str(get(sldPlaySpeed,'Value')) '.' char(10) 'Negative values: faster speed by skipping over images.' char(10) 'Positive values: slower speed by introducing pauses.'] );
        updateTooltips;
        if debugMatVis, debugMatVisFcn(2); end
    end

%Callback for menu-display toggle button
    function toggleMenuBars(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if get(tbMenuBars, 'Value') && strcmp(get(imageWin(1), 'MenuBar'), 'none')
            for ii=1:nMat
                set(imageWin(ii), 'Position', get(imageWin(ii), 'Position') + [0 0 0 -winWidthMenuBar], 'MenuBar', 'figure', 'WindowButtonMotionFcn', '');
                set(zoomWin(ii), 'Position', get(zoomWin(ii), 'Position') + [0 0 0 -winWidthMenuBar], 'MenuBar', 'figure', 'WindowButtonMotionFcn', '');
                set(plotWin(ii), 'Position', get(plotWin(ii), 'Position') + [0 0 0 -winWidthMenuBar], 'MenuBar', 'figure', 'WindowButtonMotionFcn', '');
            end
            for hh = [histWin profileTraceWin]
              if ~isempty(hh)
                set(hh, 'Position', get(hh, 'Position') + [0 0 0 -winWidthMenuBar], 'MenuBar', 'figure');
              end
            end
        elseif ~get(tbMenuBars, 'Value') && strcmp(get(imageWin(1), 'MenuBar'), 'figure')
            for ii=1:nMat
                set(imageWin(ii), 'Position', get(imageWin(ii), 'Position') + [0 0 0 winWidthMenuBar], 'MenuBar', 'none', 'WindowButtonMotionFcn', @mouseMotion);
                set(zoomWin(ii), 'Position', get(zoomWin(ii), 'Position') + [0 0 0 winWidthMenuBar], 'MenuBar', 'none', 'WindowButtonMotionFcn', @mouseMotion);
                set(plotWin(ii), 'Position', get(plotWin(ii), 'Position') + [0 0 0 winWidthMenuBar], 'MenuBar', 'none', 'WindowButtonMotionFcn', @mouseMotion);
            end
            for hh = [histWin profileTraceWin]
              if ~isempty(hh)
                set(hh, 'Position', get(hh, 'Position') + [0 0 0 winWidthMenuBar], 'MenuBar', 'none');
              end
            end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Callback for RGB toggle button
    function switchRGB(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        rgbCount = mod(rgbCount + 1  , nDim - 1);
        %RGB mode off
        if rgbCount == 0
            set(bg_colormap, 'SelectedObject', cmManual);
            set(tbSwitchRGB, 'String', '', 'CData', RGB, 'Value', 0);
            set(cmImage, 'String', 'Image','ToolTipString', '');
            set(cmZoom, 'Visible','on');
            set(cmManual, 'Visible', 'on');
            set([cmThresh popThresh], 'Visible', 'on');
            set(popLut, 'Enable', 'on');
            set(tbColorbar, 'Enable', 'on');
            set(bg_colormap,'SelectionChangeFcn', @updateColormap, 'Tag',get(bg_colormap,'UserData'));
            set([cmStretchRGBMean cmStretchRGBMax], 'Visible', 'off');
            set(projMethodPop ,'Enable', 'on');
            set([sldGamma valSldGamma strGamma], 'Visible', 'on');
            set(btMean, 'UserData', mod(get(btMean, 'UserData')  ,5));
            set(btMean, 'CData', squeeze(meanBt{get(btMean, 'UserData')+1}));
            set([sldMax_RGB sldMin_RGB valSldMax_RGB valSldMin_RGB rgbSldIm], 'Visible','off');
            set([valSldMin valSldMax], 'Callback',{@updateColormap,0});
            set([sldMin sldMax], 'Callback',@updateColormap);
            set([sld etxt btPlayAll btPlayZoom], 'Enable','on');
            drawGuiHist(1);
            drawPlots;
            updateImages;
            if ~get(tb_moveHist,'Value')
                updateGuiHistVal;
            end
            updateColormap;
        else
            %RGB mode on
            if rgbCount == 1 || any(get(bg_colormap, 'SelectedObject')==[cmManual cmImage cmZoom cmThresh]) 
                set(bg_colormap, 'SelectedObject', cmGlobal);
            end
            notXY = setdiff(1:nDim, xySel);
            rgbDim = notXY(rgbCount);
            set(tbSwitchRGB, 'String', dimNames(rgbDim), 'CData', RGB2, 'FontSize', 7, 'Value', 1);
            set(cmImage, 'String', 'Channel','ToolTipString', 'Scale each color channel to its range.');
            set(cmManual, 'Visible', 'off');
            set(cmZoom, 'Visible','off');
            set([cmThresh popThresh], 'Visible', 'off');
%             set([valSldMin valSldMax], 'Enable', 'on','Callback',{@updateImages,0});
            set([sldMin sldMax], 'Enable', 'on','Callback',@updateImages);
            set(popLut, 'Enable', 'off');
            set(tbColorbar, 'Enable', 'off', 'Value', 0);
            if rgbCount == get(projDimPop,'Value')
              projMethod =0;
              set(projMethodPop,'Enable', 'off');
            else
              projMethod = get(projMethodPop, 'Value') - 1;
              set(projMethodPop,'Enable', 'on');
            end
            %             set([sldGamma valSldGamma strGamma], 'Visible', 'off');
            set([sldMin sldMax], 'Callback',@updateImages);
            set([valSldMin_RGB valSldMax_RGB sldMin_RGB sldMax_RGB], 'Visible', 'on');
            for ii=1:nMat
                colorbar('peer', imAx(ii),'off');
            end
            set([cmStretchRGBMean cmStretchRGBMax], 'Visible', 'on');
            set(bg_colormap,'SelectionChangeFcn', @rgbDisplay,'UserData',get(bg_colormap,'Tag'),...
                'Tag','Color adjustment. Global: Scale RGB images according to slider settings. Channel: Equalize range of images before applying slider settings. MeanStretch / MaxStretch: Stretch RGB over entire dimension with color coding for position.');
            rgbDisplay;
            colormap(contrastAx,gray(255));
            if get(btMean, 'UserData') == 5
                drawPlots;
            end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Options for RGB display (Global/Channel/Image vs. Stretch)
    function rgbDisplay(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        set([sld etxt btPlayAll btPlayZoom], 'Enable','on');
        if nargin == 3
            if get(sldMin_RGB, 'Value') > get(sldMax_RGB, 'Value')
                set(sldMin_RGB, 'Value', get(sldMax_RGB, 'Value') - 0.01);
            end
        elseif any(get(bg_colormap, 'SelectedObject') ==  [cmStretchRGBMean cmStretchRGBMax])
            set([sldGamma valSldGamma strGamma], 'Visible', 'off');
            set([valSldMin_RGB valSldMax_RGB sldMin_RGB sldMax_RGB rgbSldIm], 'Visible', 'on');
            set([sld(rgbDim) etxt(rgbDim) btPlayAll(rgbDim) btPlayZoom(rgbDim)], 'Enable','off');
            drawGuiHist(dim(rgbDim));
        else
            set([sldGamma valSldGamma strGamma], 'Visible', 'on');
            set([valSldMin_RGB valSldMax_RGB sldMin_RGB sldMax_RGB rgbSldIm], 'Visible', 'off');
            if dim(rgbDim) > 2
                drawGuiHist(3);
            else
                drawGuiHist(2);
            end
        end
        updateImages('forceUpdateGuiHist');
        if get(btMean, 'UserData') == 5  %for RGB plots
            drawPlots;
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function buttonDownGuiHist(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        switch get(gui,'SelectionType')
            %Lin/log for middle click
            case 'extend'
                yl = get(histAxGui, 'YLim');
                if strcmp(get(histAxGui, 'YScale'), 'log')
                    set(histAxGui, 'YScale','lin', 'YLim',[0 yl(2)]);
                    if strcmp(get(histAxPlot(1),'Type'),'bar') 
                      for ii=1:length(histAxPlot)
                        set(histAxPlot(ii), 'BaseValue',0);
                      end
                    end
                    set(histAxBg, 'Position', [cmMinMax(currContrastSel,1) 0 cmMinMax(currContrastSel,2)-cmMinMax(currContrastSel,1) yl(2)]);
                    set(txt_linlog, 'String', 'lin');
                else
                    baseValLog = .99;
                    set(histAxBg, 'Position', [cmMinMax(currContrastSel,1) baseValLog cmMinMax(currContrastSel,2)-cmMinMax(currContrastSel,1) yl(2)]);
                    set(histAxGui, 'YScale','log', 'YLim',[baseValLog yl(2)]);
                    if length(histAxPlot) > 1
                        % Update of BaseValue does not work for some reason (in
                        % case of 3 bar plots (RGB), only one is updated). Thus
                        % the objects are created new each time :(
%                         drawGuiHist(length(histAxPlot));
%                         updateImages('forceUpdateGuiHist');
                    else
                        set(histAxPlot, 'BaseValue',baseValLog);
                    end
                    set(txt_linlog, 'String', 'log');
                end
            case 'alt' % Switch between normal x-axes and gamma-stretched x-axes
                if ~rgbCount
                    if strcmp(get(histAxGui, 'UserData'),'gamma')
                        set(histAxGui, 'UserData','');
                        histCDataGamma = applyGamma(linspace(minVal(1),maxVal(1),numel(histXDataBin)),cmMinMax(currContrastSel,:),currGamma(currContrastSel),1);
                    else
                        set(histAxGui, 'UserData','gamma');
                        histCDataGamma = linspace(minVal(1),maxVal(1),numel(histXDataBin));
                    end
                    set(contrastSldIm, 'CData',histCDataGamma);
                    updateGuiHistVal;
                end
                %In RGB mode switch order of bar plots
            case  'normal'
                if rgbCount && ~any(get(bg_colormap, 'SelectedObject') ==  [cmStretchRGBMean cmStretchRGBMax])
                    chldr{1} = findobj(histAxGui, 'Type', 'hggroup');
                    chldr{2} = setdiff(get(histAxGui, 'Children'), chldr{1});
                    set(histAxGui, 'Children', [circshift(chldr{1},1)' chldr{2}']);
                end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function meanPlots(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if varargin{1} == -1
            incr = -1;
        else
            incr = 1;
        end
        if rgbCount
            set(btMean, 'UserData', mod(get(btMean, 'UserData')+incr  ,6));
        else
            set(btMean, 'UserData', mod(get(btMean, 'UserData')+incr  ,5));
        end
        set(btMean, 'CData', squeeze(meanBt{get(btMean, 'UserData')+1}));
        drawPlots;
        if ~isempty(plotDim)
            updateObjects;
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
%% Mouse Controls
%Mouse Movement Controls
    function mouseMotion(varargin)
        if debugMatVis > 1, debugMatVisFcn(1); end
        % Update tooltips
        if get(tbTooltips, 'Value')
            updateTooltips;
        end
        if any(myGcf == imageWin) || any(myGcf == zoomWin)
            %Get current position of cursor
            p = get(myGca, 'CurrentPoint');
            p = p(1 ,1:2);
            if customDimScale
                pLoc = p;
                p = pixelLocation(p);
            end
            p = round(p([2 1]));
            p(1) = max([min([p(1),size(currIm{1},1)]),1]);            %max([min([p(1),dim(xySel(2))]),1]);
            p(2) = max([min([p(2),size(currIm{1},2)]),1]);
            %Display position / value of current point while inside Image
            %window
            if any(myGcf == imageWin)
                for ii = 1:nMat
                    if customDimScale
                        if withDimUnits
                            set(imageWin(ii), 'Name', ['Image (',varName{ii},') - Value: ',num2str(currImVal{ii}(p(1),p(2),:)) '  /  Position: (' num2str(round(pLoc(2)/pxWidth(xySel(1)))*pxWidth(xySel(1))) ' ' dimUnits{xySel(1)} ', ' num2str(round(pLoc(1)/pxWidth(xySel(2)))*pxWidth(xySel(2))) ' ' dimUnits{xySel(2)} ')  /  Pixelpos.: (', num2str(p(1)) ', ' num2str(p(2)),')']);
                        else
                            set(imageWin(ii), 'Name', ['Image (',varName{ii},') - Value: ',num2str(currImVal{ii}(p(1),p(2),:)) '  /  Position: (' num2str(pLoc(2)) ', ' num2str(pLoc(1))  ')  /  Pixelpos.: (', num2str(p(1)) ', ' num2str(p(2)),')']);
                        end
                    else
                        set(imageWin(ii), 'Name', ['Image (',varName{ii},') - Value: ',...
                            num2str(currImVal{ii}(p(1),p(2),:)) ' / Position: (', num2str(p(1)) ', ' num2str(p(2)),')']);
                    end
                end
                %Change cursor while inside / outside zoom region in
                %image window if not in roi selection mode
                if ~any(cell2mat(get(tb_newRoi, 'Value')))
                    if pointInArea(p([2 1]),zoomValXY) && (zoomValXY(3) < size(currIm{1},2) || zoomValXY(4) < size(currIm{1},1))
                        posInZoomArea = 1;
                        set(imageWin,'Pointer','custom','PointerShapeCData',hand, 'PointerShapeHotSpot',[9 9]);
                    else
                        posInZoomArea = 0;
                    end
                end
            else
                posInZoomArea = 0;
            end
            pointOnLine = zeros(1,2);
            if ~any(cell2mat(get(tb_newRoi, 'Value'))) %~get(tbRoi,'Value')
                if get(tbShowObjects, 'Value')
                    % Check whether cursor is close to position lines and change
                    % properties
                    if any(myGcf == zoomWin)
                        interv = zoomVal(xySel,2)*.01;
                    else
                        interv = dim(xySel)*.01;
                    end
                    for ii=1:2
                        if pointInArea(p(ii), [currPos(xySel(ii))-interv(ii) 2*interv(ii)])
                            pointOnLine(ii) = 1;
                        end
                    end
                end
                if any(pointOnLine)
                    if sum(pointOnLine) == 2
                        set(myGcf,'Pointer','fleur');
                    elseif pointOnLine(1)
                        set(myGcf,'Pointer','custom','PointerShapeCData',arrow_ud, 'PointerShapeHotSpot',[9 9]);
                    else
                        set(myGcf,'Pointer','custom','PointerShapeCData',arrow_lr, 'PointerShapeHotSpot',[9 9]);
                    end
                    set(myGcf, 'WindowButtonDownFcn', {@placePosLine,pointOnLine});
                elseif posInZoomArea
                    set(imageWin,'Pointer','custom','PointerShapeCData',hand, 'PointerShapeHotSpot',[9 9]);
                    set(imageWin, 'WindowButtonUpFcn', @endPanWindow);
                    set(imageWin, 'WindowButtonDownFcn', @startPanZoom);
                else
                    set(myGcf,'Pointer','arrow');
                    set(myGcf, 'WindowButtonUpFcn', '');
                    set(myGcf, 'WindowButtonDownFcn', @buttonDownCallback);
                end
            end
            %Display position / value of current point while inside Zoom
            %window
            if any(myGcf == zoomWin)
                for ii = 1:nMat
                    if customDimScale
                        if withDimUnits
                            set(zoomWin(ii), 'Name', ['Image (',varName{ii},') - Value: ',num2str(currImVal{ii}(p(1),p(2),:)) '  /  Position: (' num2str(round(pLoc(2)/pxWidth(xySel(1)))*pxWidth(xySel(1))) ' ' dimUnits{xySel(1)} ', ' num2str(round(pLoc(1)/pxWidth(xySel(2)))*pxWidth(xySel(2))) ' ' dimUnits{xySel(2)} ')  /  Pixelpos.: (', num2str(p(1)) ', ' num2str(p(2)),')']);
                        else
                            set(zoomWin(ii), 'Name', ['Image (',varName{ii},') - Value: ',num2str(currImVal{ii}(p(1),p(2),:)) '  /  Position: (' num2str(pLoc(2)) ', ' num2str(pLoc(1))  ')  /  Pixelpos.: (', num2str(p(1)) ', ' num2str(p(2)),')']);
                        end
                    else
                        set(zoomWin(ii), 'Name', ['Image (',varName{ii},') - Value: ',...
                            num2str(currImVal{ii}(p(1),p(2),:)) ' / Position: (', num2str(p(1)) ', ' num2str(p(2)),')']);
                    end
                end
            end
        end
        % Code for dragging position line with left click when mouse is close
        % to current position.
        if any(myGcf == plotWin) && ~isempty(subPlotHandles)
            p = get(myGca, 'CurrentPoint');
            p = p(1  ,1);
            if customDimScale
                p = pixelLocation(p(1),plotDim(subPlotHandles(myGcf == plotWin,:) == myGca));
            end
            inPlot = plotDim(subPlotHandles(myGcf == plotWin,:) == myGca);
            if p >= currPos(inPlot)-diff(get(subPlotHandles(subPlotHandles(myGcf == plotWin,:)==myGca),'XLim')+1)*.01 && p <= currPos(inPlot)+diff(get(subPlotHandles(subPlotHandles(myGcf == plotWin,:)==myGca),'XLim')+1)*.01
                set(myGcf,'Pointer','custom','PointerShapeCData',arrow_lr,'PointerShapeHotSpot',[8 8]);
                set(myGcf, 'WindowButtonDownFcn', @placePosLine);
            else
                set(myGcf,'Pointer','arrow');
                set(myGcf, 'WindowButtonDownFcn', @buttonDownCallback);
            end
        end
        if debugMatVis > 1, debugMatVisFcn(2); end
    end
    function px = pixelLocation(px,varargin)
        if debugMatVis > 1, debugMatVisFcn(1); end
        if nargin == 1 % Default xy case
            px = ((px-dimScale(xySel([2 1]),1)')'./(diff(dimScale(xySel([2 1]),:),1,2)).*(dim(xySel([2 1]))'-1))'+1;
        else  % Different dimension specified as second argument
            px = ((px-dimScale(varargin{1},1)')'./(diff(dimScale(varargin{1},:),1,2)).*(dim(varargin{1})'-1))'+1;
        end
        if debugMatVis > 1, debugMatVisFcn(2); end
    end
    function px = scaledLocation(px,varargin)
        if debugMatVis > 1, debugMatVisFcn(1); end
        if nargin == 1 % Default xy case
            px = ((px-[1 1])'.*diff(dimScale(xySel([2 1]),:),1,2)./(dim(xySel([2 1]))-[1 1])')'+dimScale(xySel([2 1]),1)';
        else  % Different dimension specified as second argument
            px = ((px-1)'.*diff(dimScale(varargin{1},:),1,2)./(dim(varargin{1})-1)')'+dimScale(varargin{1},1)';
        end
        if debugMatVis > 1, debugMatVisFcn(2); end
    end
    function toggleTooltipDisplay(varargin)
		if debugMatVis, debugMatVisFcn(1); end
        if strcmp(os(1:3),'MAC')
            f = macScaleFactor(2);
        else
            f = 1;
        end
        if ~get(tbTooltips, 'Value')   % Hide headings and tooltips
            set([sepLine(end) txt_titles txt_tooltips txt_plots], 'Visible','off');
            set(panel_positionControls, 'Position', get(panel_positionControls, 'Position') - f*[0 110 0 0]);
            set(panel_imageControls, 'Position', get(panel_imageControls, 'Position') - f*[0 90 0 0]);
            set(panel_windowButtons, 'Position', get(panel_windowButtons, 'Position') - f*[0 67 0 0]);
            set(sepLine(2), 'Position',get(sepLine(2), 'Position') - f*[0 90 0 0]);
            set(sepLine(1), 'Position',get(sepLine(1), 'Position') - f*[0 110 0 0]);
            set(gui, 'Position', get(gui,'Position')-f*[0 -130 0 130]);
        else % Show headings and tooltips
            set([sepLine(end) txt_titles txt_tooltips txt_plots], 'Visible','on');
            set(panel_positionControls, 'Position', get(panel_positionControls, 'Position') + f*[0 110 0 0]);
            set(panel_imageControls, 'Position', get(panel_imageControls, 'Position') + f*[0 90 0 0]);
            set(panel_windowButtons, 'Position', get(panel_windowButtons, 'Position') + f*[0 67 0 0]);
            set(sepLine(2), 'Position',get(sepLine(2), 'Position') + f*[0 90 0 0]);
            set(sepLine(1), 'Position',get(sepLine(1), 'Position') + f*[0 110 0 0]);
            set(gui, 'Position', get(gui, 'Position')+f*[0 -130 0 130]);
        end
		if debugMatVis, debugMatVisFcn(2); end
    end

%     function toggleMoveHist(varargin)
%         if get(tb_moveHist, 'Value')
%             playbackUpdateGuiHist = 1;
%         else
%             playbackUpdateGuiHist = 0;
%         end
%     end

    function updateTooltips(varargin)
        if debugMatVis > 1, debugMatVisFcn(1); end
        if get(tbTooltips, 'Value')
            rootPoint = get(0, 'PointerLocation');
            % Tooltips for the elements of the main GUI
            if pointInArea(rootPoint, get(gui, 'Position'))
                guiPos = get(gui, 'Position');
                p = rootPoint - guiPos(1:2);
                if pointInArea(p, get(panel_positionControls, 'Position'))
                    objectParent = panel_positionControls;
                elseif pointInArea(p, get(panel_imageControls, 'Position'))
                    objectParent = panel_imageControls;
                elseif pointInArea(p, get(panel_windowButtons, 'Position'))
                    objectParent = panel_windowButtons;
                else
                    objectParent = [];
                end
                if ~isempty(objectParent)
                    panelChildren = get(objectParent,'Children');
                    panelPos = get(objectParent, 'Position');
                    p = p - panelPos(1:2) + [2 2];
                    s = [];
                    ii = 0;
                    while isempty(s) && ii < numel(panelChildren)
                        ii= ii + 1;
                        if pointInArea(p, get(panelChildren(ii),'Position')) && (strcmp('on',get(panelChildren(ii),'Visible'))  || strcmp('axes',get(panelChildren(ii),'Type')))
                            s = get(panelChildren(ii), 'Tag');
                        end
                    end
                end
                if isempty(objectParent) || isempty(s)
                    s = {'Right click to bring all associated visible windows on top. Middle click to move all associated windows together. Use number keys or ctrl+number keys to go forth or back in the different dimensions, respectively. If this does not work, click on an empty space in the main Gui first.'};
                end
            end
            % Tooltips for zoom, image and plot windows
            if pointInArea(rootPoint, get(zoomWin(1), 'Position')) || pointInArea(rootPoint, get(imageWin(1), 'Position')) || pointInArea(rootPoint, get(plotWin(1), 'Position'))
                s = {'Left click: Draw new zoom area. Affects plots only if range is fixed to zoom. Middle click: Move position cross / line. Right click: Reset zoom to full view. Double click: Copy figure content to clipboard. Scroll wheel: Zoom in/out.'};
            end
            % Tooltips for the ROI Manager
            if ~isempty(roiWin) && pointInArea(rootPoint, get(roiWin, 'Position'))
                roiGuiPos = get(roiWin, 'Position');
                p = rootPoint - roiGuiPos(1:2);
                roiGuiChildren = get(roiWin, 'Children');
                s = [];
                ii = 0;
                while isempty(s) && ii < numel(roiGuiChildren)
                    ii= ii + 1;
                    if pointInArea(p, get(roiGuiChildren(ii),'Position'))
                        s = get(roiGuiChildren(ii), 'Tag');
                    end
                end
                if isempty(s)
                    s = {'ROI Manager.'};
                end
            end
            % Set tooltips string
            set (txt_tooltips, 'String', s);
        end
        if debugMatVis > 1, debugMatVisFcn(2); end
    end
    function out = pointInArea(point,rect)
        if debugMatVis > 1, debugMatVisFcn(1); end
        out = 0;
        if numel(rect) == 2
            if point(1) > rect(1) && point(1) < sum(rect)
                out = 1;
            end
        elseif point(1) > rect(1) && ...
                point(1) < rect(1) + rect(3) && ...
                point(2) > rect(2) && ...
                point(2) < rect(2) + rect(4)
            out = 1;
        end
        if debugMatVis > 1, debugMatVisFcn(2); end
    end

%Functions for dragging zoom rectangle
    function startPanZoom(varargin)
        if debugMatVis > 1, debugMatVisFcn(1); end
        switch get(myGcf,'SelectionType')
            case {'normal';'alt'}   % Control + left click: Start dragging the zoom area
                set(myGcf,'Pointer','custom','PointerShapeCData',hand2, 'PointerShapeHotSpot',[9 9]);
                pStart = get(myGca, 'CurrentPoint');
                pStart = pStart(1  ,1:2);
                if customDimScale
                    pStart = pixelLocation(pStart);
                end
                zoomStart = zoomValXY;
                set(myGcf, 'WindowButtonMotionFcn', @panWindow);
            case 'extend'   % Middle click: drag position cross (via buttonDownCallback)
                buttonDownCallback;
        end
        if debugMatVis > 1, debugMatVisFcn(2); end
    end
    function endPanWindow(varargin)
        if debugMatVis > 1, debugMatVisFcn(1); end
        set(myGcf,'Pointer','custom','PointerShapeCData',hand, 'PointerShapeHotSpot',[9 9]);
        if nargin == 3
            set(myGcf, 'WindowButtonMotionFcn', @mouseMotion);
            if debugMatVis > 1, debugMatVisFcn(2); end
            return;
        end
        currWin = myGcf;
        pCurr = (get(myGca, 'CurrentPoint'));
        pCurr = pCurr(1  ,1:2);
        if customDimScale
            pCurr = pixelLocation(pCurr);
        end
        if any(myGcf == zoomWin)
            zoomValXY = [zoomStart(1)+round(pStart(1)-pCurr(1)), zoomStart(2)+round(pStart(2)-pCurr(2)),...
                zoomStart(3), zoomStart(4)];
        else
            zoomValXY = [zoomStart(1)-pStart(1)+pCurr(1), zoomStart(2)-pStart(2)+pCurr(2),...
                zoomStart(3), zoomStart(4)];
        end
        updateZoom;
        updatePlots;
        figure(currWin);
        set(myGcf, 'WindowButtonMotionFcn', @mouseMotion);
        if debugMatVis > 1, debugMatVisFcn(2); end
    end
    function panWindow(varargin)
        if debugMatVis > 1, debugMatVisFcn(1); end
        pCurr = (get(myGca, 'CurrentPoint'));
        pCurr = pCurr(1  ,1:2);
        if customDimScale
            pCurr = pixelLocation(pCurr);
        end
        if ~any(myGcf == zoomWin)
            zoomValXY = [zoomStart(1)-pStart(1)+pCurr(1), zoomStart(2)-pStart(2)+pCurr(2),...
                zoomStart(3), zoomStart(4)];
            updateZoom;
        end
        if debugMatVis > 1, debugMatVisFcn(2); end
    end
    function buttonDownGui(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        %Check for mouse position to enable fake 'button
        %right-clicks' in main GUI
        if strcmp(get(gui,'SelectionType'),'alt')
            p1 = round(get(gui,'CurrentPoint'));
            btPosAspectRatio = get(tbAspRatio,  'Position') + get(panel_imageControls,    'Position') .* [1 1 0 0];
            btPosPlotsXLim   = get(tbPlotsXLim, 'Position') + get(panel_positionControls, 'Position') .* [1 1 0 0];
            btPosPlotsYLim   = get(tbPlotsYLim, 'Position') + get(panel_positionControls, 'Position') .* [1 1 0 0];
            btPosMean        = get(btMean,      'Position') + get(panel_positionControls, 'Position') .* [1 1 0 0];
            btPos_100pct     = get(bt_100pct,   'Position') + get(panel_imageControls, 'Position') .* [1 1 0 0];
            % Right click on play buttons: play backwards
            revPlay = NaN;
            for ii = 1:nDim
                btPosPlay = get(btPlayAll(ii), 'Position')+ get(panel_positionControls, 'Position') .* [1 1 0 0];
                btPosPlayZoom = get(btPlayZoom(ii), 'Position')+ get(panel_positionControls, 'Position') .* [1 1 0 0];
                if pointInArea(p1,btPosPlay)
                    revPlay = ii;
                elseif pointInArea(p1,btPosPlayZoom)
                    revPlay = -ii;
                end
            end
            if ~isnan(revPlay)
                if revPlay > 0
                    set(btPlayAll(revPlay), 'Value',1);
                    playCallback(0,0,revPlay,'all',0);
                else
                    set(btPlayZoom(-1*revPlay), 'Value',1);
                    playCallback(0,0,-1*revPlay,'zoom',0);
                end
            %Right click on aspect ratio button: set aspect ratio if
            %button value is 1
            elseif pointInArea(p1,btPosAspectRatio) && get(tbAspRatio, 'Value')
                gp = get(gui,'Position');
                set(tempWin, 'HandleVisibility', 'on');
                if ~isempty(tempWin) && any(get(0,'Children') == tempWin)
                    set(tempWin, 'HandleVisibility', 'on');
                    close(tempWin);
                    tempWin = [];
                end
                tempWin = figure('Position',[gp(1) gp(2)-92 150 50],...
                    'Name', 'Aspect Ratio','MenuBar', 'none', 'WindowStyle','normal', ...
                    'Resize', 'off', 'NumberTitle','off', 'HandleVisibility', 'off', 'CloseRequestFcn',@closeTempWin);
                uicontrol(tempWin, 'Style', 'Text', 'Position', [10 30 130 15], ...
                    'String','Specify Aspect Ratio','FontWeight', 'bold',...
                    'BackgroundColor', get(tempWin, 'Color'), 'HorizontalAlignment', 'left');
                % strAspRatio:
                uicontrol(tempWin, 'Style', 'Edit', 'Position', [40 10 50 15], ...
                    'String',['[',num2str(aspRatio(1)),' ',num2str(aspRatio(2)),']'],'FontWeight', 'bold',...
                    'HorizontalAlignment', 'center', 'Callback', @updateAspRatio);
                %Right click on XLim button: set plot XLim if
                %button value is 1
            elseif pointInArea(p1,btPosPlotsXLim) && get(tbPlotsXLim, 'Value')
                gp = get(gui,'Position');
                set(tempWin, 'HandleVisibility', 'on');
                if ~isempty(tempWin) && any(get(0,'Children') == tempWin)
                    close(tempWin);
                    tempWin = [];
                end
                tempWin = figure('Position',[gp(1) gp(2)-(70+nDim*20) 150 30+nDim*20],...
                    'Name', 'XLim','MenuBar', 'none', 'WindowStyle','normal', ...
                    'Resize', 'off', 'NumberTitle','off', 'HandleVisibility', 'off', 'CloseRequestFcn', @closeTempWin);
                uicontrol(tempWin, 'Style', 'Text', 'Position', [10 nDim*20+10 130 15], ...
                    'String','Specify x limits','FontWeight', 'bold',...
                    'BackgroundColor', get(tempWin, 'Color'), 'HorizontalAlignment', 'left');
                for ii=1:nDim
                    uicontrol(tempWin, 'Style', 'Text', 'Position', [10 nDim*20+20-ii*20-11 50 15], ...
                        'String',dimNames{ii} ,'FontWeight', 'bold',...
                        'BackgroundColor', get(tempWin, 'Color'), 'HorizontalAlignment', 'left');
                    % strXLim(ii):
                    uicontrol(tempWin, 'Style', 'Edit', 'Position', [50 nDim*20+20-ii*20-10 90 15], ...
                        'String',['[',num2str(plotXLim(ii,1)),' ',num2str(plotXLim(ii,2)),']'],...
                        'HorizontalAlignment', 'center', 'Callback', {@updatePlotLim,'x'});
                end
                %Right click on YLim button: set plot YLim if
                %button value is 1
            elseif pointInArea(p1, btPosPlotsYLim) && get(tbPlotsYLim, 'Value')
                gp = get(gui,'Position');
                set(tempWin, 'HandleVisibility', 'on');
                if ~isempty(tempWin) && any(get(0,'Children') == tempWin)
                    close(tempWin);
                    tempWin = [];
                end
                tempWin = figure('Position',[gp(1) gp(2)-(70+nMat*20) gp(3) 30+nMat*20],...
                    'Name', 'YLim','MenuBar', 'none', 'CloseRequestFcn', @closeTempWin,...
                    'Resize', 'off', 'NumberTitle','off', 'HandleVisibility', 'off', 'WindowStyle','normal');
                uicontrol(tempWin, 'Style', 'Text', 'Position', [10 nMat*20+10 100 15], ...
                    'String','Specify y limits','FontWeight', 'bold',...
                    'BackgroundColor', get(tempWin, 'Color'), 'HorizontalAlignment', 'left');
                for ii=1:nMat
                    uicontrol(tempWin, 'Style', 'Text', 'Position', [10 nMat*20+20-ii*20-11 150 15], ...
                        'String',varName{ii} ,'FontWeight', 'bold',...
                        'BackgroundColor', get(tempWin, 'Color'), 'HorizontalAlignment', 'left');
                    uicontrol(tempWin, 'Style', 'Edit', 'Position', [160 nMat*20+20-ii*20-10 120 15], ...
                        'String',['[',num2str(plotYLim(ii,1),'%6.2f'),'   ',num2str(plotYLim(ii,2),'%6.2f'),']'],...
                        'HorizontalAlignment', 'center', 'Callback', {@updatePlotLim,'y'});
                end
                if debugMatVis, debugMatVisFcn(2); end
                return
                %Right click on meanPlot button: Reduce setting mode by 1
            elseif pointInArea(p1, btPosMean)
                meanPlots(-1);
                %Right click on 100% button: set scaling factor
            elseif pointInArea(p1, btPos_100pct)
                gp = get(gui,'Position');
                set(tempWin, 'HandleVisibility', 'on');
                if ~isempty(tempWin) && any(get(0,'Children') == tempWin)
                    close(tempWin);
                    tempWin = [];
                end
                tempWin = figure('Position',[gp(1) gp(2)-(70+nDim*20) 150 30+nDim*20],...
                    'Name', 'Win Scale','MenuBar', 'none', 'WindowStyle','normal', ...
                    'Resize', 'off', 'NumberTitle','off', 'HandleVisibility', 'off', 'CloseRequestFcn', @closeTempWin);
                uicontrol(tempWin, 'Style', 'Text', 'Position', [10 nDim*20+10 130 15], ...
                    'String','Window scaling factor [%]','FontWeight', 'bold',...
                    'BackgroundColor', get(tempWin, 'Color'), 'HorizontalAlignment', 'left');
                uicontrol(tempWin, 'Style', 'Edit', 'Position', [50 10 40 15], ...
                    'String',[num2str(currWinScale)],'Callback',@setWinScale,...
                    'HorizontalAlignment', 'center');
            else
                %Right click on main gui: Bring all visible windows to front
                if get(tbTifPar, 'Value')
                    figure(tifParFig);
                end
                if get(tbHist,'Value') && ~withAlpha
                    figure(histWin);
                end
                if get(tbRoi,'Value')
                    figure(roiWin);
                end
                if get(tbProfile,'Value')
                    figure(profileWin);
                    figure(profileTraceWin);
                end
                if get(tbWin(1),'Value')
                    for ii=1:nMat
                        figure(imageWin(ii));
                    end
                end
                if get(tbWin(2),'Value')
                    for ii=1:nMat
                        figure(zoomWin(ii));
                    end
                end
                if get(tbWin(3),'Value')
                    for ii=1:nMat
                        figure(plotWin(ii));
                    end
                end
                if get(tbRecord, 'Value')
                    figure(movdata.gui.hvidgen);
                    figure(movdata.prev.hprev);
                end
                if get(roiBtRename, 'Value')
                    figure(tempRoiPropWin);
                end
            end
        elseif strcmp(get(gui,'SelectionType'),'extend')  % Move all windows together 
            startPosGui = get(gui, 'Position');
            p1 = startPosGui(1:2) + get(gui,'CurrentPoint');
            set(gui,'WindowButtonUpFcn',@releaseMiddleClickGui);
            set(gui,'WindowButtonMotionFcn',@moveAllWindows);
            for ii = 1:numel(data)
                startPosWins(1,ii,:) = get(imageWin(ii),'Position');
                startPosWins(2,ii,:) = get(zoomWin(ii),'Position');
                startPosWins(3,ii,:) = get(plotWin(ii),'Position');
                if ~isempty(histWin);         startPosWins(4,ii,:) = get(histWin,'Position');         end
                if ~isempty(roiWin);          startPosWins(5,ii,:) = get(roiWin,'Position');          end
                if ~isempty(profileWin);      startPosWins(6,ii,:) = get(profileWin,'Position');      end
                if ~isempty(profileTraceWin); startPosWins(7,ii,:) = get(profileTraceWin,'Position'); end
            end
        end
        if debugMatVis, debugMatVisFcn(2); end
        function moveAllWindows(varargin)
            currPosGui = get(gui, 'Position');
            p2 = currPosGui(1:2) + get(gui,'CurrentPoint');
            set(gui, 'Position', startPosGui + [p2-p1 0 0]);
            for iii = 1:numel(data)
                set(imageWin(iii),'Position', squeeze(startPosWins(1,iii,:))' + [p2-p1 0 0]);
                set(zoomWin(iii),'Position', squeeze(startPosWins(2,iii,:))' + [p2-p1 0 0]);
                set(plotWin(iii),'Position',squeeze(startPosWins(3,iii,:))' + [p2-p1 0 0]);
                if ~isempty(histWin);         set(histWin,'Position',squeeze(startPosWins(4,iii,:))' + [p2-p1 0 0]);         end
                if ~isempty(roiWin);          set(roiWin,'Position',squeeze(startPosWins(5,iii,:))' + [p2-p1 0 0]);          end
                if ~isempty(profileWin);      set(profileWin,'Position',squeeze(startPosWins(6,iii,:))' + [p2-p1 0 0]);      end
                if ~isempty(profileTraceWin); set(profileTraceWin,'Position',squeeze(startPosWins(7,iii,:))' + [p2-p1 0 0]); end
            end
        end
        function releaseMiddleClickGui(varargin)
            set(gui,'WindowButtonMotionFcn',[]);
        end
    end
% resize function for image, zoom, plot windows to apply position
% change to all corresponding figures (e.g. all plot figures)
    function resizeWin(varargin)
        if debugMatVis > 1, debugMatVisFcn(1); end
        if nMat > 1 && get(tbLinkWin, 'Value') && any(myGcf==[zoomWin imageWin plotWin])
            switch varargin{3}
                case 'zoom'
                    figHandles = zoomWin;
                case 'image'
                    figHandles = imageWin;
                case 'plot'
                    figHandles = plotWin;
            end
            if nargin == 4
                cf = varargin{4};
            else
                cf = myGcf;
            end
            newPos = get(cf, 'Position');
            newPos(2) = newPos(2)+(newPos(4)+winWidthTop+8)*(find(myGcf == figHandles)-1);
            for ii = 1:numel(figHandles)
                oldPos = get(figHandles(ii), 'Position');
                set(figHandles(ii), 'Position', [newPos(1) newPos(2)-(newPos(4)+winWidthTop+8)*(ii-1) newPos(3:4)]);
            end
        end
        if debugMatVis > 1, debugMatVisFcn(2); end
    end
% Callback function for toggle button for (un)linking figure position
% and size
    function linkWins(varargin)
        if debugMatVis > 1, debugMatVisFcn(1); end
        if get(tbLinkWin, 'Value')
            set(tbLinkWin, 'CData', icon_linkWins);
            set([imageWin zoomWin plotWin], 'ResizeFcn', '');
            resizeWin(0,0,'image', imageWin(1));
            drawnow
            resizeWin(0,0,'zoom', zoomWin(1));
            drawnow
            resizeWin(0,0,'plot', plotWin(1));
            drawnow
            set(imageWin, 'ResizeFcn', {@resizeWin, 'image'});
            set(zoomWin, 'ResizeFcn', {@resizeWin, 'zoom'});
            set(plotWin, 'ResizeFcn', {@resizeWin, 'plot'});
        else
            set(tbLinkWin, 'CData', icon_unlinkWins);
        end
        if debugMatVis > 1, debugMatVisFcn(2); end
    end
%Button Press Controls for image, zoom and plot windows
    function buttonDownCallback(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        currWin = myGcf;
        % Leave if window is empty (can happen for plotWin in ROI mode)
        if isempty(get(currWin, 'Children'))
            if debugMatVis
                display([repmat(' ',[1 debugIndent*fctLevel]) num2str(fctCount.(ST)) ': End   buttonDownCallback']);
                fctLevel = fctLevel-1;
            end
            %return
        else
          % "Normal" button clicks
          switch get(currWin,'SelectionType')
            case 'normal'  %Left click
              %Zoom by drawing zoom region
              set(myGcf, 'HandleVisibility', 'on');
              p1 = get(myGca,'CurrentPoint');
              p1 = p1(1,1:2);
              rbbox;
              p2 = get(myGca,'CurrentPoint');
              p2 = p2(1,1:2);
              if customDimScale
                if any(currWin == [imageWin zoomWin])
                  p1 = pixelLocation(p1);
                  p2 = pixelLocation(p2);
                elseif any(currWin == plotWin)
                  p1 = pixelLocation(p1(1,1),plotDim(find(subPlotHandles(plotWin==myGcf,:) == myGca)));
                  p2 = pixelLocation(p2(1,1), plotDim(find(subPlotHandles(plotWin==myGcf,:) == myGca)));
                end
              end
              p1 = ceil(p1);
              p2 = ceil(p2);
              set(myGcf, 'HandleVisibility', 'off');
              if p1 == p2
                if debugMatVis, debugMatVisFcn(2); end
                return;
              end
              if any(currWin == plotWin) && ~isempty(plotDim)
                subPlot = find(subPlotHandles(plotWin==myGcf,:) == myGca);
                if plotDim(subPlot) == xySel(2)
                  p1(1  ,2) = zoomValXY(2);
                  p2(1  ,2) = zoomValXY(2) + zoomValXY(4);
                elseif plotDim(subPlot) == xySel(1)
                  p1(1  ,2) = p1(1  ,1);
                  p2(1  ,2) = p2(1  ,1);
                  p1(1  ,1) = zoomValXY(1);
                  p2(1  ,1) = zoomValXY(1) + zoomValXY(3);
                else
                  if p1(1  ,1) ~= p2(1  ,1)
                    zoomVal(plotDim(subPlot),1) = min(p1(1  ,1),p2(1  ,1));
                    zoomVal(plotDim(subPlot),2) = max(p1(1  ,1),p2(1  ,1)) - min(p1(1  ,1),p2(1  ,1)) + 1;
                    updateZoom;
                  end
                  if debugMatVis, debugMatVisFcn(2); end
                  return;
                end
              end
              coord(1  ,:)    = round(p1(1  ,1:2));
              coord(2,:)    = round(p2(1  ,1:2));
              if coord(1  ,1) == coord(2,1) || coord(1  ,2) == coord(2,2)
                if debugMatVis, debugMatVisFcn(2); end
                return;
              end
              if coord(1  ,1) > coord(2,1)
                coord(:,1) = flipdim(coord(:,1),1);
              end
              if coord(1  ,2) > coord(2,2)
                coord(:,2) = flipdim(coord(:,2),1);
              end
              zoomValXY(1:2) = coord(1  ,:);
              zoomValXY(3:4) = coord(2,:) - coord(1  ,:);
              updateZoom;
              updatePlots;
            case 'extend'    %Middle click
              placePosLine(0,0,[1 1]);
            case 'alt'  %Right click
              %Unzoom
              if any(currWin == [zoomWin imageWin])
                zoomValXY = [1  ,1  ,size(currIm{1},2),size(currIm{1},1)]; %dim(xySel(2)),dim(xySel(1))];
              else
                subPlot = find(subPlotHandles(plotWin==myGcf,:) == myGca);
                zoomVal(plotDim(subPlot),:) = [1 dim(plotDim(subPlot))];
                % Update zoomValXY here, otherwise it'll overwrite zoomVal
                % in updateZoom
                zoomValXY([1,3]) = zoomVal(xySel(2),:);
                zoomValXY([2,4]) = zoomVal(xySel(1),:);
              end
              updateZoom;
              updateObjects;
            case 'open'   % Double click: Copy content of  current figure to clipboard (only MS Windows)
              if strcmp(os(1:5),'PCWIN')
                  a = listdlg('PromptString','Copy/save/export figure content:',...
                              'SelectionMode','single',...
                              'ListString',{'Clipboard: Vector graphics';'Clipboard: Bitmap';'Save to file';'Export to workspace'});
%                 a = questdlg('Copy as vector graphics (enhanced meta-file) or 8-bit bitmap (bmp) or save as image file?','Copy figure content to clipboard or save to file','Vector graphics','Bitmap','Save','Vector graphics');

                switch a
                    case 1
                        print(myGcf, '-dmeta');  %#ok
                    case 2
                        print(myGcf, '-dbitmap');    %#ok
                    case 3
                        [fNameFig,pNameFig] = uiputfile({'*.jpg';'*.png';'*.tif';'*.pdf'},'Select file location, name and file type!');
                        saveas(myGcf, [pNameFig  fNameFig]);
                    case 4
                        dataSetIdx = mod(find(myGcf == [zoomWin imageWin plotWin]),numel(data))+1;
                        if isempty(varName{dataSetIdx})
                            wsExportName = 'matVisDataExport';
                        else
                            wsExportName = [varName{dataSetIdx} '_export'];
                        end
                        if any(myGcf == zoomWin)
                            assignin('base','test',myGcf);
                            assignin('base',wsExportName, currIm{dataSetIdx}(zoomVal(xySel(1),1):sum(zoomVal(xySel(1),:))-1,zoomVal(xySel(2),1):sum(zoomVal(xySel(2),:))-1,:));
                        else
                            assignin('base',wsExportName, currIm{dataSetIdx});
                        end
                end
              end
          end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

    function placePosLine(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        currWin = myGcf;
        dragging = 1;
        if nargin == 3
            xy = varargin{3};
            set(currWin, 'WindowButtonMotionFcn', {@mousePosition,find(xy)});
        else
            set(currWin, 'WindowButtonMotionFcn', @mousePosition);
        end
        set(currWin, 'WindowButtonUpFcn', @endMousePosition);
        pt = get(myGca, 'CurrentPoint');
        pt = pt(1 ,[1 2]);
        if any(currWin == [zoomWin imageWin])
            if customDimScale
                pt = pixelLocation(pt);
            end
            pt = round(pt([2 1]));
            if sum(xy) == 2
                set(myGcf,'Pointer','custom','PointerShapeCData',nan(16,16));
            end
            set([rectPlotArea_Zoom rectPlotArea_Im], 'Visible','on');
        end
        %In plot window set current position of respective dimension
        if any(currWin == plotWin) && ~isempty(plotDim) && ~isempty(subPlotPlots)
            if customDimScale
                pt = pixelLocation(pt(1),plotDim(find(subPlotHandles(myGcf == plotWin,:) == myGca)));
            end
            pt = round(pt(1));
            subPlot = find(subPlotHandles(myGcf == plotWin,:) == myGca);
            currPos(plotDim(subPlot)) = pt;
            updateSelection(plotDim(subPlot));
        else
            for ii=1:2
                if xy(ii)
                    currPos(xySel(ii)) = min([dim(xySel(ii)),max([1  ,pt(ii)])]);
                end
            end
            for ii = 1:nMat
                set(imageWin(ii), 'Name', ['Image (',varName{ii},') - Pos: (', num2str(currPos(xySel)),')  Val: ',...
                    num2str(currImVal{ii}(currPos(xySel(1)),currPos(xySel(2)),:))]);
                set(zoomWin(ii), 'Name', ['Zoom (',varName{ii},') - Pos: (', num2str(currPos(xySel)),') Val: ',...
                    num2str(currImVal{ii}(currPos(xySel(1)),currPos(xySel(2)),:)), '  Zoom: (', num2str(zoomValXY), ')']);
            end
            if get(btMean, 'UserData') == 4 && any(currWin == imageWin)
                zoomValXY = [currPos(xySel(2))-round(zoomValXY(3)/2) currPos(xySel(1))-round(zoomValXY(4)/2) zoomValXY(3) zoomValXY(4)];
                updateZoom;
            end
            updateSelection(find(xy));
            updateObjects;
            updatePlots;
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

% Enable panning in zoom windows while control button is pressed
    function zoomKeyPressFcn(src,evnt)  %#ok
        if debugMatVis, debugMatVisFcn(1); end
        if strcmp(evnt.Modifier,'control')
            oldButtonDownFcn = get(src, 'WindowButtonDownFcn');
            set(src, 'KeyReleaseFcn', {@zoomKeyReleaseFcn, oldButtonDownFcn});
            set(src, 'KeyPressFcn', '');
            set(src, 'WindowButtonDownFcn', @startPanZoom);
            set(src,'Pointer','custom','PointerShapeCData',hand, 'PointerShapeHotSpot',[9 9]);
            set(src, 'WindowButtonUpFcn', @endPanWindow);
            set(src, 'WindowButtonDownFcn', @startPanZoom);
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

    function zoomKeyReleaseFcn(src,evt,fcnHndl) %#ok
        if debugMatVis, debugMatVisFcn(1); end
        if strcmp(func2str(get(zoomWin(1), 'WindowButtonMotionFcn')), 'matVis/panWindow')
            if debugMatVis, debugMatVisFcn(2); end
            return
        end
        %         set(src, 'WindowButtonDownFcn', @buttonDownCallback);
        set(src,'Pointer','arrow');
        set(src, 'WindowButtonUpFcn', '');
        %         set(src, 'WindowButtonDownFcn', @buttonDownCallback);
        set(src, 'KeyPressFcn', @zoomKeyPressFcn);
        if nargin==3
            set(src, 'WindowButtonDownFcn', fcnHndl);
        else
            set(src, 'WindowButtonDownFcn', @buttonDownCallback);
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function scrollWheelCallbackGui(varargin)
        if debugMatVis > 1, debugMatVisFcn(1); end
        % When the scroll wheel is used over the main gui, the dimensions can be scrolled through if mouse is above one of the dimension sliders
        % Look for mouse position to find correct slider
        mousePos = round(get(gui, 'CurrentPoint'));
        % Check for pressed modifier keys to allow faster scrolling ("frame skipping")
        mdf  = get(gui,'currentModifier');
        % If two modifier keys are pressed 
        if isempty(mdf)
            mdf = '';
        elseif numel(mdf) == 2
            mdf = cell2mat(mdf);
        else 
            mdf = mdf{1};
        end
        switch mdf
            case 'control'
                f = 10;
            case 'shift'
                f = 100;
            case 'shiftcontrol'
                f = 1000;
            otherwise
                f = 1;
        end
        for ii = 1:nDim
            sldPos(:,ii) = get(sld(ii),'Position') + (get(panel_positionControls,'Position').*[1 1 0 0]);
        end
        sldSel = NaN;
        for ii = 1:nDim
            if pointInArea(mousePos,sldPos(:,ii))
                sldSel = ii;
            end
        end
        if ~isnan(sldSel)
            currPos(sldSel) = max(1,min(dim(sldSel),currPos(sldSel) - f*varargin{2}.VerticalScrollCount));
            updateSelection(sldSel);
        end
        if debugMatVis > 1, debugMatVisFcn(2); end
    end
    function scrollWheelCallback(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        z = varargin{2}.VerticalScrollCount;
        p = get(myGca, 'CurrentPoint');
        p = p(1  ,1:2);
        if customDimScale
            p = pixelLocation(p);
        end
        if any(myGcf == zoomWin)
            if any(p < zoomValXY(1:2)) || any(p > zoomValXY(1:2) + zoomValXY(3:4))
              if debugMatVis, debugMatVisFcn(2); end
              return
            end
            mousePos = p - zoomValXY(1:2);
            relMousePos = mousePos ./ zoomValXY(3:4);
            figPos = get(myGcf, 'Position');
            if z<0
                if zoomValXY(3)/zoomValXY(4) > figPos(3)/figPos(4)
                    zoomValXY(3) = 1.1^z * zoomValXY(3);
                else
                    zoomValXY(4) = 1.1^z * zoomValXY(4);
                end
            else
                zoomValXY(3:4) = 1.1^z * zoomValXY(3:4);
            end
            zoomValXY(1:2) = p - zoomValXY(3:4) .* relMousePos;
        else
            if any(p < [1 1]) || any(p > [size(currIm{1},2) size(currIm{1},1)])
                if debugMatVis, debugMatVisFcn(2); end
                return
            end
            zoomValXY(3:4) = 1.1^z * zoomValXY(3:4);
            zoomValXY(1:2) = p - zoomValXY(3:4) / 2;
        end
        updateZoom;
        if get(tbWin(3),'Value') && get(btMean, 'UserData') == 4 % Plot average over zoom area
            updatePlots;
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Drag Position Cross
    function mousePosition(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        currWin = myGcf;
        p = get(myGca, 'CurrentPoint');
        p = p(1  ,1:2);
        if any(currWin == plotWin) && ~isempty(plotDim)
            if customDimScale
                p = pixelLocation(p(1),plotDim(find(subPlotHandles(myGcf == plotWin,:) == myGca)));
            end
            p = round(p);
            dimNum = plotDim(subPlotHandles(myGcf == plotWin,:) == myGca);
            currPos(dimNum) = min([dim(dimNum),max([1  ,p(1)])]);
        else
            if customDimScale
                p = pixelLocation(p);
            end
            p = round(p);
            dimNum = varargin{3};
            for ii=1:2
                if any(dimNum == ii)
                    currPos(xySel(ii)) = min([dim(xySel(ii)),max([1  ,p(3-ii)])]);
                end
            end
            for ii = 1:nMat
                set(imageWin(ii), 'Name', ['Image (',varName{ii},') - Pos: (', num2str(currPos(xySel)),')  Val: ',...
                    num2str(currImVal{ii}(currPos(xySel(1)),currPos(xySel(2)),:))]);
                set(zoomWin(ii), 'Name', ['Zoom (',varName{ii},') - Pos: (', num2str(currPos(xySel)),') Val: ',...
                    num2str(currImVal{ii}(currPos(xySel(1)),currPos(xySel(2)),:)), '  Zoom: (', num2str(zoomValXY), ')']);
            end
            if get(btMean, 'UserData') == 4 && any(currWin == imageWin)
                zoomValXY = [currPos(xySel(2))-round(zoomValXY(3)/2) currPos(xySel(1))-round(zoomValXY(4)/2) zoomValXY(3) zoomValXY(4)];
                updateZoom;
            end
        end
        updateSelection(dimNum);
        if debugMatVis, debugMatVisFcn(2); end
    end

%Release Mouse Button after Dragging
    function endMousePosition(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        dragging = 0;
        set(myGcf, 'WindowButtonMotionFcn', @mouseMotion,'Pointer','arrow');
        mouseMotion; % Call to update pointer if mouse is on position line
        set([rectPlotArea_Zoom rectPlotArea_Im], 'Visible','off');
        %         updateSelection(0);
        % Update GUI hist
        if ~get(tb_moveHist,'Value')
            updateGuiHistVal
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
%% Update Functions
%Update current selection (called during change of current position)
    function updateSelection(dimNum)
        if debugMatVis, debugMatVisFcn(1); end
        %         if nargin == 3  % ???
        %             set(sld(dimNum), 'Value', currPos(dimNum));
        %             set(etxt(dimNum), 'String', num2str(currPos(dimNum)));
        %         else
        if nargin == 0
            dimNum = [];
        end
        for ii = dimNum
            set(sld(ii), 'Value', currPos(ii));
            set(etxt(ii), 'String', num2str(currPos(ii)));
        end
        %         end
        % Call from position buttons
        if nargin == 0
            updateObjects;
            updateImages;
            if (get(tbProfile,'Value') && nProfiles)
               updateProfileSelection; 
            end
        else
            % Other calls
            imageUpdated = 0;
            if numel(unique([xySel dimNum]))>2 %(numel(dimNum) > 1) ||(numel(dimNum) == 1 && ~any(xySel == dimNum))
                updateImages;%(dimNum)
                imageUpdated = 1;
                if ~isempty(globalHist) && strcmp(get(tbHist, 'Enable'),'on') && ~calculatingGlobalHist
                    updateHist;
                end
            end
            if numel(dimNum) > 1 || any(plotDim ~= dimNum)
                % Not clear why this has to be updated here. Reintroduce if
                % necessary. However note that (eg) when projection is
                % selected and the cross hair in the image should be moved
                % - with the mouse - the projection will be recomputed
                % which might take some time!
%                 if ~imageUpdated
%                     updateCurrIm;
%                 end
                updatePlots;
            end
            updateObjects;
            if (get(tbProfile,'Value') && nProfiles)
               updateProfileProperties; 
            end
            if get(tbRoi,'Value') & nRois & rgbCount~=dimNum
              updateRoiSelection(get(roiListbox, 'Value')); % updateRoiProperties(0);
            end
        end
        if any(get(bg_colormap, 'SelectedObject') == [cmImage cmZoom cmThresh]) && ~get(tbSwitchRGB, 'Value')
            updateColormap;
        end
        if movdata.rec
            vidWriter('append');
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Function called if zoomProj Button is pressed (toggle projection modus
%from all data along projDim to data inside zoom interval)
    function  zoomProj(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if get(bt_zoomProj, 'Value')
            set(bt_zoomProj, 'CData', lockClosed);
        else
            set(bt_zoomProj, 'CData', lockOpen);
        end
        updateImages;
        if debugMatVis, debugMatVisFcn(2); end
    end

%Function for rgb map for stretched RGB mode
    function rgbV = rgbValue(col, exp)
        if col<0
            rgbV=[ 0 0 1];
        elseif col<1/4
            rgbV=[ 0 (col*4)^exp 1];
        elseif col<2/4
            rgbV=[ 0 1 (2-col*4)^exp];
        elseif col<3/4
            rgbV=[(col*4-2)^exp 1 0];
        elseif col<1
            rgbV=[1 (4-col*4)^exp 0];
        else
            rgbV=[1 0 0];
        end
    end
    function updateContrastSel(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        currContrastSel = get(popContrastSel, 'Value');
        set(sldMin, 'Value',cmMinMax(currContrastSel,1));
        set(sldMax, 'Value',cmMinMax(currContrastSel,2));
        set(sldGamma, 'Value',currGamma(currContrastSel));
        if isinteger(data{currContrastSel})   % (currContrastSel<=nMat && isinteger(data{currContrastSel})) || (currContrastSel>nMat && isinteger(alphaMap{currContrastSel-nMat}))
            cmMinMax = round(cmMinMax);
            set(valSldMin, 'String', num2str(cmMinMax(currContrastSel  ,1)));
            set(valSldMax, 'String', num2str(cmMinMax(currContrastSel  ,2)));
        else
            set(valSldMin, 'String', num2str(cmMinMax(currContrastSel  ,1),'%6.3f'));
            set(valSldMax, 'String', num2str(cmMinMax(currContrastSel  ,2),'%6.3f'));
        end
        if withAlpha
            set(edt_alphaMin, 'String', num2str(cmMinMax(nMat+currContrastSel,1),'%6.2f'));
            set(edt_alphaMax, 'String', num2str(cmMinMax(nMat+currContrastSel,2),'%6.2f'));
        end
        updateGuiHistVal;
        updateSldLim;
        if withAlpha
            updateAlphaContrast;
        end
        updateGamma(0,0,'sld','data');
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateLinkContrast(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        linkContrastSettings = get(tb_linkContrastAdjustments, 'Value');
        if linkContrastSettings
            set(tb_linkContrastAdjustments, 'CData', lockClosed);
            currGamma(1:nMat) = ones(1,nMat)*currGamma(currContrastSel);
            if withAlpha
                currGamma(nMat+1:2*nMat) = currGamma(currContrastSel+nMat);
                cmMinMax(nMat+1:2*nMat,:) = repmat(cmMinMax(currContrastSel+nMat,:),[nMat 1]);
                updateAlphaContrast;
            end
            updateImages;
            updateColormap;
        else
            set(tb_linkContrastAdjustments, 'CData', lockOpen);
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateFilter(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        % Update selection of filter type (callback from filter popup menu)
        sel = get(popFilter, 'String');
        switch sel{get(popFilter, 'Value')}
            case 'gauss'
                set(txtFilterPar, 'String', 'Sigma');
            case {'median', 'mean', 'unsharp-dip'}
                set(txtFilterPar, 'String', 'Size');
            case 'unsharp'
                set(txtFilterPar, 'String', 'Alpha');
            case 'none'
                set(txtFilterPar, 'String', 'Par');
        end
        updateImages;
        drawPlots;
        if any(get(bg_colormap, 'SelectedObject') == [cmImage cmZoom cmThresh]) && ~get(tbSwitchRGB, 'Value')
            updateColormap;
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function out = filterImage(in)
        if debugMatVis, debugMatVisFcn(1); end
        % Filter images (called from updateCurrIm, individually for normal,
        % projection and RGB mode)
        filtSz = zeros(1,nDim);
        % Evaluate string in Par-text field (entering of vectors allowed
        % for anistotropic filtering)
        fs = eval(get(editFilterPar, 'String')); % use eval to allow vector inputs (for anisotropic filtering)
        % Flag for use on whole data set if dimensions other than
        % x/y-dimensions are to be included
        wholeDataSet = 0;
        % If in RGB-mode only filter in x/y
        if rgbCount
            fs = fs(1:min(2,numel(fs)));
        end
        % Specifiy par-vector depending on input
        if numel(fs) == 1 % Single number: isotropic filtering in x/y
            filtSz = fs*[1 1];
        elseif numel(fs) == 2 % Two numbers: anisotropic filtering in x/y
            filtSz = fs([2 1]);
        else % More than two numbers: anisotropic filtering of whole data set
            filtSz(1:numel(fs)) = fs;
            wholeDataSet = 1;
            % Create index for current image
            for ii = 1:nDim
                imIndex{ii} = currPos(ii);  %#ok
            end
            imIndex{xySel(1)} = ':';
            imIndex{xySel(2)} = ':';
        end
        if any(filtSz)
            sel = get(popFilter, 'String');
            for ii = 1:nMat
                if withDipimage  % Use DIPImage functions
                    if wholeDataSet
                        in = data{ii};  % replace in matrix (i.e. curr image) by whole data set
                    end
                    switch sel{get(popFilter, 'Value')}
                        case 'none'
                            out = in;
                        case 'gauss'
                            out = double(gaussf(double(in),filtSz));
                        case 'median'
                            out = double(medif(double(in),filtSz,'elliptic'));
                        case 'mean'
                            out = double(unif(double(in),filtSz,'elliptic'    ));
                        case 'unsharp-ml'  % Image processing toolbox implementation
                            if wholeDataSet
                                wholeDataSet = 0;
                                in = squeeze(in(imIndex{:}));
                                if xySel(1) > xySel(2)
                                    in = in';
                                end
                            end
                            if exist('imsharpen', 'file')  % New improved version, starting from Matlab 2013a
                                if numel(filtSz == 2)
                                    out = imsharpen(in, 'Radius', filtSz(1), 'Amount', filtSz(2));
                                else
                                    out = imsharpen(in, 'Radius', filtSz(1));
                                end
                            else
                                h = fspecial('unsharp', min([1, filtSz(1)/10]));
                                out = imfilter(in, h, 'replicate');
                            end
                        case 'unsharp-dip'  % DIPImage implementation
                            out = double(double(in) - laplace(double(in), filtSz));
                    end
                    if wholeDataSet
                        out = squeeze(out(imIndex{:}));  % Extract filtered image of current location
                        if xySel(1) > xySel(2)  % Invert matrix if necessary
                            out = out';
                        end
                    end
                else  % Use Image Processing Toolbox
                    switch sel{get(popFilter, 'Value')}
                        case 'none'
                            out = in;
                        case 'gauss'
                            h = fspecial('gaussian', filtSz(1:2)*5, filtSz(1));
                            out = imfilter(in, h, 'replicate');
                        case 'median'
                            out = medfilt2(in, filtSz(1:2));
                        case 'mean'
                            h = fspecial('disk', filtSz(1));
                            out = imfilter(in, h, 'replicate');
                        case 'unsharp-ml'
                            h = fspecial('unsharp', min([1, filtSz(1)/10]));
                            out = imfilter(in, h, 'replicate');
                    end
                end
            end
        else
            out = in;
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Update current image values
    function updateCurrIm(varargin)   %#ok
        if debugMatVis, debugMatVisFcn(1); end
        if nargin > 0 && strcmp(varargin{1}, 'forceUpdateGuiHist')
            forceUpdateGuiHist = 1;
        else
            forceUpdateGuiHist = 0;
        end
        % Create index for current image
        for ii = 1:nDim
            imIndex{ii} = currPos(ii);  %#ok
        end
        imIndex{xySel(1)} = ':';
        imIndex{xySel(2)} = ':';
        if withAlpha
            for ii = 1:nMat
                cA = squeeze(alphaMap{ii}(imIndex{:}));
                if xySel(1) > xySel(2)
                    currAlphaMap{ii} = cA';
                else
                    currAlphaMap{ii} = cA;
                end
            end
        elseif get(tbSwitchRGB, 'Value')
            if (~get(cmStretchRGBMean, 'Value') && ~get(cmStretchRGBMax, 'Value'))
                imIndex{rgbDim} = mod((currPos(rgbDim)-2:currPos(rgbDim)), dim(rgbDim))+1;  
            else
              if get(tb_lockPlots2Zoom(rgbDim), 'Value')
                imIndex{rgbDim} = zoomVal(rgbDim,1):sum(zoomVal(rgbDim,:))-1;
              else
                imIndex{rgbDim} = ':';
              end
            end
        end
        if with2DHist
            updateAlphaContrast;
        end
        % No projection selected
        if projMethod == 0
                for ii = 1:nMat
                    c = squeeze(data{ii}(imIndex{:}));
                    if get(tbSwitchRGB, 'Value')
                      [s,ind] = sort([xySel,rgbDim]);
                      currIm{ii} = ipermute(c, ind);
                    else
                      if xySel(1) > xySel(2)
                        currIm{ii} = permute(c,[2 1 3]);
                      else
                        currIm{ii} = c;
                      end
                    end
                end
                %Projection of data if selected
            else
                busy(1);
                if get(bt_zoomProj, 'Value')
                  imIndex{projDim} = zoomVal(projDim,1):sum(zoomVal(projDim,:))-1;
                else
                  imIndex{projDim} = ':';
                end
                % Find number of dimension of xySel(1), xySel(2) and projDim with
                % respect to extracted 3D data volume
                if get(tbSwitchRGB, 'Value') == 0
                  xx  = find(xySel(1) == sort([xySel projDim]));
                  yy  = find(xySel(2) == sort([xySel projDim]));
                  p   = find(projDim == sort([xySel projDim]));
                else
                  xx  = find(xySel(1) == sort([xySel projDim rgbDim]));
                  yy  = find(xySel(2) == sort([xySel projDim rgbDim]));
                  p   = find(projDim  == sort([xySel projDim rgbDim]));
                  rgb = find(rgbDim   == sort([xySel projDim rgbDim]));
                end
                for ii = 1:nMat
                    % Sort dimension to [xySel(1) xySel(2) projDim]
                    if get(tbSwitchRGB, 'Value') == 0
                      c = squeeze(permute(squeeze(data{ii}(imIndex{:})),[xx yy p]));
                      if withAlpha
                        cA = squeeze(permute(squeeze(alphaMap{ii}(imIndex{:})),[xx yy p]));
                        sz1 = size(cA);
                      end
                    else
                      c = squeeze(permute(squeeze(data{ii}(imIndex{:})),[xx yy p rgb]));
                      sz1 = size(c);
                    end
                    switch projMethod
                        case 1      %maximum projection
                          if withAlpha
                            [currAlphaMap{ii},cAInd] = max(cA, [], 3); %squeeze() % used to be nanmax
                            cAInd  = (1:prod(sz1(1:2)))' + prod(sz1(1:2)) * (cAInd(:)-1);
                            currIm{ii}  = reshape(c(cAInd), sz1(1:2));
                          elseif get(tbSwitchRGB, 'Value')
                            [~,cInd] = max(sum(c,4), [], 3); %squeeze() % used to be nanmax
                            cInd  = (1:prod(sz1(1:2)))' + prod(sz1(1:2)) * (cInd(:)-1);
                            c = reshape(c, [prod(sz1(1:3)) sz1(4)]);
                            currIm{ii}  = reshape(c(cInd,:), sz1([1 2 4]));
                          else
                            currIm{ii} = max(c, [], 3); % used to be nanmax
                          end
                        case 2      %minimum projection
                          if withAlpha
                            warning(sprintf('MIN PROJECTION not understood in combination with alphaMap\nmin values of datamatrix is shown instead'))
                          end
                          if get(tbSwitchRGB, 'Value')
                            [~,cInd] = min(sum(c,4), [], 3); %squeeze() % used to be nanmax
                            cInd  = (1:prod(sz1(1:2)))' + prod(sz1(1:2)) * (cInd(:)-1);
                            c = reshape(c, [prod(sz1(1:3)) sz1(4)]);
                            currIm{ii}  = reshape(c(cInd,:), sz1([1 2 4]));
                          else
                            currIm{ii} = min(c, [], 3);    % Used to be nanmin
                          end
                        case 3      %mean projection
                          if withAlpha
                            currIm{ii}       = sum(c.*cA,3, 'omitnan')./sum(cA.*~isnan(c),3, 'omitnan'); % weighted mean ratio
                            currAlphaMap{ii} = sum(cA.^2,3, 'omitnan')./sum(cA,           3, 'omitnan'); % mean Intensity
                          else
                            currIm{ii} = squeeze(mean(c,  3, 'omitnan')); % used to be nanmean
                          end
                        case 4      %standard deviation projection
                          if withAlpha
                            % standard error:  SE = sum( (x - <x>)^2.*w, 3) ./ sum(w, 3)
                            currIm{ii}       = sqrt( sum( (c - repmat(sum(c.*cA,3, 'omitnan')./sum(cA.*~isnan(c),3, 'omitnan'),[1 1 sz1(3)]) ).^2 .* cA,3, 'omitnan')./sum(cA.*~isnan(c),3, 'omitnan'));  % standard error of currIm{ii} -> SEM ././sqrt(sum(~isnan(A_rr),3))
                            % standard error of mean:  SEM = sum( (x - <x>)^2.*w, 3) ./ sum(w, 3) ./ sz(3)
                            %currIm{ii}       = sqrt( sum( (c - repmat(sum(c.*cA,3, 'omitnan')./sum(cA.*~isnan(c),3, 'omitnan'),[1 1 sz1(3)]) ).^2 .* cA,3, 'omitnan')./sum(cA.*~isnan(c),3, 'omitnan'))./sqrt(sum(~isnan(c),3));
                            currAlphaMap{ii} = sum(cA.^2,3, 'omitnan')  ./sum(cA,3, 'omitnan');           % mean Intensity
                          elseif get(tbSwitchRGB, 'Value')
                            [~,cInd] = std(sum(c,4), [], 3, 'omitnan'); %squeeze() % used to be nanmax
                            cInd  = (1:prod(sz1(1:2)))' + prod(sz1(1:2)) * (cInd(:)-1);
                            c = reshape(c, [prod(sz1(1:3)) sz1(4)]);
                            currIm{ii}  = reshape(c(cInd,:), sz1([1 2 4]));
                          else
                            currIm{ii} = std(double(c), [], 3, 'omitnan');
                          end
                        case 5      %variance projection
                          if withAlpha
                            currIm{ii}       = sum( (c - repmat(sum(c.*cA,3, 'omitnan')./sum(cA.*~isnan(c),3, 'omitnan'),[1 1 sz1(3)]) ).^2 .* cA,3, 'omitnan')...
                              ./sum(cA.*~isnan(c),3, 'omitnan');  % standard error of currIm{ii} -> SEM ././sqrt(sum(~isnan(A_rr),3))
                            currAlphaMap{ii} = sum(cA.^2,3, 'omitnan')  ./sum(cA,3, 'omitnan');          % mean Intensity
                          elseif get(tbSwitchRGB, 'Value')
                            [~,cInd] = var(sum(c,4), [], 3, 'omitnan'); 
                            cInd  = (1:prod(sz1(1:2)))' + prod(sz1(1:2)) * (cInd(:)-1);
                            c = reshape(c, [prod(sz1(1:3)) sz1(4)]);
                            currIm{ii}  = reshape(c(cInd,:), sz1([1 2 4]));
                          else
                            currIm{ii} = var(double(c), [], 3, 'omitnan');
                          end
                        case 6    %Tile images ('montage')
                            nCols = max(1,round(sqrt(dim(xySel(1))/dim(xySel(2))*size(c,3)/aspRatio(2) * aspRatio(1))));
                            nRows = ceil(size(c,3) / nCols);
                            if nRows > 1 && nCols > 1
                                nBlack = nRows*nCols - size(c,3);
                            else
                                nBlack = 0;
                                if nRows == 1
                                    nCols = size(c,3);
                                else
                                    nRows = size(c,3);
                                end
                            end
                            if nBlack > 0
                                c(:,:,size(c,3)+1:nRows*nCols) = minVal(1);
                                if withAlpha
                                  cA(:,:,size(c,3)+1:nRows*nCols) = 0;
                                end
                            end
                            if islogical(c)
                                currIm{ii} = false(dim(xySel(1))*nRows, dim(xySel(2))*nCols);
                            else
                                currIm{ii} = ones(dim(xySel(1))*nRows, dim(xySel(2))*nCols, class(c));
                            end
                            if withAlpha
                              currAlphaMap{ii} = zeros(dim(xySel(1))*nRows, dim(xySel(2))*nCols, class(c));
                            end
                            for r=1:nRows
                                currIm{ii}((r-1)*dim(xySel(1))+1:r*dim(xySel(1)),:) = c(:,(r-1)*nCols*dim(xySel(2))+1:r*nCols*dim(xySel(2)));
                                if withAlpha
                                  currAlphaMap{ii}((r-1)*dim(xySel(1))+1:r*dim(xySel(1)),:) = cA(:,(r-1)*nCols*dim(xySel(2))+1:r*nCols*dim(xySel(2)));
                                end
                            end
                    end
                    clear c cA sz1;
                end
                busy(0);
            end
        %Non-RGB mode
        if get(tbSwitchRGB, 'Value') == 0
            currImVal = currIm; % Remember "original" values
            if withAlpha
                currAlphaMapVal = currAlphaMap;
            end
            % Apply filter
            if withFilter && get(popFilter, 'Value')>1
                for ii=1:nMat
                    currIm{ii} = filterImage(currIm{ii});
                end
            end
            %Adjust for gamma value
            for ii = 1:nMat
                if currGamma(ii) ~= 1
                    currIm{ii} = applyGamma(currIm{ii},cmMinMax(ii,:),currGamma(ii),1);
                end
                if withAlpha && currGamma(ii+nMat) ~= 1
                    currAlphaMap{ii} = applyGamma(currAlphaMap{ii},cmMinMax(ii+nMat,:),currGamma(ii+nMat),1);
                end
            end
            % Set colorbar in between contrast sliders
            if ~strcmp(get(histAxGui, 'UserData'),'gamma')
                if currGamma(currContrastSel) ~= 1
                    histCDataGamma = applyGamma(linspace(histXDataBin(1),histXDataBin(end),numel(histXDataBin)),cmMinMax(currContrastSel,:),currGamma(currContrastSel),1);
                else
                    histCDataGamma = linspace(histXDataBin(1),histXDataBin(end),numel(histXDataBin));
                end
                set(contrastSldIm, 'CData',histCDataGamma);
            end
            % Set correct y-ticks in colorbar in image window
            if get(tbColorbar, 'Value') == 1
                for ii=1:nMat
                    cb = findobj(imageWin(ii), 'Tag','Colorbar');
                    ytl = applyGamma([cmMinMax(ii,1) get(cb,'YTick') cmMinMax(ii,2)],cmMinMax(ii,:),currGamma(ii),-1);
                    set(cb, 'YTickLabel', ytl(2:end-1));
                end
            end
        else
            % RGB mode
            % Update slider for RGB mode (update otherwise by updateColormap function)
            if get(cmStretchRGBMean, 'Value') || get(cmStretchRGBMax, 'Value')
                busy(1);
            end
            if nargin == 3
                set(sldMin, 'Value', max(get(sldMin,'Min'),str2double(get(valSldMin, 'String'))));
                set(sldMax, 'Value', min(get(sldMax,'Max'),str2double(get(valSldMax, 'String'))));
                set(sldMin_RGB, 'Value', str2double(get(valSldMin_RGB, 'String')));
                set(sldMax_RGB, 'Value', str2double(get(valSldMax_RGB, 'String')));
            end
            if isinteger(data{currContrastSel})
                cmMinMax(currContrastSel  ,:) = round([get(sldMin, 'Value') get(sldMax, 'Value')]);
            else
                cmMinMax(currContrastSel  ,:) = [get(sldMin, 'Value') get(sldMax, 'Value')];
            end
            if ~(cmMinMax(currContrastSel  ,1) < cmMinMax(currContrastSel  ,2))
                w = get(sldMax, 'SliderStep');
                cmMinMax(currContrastSel  ,2) = cmMinMax(currContrastSel  ,1) + w(1)*(maxVal(1)-minVal(1));
                set(sldMin, 'Value', cmMinMax(currContrastSel  ,1));
                set(sldMax, 'Value', cmMinMax(currContrastSel  ,2));
            end
            if isinteger(data{currContrastSel})
                cmMinMax = round(cmMinMax);
                set(valSldMin, 'String', num2str(cmMinMax(currContrastSel  ,1)));
                set(valSldMax, 'String', num2str(cmMinMax(currContrastSel  ,2)));
            else
                set(valSldMin, 'String', num2str(cmMinMax(currContrastSel  ,1),'%6.3f'));
                set(valSldMax, 'String', num2str(cmMinMax(currContrastSel  ,2),'%6.3f'));
            end
            for ii = 2:nMat
                cmMinMax(ii,:) = double(cmMinMax(currContrastSel  ,:))/double(maxVal(1))*double(maxVal(ii));
                cmMinMax(isnan(cmMinMax)) = 0;
            end
            if (~get(cmStretchRGBMean, 'Value') && ~get(cmStretchRGBMax, 'Value'))
                %"Normal" RGB mode (Global, Channel or Image)
                for ii = 1:nMat
                    sz = size(currIm{ii});
                    if dim(rgbDim) == 2
                        currIm{ii}(:,:,3) = min(currIm{ii}(:));
                        colDim = 2;
                    else
                        colDim = 3;
                    end
                    [s,ind] = sort([xySel,rgbDim]);
                    %currIm{ii} = ipermute(currIm{ii}, ind);
                    % Apply filter
                    if withFilter && get(popFilter, 'Value')>1
                      for jj = 1:size(currIm{ii},3)
                        currIm{ii}(:,:,jj) = filterImage(currIm{ii}(:,:,jj));
                      end
                    end
                    currImVal = currIm; % Remember "original" values
                    switch get(bg_colormap, 'SelectedObject')
                        %Adjust RGB channels according to selected display mode
                        case cmGlobal       %Global (all data points)
                            % Get values for histogram
                            if (forceUpdateGuiHist || updateGuiHistState) && ii==1  % Disable update of histogram when necessary
                              X = single(currIm{ii}(:,:,1:min(dim(rgbDim),sz(3))));
                              histValCurrIm = [];
                              for nX = 1: size(X,3)
                                histValCurrIm(nX,:) = histcounts(X(:,:,nX),histXData);
                              end
                            end
                        case cmImage        %Image (each of the three images individual).  Before: Channel (by selected rgb-dimension) (see below)
                            currIm{ii} = double(currIm{ii});
                            for kk = 1:colDim
                                cc = currIm{ii}(:,:,kk);
                                % Scale each channel independently to
                                % [minVal(ii) maxVal(ii)]
                                minCC = min(cc(:)); maxCC = max(cc(:));
                                if ~(minCC == maxCC)
                                    %cc = (cc-min(cc(:))) / (max(cc(:))-min(cc(:))) * double(maxVal(ii)-minVal(ii)) + double(minVal(ii));
                                    cc = scaleMinMax(double(cc), minVal(ii), maxVal(ii), minCC,maxCC);
                                end
                                currIm{ii}(:,:,kk) = cc;
                                % Get values for histogram
                            end
                            % Get values for histogram
                            if ii==1 %(forceUpdateGuiHist || updateGuiHistState) && ii==1
                              X = single(currIm{ii}(:,:,1:min(dim(rgbDim),sz(3))));
                              histValCurrIm = [];
                              for nX = 1: size(X,3)
                                histValCurrIm(nX,:) = histcounts(X(:,:,nX),histXData);
                              end
                            end
                            if colDim == 2      %otherwise channel 3 remains constant
                                currIm{ii}(:,:,3) = 0;
                            end
                    end
                    %Adjust for gamma value
                    currImVal = currIm; % Remember "original" values
                    if currGamma(ii) ~= 1
                        for jj = 1:min(size(currIm{ii},3),dim(rgbDim))   
                            currIm{ii}(:,:,jj) = applyGamma(currIm{ii}(:,:,jj),cmMinMax(ii,:),currGamma(ii),1);
                        end
                    end
                    if withAlpha
                        currAlphaMapVal = currAlphaMap;
                    end
                    % Not really clear how / if to implement alpha maps for
                    % RGB images. For now, use alpha map of current green channel
                    if withAlpha && currGamma(ii+nMat) ~= 1
                        currAlphaMap{ii} = applyGamma(currAlphaMap{ii},cmMinMax(ii+nMat,:),currGamma(ii+nMat),1);
                    end
                    % Set colorbar in between contrast sliders
                    if strcmp(get(histAxGui, 'UserData'),'gamma')
                        if currGamma(currContrastSel)
                            histCDataGamma = applyGamma(linspace(histXDataBin(1),histXDataBin(end),numel(histXDataBin)),cmMinMax(currContrastSel,:),currGamma(currContrastSel),1);
                        else
                            histCDataGamma = linspace(histXDataBin(1),histXDataBin(end),numel(histXDataBin));
                        end
                        set(contrastSldIm, 'CData',histCDataGamma);
                    end
                    % Set correct y-ticks in colorbar in image window
                    if get(tbColorbar, 'Value') == 1
                        cb = findobj(imageWin(ii), 'Tag','Colorbar');
                        ytl = applyGamma([cmMinMax(ii,1) get(cb,'YTick') cmMinMax(ii,2)],cmMinMax(ii,:),currGamma(ii),-1);
                        set(cb, 'YTickLabel', ytl(2:end-1));
                    end
                    % Apply contrast adjustment
                    min_currStack_rel = (single(cmMinMax(ii,1))-minVal(ii)) / (maxVal(ii) - minVal(ii)); % contrast minimum scaled to interval [0,1]
                    max_currStack_rel = (single(cmMinMax(ii,2))-minVal(ii)) / (maxVal(ii) - minVal(ii)); % contrast maximum scaled to interval [0,1]
                    currIm{ii} = rgbContrastAdjust(single(currIm{ii}),min_currStack_rel, max_currStack_rel);
                end
                if get(tbInvert, 'Value')
                    currIm{ii} = 1 - currIm{ii};
                end
                if get(tbFlip, 'Value') 
                  currIm{ii} = flipdim(currIm{ii},3); 
                end
                %                 if (forceUpdateGuiHist || updateGuiHistState)
                updateGuiHist(histValCurrIm);
                %                 end
            else
                %Stretch RGB mode
%                 rgbStretchSldVal = [get(sldMin_RGB, 'Value') get(sldMax_RGB, 'Value')];
%                 set(valSldMin_RGB, 'String', num2str(rgbStretchSldVal(1),'%6.3f'));
%                 set(valSldMax_RGB, 'String', num2str(rgbStretchSldVal(2),'%6.3f'));
                %stackIndex = imIndex;
                %stackIndex{rgbDim} = ':';
                for ii = 1:nMat
                    currStack = currIm{ii};
                    minVal(ii) = min(currStack(:));
                    maxVal(ii) = max(currStack(:));
                    sz = size(currStack);
                    if length(sz) == 2
                        sz(3) = 1;
                    end
                    if (forceUpdateGuiHist || updateGuiHistState) && ii==1 
                      if get(cb_showHist,'Value')
                        histValCurrIm = [];
                        for nX = 1: size(currStack,3)
                          histValCurrIm(nX,:) = histcounts(single(currStack(:,:,nX)),histXData);
                        end
                        %histValCurrIm = histcounts(single(reshape(currStack,[prod(sz(1:2)) sz(3)])),histXData)'; %
                      else
                        histValCurrIm = ones(sz(3), numel(histXDataBin));
                      end
                    end
                    rgbVal = [];
                    for kk = 1:sz(3)
                        rgbVal(:,kk)=rgbValue(get(sldMin_RGB,'Value')+(get(sldMax_RGB,'Value')-get(sldMin_RGB,'Value'))*(kk-1)/(sz(3)-1),1);
                    end
                    rgbValPlots = rgbVal;
                    currStack = reshape(double(currStack),[sz(1)*sz(2) sz(3)]);
                    % Mean RGB Stretch
                    if get(cmStretchRGBMean, 'Value')
                        if sz(3) > 2
                            % More than two "color-channels"
                            s = sum(rgbVal,2);
                            rgbVal = rgbVal / max(s);  % rgbVal = rgbVal./repmat(s,[1 sz(3)]);
                            currStack = currStack * rgbVal';
                        else % Stretch 2 channels across complete RGB-rainbow (instead of using only blends of red and blue)
                            % Averaged and normalized values used for
                            % brightness of image:
                            d_brightness = sum(currStack,2, 'omitnan');
                            d_brightness = d_brightness./max(d_brightness(:));
                            % RGB color balance
                            currStack = currStack - min(currStack(:));
                            d_color = currStack(:,1)./sum(currStack,2, 'omitnan'); % Weighted color balance
                            d_color(isnan(d_color)) = 0;                % avoid NaNs
                            rgbVal2 = zeros(256,3);
                            for kk = 1:256
                                rgbVal2(kk,:)=rgbValue(get(sldMin_RGB,'Value')+(get(sldMax_RGB,'Value')-get(sldMin_RGB,'Value'))*(kk-1)/255,1);
                            end
                            d_color = round(255*d_color)+1;     % index for rgb_val
                            d_color = [d_color+512; d_color+256; d_color];
                            currStack = rgbVal2(d_color(:));
                            currStack = currStack.*repmat(d_brightness,[3 1]);
                        end
                    else % Max RGB Stretch
                        [currStack ind] = max(currStack,[],2);
                        currStack = repmat(double(currStack(:)),[1 3]) .* rgbVal(:,ind(:))';
                    end
                    
                    % Apply contrast adjustment
                    min_currStack_rel = (double(cmMinMax(ii,1))-minVal(ii)) / (maxVal(ii) - minVal(ii)); % contrast minimum scaled to interval [0,1]
                    max_currStack_rel = (double(cmMinMax(ii,2))-minVal(ii)) / (maxVal(ii) - minVal(ii)); % contrast maximum scaled to interval [0,1]
                    %                     currIm{ii} = reshape(currStack,[sz(1:2) 3]);
                    currIm{ii} = reshape(currStack,[sz(1:2) 3]);
                    % Apply filter
                    if withFilter && get(popFilter, 'Value')>1 
                      for jj = 1:size(currIm{ii},3)
                        currIm{ii}(:,:,jj) = filterImage(currIm{ii}(:,:,jj));
                      end
                    end
                    currIm{ii} = rgbContrastAdjust(currIm{ii},min_currStack_rel, max_currStack_rel);
                end
                if get(tbInvert, 'Value')
                    currIm{ii} = 1 - currIm{ii};
                end
                if max(histValCurrIm(:))>1
                  updateGuiHist(histValCurrIm);
                end
                set(contrastAx, 'CLim',[cmMinMax(currContrastSel,1) cmMinMax(currContrastSel,2)]);
                currImVal = currIm;
            end
            busy(0);
        end
        if with2DHist
            update2DHist;
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function out = rgbContrastAdjust(in, relMin, relMax)
        if debugMatVis, debugMatVisFcn(1); end
        % Contrast adjustment: simply scaling and cutting the
        % values larger than 1 can generate shifts in hue when
        % adjusting the contrast. Thus, the values are scaled
        % to the maximum of any channel, thereby preserving hue
        % and only changing brightness:
        sz = size(in);
        % Reshape and avoid negative values by shifting data
        in = single(reshape(in-min(in(:)), [], sz(3)));
        % Determine the brightness of each pixel, given by its norm
        brightness = single(sqrt(sum(in.^2,2, 'omitnan')));
        % Determine max vector length for the RGB triple of each pixel
        % (e.g. pixels with a hue along [1 1 1] can have a larger RGB
        % vector than pixels along the RGB triple [1 0 0])
        maxBright = single(in ./ repmat(max(in,[],2), [1 sz(3)]));
        maxBright(max(in,[],2)==0,:) = 0; % sets NaNs in maxBright obtained from division by zero max(in,[],2) zo zero
        maxBright = single(sqrt(sum(maxBright.^2,2)));
        % Determine the hue of each pixel, given by its normalized RGB
        % vector
        hue = single(in ./ repmat(brightness, [1 sz(3)]));
        hue(brightness==0,:) = 0; % sets NaNs in hue obtained from division by zero brightness zo zero
        clear in;
        % Normalize intensities according to max intensities
        brightness = single(brightness ./ maxBright);
        brightness(isnan(brightness))=0;
        % Apply contrast adjustments by linearly scaling the brightness
        % values according to the slider settings
        brightness = single(scaleMinMax(brightness)); %, 0,1, relMin, relMax %Scale to interval [0,1] so that the relativ min/max values can be used
        brightness = single((brightness - relMin) / (relMax - relMin)); % Scale to relative min/max values
        brightness(brightness < 0) = 0;  % cut off values below black point
        brightness(brightness > 1) = 1;  % cut off values above white point
        brightness = single(brightness .* maxBright); % Rescale to max. intensities for RGB triples
        clear maxBright
        out = single(hue .* repmat(brightness, [1 sz(3)]));
        clear hue brightness
        out(out>1) = 1;
        out = reshape(out, sz);
        if debugMatVis, debugMatVisFcn(2); end
    end
    function out = applyGamma(in,contrastMinMax, gammaVal, gammaSign)
        if debugMatVis, debugMatVisFcn(1); end
        in = double(in);
        % Stretch data inside the contrast range to the interval [0 1]
        out = scaleMinMax(in,0, 1, contrastMinMax(1), contrastMinMax(2), 1);
        % Apply gamma
        out = out.^(gammaVal^gammaSign);
        % Stretch to original range
        out = scaleMinMax(out,contrastMinMax(1), contrastMinMax(2), 0, 1, 1);
        if debugMatVis, debugMatVisFcn(2); end
    end

% Update values for small histogram in Main Gui.
    function updateGuiHistVal(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if rgbCount == 0
            if strcmp(get(histAxGui, 'USerData'),'gamma')
                guiHistVal = histcounts(single(currIm{currContrastSel}(:)),histXData);
            else
                guiHistVal = histcounts(single(currImVal{currContrastSel}(:)),histXData);
            end
            set(histAxPlot, 'YData', guiHistVal, 'XData',histXDataBin);
            updateGuiHist(guiHistVal);
        else
            % For RGB modes, histogram values are determined in updateCurrIm
            % and passed to updateCurrIm as input h
            if nargin == 3 && strcmp(varargin{3}, 'forceUpdateGuiHist')
                updateCurrIm(varargin{3});
            else
                updateCurrIm; % Calls updateGuiHist
            end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateGuiHist(h,varargin)
        if debugMatVis, debugMatVisFcn(1); end
        % if get(bt_zoomProj, 'Value') size(h,1) might be smaller than length(histAxPlot)
        n = size(h,1);
        if n >1
          pCol = plotColors(n);
          if rgbCount && ~any(get(bg_colormap, 'SelectedObject') ==  [cmStretchRGBMean cmStretchRGBMax])
            % Channel RGB-mode requires switching of RGB
            pCol = flipdim(pCol,2);
          end
          if get(tbFlip, 'Value') % flip colormap
            pCol = flipdim(pCol,2);
          end
          for ii = 1:n
            set(histAxPlot(ii), 'YData', h(ii,:), 'XData',histXDataBin,'Color', pCol(ii,:),'Visible','on');
          end
        else
          set(histAxPlot, 'YData', h, 'XData',histXDataBin,'Visible','on');
        end
        for ii = n+1:length(histAxPlot)
          set(histAxPlot(ii), 'YData', ones(1,length(histXDataBin)), 'XData',histXDataBin,'Visible','off');
        end

        y0 = strcmp(get(histAxGui, 'YScale'),'log');
        try
            set(histAxGui, 'YLim', [0 1.05*max(h(:))]);
            set(histAxBg, 'Position', [cmMinMax(currContrastSel,1) y0 cmMinMax(currContrastSel,2)-cmMinMax(currContrastSel,1) 1.1*max(max(h(:,2:end-1)))]);  % [cmMinMax(currContrastSel,1) 1 cmMinMax(currContrastSel,2)-cmMinMax(currContrastSel,1) max(max(h(:,2:end-1)))]
        catch    %#ok
            set(histAxGui, 'YLim', [0 1.05*max(h(:))]);
            if cmMinMax(currContrastSel,1) ~= cmMinMax(currContrastSel,2)
                set(histAxBg, 'Position', [cmMinMax(currContrastSel,1) y0 cmMinMax(currContrastSel,2)-cmMinMax(currContrastSel,1) 1.1*max(h(:))]);
            end
        end
        if cmMinMax(currContrastSel,1) ~= cmMinMax(currContrastSel,2)
            set(contrastAx, 'CLim',[cmMinMax(currContrastSel,1) cmMinMax(currContrastSel,2)]);
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function state = updateGuiHistState
        state = (isPlaying && get(tb_playHist, 'Value')) || (dragging && get(tb_moveHist,'Value')) || (~isPlaying && ~dragging);
    end
% Redraw small histogram in Main Gui. Necessary when toggling between
% RGB / non-RGB modes and when switching between lin/log in RGB mode
% (see buttonDownGuiHist for explanation). In any case the function
% that calls drawGuiHist also calls (indirectly) updateGuiHist, so the
% correct YData values don't have to be set here.
    function drawGuiHist(n)
        if debugMatVis, debugMatVisFcn(1); end
        delete(histAxPlot);
        delete(findobj(histAxGui,'Type','text','userdata','histText'))
        histAxPlot = [];
        if rgbCount && any(get(bg_colormap, 'SelectedObject') ==  [cmStretchRGBMean cmStretchRGBMax])
            set(cb_showHist,'Visible','on')
            hold(histAxGui,'on');
            pCol = plotColors(n);
            for ii=1:n
                histAxPlot(ii) = plot(histXDataBin,zeros(1,numel(histXDataBin)),'Parent',histAxGui,'Color', pCol(ii,:));
            end
            set(histAxPlot, 'ButtonDownFcn', @buttonDownGuiHist);
            
%             text(1,1,1,'No histogram available.','parent',histAxGui,'units','normalized','Position',[.5 .5 0],...
%                 'HorizontalAlignment','center','Userdata','histText');
        else
            set(cb_showHist,'Visible','off')
            hold(histAxGui,'on');
            if n == 1
              ii = 1;
              histAxPlot(ii) = bar(histXDataBin,zeros(1,numel(histXDataBin)),'Parent',histAxGui,...
                'FaceColor', round(n/3) * [ii==1 ii==2 ii==3],...
                'EdgeColor', round(n/3) * [ii==1 ii==2 ii==3],...
                'BarWidth',1,'BaseValue',0.9*strcmp(get(histAxGui, 'YScale'), 'log'));
                set(get(histAxPlot,'BaseLine'),'LineStyle','none');
            else
              pCol = [1 0 0;0 1 0;0 0 1];
              for ii=1:n
                histAxPlot(ii) = plot(histXDataBin,ones(1,numel(histXDataBin)),'Parent',histAxGui,'Color', pCol(ii,:));
                %histAxPlot(ii) = bar(histXDataBin,zeros(1,numel(histXDataBin)),'Parent',histAxGui,'FaceColor', round(n/3) * [ii==1 ii==2 ii==3],'EdgeColor',round(n/3) * [ii==1 ii==2 ii==3],'BarWidth',1,'BaseValue',0.9*strcmp(get(histAxGui, 'YScale'), 'log'));
              end
                %hap_bl = get(histAxPlot,'BaseLine');
                %set([hap_bl{:}],'LineStyle','none');
            end
            set(histAxPlot, 'ButtonDownFcn', @buttonDownGuiHist);
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Draw Image and Zoom window axes
    function drawImages(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if get(tbWin(1), 'Value') == 1 ||   get(tbWin(2), 'Value') == 1    %Image window(s)
            updateCurrIm;
            if get(tbWin(1), 'Value') == 1
                for ii = 1:nMat
                    set(0,'CurrentFigure',imageWin(ii));
                    set(imageWin(ii), 'HandleVisibility', 'on');
                    if customDimScale
                        imHandle(ii) = imagesc(dimScale(xySel(2),:),dimScale(xySel(1),:),currIm{ii}, cmMinMax(ii,:));
                    else
                        imHandle(ii) = imagesc(currIm{ii}, cmMinMax(ii,:));
                    end
                    imAx(ii) = myGca;
                    set(imageWin(ii), 'HandleVisibility', 'off');
                    if ~withDimUnits
                        xlabel(imAx(ii), [dimNames{xySel(2)}]);
                        ylabel(imAx(ii), [dimNames{xySel(1)}]);
                    else
                        xlabel(imAx(ii), [dimNames{xySel(2)} ' [' dimUnits{xySel(2)} ']']);
                        ylabel(imAx(ii), [dimNames{xySel(1)} ' [' dimUnits{xySel(1)} ']']);
                    end
                    if withAlpha
                        set(imHandle(ii), 'AlphaData', currAlphaMap{ii});
                        set(imAx(ii), 'XColor','w','YColor','w');
                    end
                end
                if customDimScale
                    axis(imAx, [dimScale(xySel(2),1) dimScale(xySel(2),2), dimScale(xySel(1),1), dimScale(xySel(1),2)]);
                end
                showColorbar;
            end
            if get(tbWin(2), 'Value') == 1          %Zoom window(s)
                for ii = 1:nMat
                    set(0,'CurrentFigure',zoomWin(ii));
                    set(zoomWin(ii), 'HandleVisibility', 'on');
                    if customDimScale
                        zoomHandle(ii) = imagesc(dimScale(xySel(2),:),dimScale(xySel(1),:), currIm{ii}, cmMinMax(ii,:));
                    else
                        zoomHandle(ii) = imagesc(currIm{ii}, cmMinMax(ii,:));
                    end
                    zoomAx(ii) = myGca;
                    %                     colormap(zoomAx(ii), cmap);
                    set(zoomWin(ii), 'HandleVisibility', 'off');
                    set(myGca, 'Position', [0 0 1 1]);
                    if withAlpha
                        set(zoomHandle(ii), 'AlphaData', currAlphaMap{ii});
                    end
                end
                if customDimScale
                    axis(zoomAx, [dimScale(xySel(2),1) dimScale(xySel(2),2), dimScale(xySel(1),1), dimScale(xySel(1),2)]);
                else
                    axis(zoomAx, [zoomValXY(1)-0.5, zoomValXY(1)+zoomValXY(3)-.5, zoomValXY(2)-0.5, zoomValXY(2)+zoomValXY(4)-.5]);
                end
            end
            if withAlpha
                set([imAx zoomAx], 'Color', 'k');
                set([imHandle zoomHandle],'AlphaDataMapping','scaled');
                updateAlphaContrast;
            end
            if any(get(bg_colormap, 'SelectedObject') == [cmImage cmThresh])
                updateColormap;
            else
                for ii=1:nMat
                    colormap(imAx(ii), cmap);
                    colormap(zoomAx(ii), cmap);
                end
            end
        end
        if strcmp(get(bt_XDir,'String'),'<')
            set([zoomAx, imAx], 'XDir','reverse');
        end
        if strcmp(get(bt_YDir,'String'),'^')
            set([zoomAx, imAx], 'YDir','normal');
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function setAxisDir(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        % Set directoin of image/zoom axis (Callback of buttons bt_XDir and
        % bt_YDir)
        switch varargin{3}
            case 'x'
                switch get(zoomAx(1), 'XDir')
                    case 'normal'
                        set(bt_XDir,'String','<');
                        set([zoomAx imAx], 'XDir','reverse');
                    case 'reverse'
                        set(bt_XDir,'String','>');
                        set([zoomAx imAx], 'XDir','normal');
                end
            case 'y'
                switch get(zoomAx(1), 'YDir')
                    case 'normal'
                        set(bt_YDir,'String','v');
                        set([zoomAx imAx], 'YDir','reverse');
                    case 'reverse'
                        set(bt_YDir,'String','^');
                        set([zoomAx imAx], 'YDir','normal');
                end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateRgbStretchPar(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        switch varargin{1}
            case valSldMin_RGB,
                rgbStretchSldVal(1) = str2double(get(valSldMin_RGB, 'String'));
                set(sldMin_RGB, 'Value', rgbStretchSldVal(1));
            case valSldMax_RGB
                rgbStretchSldVal(2) = str2double(get(valSldMax_RGB, 'String'));
                set(sldMax_RGB, 'Value', rgbStretchSldVal(2));
            case sldMin_RGB
                rgbStretchSldVal(1) = get(sldMin_RGB, 'Value');
            case sldMax_RGB
                rgbStretchSldVal(2) = get(sldMax_RGB, 'Value');
        end 
        set(valSldMin_RGB, 'String', num2str(rgbStretchSldVal(1),'%6.3f'));
        set(valSldMax_RGB, 'String', num2str(rgbStretchSldVal(2),'%6.3f'));
        updateImages;
        if debugMatVis, debugMatVisFcn(2); end
    end
%Update axes in Image and zoom window
    function updateImages(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        %if nargin; dimNum = varargin{1}; end
        currGamma(currContrastSel) = get(sldGamma(1), 'Value');
        if (~rgbCount || (rgbCount && ~updateGuiHistState)) && (get(tbWin(1), 'Value') == 1 || get(tbWin(2), 'Value') == 1 || get(tbHist,'Value'))  % in RGB mode, updateCurrIm is called from updateGuiHist
            updateCurrIm(varargin{:});
        end
        if updateGuiHistState %&& ~rgbCount % AZ: added '&& ~rgbCount' because 'updateGuiHistVal' was called twice in RGB mode, but this was wrong
            updateGuiHistVal;
        end
        if get(tbWin(1), 'Value') == 1
            for ii = 1:nMat
                set(imHandle(ii), 'CData', currIm{ii});
                if withAlpha
                    set(imHandle(ii), 'AlphaData', currAlphaMap{ii});
                end
            end
        end
        if get(tbWin(2), 'Value') == 1
            for ii = 1:nMat
                set(zoomHandle(ii), 'CData', currIm{ii});
                if withAlpha
                    set(zoomHandle(ii), 'AlphaData', currAlphaMap{ii});
                end
            end
        end
        % If called from filter items and image or zoom contrast mode
        % selected, colormap has to be updated
        if withFilter && nargin && isnumeric(varargin{1}) && any(varargin{1} == [editFilterPar popFilter]) && any(get(bg_colormap, 'SelectedObject') == [cmImage cmZoom])
            updateColormap;
        end
        %         if
        if get(tbRoi, 'Value')
          % after optimizing ROI properties I switched
          % "updateRoiProperties(1)" off, since it was called twice
%          updateRoiProperties(1); % I do not know, which dim was called + changed
%           if ~isempty(plotDim) && exist('dimNum','var')
%             updateRoiProperties(length(plotDim)>1 || any(plotDim ~= dimNum)); % Update properties of ROIs, with plot updates only if necessary
%           else
%             updateRoiProperties(0);
%           end
        end

        if debugMatVis, debugMatVisFcn(2); end
    end

%Draw Objects (Position lines, zoom rectangle)
    function drawObjects(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if get(tbWin(1), 'Value') == 1 || projMethod == 6
          if oldMATLAB
            rectPlotArea_Im = []; %#ok
          else
            rectPlotArea_Im = matlab.graphics.primitive.Rectangle.empty;
          end
          for ii = 1:nMat
            set(0,'CurrentFigure',imageWin(ii)); %axes(imAx);
            set(imageWin(ii), 'HandleVisibility', 'on');
            if customDimScale
              zoomReg(ii) = rectangle('Position',zoomPxXY-[pxWidth(xySel([2 1])) 0 0]/2, 'EdgeColor', 'r');
              currLoc = scaledLocation(currPos(xySel));
              lineHorIm(ii) = line(currLoc(2)*[1 1],get(myGca, 'YLim'), 'Color', 'w');
              lineVertIm(ii) = line(get(myGca, 'XLim'), currLoc(1)*[1 1], 'Color', 'w');
            else
              zoomReg(ii) = rectangle('Position',zoomValXY-[.5 .5 0 0], 'EdgeColor', 'r');
              lineHorIm(ii) = line([currPos(xySel(2)) currPos(xySel(2))],get(myGca, 'YLim'), 'Color', 'w');
              lineVertIm(ii) = line(get(myGca, 'XLim'), [currPos(xySel(1)) currPos(xySel(1))], 'Color', 'w');
            end
            for jj=3:-1:1
              if customDimScale
                rectPlotArea_Im(ii,jj) = rectangle('Position',[currLoc(2)-jj*pxWidth(xySel(2)) currLoc(1)-jj*pxWidth(xySel(1)) 2*jj*pxWidth(xySel(2)) 2*jj*pxWidth(xySel(1))],'EdgeColor',[jj==1 jj==2 jj==3],'Parent',imAx(ii),'Visible','off');
              else
                rectPlotArea_Im(ii,jj) = rectangle('Position',[currPos(xySel(2))-jj currPos(xySel(1))-jj currPos(xySel(2))+jj 2*jj],'EdgeColor',[jj==1 jj==2 jj==3],'Parent',imAx(ii),'Visible','off');
              end
            end
            set(imageWin(ii), 'HandleVisibility', 'off');
          end
        end
        if get(tbWin(2), 'Value') == 1
            if oldMATLAB
                rectPlotArea_Zoom = []; %#ok
            else
                rectPlotArea_Zoom = matlab.graphics.primitive.Rectangle.empty;
            end
            for ii = 1:nMat
                set(0,'CurrentFigure',zoomWin(ii)); %axes(zoomAx);
                set(zoomWin(ii), 'HandleVisibility', 'on');
                if customDimScale
                    currLoc = scaledLocation(currPos(xySel));
                    lineHorZoom(ii) = line(currLoc(2)*[1 1],get(myGca, 'YLim'), 'Color', 'w');
                    lineVertZoom(ii) = line(get(myGca, 'XLim'), currLoc(1)*[1 1], 'Color', 'w');
                else
                    lineHorZoom(ii) = line([currPos(xySel(2)) currPos(xySel(2))],get(myGca, 'YLim'), 'Color', 'w');
                    lineVertZoom(ii) = line(get(myGca, 'XLim'), [currPos(xySel(1)) currPos(xySel(1))], 'Color', 'w');
                end
                for jj=3:-1:1
                    if customDimScale
                        rectPlotArea_Zoom(ii,jj) = rectangle('Position',[currLoc(2)-jj*pxWidth(xySel(2)) currLoc(1)-jj*pxWidth(xySel(1)) 2*pxWidth(xySel(2)) 2*jj*pxWidth(xySel(1))],'EdgeColor',[jj==1 jj==2 jj==3],'Parent',zoomAx(ii),'Visible','off');
                    else
                        rectPlotArea_Zoom(ii,jj) = rectangle('Position',[currPos(xySel(2))-jj currPos(xySel(1))-jj currPos(xySel(2))+jj 2*jj],'EdgeColor',[jj==1 jj==2 jj==3],'Parent',zoomAx(ii),'Visible','off');
                    end
                end
                set(zoomWin(ii), 'HandleVisibility', 'off');
            end
        end
        updateObjects;
        toggleShowObjects;
        updateAspRatio;
        if debugMatVis, debugMatVisFcn(2); end
    end

%Update Objects (i.e. position and zoom indicators)
    function updateObjects(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if get(tbWin(1), 'Value') || get(tbWin(2), 'Value')
            meanPlotVal = 2*mod(get(btMean, 'UserData'),5)+1;
            if customDimScale
                set(zoomReg, 'Position',zoomPxXY-[pxWidth(xySel([2 1])) 0 0]/2);
                currLoc = scaledLocation(currPos(xySel([2 1])));
                set(lineHorIm, 'XData',currLoc(1)*[1 1],'YData', get(imAx(1), 'YLim'));
                set(lineVertIm, 'XData', get(imAx(1), 'XLim'), 'YData', currLoc(2)*[1 1]);
                set(lineHorZoom, 'XData',currLoc(1)*[1 1], 'YData',  get(zoomAx(1), 'YLim'));
                set(lineVertZoom, 'XData', get(zoomAx(1), 'XLim'), 'YData', currLoc(2)*[1 1]);
                if meanPlotVal < 7
                    set([rectPlotArea_Zoom rectPlotArea_Im], 'Position',[currLoc(1)-meanPlotVal*pxWidth(xySel(2))/2 currLoc(2)-meanPlotVal*pxWidth(xySel(1))/2 meanPlotVal*pxWidth(xySel(2)) meanPlotVal*pxWidth(xySel(1))]);
                elseif meanPlotVal == 7
                    for ii = 1:3
                        set([rectPlotArea_Zoom(:,ii) rectPlotArea_Im(:,ii)], 'Position',[currLoc(1)-(2*ii-1)*pxWidth(xySel(2))/2 currLoc(2)-(2*ii-1)*pxWidth(xySel(1))/2 (2*ii-1)*pxWidth(xySel(2)) (2*ii-1)*pxWidth(xySel(1))]);
                    end
                else
                    set([rectPlotArea_Zoom rectPlotArea_Im],'Position',zoomValXY-[.5 .5 0 0]);
                end
            else
                set(lineHorIm, 'XData', [currPos(xySel(2)) currPos(xySel(2))],'YData', get(imAx(1), 'YLim'));
                set(lineVertIm, 'XData', get(imAx(1), 'XLim'), 'YData', [currPos(xySel(1)) currPos(xySel(1))]);
                set(zoomReg, 'Position',zoomValXY-[.5 .5 0 0]);
                set(lineHorZoom, 'XData',[currPos(xySel(2)) currPos(xySel(2))], 'YData',  get(zoomAx(1), 'YLim'));
                set(lineVertZoom, 'XData', get(zoomAx(1), 'XLim'), 'YData', [currPos(xySel(1)) currPos(xySel(1))]);
                if meanPlotVal < 7
                    set([rectPlotArea_Zoom rectPlotArea_Im], 'Position',[currPos(xySel(2))-meanPlotVal/2 currPos(xySel(1))-meanPlotVal/2 meanPlotVal meanPlotVal]);
                elseif meanPlotVal == 7
                  if ~isempty(rectPlotArea_Zoom)
                    for ii = 1:3
                        set([rectPlotArea_Zoom(:,ii)], 'Position',[currPos(xySel(2))-(2*ii-1)/2 currPos(xySel(1))-(2*ii-1)/2 2*ii-1 2*ii-1]);
                    end
                  end
                  if ~isempty(rectPlotArea_Im)
                    for ii = 1:3
                        set([rectPlotArea_Im(:,ii)], 'Position',[currPos(xySel(2))-(2*ii-1)/2 currPos(xySel(1))-(2*ii-1)/2 2*ii-1 2*ii-1]);
                    end
                  end
                else
                    set([rectPlotArea_Zoom rectPlotArea_Im],'Position',zoomValXY-[.5 .5 0 0]);
                end
            end
        end
        if get(tbWin(3), 'Value')
            for jj = 1:nMat                                  %Number of data matrices
                for ii = 1 : nPlots                          %Number plot dimensions (subplots)
                    if customDimScale
                        currLoc = scaledLocation(currPos(plotDim(ii)), plotDim(ii));
                        set(posLine(jj,ii), 'XData', currLoc*[1 1],...
                            'YData', get(subPlotHandles(jj,ii), 'YLim'));
                    else
                        set(posLine(jj,ii), 'XData', [currPos(plotDim(ii)) currPos(plotDim(ii))],...
                            'YData', get(subPlotHandles(jj,ii), 'YLim'));
                    end
                end
            end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Update current plot values
    function updatePlotValues(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        for ii = 1:nDim
            imIndex{ii} = currPos(ii);   %#ok
        end
        %plotValues = [];
        for jj = 1:nMat
            for ii = 1:nPlots
                plotIndex = imIndex;
                if ~get(tbRoi, 'Value')
                    w = get(btMean, 'UserData');  %0: no averaging, 1: 3x3 average, 2: 5x5 average, 3:1x1,3x3,5x5 as red,green and blue plots, 4: across zoom region
                    % Averaging of area [x+w:x-w, y+w:y-w]
                    if w < 3
                        w(2) = w;  %for x and y
                    elseif w == 3
                        w = [0 0;1 1;2 2];  % for triple-plot
                    elseif w == 5  % RGB plot: no averaging
                        w = [0 0];
                    end
                    for kk = 1:size(w,1)
                        if w == 4  %Average over zoom region
                            plotIndex{xySel(1)} = round(zoomValXY(2)):zoomValXY(2)+zoomValXY(4)-1;
                            plotIndex{xySel(2)} = round(zoomValXY(1)):zoomValXY(1)+zoomValXY(3)-1;
                        else       %Average over 'w-interval'
                            plotIndex{xySel(1)} = max(1  ,currPos(xySel(1))-w(kk,2)) : min(currPos(xySel(1))+w(kk,2),dim(xySel(1)));
                            plotIndex{xySel(2)} = max(1  ,currPos(xySel(2))-w(kk,1)) : min(currPos(xySel(2))+w(kk,1),dim(xySel(2)));
                        end
                        plotIndex{plotDim(ii)} = ':';
                        % Plot along x direction
                        if plotDim(ii) == xySel(1)
                            if get(btMean, 'UserData')==5   %RGB plot
                                if get(cmStretchRGBMean, 'Value') || get(cmStretchRGBMax, 'Value') % Stretch RGB mode
                                    for ll = 1:dim(rgbDim)
                                        plotIndex{rgbDim} = ll;
                                        plotValues{jj,ii,ll} = data{jj}(plotIndex{:});
                                    end
                                else   %Normal RGB mode
                                    for ll=1:3
                                        plotValues{jj,ii,ll} = currImVal{jj}(:,currPos(xySel(2)),ll);
                                    end
                                end
                            else
                                if withFilter && get(popFilter, 'Value')>1  % If filter is selected, neglect "averaging button" and display filtered and unfiltered data instead.
                                    plotValues{jj,ii,1} = currImVal{jj}(:,currPos(xySel(2)));
                                    plotValues{jj,ii,2} = currIm{jj}(:,currPos(xySel(2)));
                                else
                                  if withAlpha
                                    plotValues{jj,ii,kk} =  sum(currImVal{jj}(:,plotIndex{xySel(2)}).*currAlphaMapVal{jj}(:,plotIndex{xySel(2)}),2, 'omitnan')./...
                                                            sum(currAlphaMapVal{jj}(:,plotIndex{xySel(2)}).*~isnan(currImVal{jj}(:,plotIndex{xySel(2)})),2, 'omitnan');
                                    plotValuesAlpha{jj,ii,kk} =  mean(currAlphaMapVal{jj}(:,plotIndex{xySel(2)}),    2, 'omitnan');
                                  else
                                    plotValues{jj,ii,kk} = mean(currImVal{jj}(:,plotIndex{xySel(2)}),2, 'omitnan');
                                  end
                                end
                            end
                            % Plot along y direction
                        elseif plotDim(ii) == xySel(2)
                            if get(btMean, 'UserData')==5   %RGB plot
                                if get(cmStretchRGBMean, 'Value') || get(cmStretchRGBMax, 'Value')  % Stretch RGB mode
                                    for ll = 1:dim(rgbDim)
                                        plotIndex{rgbDim} = ll;
                                        plotValues{jj,ii,ll} = data{jj}(plotIndex{:});
                                    end
                                else   %Normal RGB mode
                                    for ll=1:3
                                        plotValues{jj,ii,ll} = currImVal{jj}(currPos(xySel(1)),:,ll);
                                    end
                                end
                            else
                                if withFilter && get(popFilter, 'Value')>1  % If filter is selected, neglect "averaging button" and display filtered and unfiltered data instead.
                                    plotValues{jj,ii,1} = currImVal{jj}(currPos(xySel(1)),:);
                                    plotValues{jj,ii,2} = currIm{jj}(currPos(xySel(1)),:);
                                else
                                  if withAlpha
                                    plotValues{jj,ii,kk} =  sum(currImVal{jj}(plotIndex{xySel(1)},:).*currAlphaMapVal{jj}(plotIndex{xySel(1)},:),1, 'omitnan')./...
                                                            sum(currAlphaMapVal{jj}(plotIndex{xySel(1)},:).*~isnan(currImVal{jj}(plotIndex{xySel(1)},:)),1, 'omitnan');
                                    plotValuesAlpha{jj,ii,kk} =  mean(currAlphaMapVal{jj}(plotIndex{xySel(1)},:),    1, 'omitnan');
                                  else
                                    plotValues{jj,ii,kk} = mean(currImVal{jj}(plotIndex{xySel(1)},:),1, 'omitnan');
                                  end
                                end
                            end
                            % Plot along any other direction
                        elseif get(btMean, 'UserData')==5 && rgbDim ~= plotDim(ii)  %RGB plot - not possible if plot dimension and RGB dimension are identical
                            if get(tbSwitchRGB, 'Value') && (get(cmStretchRGBMean, 'Value') || get(cmStretchRGBMean, 'Value'))% Stretch RGB mode
                                for ll = 1:dim(rgbDim)
                                    plotIndex{rgbDim} = ll;
                                    plotValues{jj,ii,ll} = data{jj}(plotIndex{:});
                                end
                            else
                                for ll = -1:1      %Normal RGB mode
                                    plotIndex{rgbDim} = min(max(1,currPos(rgbDim)-ll),dim(rgbDim));
                                    plotValues{jj,ii,ll+2} = data{jj}(plotIndex{:});
                                end
                            end
                        else  % plot along potential RGB dimension
                          if withAlpha
                            %sum(c.*cA,3, 'omitnan')./sum(cA.*~isnan(c),3, 'omitnan'); % weighted mean ratio
                            plotValues{jj,ii,kk} =  sum(sum(data{jj}(plotIndex{:})    .*alphaMap{jj}(plotIndex{:}),    xySel(1), 'omitnan'),xySel(2), 'omitnan')./...
                                                    sum(sum(alphaMap{jj}(plotIndex{:}).*~isnan(data{jj}(plotIndex{:})),xySel(1), 'omitnan'),xySel(2), 'omitnan');
                            plotValuesAlpha{jj,ii,kk} =  mean(mean(alphaMap{jj}(plotIndex{:}),    xySel(1), 'omitnan'),xySel(2), 'omitnan');
                          else
                            plotValues{jj,ii,kk} = mean(mean(data{jj}(plotIndex{:}),xySel(1), 'omitnan'),xySel(2), 'omitnan');
                          end
                        end
                    end
                elseif nRois>0
                    %Get plot values for fs
                    deltaIndex = prod(dim(1:plotDim(ii)-1));                             %Difference value of indices along plotted dimension
                    selRois = get(roiListbox, 'Value');                                 %Get numbers of selected rois
                    for kk=1:size(selRois,2)                                             %Iterate over number of selected rois
                        plotIndex = imIndex;
                        for m=1:nDim                                                    %Fill indices of all dimension with current-position-vectors of length of roi-size
                            plotIndex{m} = plotIndex{m} * ones(size(roiList(selRois(kk)).index.x,1),1);
                        end
                        plotIndex{xySel(1)} = roiList(selRois(kk)).index.x;                       %Fill xySel dimension indices with roi indices
                        plotIndex{xySel(2)} = roiList(selRois(kk)).index.y;
                        plotIndex{plotDim(ii)} = ones(size(roiList(selRois(kk)).index.x,1),1);        %Fill plot-dimension with ones (for first point)
                        plotIndex = sub2ind(dim, plotIndex{:});                         %Determine linear index of roi pixels for first point
                        plotIndex = repmat(plotIndex, [1 dim(plotDim(ii))]);                %Replicate linear index
                        plotIndex = plotIndex + repmat(deltaIndex * (0:dim(plotDim(ii))-1),[size(roiList(selRois(kk)).index.x,1) 1]);   %Extend to all other points by adding deltaInd for each step
                        if withAlpha
                          plotValues{jj,ii,kk} =  sum(data{jj}(plotIndex)    .*alphaMap{jj}(plotIndex),    1, 'omitnan')./...
                                                  sum(alphaMap{jj}(plotIndex).*~isnan(data{jj}(plotIndex)),1, 'omitnan');
                          plotValuesAlpha{jj,ii,kk} =  mean(alphaMap{jj}(plotIndex),    1, 'omitnan');
                        else
                          plotValues{jj,ii,kk} = mean(data{jj}(plotIndex),1, 'omitnan');
                        end
                    end
                end
            end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Draw plots
    function drawPlots(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if get(tbWin(3), 'Value') == 1
            %Get number of subplots
            nPlots = sum(plotSel);
            plotDim = find(plotSel == 1);
            nPlotRow = ceil(nPlots/2);
            if nPlots == 1
                nPlotCol = 1;
            else
                nPlotCol = 2;
            end
            plotValues = [];plotValuesAlpha = [];
            updatePlotValues;
            if oldMATLAB
                subPlotHandles = [];
                subPlotPlots = [];
                if withAlpha
                  subPlotPlotsAlpha = [];
                end
            else
                subPlotHandles = matlab.graphics.chart.primitive.Line.empty;                     %Handles to subplot axes in Plots Window
                subPlotPlots = matlab.graphics.chart.primitive.Line.empty;                       %Handles to data displayed in Plots Window
                if withAlpha
                  subPlotPlotsAlpha = matlab.graphics.chart.primitive.Line.empty;                       %Handles to data displayed in Plots Window
                end
            end
            if get(tbRoi, 'Value')
                roiSel = get(roiListbox, 'Value');
                l = '';
            end
            for jj = 1:nMat
                set(0,'CurrentFigure',plotWin(jj));
                set(plotWin(jj), 'HandleVisibility', 'on');
                clf;
                colorcodePlots = plotColors(size(plotValues,3));
                for ii = 1 : nPlots
                    for kk = 1:size(plotValues,3)
                        subPlotHandles(jj,ii) = subplot(nPlotCol, nPlotRow, ii);
                        if ~get(tbRoi, 'Value')
                            if get(btMean, 'UserData') == 5  %RGB plot mode
                                if get(cmStretchRGBMean,'Value') || get(cmStretchRGBMax, 'Value')  %RGB stretch mode
                                    try   %for dimensions with less plots, e.g. rgbDim in RGB plot mode
                                        if customDimScale
                                            subPlotPlots{jj,ii,kk} = plot(subPlotHandles(jj,ii), linspace(dimScale(plotDim(ii),1),dimScale(plotDim(ii),2),dim(plotDim(ii))),squeeze(plotValues{jj,ii,kk}),...
                                                'Color', rgbValPlots(:,kk));
                                        else
                                            subPlotPlots{jj,ii,kk} = plot(subPlotHandles(jj,ii), squeeze(plotValues{jj,ii,kk}),...
                                                'Color', rgbValPlots(:,kk));
                                        end
                                        
                                    catch    %#ok
                                    end
                                else           % normal RGB mode
                                    try   %for dimensions with less plots, e.g. rgbDim in RGB plot mode
                                        if customDimScale
                                            subPlotPlots{jj,ii,kk} = plot(subPlotHandles(jj,ii),  linspace(dimScale(plotDim(ii),1),dimScale(plotDim(ii),2),dim(plotDim(ii))),squeeze(plotValues{jj,ii,kk}),...
                                                'Color', [kk==1 kk==2 kk==3]);
                                        else
                                            subPlotPlots{jj,ii,kk} = plot(subPlotHandles(jj,ii), squeeze(plotValues{jj,ii,kk}),...
                                                'Color', [kk==1 kk==2 kk==3]);
                                        end
                                    catch    %#ok
                                    end
                                end
                            else
                                if ~isempty(plotValues{jj,ii,kk})
                                    if customDimScale
                                        subPlotPlots{jj,ii,kk} = plot(subPlotHandles(jj,ii), linspace(dimScale(plotDim(ii),1),dimScale(plotDim(ii),2),dim(plotDim(ii))),squeeze(plotValues{jj,ii,kk}),...
                                            'Color', [kk==1 kk==2 kk==3], 'LineWidth', kk/2,'LineStyle','-');
                                        if withAlpha && ~verLessThan('matlab','9.0')
                                          yyaxis(subPlotHandles(jj,ii),'right')
                                          subPlotPlotsAlpha{jj,ii,kk} = plot(subPlotHandles(jj,ii), linspace(dimScale(plotDim(ii),1),dimScale(plotDim(ii),2),dim(plotDim(ii))),squeeze(plotValuesAlpha{jj,ii,kk}),...
                                            'Color', [kk==1 kk==2 kk==3], 'LineWidth', kk/2,'LineStyle',':','Marker','none');
                                          yyaxis(subPlotHandles(jj,ii),'left')
                                        end
                                    else
                                        subPlotPlots{jj,ii,kk} = plot(subPlotHandles(jj,ii), squeeze(plotValues{jj,ii,kk}),...
                                            'Color', [kk==1 kk==2 kk==3], 'LineWidth', kk/2,'LineStyle','-');
                                          if withAlpha && ~verLessThan('matlab','9.0')
                                            yyaxis(subPlotHandles(jj,ii),'right')
                                            subPlotPlotsAlpha{jj,ii,kk} = plot(subPlotHandles(jj,ii),squeeze(plotValuesAlpha{jj,ii,kk}),...
                                              'Color', [kk==1 kk==2 kk==3], 'LineWidth', kk/2,'LineStyle',':','Marker','none');
                                            yyaxis(subPlotHandles(jj,ii),'left')
                                          end
                                    end
                                end
                            end
                        else
                            if ~(isequal(roiSel,0) || isempty(roiSel))
                                if customDimScale
                                    subPlotPlots{jj,ii,kk} = plot(subPlotHandles(jj,ii), linspace(dimScale(plotDim(ii),1),dimScale(plotDim(ii),2),dim(plotDim(ii))),squeeze(plotValues{jj,ii,kk}),...
                                      'Color',colorcodePlots(kk,:),'LineStyle','-');
                                    if withAlpha && ~verLessThan('matlab','9.0')
                                      yyaxis(subPlotHandles(jj,ii),'right')
                                      subPlotPlotsAlpha{jj,ii,kk} = plot(subPlotHandles(jj,ii), linspace(dimScale(plotDim(ii),1),dimScale(plotDim(ii),2),dim(plotDim(ii))),squeeze(plotValuesAlpha{jj,ii,kk}),...
                                        'Color', colorcodePlots(kk,:),'LineStyle',':','Marker','none');%, 'LineWidth', kk/2
                                      yyaxis(subPlotHandles(jj,ii),'left')
                                    end
                                else
                                    subPlotPlots{jj,ii,kk} = plot(subPlotHandles(jj,ii), squeeze(plotValues{jj,ii,kk}),'Color',colorcodePlots(kk,:),'LineStyle','-');
                                    if withAlpha && ~verLessThan('matlab','9.0')
                                      yyaxis(subPlotHandles(jj,ii),'right')
                                      subPlotPlotsAlpha{jj,ii,kk} = plot(subPlotHandles(jj,ii), squeeze(plotValuesAlpha{jj,ii,kk}),...
                                        'Color', colorcodePlots(kk,:),'LineStyle',':','Marker','none');%, 'LineWidth', kk/2
                                      yyaxis(subPlotHandles(jj,ii),'left')
                                    end
                                end
                                l{kk} = roiList(roiSel(kk)).name;  %Legend strings for roi plots
                            end
                        end
                        hold all;
                    end
                    if get(tbRoi, 'Value') && ~isempty(subPlotPlots)
                        legend(l);
                    else
                        legend('off');
                    end
                    hold off;
                    if get(tbPlotsXLim, 'Value') == 1 && sum(plotDim(ii) == xySel) > 0
                        kk = find(xySel == plotDim(ii));
                        set(subPlotHandles(jj,ii), 'XLim', [zoomValXY(kk) zoomValXY(kk)+zoomValXY(kk+2)]);
                    else
                        axis(subPlotHandles(jj,ii),'tight');
                    end
                    posLine(jj,ii) = line([currPos(plotDim(ii)), currPos(plotDim(ii))], get(subPlotHandles(jj,ii), 'YLim'));
                    if customDimScale
                        if withDimUnits
                            xlabel([dimNames{plotDim(ii)} ' [' dimUnits{plotDim(ii)} ']']);
                        else
                            xlabel(dimNames{plotDim(ii)});
                        end
                    else
                        xlabel(dimNames(plotDim(ii)));
                    end
                end
                set(plotWin(jj), 'HandleVisibility', 'off');
            end
        end
        updatePlots('drawPlots');  %Necessary to make sure all axis limits etc. are correctly set
        if ~isempty(plotDim)
            updateObjects;  %Necessary to make sure the position lines have correct length
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function  colors = plotColors(nColors)
        if debugMatVis, debugMatVisFcn(1); end
%         if exist('color_map.m','file')
            colors = color_map(nColors,3,1);
%         else
%             if nColors == 1
%                 colors = [1 0 0];
%             elseif nColors < 8
%                 colors = jet(7);
%             else
%                 colors = jet(nColors);
%             end
%         end
        if debugMatVis, debugMatVisFcn(2); end
   end

%Update plots
    function updatePlots(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if get(tbRoi,'Value') && nRois == 0
            %return
        elseif get(tbWin(3), 'Value')
            if ~(nargin > 0 && strcmp(varargin{1}, 'drawPlots'))
                updatePlotValues;
            end
            for jj = 1:nMat                                  %Number of data matrices
                for ii = 1 : nPlots                          %Number plot dimensions (subplots)
                    for kk = 1:size(plotValues,3)            %Number of plots in one axes (more than 1 for mean plots and rois)
                        try   %For different number of plots in different plot axes
                            set(subPlotPlots{jj,ii,kk}, 'YData', squeeze(plotValues{jj,ii,kk}));
                            if withAlpha
                              set(subPlotPlotsAlpha{jj,ii,kk}, 'YData', squeeze(plotValuesAlpha{jj,ii,kk}));
                            end
                        catch    %#ok
                           fprintf('''set(subPlotPlots{jj,ii,kk}, ''YData'', squeeze(plotValues{jj,ii,kk}));'' failed\n')
                        end
                        try                                   %for dimensions with less plot entries, e.g. rgbDim for RGB plot mode
                            pV(kk,:) = squeeze(plotValues{jj,ii,kk});  %#ok
                            if withAlpha
                              pVA(kk,:) = squeeze(plotValuesAlpha{jj,ii,kk});  %#ok
                            end
                        catch    %#ok
                        end
                    end
                    %Set XLim
                    if get(tbPlotsXLim, 'Value')
                        try
                            if customDimScale
                                set(subPlotHandles(jj,ii), 'XLim', plotXLimScale(plotDim(ii),:));
                            else
                                set(subPlotHandles(jj,ii), 'XLim', plotXLim(plotDim(ii),:));
                            end
                        catch   %#ok
                        end
                    else
                        try
                            if customDimScale
                                set(subPlotHandles(jj,ii), 'XLim', [dimScale(plotDim(ii),1) dimScale(plotDim(ii),2)]);
                            else
                                set(subPlotHandles(jj,ii), 'XLim', [1 dim(plotDim(ii))]);
                            end
                        catch   %#ok
                        end
                    end
                    %Set YLim
                    if get(tbPlotsYLim, 'Value')
                        set(subPlotHandles(jj,ii), 'YLim', plotYLim(jj,:));
                    else
                        minValY = min(min(pV(:,plotXLim(plotDim(ii),1):plotXLim(plotDim(ii),2))));
                        maxValY = max(max(pV(:,plotXLim(plotDim(ii),1):plotXLim(plotDim(ii),2))));
                        if  ~isnan(minValY) && minValY ~= maxValY
                            set(subPlotHandles(jj,ii), 'YLim', [minValY-(maxValY-minValY)/25 maxValY+(maxValY-minValY)/25]);
                        end
                        if withAlpha && ~verLessThan('matlab','9.0')
                          minValYA = min(min(pVA(:,plotXLim(plotDim(ii),1):plotXLim(plotDim(ii),2))));
                          maxValYA = max(max(pVA(:,plotXLim(plotDim(ii),1):plotXLim(plotDim(ii),2))));
                          if  ~isnan(minValYA) && minValYA ~= maxValYA
                            yyaxis(subPlotHandles(jj,ii),'right')
                            set(subPlotHandles(jj,ii), 'YLim', [minValYA-(maxValYA-minValYA)/25 maxValYA+(maxValYA-minValYA)/25]);
                            yyaxis(subPlotHandles(jj,ii),'left')
                          end
                        end
                    end
                    %                     %Apply zoom intervals to axes if "Zoom axes" is selected
                    %                     if get(tbPlotsXLim, 'Value') == 1 && sum(plotDim(ii) == xySel) > 0
                    %                         kk = abs(find(xySel == plotDim(ii))-3);
                    %                         set(subPlotHandles(jj,ii), 'XLim', [zoomValXY(kk) zoomValXY(kk)+zoomValXY(kk+2)]);
                    %                     else
                    %                         %axis(subPlotHandles(jj,ii),'tight');
                    %                     end
                    %                     %Scale plots in y-direction if selected
                    %                     if get(tbPlotsYLim, 'Value') == 1
                    %                         try set(subPlotHandles(jj,ii), ...
                    %                                 'YLim',[min([plotValues{jj,ii,:}]) max([plotValues{jj,ii,:}])]); catch    %#ok end
                    %                     else
                    %                         set(subPlotHandles(jj,ii),...
                    %                             'YLim', [minVal(jj) maxVal(jj)]);
                    %                     end
                    clear pV pVA
                end
                if ~isempty(subPlotPlots)
                    if get(tbMarker, 'Value') == 1
                        set([subPlotPlots{:}], 'Marker', 'd', 'MarkerSize', 3);
                    else
                        set([subPlotPlots{:}],'Marker', 'none');
                    end
                end
            end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function setColormapFromEdit(varargin)
        % Called from edit fields for contrast values
            if str2double(get(valSldMin, 'String')) < get(sldMin, 'Min')
                set(valSldMin, 'String', num2str(get(sldMin, 'Min')));
            end
            if str2double(get(valSldMax, 'String')) > get(sldMax, 'Max')
                set(valSldMax, 'String', num2str(get(sldMax, 'Max')));
            end
            set(sldMin, 'Value', max(get(sldMin, 'Min'), str2double(get(valSldMin, 'String'))));
            set(sldMax, 'Value', min(get(sldMin, 'Max'), str2double(get(valSldMax, 'String'))));
            if rgbCount
                updateImages;
            else
                updateColormap;
            end
    end
%Adjust Colormap
    function updateColormap(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        %Get Slider Values in Manual Mode
        if (~(any(get(bg_colormap, 'SelectedObject') == [cmManual cmThresh])) && ~rgbCount)
            set(sldMin, 'Enable', 'off');
            set(sldMax, 'Enable', 'off');
            set(valSldMin, 'Enable', 'off');
            set(valSldMax, 'Enable', 'off');
        elseif get(bg_colormap, 'SelectedObject') == cmThresh
            set(sldMin, 'Enable', 'off');
            set(sldMax, 'Enable', 'on');
            set(valSldMin, 'Enable', 'off');
            set(valSldMax, 'Enable', 'on');
        else
            set(sldMin, 'Enable', 'on');
            set(sldMax, 'Enable', 'on');
            set(valSldMin, 'Enable', 'on');
            set(valSldMax, 'Enable', 'on');
        end
        switch get(bg_colormap, 'SelectedObject')
            case cmGlobal       %Global Min / Max Values (each data set seperateley)
                cmMinMax = cat(1  ,minVal, maxVal)';
            case cmImage        %Scaled to Min / Max Values of currently displayed image(s)
                for ii = 1:nMat
                    cmMinMax(ii,:) = [min(currImVal{ii}(:)) max(currImVal{ii}(:))]';
                end
            case cmZoom
                for ii = 1:nMat
                    w = currImVal{ii}(zoomValXY(2):zoomValXY(2)+zoomValXY(4)-1,zoomValXY(1):zoomValXY(1)+zoomValXY(3)-1);
                    cmMinMax(ii,:) = [min(w(:)) max(w(:))]';
                end
            case cmManual       %Values from Min / Max sliders
                if isinteger(data{currContrastSel})
                    cmMinMax(currContrastSel  ,:) = round([get(sldMin, 'Value') get(sldMax, 'Value')]);
                else
                    cmMinMax(currContrastSel  ,:) = [get(sldMin, 'Value') get(sldMax, 'Value')];
                end
                if linkContrastSettings
                    if ~(cmMinMax(currContrastSel  ,1) < cmMinMax(currContrastSel  ,2))
                        w = get(sldMax, 'SliderStep');
                        cmMinMax(currContrastSel  ,2) = cmMinMax(currContrastSel  ,1) + w(1)*(maxVal(1)-minVal(1));
                        set(sldMin, 'Value', cmMinMax(currContrastSel  ,1));
                        set(sldMax, 'Value', cmMinMax(currContrastSel  ,2));
                    end
                    idx = 1:nMat;
                    idx(currContrastSel) = [];
                    for ii = idx
                        % Min val
                        cmMinMax(ii,1) = minVal(ii)+(cmMinMax(currContrastSel,1)-minVal(currContrastSel))*(maxVal(ii)-minVal(ii))/(maxVal(currContrastSel)-minVal(currContrastSel));
                        % Max val
                        cmMinMax(ii,2) = minVal(ii)+(cmMinMax(currContrastSel,2)-minVal(currContrastSel))*(maxVal(ii)-minVal(ii))/(maxVal(currContrastSel)-minVal(currContrastSel));
                        cmMinMax(isnan(cmMinMax)) = 0;
                    end
                end
            case cmThresh
                if isinteger(data{currContrastSel})
                    cmMinMax(currContrastSel  ,:) = round([get(sldMin, 'Value') get(sldMax, 'Value')]);
                else
                    cmMinMax(currContrastSel  ,:) = [get(sldMin, 'Value') get(sldMax, 'Value')];
                end
                thAlg = get(popThresh, 'String');
                thAlg = thAlg{get(popThresh, 'Value')};
                switch thAlg
                    case 'triangle'  % Method from Zack GW, Rogers WE, Latt SA (1977)
                        % Dipimage implementation
                        if withDipimage
                            [~, cmMinMax(currContrastSel,1)] = threshold(currImVal{1},'triangle');
                        else
                            % Implementation adopted from file exchange
                            % submission by Dr B. Panneton (triangle_th.m):
                            % Might produce slightly different results than
                            % dipimage implementation!
                            %   Find maximum of histogram and its location along the x axis
                            guiHistVal = histcounts(currImVal{1}(:), histXData);
                            [~,xmax] = max(guiHistVal);
                            num_bins = numel(histXDataBin);
                            xmax = round(mean(xmax));   %can have more than a single value!
                            h = guiHistVal(xmax);
                            %   Find location of first and last non-zero values.
                            %   Values<h/10000 are considered zeros.
                            indi = find(guiHistVal>h/10000);
                            fnz  = indi(1);
                            lnz  = indi(end);
                            %   Pick side as side with longer tail. Assume one tail is longer.
                            lspan=xmax-fnz;
                            rspan=lnz-xmax;
                            if rspan>lspan  % then flip guiHistVal
                                ghv=fliplr(guiHistVal);
                                a=num_bins-lnz+1;
                                b=num_bins-xmax+1;
                                isflip=1;
                            else
                                ghv=guiHistVal;
                                isflip=0;
                                a=fnz;
                                b=xmax;
                            end
                            %   Compute parameters of the straight line from first non-zero to peak
                            %   To simplify, shift x axis by a (bin number axis)
                            m=h/(b-a);
                            %   Compute distances
                            x1=0:(b-a);
                            y1=ghv(x1+a);
                            bb=y1+x1/m;
                            x2=bb/(m+1/m);
                            y2=m*x2;
                            L=((y2-y1).^2+(x2-x1).^2).^0.5;
                            %   Obtain threshold as the location of maximum L.
                            thresh=find(max(L)==L);
                            thresh=a+mean(thresh);
                            %   Flip back if necessary
                            if isflip
                                thresh=num_bins-thresh+1;
                            end
                            thresh = thresh/num_bins;
                            cmMinMax(currContrastSel,1)= histXDataBin(1)+(histXDataBin(end)-histXDataBin(1))*thresh;
                        end
                    case 'otsu'
                        % Dipimage implementation
                        if withDipimage
                            [~, cmMinMax(currContrastSel,1)] = threshold(currImVal{1},'otsu');
                        else
                            % Matlab implementation
                            thresh = graythresh(scaleMinMax(currImVal{1})); % Matlab implementation of Otsu's algorithm
                            range = [min(currImVal{1}(:)) max(currImVal{1}(:))];
                            cmMinMax(currContrastSel,1) = range(1)+(range(2)-range(1))*thresh;
                        end
                    case 'max. entropy' % Maximum entropy algorithm adopted from file exchange submission maxentropy.m from F.Gargouri
                        n = dim(xySel(1));
                        m = dim(xySel(2));
                        ghv = histcounts(uint8(255*scaleMinMax(double(currImVal{1}(:)))),0:255);
                        %normalize the histogram ==>  hn(k)=h(k)/(n*m) ==> k  in [1 256]
                        hn=ghv/(n*m);
                        %Cumulative distribution function
                        cdfct(1) = hn(1);
                        for ll = 2:256
                            cdfct(ll)=cdfct(ll-1)+hn(ll);
                        end
                        hl = zeros(1,256);
                        hh = zeros(1,256);
                        for tt=1:256
                            %low range entropy
                            cl=double(cdfct(tt));
                            if cl>0
                                for ii=1:tt
                                    if hn(ii)>0
                                        hl(tt) = hl(tt)- (hn(ii)/cl)*log(hn(ii)/cl);
                                    end
                                end
                            end
                            %high range entropy
                            ch=double(1.0-cl);  %constraint cl+ch=1
                            if ch>0
                                for ii=tt+1:256
                                    if hn(ii)>0
                                        hh(tt) = hh(tt)- (hn(ii)/ch)*log(hn(ii)/ch);
                                    end
                                end
                            end
                        end
                        % choose best threshold
                        h_max =hl(1)+hh(1);
                        thresh = 0;
                        entropie(1)=h_max;
                        for tt=2:256
                            entropie(tt)=hl(tt)+hh(tt);
                            if entropie(tt)>h_max
                                h_max=entropie(tt);
                                thresh=tt-1;
                            end
                        end
                        cmMinMax(currContrastSel,1)  = (double(thresh)/255)*(maxVal(1)-minVal(1))+minVal(1);
                    case 'isodata'
                        % Dipimage implementation
                        if withDipimage
                            [~, cmMinMax(currContrastSel,1)] = threshold(currImVal{1},'isodata');
                        else
                            % Isodata algorithm from  Ridler and Calvard;
                            % implementation adopted from Jing Tian's file
                            % exchange function func_threshold. Might
                            % produce slightly different results than
                            % Dipimage implementation
                            % STEP 1: Compute mean intensity of image from histogram, set T(1) = mean(I)
                            [counts, N] = histcounts(currImVal{1}(:),histXData);
                            ct = 1;
                            mu1 = cumsum(counts);
                            T(ct) = (sum(N.*counts)) / mu1(end);
                            
                            % STEP 2: compute the sample mean of the data classified by the above threshold
                            mu2 = cumsum(counts(N<=T(ct)));
                            MBT = sum(N(N<=T(ct)).*counts(N<=T(ct)))/mu2(end);
                            
                            mu3 = cumsum(counts(N>T(ct)));
                            MAT = sum(N(N>T(ct)).*counts(N>T(ct)))/mu3(end);
                            ct=ct+1;
                            
                            % new T = (MAT+MBT)/2
                            T(ct) = (MAT+MBT)/2;
                            
                            % STEP 3 repeat step 2 while T(i)~=T(i-1)
                            thre = T(ct);
                            
                            while abs(T(ct)-T(ct-1))>=1
                                mu2 = cumsum(counts(N<=T(ct)));
                                MBT = sum(N(N<=T(ct)).*counts(N<=T(ct)))/mu2(end);
                                
                                mu3 = cumsum(counts(N>T(ct)));
                                MAT = sum(N(N>T(ct)).*counts(N>T(ct)))/mu3(end);
                                
                                ct=ct+1;
                                T(ct) = (MAT+MBT)/2;
                                thre = T(ct);
                            end
                            cmMinMax(currContrastSel,1) = thre;
                        end
                    case 'minerror'
                        % Dipimage implementation
                        [~,cmMinMax(currContrastSel,1)] = threshold(currImVal{1},'minerror');
                    case 'background'
                        % Dipimage implementation
                        [~,cmMinMax(currContrastSel,1)] = threshold(currImVal{1},'background');
                end
                if cmMinMax(currContrastSel,2) <= cmMinMax(currContrastSel,1)
                    if isinteger(data{currContrastSel})
                        cmMinMax(currContrastSel,2) = cmMinMax(currContrastSel,1)+1;
                    else
                        cmMinMax(currContrastSel,2) = cmMinMax(currContrastSel,1)*((maxVal(1)-minVal(1))/256)^sign(cmMinMax(currContrastSel,1));
                    end
                    set(sldMax, 'Value', cmMinMax(currContrastSel  ,2));
                end
                if isinteger(data{currContrastSel})
                    cmMinMax(currContrastSel,1) = round(cmMinMax(currContrastSel,1));
                end
        end
        %Set displayed Min / Max Values
        if isinteger(data{currContrastSel})
            cmMinMax(currContrastSel,:) = round(cmMinMax(currContrastSel,:));
            set(valSldMin, 'String', num2str(cmMinMax(currContrastSel  ,1)));
            set(valSldMax, 'String', num2str(cmMinMax(currContrastSel  ,2)));
        else
            set(valSldMin, 'String', num2str(cmMinMax(currContrastSel  ,1),'%6.3f'));
            set(valSldMax, 'String', num2str(cmMinMax(currContrastSel  ,2),'%6.3f'));
        end
        %Set Colortable (Look Up Table)
        lutStr = get(popLut,'String');
        switch lutStr{get(popLut, 'Value')}
            case 'Gray'
                cmap = gray(255);
            case 'Gray (Range)'
                cmap = color_map(255, 5);
            case '4 x Gray'
                cmap = repmat(gray(64),[4 1]);
            case 'Parula'
                cmap = parula(255);
            case 'Jet'
                cmap = jet(255);
            case 'HSV'
                cmap = hsv(255);
            case 'Hot'
                cmap = hot(255);
            case 'Cool'
                cmap = cool(255);
            case 'Red 1'  %Black - Red
                cmap = linspace(0,1,255)';
                cmap(:,2:3) = 0;
            case 'Red 2'  %Black - Red - White
                cmap = linspace(0,1,128)';
                cmap(129:255) = 1;
                cmap(:,2:3) = 0;
                cmap(129:255,2:3) = repmat(linspace(0,1,127)',[1,2]);
            case 'Green 1'  %Black - Green
                cmap = linspace(0,1,255)';
                cmap(:,2:3) = 0;
                cmap = circshift(cmap,[0 1]);
            case 'Green 2'  %Black - Green - White
                cmap = linspace(0,1,128)';
                cmap(129:255) = 1;
                cmap(:,2:3) = 0;
                cmap(129:255,2:3) = repmat(linspace(0,1,127)',[1,2]);
                cmap = circshift(cmap,[0 1]);
            case 'Blue 1'  %Black - Blue
                cmap = linspace(0,1,255)';
                cmap(:,2:3) = 0;
                cmap = circshift(cmap,[0 2]);
            case 'Blue 2'  %Black - Blue - White
                cmap = linspace(0,1,128)';
                cmap(129:255) = 1;
                cmap(:,2:3) = 0;
                cmap(129:255,2:3) = repmat(linspace(0,1,127)',[1,2]);
                cmap = circshift(cmap,[0 2]);
            case {'Rainbow1';'Rainbow2';'Rainbow3';'Rainbow4'}
                cmap = color_map(255, str2double(lutStr{get(popLut, 'Value')}(8)));
            case {'Blue-Gray-Red (0 centered)';'Blue-Gray-Yellow (0 centered)';'Magenta-Gray-Green (0 centered)'}
                mn = cmMinMax(currContrastSel,1); %minScale
                mx = cmMinMax(currContrastSel,2); %maxScale
                myGV = .85; % GrayValue
                if sign(mn)*sign(mx) == 1 % Use as "normal" colormap if all values are either positive or negative
                    cmap(:,1) = [linspace(0,myGV,128) linspace(myGV,1,127)];
                    cmap(:,2) = [linspace(0,myGV,128) linspace(myGV,0,127)];
                    cmap(:,3) = [linspace(1,myGV,128) linspace(myGV,0,127)];
                else
                    % Set 0 to gray and negative values blue, positive values yellow.
                    % The color gradient is identical in  negative and
                    % positive directions, independent of the distance
                    % of mn and mx from zero.
                    mnLength = round(255*abs(mn)/(mx-mn));
                    mxLength = 255-mnLength;
                    cmap = myGV * ones(255,3);
                    mxCorrL = 0; mnCorrL = 0;
                    mxCorrH = 1; mnCorrH = 1;
                    if abs(mn) > mx
                        mxCorrL = (1-mx/abs(mn))*myGV;
                        mxCorrH = mx/abs(mn)*(1-myGV)+myGV;
                    elseif abs(mn) < mx
                        mnCorrL = (1-abs(mn)/mx)*myGV;
                        mnCorrH = abs(mn)/mx*(1-myGV)+myGV;
                    else
                    end
                    cmap(1:mnLength,1) = linspace(mnCorrL,myGV,mnLength);
                    cmap(1:mnLength,2) = linspace(mnCorrL,myGV,mnLength);
                    cmap(1:mnLength,3) = linspace(mnCorrH,myGV,mnLength);
                    cmap(mnLength+1:end,1) = linspace(myGV,mxCorrH,mxLength);
                    cmap(mnLength+1:end,2) = linspace(myGV,mxCorrL,mxLength);
                    cmap(mnLength+1:end,3) = linspace(myGV,mxCorrL,mxLength);
                end
                if ~strcmp(lutStr{get(popLut, 'Value')},'Blue-Gray-Red (0 centered)')
                    cmap(:,2)=cmap(:,1); %LUT for 'Blue-Gray-Yellow (0 centered)'
                    if ~strcmp(lutStr{get(popLut, 'Value')},'Blue-Gray-Yellow (0 centered)')
                        cmap(:,1)=cmap(:,3); %LUT for 'Magenta-Gray-Green (0 centered)'
                    end
                end
                % set gamma value to 1
                if any(currGamma(1) ~= 1)
                    set(sldGamma(1), 'Value',1);
                    set(valSldGamma(1),'String','1.000');
                    updateImages;
                end

%             case 'Blue-Gray-Yellow (0 centered)'
%                 mn = cmMinMax(currContrastSel,1);
%                 mx = cmMinMax(currContrastSel,2);
%                 if sign(mn)*sign(mx) == 1 % Use as "normal" colormap if all values are either positive or negative
%                     cmap(:,1) = linspace(0,1,255);
%                     cmap(:,2) = linspace(0,1,255);
%                     cmap(:,3) = linspace(1,0,255);
%                 else
%                     % Set 0 to gray and negative values blue, positive values yellow.
%                     % The color gradient is identical in  negative and
%                     % positive directions, independent of the distance
%                     % of mn and mx from zero.
%                     mnLength = round(255*abs(mn)/(mx-mn));
%                     mxLength = 255-mnLength;
%                     cmap = 0.5 * ones(255,3);
%                     if abs(mn) > mx
%                         mnCorr = 0.5;
%                         mxCorr = mx/abs(mn)/2;
%                     elseif abs(mn) < mx
%                         mnCorr = abs(mn)/mx/2;
%                         mxCorr = 0.5;
%                     else
%                         mnCorr = 0.5;
%                         mxCorr = 0.5;
%                     end
%                     cmap(1:mnLength,1) = linspace(0.5-mnCorr,.5,mnLength);
%                     cmap(1:mnLength,2) = linspace(0.5-mnCorr,.5,mnLength);
%                     cmap(1:mnLength,3) = linspace(0.5+mnCorr,.5,mnLength);
%                     cmap(mnLength+1:end,1) = linspace(.5,0.5+mxCorr,mxLength);
%                     cmap(mnLength+1:end,2) = linspace(.5,0.5+mxCorr,mxLength);
%                     cmap(mnLength+1:end,3) = linspace(.5,0.5-mxCorr,mxLength);
%                 end
%             case 'Magenta-Gray-Green (0 centered)'
%                 mn = cmMinMax(currContrastSel,1);
%                 mx = cmMinMax(currContrastSel,2);
%                 if sign(mn)*sign(mx) == 1 % Use as "normal" colormap if all values are either positive or negative
%                     cmap(:,1) = linspace(1,0,255);
%                     cmap(:,2) = linspace(0,1,255);
%                     cmap(:,3) = linspace(1,0,255);
%                 else
%                     % Set 0 to gray and negative values blue, positive values yellow.
%                     % The color gradient is identical in  negative and
%                     % positive directions, independent of the distance
%                     % of mn and mx from zero.
%                     mnLength = round(255*abs(mn)/(mx-mn));
%                     mxLength = 255-mnLength;
%                     cmap = 0.5 * ones(255,3);
%                     if abs(mn) > mx
%                         mnCorr = 0.5;
%                         mxCorr = mx/abs(mn)/2;
%                     elseif abs(mn) < mx
%                         mnCorr = abs(mn)/mx/2;
%                         mxCorr = 0.5;
%                     else
%                         mnCorr = 0.5;
%                         mxCorr = 0.5;
%                     end
%                     cmap(1:mnLength,1) = linspace(0.5+mnCorr,.5,mnLength);
%                     cmap(1:mnLength,2) = linspace(0.5-mnCorr,.5,mnLength);
%                     cmap(1:mnLength,3) = linspace(0.5+mnCorr,.5,mnLength);
%                     cmap(mnLength+1:end,1) = linspace(.5,0.5-mxCorr,mxLength);
%                     cmap(mnLength+1:end,2) = linspace(.5,0.5+mxCorr,mxLength);
%                     cmap(mnLength+1:end,3) = linspace(.5,0.5-mxCorr,mxLength);
%                 end
            otherwise
                cmap = defaultColormap{1};
                if get(bg_colormap, 'SelectedObject') == cmGlobal
                    cmMinMax = repmat([0 255], [nMat 1]);
                end
        end
        % Old version to apply gamma correction via colormap. Problem: Colormap is always stretched
        % with equi-distance over data and can not have more than 255 entries. Now by
        % re-calculating the data (see updateCurrIm).
        %         if currGamma ~= 1
        %             interpVal = 1;  % Interpolation value to avoid 'mosaicing'
        %             if currGamma >= 1
        %                 cmap_gamma = round(scaleMinMax(brighten(1:255*interpVal, 1/currGamma - 1))*(255*interpVal-1)+1);
        %             else
        %                 cmap_gamma = round(scaleMinMax(brighten(1:255*interpVal, 1 - currGamma))*(255*interpVal-1)+1);
        %             end
        %             cmap_old = cmap;
        %             cmap = zeros(interpVal*255,3);
        %             for ii = 1:3
        % %                 cmap_tmp = interp(cmap_old(:,ii),interpVal);
        %                 cmap(:,ii) = cmap_old(cmap_gamma,ii);
        %             end
        %             cmap(cmap>1)=1;
        %             cmap(cmap<0)=0;
        %         end
        if any(currGamma ~= 1)
            updateImages;
        end
        if get(tbInvert, 'Value')
            cmap = 1 - cmap;
        end
        if get(tbFlip, 'Value')
            cmap = flipdim(cmap,1);
        end
        if get(tbWin(1), 'Value')
            for ii = 1:nMat
                if cmMinMax(ii,1) ~= cmMinMax(ii,2)
                    set(imAx(ii), 'CLim', cmMinMax(ii,:));
                    colormap(imAx(ii), cmap);
                end
            end
        end
        if get(tbWin(2), 'Value')
            for ii = 1:nMat
                if cmMinMax(ii,1) ~= cmMinMax(ii,2)
                    set(zoomAx(ii), 'CLim', cmMinMax(ii,:));
                    colormap(zoomAx(ii), cmap);
                end
            end
        end
        if get(tbProfile, 'Value')
            if cmMinMax(ii,1) ~= cmMinMax(ii,2)
                set(profileStruct.gui.traceim.ax, 'CLim', cmMinMax(ii,:));
                colormap(profileStruct.gui.traceim.ax, cmap);
            end
        end
        if get(tbRoi, 'Value')
            if cmMinMax(ii,1) ~= cmMinMax(ii,2)
                set(roiAxes, 'CLim', cmMinMax(ii,:));
                colormap(roiAxes, cmap);
            end
        end
        if get(tbColorbar, 'Value') == 1 && all(currGamma == 1)
            showColorbar;
%             for ii = 1:nMat
%                 colorbar('peer', imAx(ii));
%             end
        end
        updateHistObjects;
        if rgbCount
            colormap(contrastAx, gray(255));
        else
            colormap(contrastAx,cmap);
        end
        w = get(histAxGui, 'Ylim');
        if cmMinMax(currContrastSel,1)~=cmMinMax(currContrastSel,2)
            set(contrastAx,'CLim', cmMinMax(currContrastSel,:));
            if strcmp(get(histAxGui, 'YScale'),'log')
                set(histAxBg, 'Position', [cmMinMax(currContrastSel,1) 1 cmMinMax(currContrastSel,2)-cmMinMax(currContrastSel,1) 1.1*w(2)]);
            else
                set(histAxBg, 'Position', [cmMinMax(currContrastSel,1) 0 cmMinMax(currContrastSel,2)-cmMinMax(currContrastSel,1) 1.1*w(2)]);
            end
        end
        % Update contrast setting in 2D hist
        % if contrast values are changed in main matVis
        if with2DHist 
            [mv zv(1,1)] = min(abs(x_2DHist-cmMinMax(1,1)));
            [mv zv(1,2)] = min(abs(x_2DHist-cmMinMax(1,2)));
            zv(1,2) =  min([numel(x_2DHist),zv(1,2)])-zv(1,1)+1;
            [mv zv(2,1)] = min(abs(y_2DHist-cmMinMax(nMat+1,1)));
            zv(2,2) = min([numel(y_2DHist),find(y_2DHist>=cmMinMax(nMat+1,2),1)])-find(y_2DHist>=cmMinMax(nMat+1,1),1)+1;
            matVis2DHist.fctHandles.loadNewSettings('zoomVal',zv);
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateAlphaContrast(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if nargin && strcmp(varargin{3},'etxt')
            cmMinMax(nMat+currContrastSel,1) = str2double(get(edt_alphaMin, 'String'));
            cmMinMax(nMat+currContrastSel,2) = str2double(get(edt_alphaMax, 'String'));
        end
        if cmMinMax(nMat+currContrastSel,2) <= cmMinMax(nMat+currContrastSel,1)
            cmMinMax(nMat+currContrastSel,2) = cmMinMax(nMat+currContrastSel,1)+0.01;
        end
        set(edt_alphaMin, 'String',num2str(cmMinMax(nMat+currContrastSel,1),'%6.2f'));
        set(edt_alphaMax, 'String',num2str(cmMinMax(nMat+currContrastSel,2),'%6.2f'));
        for ii=1:nMat
            if cmMinMax(nMat+ii,1)==cmMinMax(nMat+ii,2)
                set([imAx(ii) zoomAx(ii)], 'ALim', [-eps 0]+cmMinMax(nMat+ii,:));
            else
                set([imAx(ii) zoomAx(ii)], 'ALim', cmMinMax(nMat+ii,:));
            end
        end
        % Update alpha contrast setting in 2D hist
        % if alpha contrast values are changed in main matVis
        if with2DHist
            zv(1,1) = find(x_2DHist>=cmMinMax(1,1),1) ;
            zv(1,2) =  min([numel(x_2DHist),find(x_2DHist>=cmMinMax(1,2),1)])-find(x_2DHist>=cmMinMax(1,1),1)+1;
            zv(2,1) = find(y_2DHist>=cmMinMax(nMat+1,1),1) ;
            zv(2,2) = min([numel(y_2DHist),find(y_2DHist>=cmMinMax(nMat+1,2),1)])-find(y_2DHist>=cmMinMax(nMat+1,1),1)+1;
            matVis2DHist.fctHandles.loadNewSettings('zoomVal',zv);
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateInvertMode(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        updateColormap;
        updateImages;
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateFlipMode(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        updateColormap;
        updateImages;
        if debugMatVis, debugMatVisFcn(2); end
    end
    function res = color_map(map_size, fun_type, exp)
        if debugMatVis, debugMatVisFcn(1); end
        % Function kindly provided by André Zeug (last modified 21.03.2006).
        % This function creates a colormap (n x 3 array of RGB triples between 0
        % and 1). It can be used analogously to 'gray', 'jet', 'hsv' etc. in
        % colormap( color_map(map_size, fun_type, exp) ).
        
        if nargin < 3
            exp=1;
        end
        
        switch fun_type
            case 1
                fun=@RGB_value1;
            case 2
                fun=@RGB_value2;
            case 3
                fun=@RGB_value3;
            case 4
                fun=@RGB_value4;
            case 5
                fun=@RGB_value5;
        end
        
        if map_size>1
            ms=ceil(map_size);
            for kk=1:ms
                res(kk,:)=fun((kk-1)/(ms-1), exp); %#ok
            end
        else
            res=fun(map_size, exp);
        end
        % black, blue, cyan, green, yellow, red, white
        function rgb_val=RGB_value1(col, exp)
            if col<1/6
                rgb_val=[ 0 0 (col*6)^exp];
            elseif col<2/6
                rgb_val=[ 0 (col*6-1)^exp 1];
            elseif col<3/6
                rgb_val=[ 0 1 (3-col*6)^exp];
            elseif col<4/6
                rgb_val=[(col*6-3)^exp 1 0];
            elseif col<5/6
                rgb_val=[1 (5-col*6)^exp 0];
            else
                rgb_val=[1 (col*6-5)^exp (col*6-5)^exp];
            end
        end
        % black, blue, cyan, green, yellow, red, magenta, white
        function rgb_val=RGB_value2(col, exp)
            if col<1/7
                rgb_val=[ 0 0 (col*7)^exp];
            elseif col<2/7
                rgb_val=[ 0 (col*7-1)^exp 1];
            elseif col<3/7
                rgb_val=[ 0 1 (3-col*7)^exp];
            elseif col<4/7
                rgb_val=[(col*7-3)^exp 1 0];
            elseif col<5/7
                rgb_val=[1 (5-col*7)^exp 0];
            elseif col<6/7
                rgb_val=[1 0 (col*7-5)^exp];
            else
                rgb_val=[1 (col*7-6)^exp 1];
            end
        end
        % blue, cyan, green, yellow, red
        function rgb_val=RGB_value3(col, exp)
            if col<1/4
                rgb_val=[ 0 (col*4)^exp 1];
            elseif col<2/4
                rgb_val=[ 0 1 (2-col*4)^exp];
            elseif col<3/4
                rgb_val=[(col*4-2)^exp 1 0];
            else
                rgb_val=[1 (4-col*4)^exp 0];
            end
            %rgb_val=rgb_val./repmat(sum(rgb_val,2),[1 3]);
        end
        % blue, cyan, green, yellow, red, magenta
        function rgb_val=RGB_value4(col, exp)
            if col<1/5
                rgb_val=[ 0 (col*5)^exp 1];
            elseif col<2/5
                rgb_val=[ 0 1 (2-col*5)^exp];
            elseif col<3/5
                rgb_val=[(col*5-2)^exp 1 0];
            elseif col<4/5
                rgb_val=[1 (4-col*5)^exp 0];
            else
                rgb_val=[1 0 (col*5-4)^exp];
            end
        end
        % range indicator with 0: blue, 1: red
        function rgb_val=RGB_value5(col, exp)
            if col==0
                rgb_val=[0 0 1];
            elseif col==1
                rgb_val=[1 0 0];
            else
                rgb=col^exp;
                rgb_val=[rgb rgb rgb];
            end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateGamma(h,e,source,type)    %#ok
        if debugMatVis, debugMatVisFcn(1); end
        if strcmp(type,'data')
            idx = 1;
        else
            idx = 2;
        end
        switch source
            case 'etxt'
                currGamma(currContrastSel+(idx-1)*nMat) = min(5,max(0,str2double(get(valSldGamma(idx), 'String'))));
                set(sldGamma(idx), 'Value', currGamma(currContrastSel+(idx-1)*nMat));
            case 'sld'
                currGamma(currContrastSel+(idx-1)*nMat) = get(sldGamma(idx), 'Value');
            case 'value'
                set(sldGamma(idx), 'Value', currGamma(currContrastSel+(idx-1)*nMat));
        end
        if linkContrastSettings
            currGamma((1:nMat)+(idx-1)*nMat) = ones(1,nMat)*currGamma(currContrastSel+(idx-1)*nMat);
        end
        set(valSldGamma(idx), 'String', num2str(currGamma(currContrastSel+(idx-1)*nMat),'%6.3f'));
        updateImages;
        if idx == 1
            if ~get(tb_moveHist,'Value') && strcmp(get(histAxGui, 'UserData'),'gamma')
                updateGuiHistVal;
            end
            % If currGamma == 1, the colorbar in between the sliders and the
            % yticks of the colorbar in the image window are NOT updated in
            % updateCurrIm, so they have to be updated here.
            if currGamma(currContrastSel) == 1
                set(contrastSldIm, 'CData',histXDataBin);
                if get(tbColorbar, 'Value') == 1
                    for ii=1:nMat
                        cb = findobj(imageWin(ii), 'Tag','Colorbar');
                        set(cb, 'YTickLabel', get(cb,'YTick'));
                    end
                end
            end
            updateHistObjects;
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Update Zoomregion
    function updateZoom(hObject, event, dimNum, source)   %#ok
        if debugMatVis, debugMatVisFcn(1); end
        % If called from etxt or slider, update zoom values first
        if nargin > 2
            switch source
                case 'etxt'
                    zv(1) = round(str2double(get(etxt_down(dimNum), 'String')));
                    zv(2) = round(str2double(get(etxt_up(dimNum), 'String')));
                    if zv(1) > zv(2)
                        zv = zv([2 1]);
                    end
                    zoomVal(dimNum,1) = max([1,zv(1)]);
                    zoomVal(dimNum,2) = min([dim(dimNum),zv(2) - zv(1) + 1]);
                    set(etxt_down(dimNum), 'String', num2str(zoomVal(dimNum,1)));
                    set(etxt_up(dimNum), 'String', num2str(sum(zoomVal(dimNum,:))-1));
                    set(sld_down(dimNum), 'Value', zoomVal(dimNum,1));
                    set(sld_up(dimNum), 'Value', sum(zoomVal(dimNum,:))-1);
                case 'sld'
                    zv(1) = round(get(sld_down(dimNum), 'Value'));
                    zv(2) = round(get(sld_up(dimNum), 'Value'));
                    if zv(1) > zv(2)
                        zv = zv([2 1]);
                    end
                    zoomVal(dimNum,1) = zv(1);
                    zoomVal(dimNum,2) = zv(2) - zv(1) + 1;
                    set(etxt_down(dimNum), 'String', num2str(zoomVal(dimNum,1)));
                    set(etxt_up(dimNum), 'String', num2str(sum(zoomVal(dimNum,:))-1));
                    set(sld_down(dimNum), 'Value', zoomVal(dimNum,1));
                    set(sld_up(dimNum), 'Value', sum(zoomVal(dimNum,:))-1);
                case 'lockBt'
                    if get(tb_lockPlots2Zoom(dimNum), 'Value')
                        set(tb_lockPlots2Zoom(dimNum), 'CData', lockClosed);
                        set(tbPlotsXLim, 'Value',1); % "Push" button as selecting the zoom lock will only have an effect with the XLim button pressed
                    else
                        set(tb_lockPlots2Zoom(dimNum), 'CData', lockOpen);
                        plotXLim(dimNum,:) = [1 dim(dimNum)];
                        if customDimScale
                            plotXLimScale(dimNum,:) = dimScale(dimNum,:);
                        end
                    end
              case 'ROIupdate'
                zoomValXY([1,3]) = zoomVal(xySel(2),:);
                zoomValXY([2,4]) = zoomVal(xySel(1),:);
            end
            for ii = dimNum
              set(txt_zoomWidth(ii), 'String', [' = ' num2str(zoomVal(ii,2))]);
            end
            if length(dimNum) < 3 && (any(xySel == dimNum) || (projMethod && any(projDim == dimNum)))
                zoomValXY([1,3]) = zoomVal(xySel(2),:);
                zoomValXY([2,4]) = zoomVal(xySel(1),:);
            end
        end
        %Border values for zoom area
        zoomBorders(1) = 1;
        zoomBorders(2) = 1;
        zoomBorders(3) = size(currIm{1},1);
        zoomBorders(4) = size(currIm{1},2);
        %Correct values outside borders
        zoomValXY(3) = min([zoomValXY(3) zoomBorders(4)]);
        zoomValXY(4) = min([zoomValXY(4) zoomBorders(3)]);
        zoomValXY(1) = max([zoomValXY(1) zoomBorders(2)]);
        zoomValXY(1) = min([zoomValXY(1) (zoomBorders(4) - zoomValXY(3) + 1)]);
        zoomValXY(2) = max([zoomValXY(2) zoomBorders(1)]);
        zoomValXY(2) = min([zoomValXY(2) (zoomBorders(3) - zoomValXY(4) + 1)]);
        zoomValXY = round(zoomValXY);
        %Adjust for display of complete pixels at borders
        zoomVal(xySel(2),:) = zoomValXY([1,3]);
        zoomVal(xySel(1),:) = zoomValXY([2,4]);
        %Set zoom axes
        if customDimScale
            zoomPxXY([1 2]) = scaledLocation(zoomValXY([1 2]));
            zoomPxXY([3 4]) = zoomValXY([3 4]) .* diff(dimScale(xySel([2 1]),:),1,2)' ./ (dim(xySel([2 1]))-1);
            axis(zoomAx, [zoomPxXY(1)-pxWidth(xySel(2))/2, zoomPxXY(1)+zoomPxXY(3)-pxWidth(xySel(2))/2, zoomPxXY(2)-pxWidth(xySel(1))/2, zoomPxXY(2)+zoomPxXY(4)-pxWidth(xySel(1))/2]);
        else
            axis(zoomAx, [zoomValXY(1)-0.5, zoomValXY(1)+zoomValXY(3)-.5, zoomValXY(2)-0.5, zoomValXY(2)+zoomValXY(4)-.5]);
        end
        %Set zoom region in image axes
        if get(tbWin(1), 'Value') == 1
            if customDimScale
                set(zoomReg, 'Position',zoomPxXY-[pxWidth(xySel([2 1])) 0 0]/2);
            else
                set(zoomReg, 'Position',zoomValXY-[.5 .5 0 0]);
            end
        end
        % Update slider and etxt
        for ii=1:length(dim)
            set(etxt_down(ii), 'String', num2str(zoomVal(ii,1)));
            set(etxt_up(ii), 'String', num2str(zoomVal(ii,2)+zoomVal(ii,1)-1));
            set(sld_down(ii), 'Value', zoomVal(ii,1));
            set(sld_up(ii), 'Value', zoomVal(ii,2)+zoomVal(ii,1)-1);
            set(txt_zoomWidth(ii), 'String', [' = ' num2str(zoomVal(ii,2))]);
            if customDimScale
                set(txt_zoomWidthScale(ii), 'String', [' = ' num2str(pxWidth(ii)*zoomVal(ii,2),'%6.2f') ' ' dimUnits{ii}]);
                pxSc = scaledLocation(zoomVal(ii,1),ii);
                set(txt_zoomDownScaled(ii), 'String',  num2str(pxSc,'%6.1f'));
                pxSc = scaledLocation(zoomVal(ii,2)+zoomVal(ii,1)-1,ii);
                set(txt_zoomUpScaled(ii), 'String',  num2str(pxSc,'%6.1f'));
            end
        end
        % Update plotXLim Values
        for ii=1:nDim
            if get(tb_lockPlots2Zoom(ii), 'Value')
                if customDimScale
                    pl(1) = scaledLocation(zoomVal(ii,1),ii);
                    pl(2) = scaledLocation(zoomVal(ii,1)+zoomVal(ii,2)-1,ii);
                    plotXLimScale(ii,:) = [pl(1) pl(2)];%[pl(1);pl(1)+pl(2)-pxWidth(ii)];
                end
                plotXLim(ii,:) = [zoomVal(ii,1);(zoomVal(ii,1)+zoomVal(ii,2)-1)];
            end
        end
        if get(tbHist, 'Value')
            updateHist;
        end
        if projMethod && get(bt_zoomProj, 'Value') && nargin > 2 && any(projDim == dimNum)
            if projMethod == 6
                drawImages;
                if get(tbRoi, 'Value')
                    drawRois;
                end
                drawObjects
            else
                updateImages;
            end
        else
            %             updateColormap;
            updateObjects; % Adjust length of position lines in zoom window to new zoom setting
        end
        if get(tbSwitchRGB, 'Value') && (get(cmStretchRGBMean, 'Value') || get(cmStretchRGBMax, 'Value')) && get(tb_lockPlots2Zoom(rgbDim), 'Value')
            updateImages;
        end
        if get(tbPlotsXLim, 'Value') && get(tbWin(3),'Value')
            updatePlots;
        end
        if strcmp('Zoom', get(get(bg_colormap, 'SelectedObject'),'String'))
            updateColormap;
        end
        % If matVis is called to display 2D Hist of another matVis
        % instance, send zoom settings to this matVis (zoom settings here
        % correspond to colormap settings in 'parent'-matVis)
        if isMatVis2DHist && ~calledFrom2DHist
            calledFrom2DHist = 1;
            if customDimScale
                sendFrom2DHist('cmMinMax', [zoomPxXY(2) (zoomPxXY(2)+zoomPxXY(4)-pxWidth(1));
                    zoomPxXY(1) (zoomPxXY(1)+zoomPxXY(3)-pxWidth(2))]);
            else
                sendFrom2DHist('cmMinMax', [(zoomValXY(1)-1)/256*(xVal2DHIst(2)-xVal2DHIst(1))+xVal2DHIst(1) ...
                    (zoomValXY(1)+zoomValXY(3)-1)/256*(xVal2DHIst(2)-xVal2DHIst(1))+xVal2DHIst(1);...
                    (zoomValXY(2)-1)/256*(yVal2DHIst(2)-yVal2DHIst(1))+yVal2DHIst(1) ...
                    (zoomValXY(2)+zoomValXY(4)-1)/256*(yVal2DHIst(2)-yVal2DHIst(1))+yVal2DHIst(1)]);
            end
        else
            calledFrom2DHist = 0;
        end
        if movdata.rec
            vidWriter('append');
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function setWinScale(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        % Set the factor used to scale the data pixels to screen pixels
        currWinScale = str2double(get(varargin{1}, 'String'));
        set(bt_100pct, 'String',[num2str(currWinScale) '%']);
        set100Pct;
        if debugMatVis, debugMatVisFcn(2); end
    end
    function adjustGuiSize(guiHandle,pwr)
        if debugMatVis, debugMatVisFcn(1); end
        switch os(1:2)
            case 'MA'  %mac
                % Reset size of uicontrols on Macs
                    scFct = macScaleFactor;
            case 'GL'  % Linux
                % Rescale fonts under linux
                w = findobj(gui,'-property','String');
                set(w,'FontSize',8);
                set(txt_titles, 'FontSize', 12);
                set(sepLine, 'FontSize', 6);
                scFct = 1;
            case 'PC' % PC
                % Rescale GUI size under windows if  'ScreenPixelsPerInch' varies from default
                % (96)
                scFct =  get(0, 'ScreenPixelsPerInch')/96*[1 1];
        end
        if ~any(scFct == 1)
            if nargin == 2
                scFct = scFct.^pwr;
            end
            handleList = findobj(guiHandle,'-property','Position');
            for ii=1:numel(handleList)
                p = get(handleList(ii),'Position');
                if numel(p) == 3
                    set(handleList(ii), 'Position',p(1:2).*scFct);
                else
                    set(handleList(ii), 'Position',p.*[scFct scFct]);
                end
            end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function set100Pct(varargin) %#ok
        if debugMatVis, debugMatVisFcn(1); end
        set(zoomWin, 'ResizeFcn','');
        for ii=1:nMat
            winPos = get(zoomWin(ii), 'Position');
            if winPos(2)+round(zoomValXY(4)*currWinScale/100) > scnSize(4)
                winPos(2) = scnSize(4)-round(zoomValXY(4)*currWinScale/100)-30;
            end
            set(zoomWin(ii), 'Position', [winPos(1) winPos(2) round(screenSizeScaling*currWinScale/100*zoomValXY(3)) round(screenSizeScaling*currWinScale/100*zoomValXY(4))]);
        end
        drawnow
        set(zoomWin, 'ResizeFcn',{@resizeWin, 'zoom'});
        figure(zoomWin(1));
        resizeWin(0,0,'zoom');
        if debugMatVis, debugMatVisFcn(2); end
    end
    function scalingFct = screenSizeScaling
        % Determine and return scaling due to OS settings (e.g. Windows
        % font scaling affects the way Matlab reports screen resolution)
        % So far only tested on Windows!
        if debugMatVis, debugMatVisFcn(1); end
        ge = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment;
        gd = ge.getDefaultScreenDevice;
        actualScreensize = [gd.getDisplayMode.getWidth gd.getDisplayMode.getHeight];
        matlabScreensize = get(0,'Screensize');
        scalingFct = matlabScreensize(3)/actualScreensize(1);
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateAspRatio(hObject, event)    %#ok
        if debugMatVis, debugMatVisFcn(1); end
        if myGcf == tempWin  %Called from aspect ratio window
            aspRatio = eval(get(hObject,'String'));
            set([imAx zoomAx], 'DataAspectRatio', [aspRatio 1]);
            set(tbAspRatio, 'String', [num2str(aspRatio(1)),':',num2str(aspRatio(2))]);
            close(tempWin);
            tempWin = [];
            s = get(projMethodPop, 'String');
            if strcmp(s{get(projMethodPop, 'Value')}, 'Tile')
                drawImages;
                drawObjects;
            end
            %Called by left click on button: switch between aspect ratio and
            %fill axes
        elseif get(tbAspRatio, 'Value') == 1
            axis(zoomAx, 'equal',  'off',  'ij');
            axis(imAx, 'equal', 'tight',  'ij');
            set([imAx zoomAx], 'DataAspectRatio', [aspRatio 1]);
            set(tbAspRatio, 'String', [num2str(aspRatio(1)),':',num2str(aspRatio(2))]);
        else
            axis(zoomAx,  'normal', 'off',  'ij');
            axis(imAx,  'normal', 'tight', 'ij');
            set(tbAspRatio, 'String', 'fill');
        end
        updateZoom
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updatePlotLim(hObject, event, xy)    %#ok
        if debugMatVis, debugMatVisFcn(1); end
        set(tempWin, 'HandleVisibility', 'on');
        switch xy
            case 'x'
                etxtFields = findobj('Parent', tempWin, 'Style', 'Edit');
                etxtFields = flipdim(etxtFields, 1);
                for ii=1:nDim
                    plotXLim(ii,:) = eval(get(etxtFields(ii), 'String'));
                end
            case 'y'
                etxtFields = findobj('Parent', tempWin, 'Style', 'Edit');
                etxtFields = flipdim(etxtFields, 1);
                for ii=1:nMat
                    plotYLim(ii,:) = eval(get(etxtFields(ii), 'String'));
                end
        end
        set(tempWin, 'HandleVisibility', 'off');
        if customDimScale
            for ii=1:nDim
                pl(1) = scaledLocation(plotXLim(ii,1),ii);
                pl(2) = scaledLocation(plotXLim(ii,2),ii);
                plotXLimScale(ii,:) = [pl(1);pl(2)];
            end
        end
        updatePlots;
        if debugMatVis, debugMatVisFcn(2); end
    end

%Show / Hide Windows
    function windowVisibility(hObject, event, numWin)    %#ok
        if debugMatVis, debugMatVisFcn(1); end
        if nargin == 0
            if get(tbWin(1), 'Value') && strcmp(get(imageWin(1), 'Visible'), 'off')
                set(imageWin, 'Visible', 'on');
            elseif ~get(tbWin(1), 'Value')
                set(imageWin, 'Visible', 'off');
            end
            if get(tbWin(2), 'Value') && strcmp(get(zoomWin(1), 'Visible'), 'off')
                set(zoomWin, 'Visible', 'on');
            elseif ~get(tbWin(2), 'Value')
                set(zoomWin, 'Visible', 'off');
            end
            if get(tbWin(3), 'Value') && strcmp(get(plotWin(1), 'Visible'), 'off')
                set(plotWin, 'Visible', 'on');
            elseif ~get(tbWin(3), 'Value')
                set(plotWin, 'Visible', 'off');
            end
        elseif numWin == 0
            set(myGcf, 'Visible', 'off');
            if numel(data) == 1
                if any(imageWin == myGcf)
                    set(tbWin(1), 'Value', 0);
                    set(imageWin, 'Visible', 'off');
                elseif any(zoomWin == myGcf)
                    set(tbWin(2), 'Value', 0);
                    set(zoomWin, 'Visible', 'off');
                else
                    set(tbWin(3), 'Value', 0);
                    set(plotWin, 'Visible', 'off');
                end
            end
        else
            onOff = get(tbWin(numWin), 'Value');
            if onOff == 1
                if numWin == 1
                    drawImages;
                    drawObjects;
                    if get(tbRoi, 'Value')
                        drawRois;
                    end
                    set(imageWin, 'Visible', 'on');
                elseif numWin == 2
                    drawImages;
                    drawObjects;
                    if get(tbRoi, 'Value')
                        drawRois;
                    end
                    set(zoomWin, 'Visible', 'on');
                else
                    drawPlots;
                    set(plotWin, 'Visible', 'on');
                end
            else
                if numWin == 1
                    set(imageWin, 'Visible', 'off');
                elseif numWin == 2
                    set(zoomWin, 'Visible', 'off');
                else
                    set(plotWin, 'Visible', 'off');
                end
            end
        end
        figure(gui);
        if debugMatVis, debugMatVisFcn(2); end
    end

%Show / hide objects (lines and rectangles)
    function toggleShowObjects(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if get(tbShowObjects, 'Value')
            set(zoomReg, 'Visible', 'on');
            set(lineHorIm, 'Visible', 'on');
            set(lineVertIm, 'Visible', 'on');
            set(lineHorZoom, 'Visible', 'on');
            set(lineVertZoom, 'Visible', 'on');
            %drawObjects;
        else
            set(zoomReg, 'Visible', 'off');
            set(lineHorIm, 'Visible', 'off');
            set(lineVertIm, 'Visible', 'off');
            set(lineHorZoom, 'Visible', 'off');
            set(lineVertZoom, 'Visible', 'off');
        end
        if projMethod == 6
            set(zoomReg, 'Visible', 'on');
            set([lineHorIm lineVertIm lineHorZoom lineVertZoom], 'Visible', 'off');
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
%% Features
%Show Window displaying parameter of customTif files
    function showTifPar(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        function hideTifPar(varargin)
            set(tifParFig, 'Visible', 'off');
            set(tbTifPar, 'Value', 0);
        end
        if isempty(tifParFig)
            tifParFig = figure('name', ['CustomTif Parameter',' (',varName{1},')'], 'Number', 'off',...
                'CloseRequestFcn', @hideTifPar, 'HandleVisibility', 'off');
            set(tifParFig, 'HandleVisibility', 'on');
            axis off;
            %Dimension
            v = struct2cell(par{1}.dim);
            n = fieldnames(par{1}.dim);
            text(0,1  ,'Dimension', 'FontWeight', 'Bold');
            for ii = 1:5
                text(0,1-0.05*ii,[n{ii},':']);
                text(0.15,1-0.05*ii,num2str(v{ii}), 'HorizontalAlignment', 'right');
            end
            %Offset
            v = struct2cell(par{1}.offset);
            n = fieldnames(par{1}.offset);
            text(0.25,1  ,'Offset', 'FontWeight', 'Bold');
            for ii = 1:3
                text(0.25,1-0.05*ii,[n{ii},':']);
                text(0.4,1-0.05*ii,num2str(v{ii}), 'HorizontalAlignment', 'right');
            end
            %Scale
            v = struct2cell(par{1}.scale);
            n = fieldnames(par{1}.scale);
            text(0.5,1  ,'Scale', 'FontWeight', 'Bold');
            for ii = 1:3
                text(0.5,1-0.05*ii,[n{ii},':']);
                text(0.65,1-0.05*ii,num2str(v{ii}), 'HorizontalAlignment', 'right');
            end
            %Excitation
            v = struct2cell(par{1}.config.ex);
            n = fieldnames(par{1}.config.ex);
            text(0,0.65,'Excitation', 'FontWeight', 'Bold');
            for ii = 1:3
                text(0,0.65-0.05*ii,[n{ii},':']);
                text(0.15,0.65-0.05*ii,num2str(v{ii}), 'HorizontalAlignment', 'right');
            end
            %Emission
            v = struct2cell(par{1}.config.em);
            n = fieldnames(par{1}.config.em);
            text(0.25,0.65,'Emission', 'FontWeight', 'Bold');
            for ii = 1:3
                text(0.25,0.65-0.05*ii,[n{ii},':']);
                text(0.4,0.65-0.05*ii,num2str(v{ii}), 'HorizontalAlignment', 'right');
            end
            %Aquisition
            v = struct2cell(par{1}.config.aq);
            n = fieldnames(par{1}.config.aq);
            text(0.5,0.65,'Aquisition', 'FontWeight', 'Bold');
            for ii = 1:5
                text(0.5, 0.65-0.05*ii,[n{ii},':']);
                text(0.7,0.65-0.05*ii,num2str(v{ii}), 'HorizontalAlignment', 'right');
            end
            for ii = 6:10
                text(0.8, 0.65-0.05*(ii-5),[n{ii},':']);
                text(1  ,0.65-0.05*(ii-5),num2str(v{ii}), 'HorizontalAlignment', 'right');
            end
            %Temperature
            text(0, 0.45, 'Temperature', 'FontWeight', 'Bold');
            text(0, 0.4, 'Temp:');
            text(0.15, 0.4, num2str(par{1}.config.temp), 'HorizontalAlignment', 'right');
            %Objective
            text(0.25, 0.45, 'Objective', 'FontWeight', 'Bold');
            text(0.25, 0.4, 'obj:');
            text(0.4, 0.4, num2str(par{1}.config.objective), 'HorizontalAlignment', 'right');
            %Description
            text(0, 0.3, 'Description', 'FontWeight', 'Bold', 'VerticalAlignment', 'top');
            figSize = get(tifParFig, 'Position');
            txtDescr = uicontrol(tifParFig, 'Position', [0.125*figSize(3) 10 0.7*figSize(3) figSize(4)*0.3-10], 'Style','edit','String',par{1}.description,...
                'HorizontalAlignment', 'left');
            set(tifParFig, 'HandleVisibility', 'off');
            % Save button
            uicontrol(tifParFig, 'Position', [0.3*figSize(3) figSize(4)*0.3+7 100 20], 'Style', 'pushbutton', 'Callback', @saveCTif, 'String','Save description');
        end
        if get(tbTifPar, 'Value')
            set(tifParFig, 'Visible', 'on');
        else
            set(tifParFig, 'Visible', 'off');
        end
        if debugMatVis
            display([repmat(' ',[1 debugIndent*fctLevel]) num2str(fctCount.(ST)) ': End   showTifPar']);
            fctLevel = fctLevel-1;
        end
        function saveCTif(varargin)
            par{1}.description = get(txtDescr, 'String');
            writeCustomTif(data{1}, par{1},  [filePath varName{1}]);
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Display and update histograms in histogramm window
    function updateHist(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if get(tbHist, 'Value') && ~withAlpha
            if nargin == 3 % Called as button down callback
                if strcmp(get(histAx(1), 'YScale'), 'log')
                    set(histAx(1), 'YScale','lin');
                else
                    set(histAx(1), 'YScale','log');
                end
            else
                if isinteger(data{1})
                    cI = single(currIm{1});
                else
                    cI = currIm{1};
                end
                cZ = cI(round(zoomValXY(2)):zoomValXY(2)+zoomValXY(4)-1,round(zoomValXY(1)):zoomValXY(1)+zoomValXY(3)-1);
                histVal(1  ,:) = histcounts(cI(:),histXData) * prod(dim) / dim(xySel(1)) / dim(xySel(2));
                histVal(2,:) = histcounts(cZ(:), histXData) /zoomValXY(3)/zoomValXY(4)*prod(dim);
                set(histAx, 'YLim',[0 max(histVal(:))]);
                set(histObj([1 3]), 'YData', [1 max(histVal(:))]);
                set(histPlots(1), 'YData', histVal(1  ,:));
                set(histPlots(2), 'YData', histVal(2,:));
            end
            updateHistObjects;
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

% Update line objects in histogram window
    function updateHistObjects
        if debugMatVis, debugMatVisFcn(1); end
        if get(tbHist, 'Value') && ~withAlpha
            set(histObj(1), 'XData', [cmMinMax(currContrastSel  ,1) cmMinMax(currContrastSel  ,1)]);
            set(histObj(3), 'XData', [cmMinMax(currContrastSel  ,2) cmMinMax(currContrastSel  ,2)]);
            transferVal(1  ,:) = linspace(single(cmMinMax(currContrastSel,1)),single(cmMinMax(currContrastSel,2)),255);
            transferVal(2,:) = ((transferVal(1  ,:)-single(cmMinMax(currContrastSel,1)))/(single(cmMinMax(currContrastSel,2)-cmMinMax(currContrastSel,1)))).^currGamma(currContrastSel);
            set(histObj(2), 'XData', transferVal(1  ,:), 'YData', transferVal(2,:)/max(transferVal(2,:))*max(histVal(:)));
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function showHist(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if isempty(globalHist)
            updateGlobalHist;
        end
        if nargin == 3 || ~get(tbHist, 'Value')
            set(histWin, 'Visible', 'off');
            set(tbHist, 'Value', 0);
        else
            set(histWin, 'Visible', 'on');
            set(tbHist, 'Value', 1);
            updateHist;
        end
        calculatingGlobalHist = 0;
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateSldLim(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if nargin && any(varargin{1} == [sldLimMin sldLimMax])
            sldLimVal(currContrastSel,:) = [str2double(get(sldLimMin, 'String')) str2double(get(sldLimMax,'String'))];
        else
            set(sldLimMin, 'String', num2str(sldLimVal(currContrastSel,1),'%6.3f'));
            set(sldLimMax, 'String', num2str(sldLimVal(currContrastSel,2),'%6.3f'));
        end
        histXData = linspace(sldLimVal(currContrastSel,1)  ,sldLimVal(currContrastSel,2),numel(histXData));
        if get(sldMin, 'Value') < histXData(1)
            set(sldMin, 'Value', histXData(1));
        end
        if get(sldMax, 'Value') > histXData(end)
            set(sldMax, 'Value', histXData(end));
        end
        set([sldMin sldMax], 'Min', histXData(1), 'Max', histXData(end));
        set([histAxPlot contrastSldIm], 'XData', histXDataBin);
        set([contrastAx histAxGui], 'XLim', [histXData(1) histXData(end)]);
        set(contrastSldIm, 'CData', linspace(histXData(1),histXData(end),255));
        updateColormap;
        updateGuiHistVal;
        %        start(timerHist);
        if get(tbHist, 'Value')
            updateGlobalHist;
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateGlobalHist(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        set(tbHist, 'Enable','off');
        % Delete timer object
        %         stop(timerHist);
        % Calculate global hist
        calculatingGlobalHist = 1;
        set(histWin, 'Name','Calculating Histogram ...');
        drawnow;
        cla(histAx(1));
        cla(histAx(2));
        %xx = histXData(:)';
        %binwidth = [diff(xx) 0];
        %xx = [xx(1)-binwidth(1)/2 xx+binwidth/2];
        % Shift bins so the interval is ( ] instead of [ ).
        %bins = xx + eps(xx);
        ii = 1;
        dataInt = 500000;
        globalHist = zeros(ceil(prod(dim)/dataInt),numel(histXData)-1);
        while  ((ii-1)*dataInt <=  prod(dim))
            clear w;
            set(tbHist, 'CData',[],'String', [num2str(min(100,round(100*ii*dataInt/prod(dim)))) '%']);drawnow;
            globalHist(ii,:) = histcounts(data{1}((ii-1)*dataInt+1:min(prod(dim),ii*dataInt))',histXData);
            ii = ii + 1;
            if shuttingDown  % Leave loop when matVis is being closed
                if debugMatVis, debugMatVisFcn(2); end
                return
            end
        end
        globalHist = (sum(globalHist,1))';
        % Combine first bin with 2nd bin and last bin with next to last bin
        %globalHist(2,:) = globalHist(2,:)+globalHist(1,:);
        %globalHist(end-1,:) = globalHist(end-1,:)+globalHist(end,:);
        %globalHist = globalHist(2:end-1,:);
        if isinteger(data{1})
            cI = single(currIm{1});
        else
            cI = currIm{1};
        end
        cZ = cI(round(zoomValXY(2)):floor(zoomValXY(2)+zoomValXY(4)-1),round(zoomValXY(1)):floor(zoomValXY(1)+zoomValXY(3)-1));
        histVal(1  ,:) = histcounts(cI(:),histXData) * prod(dim) / dim(xySel(1)) / dim(xySel(2));
        histVal(2,:) = histcounts(cZ(:),histXData)/zoomValXY(3)/zoomValXY(4)*prod(dim);
        histVal(3,:) = globalHist;
        histPlots = stairs(histXDataBin, histVal','Parent',histAx(1));
        set(histPlots(1), 'Color','c');
        set(histPlots(2), 'Color','m');
        set(histPlots(3), 'Color','k');
        %         transferVal(1  ,:) = histXDataBin;
        %         transferVal(2,:) = histXData.^currGamma;
        %         transferVal(2,:) = (transferVal(2,:)-min(transferVal(2,:)))/(max(transferVal(2,:))-min(transferVal(2,:)))*max(histVal(:));
        %         blackPt = find(transferVal(1  ,:)>histXData(1),1);
        %         whitePt = find(transferVal(1  ,:)>=histXData(end));
        %         transferVal(2,1:blackPt) = 0;
        %         transferVal(2,whitePt+1:end) = 0;
        %         histObj(1) = line(get(sldMin,'Value')*[1 1],[1 max(histVal(:))], 'Color','b','Parent',histAx(2));
        %         histObj(2) = line(transferVal(1,blackPt:whitePt),transferVal(2,blackPt:whitePt), 'Color','g','Parent',histAx(2));
        %         histObj(3) = line(get(sldMax,'Value')*[1 1],[1 max(histVal(:))], 'Color','r','Parent',histAx(2));
        histObj(1) = line([1 1],[1 max(histVal(:))], 'Color','b','Parent',histAx(2));
        histObj(2) = line(1,1, 'Color','g','Parent',histAx(2));
        histObj(3) = line([1 1],[1 max(histVal(:))],'Color','r','Parent',histAx(2));
        set(histAx, 'YLim',[0 max(histVal(:))], 'XLim',[histXData(1) histXData(end)]);
        legend(histAx(1),{'Current image (normalized)';'Zoom region (normalized)';'Complete data set'});
        %         title(histAx(1),{'Current image (normalized, cyan)';'Zoom region
        %         (normalized, magenta)';'Complete data set (black)'});
        set(histAx(2), 'Position',get(histAx(1),'Position'));
        set(histWin, 'HandleVisibility', 'off','Name','Histogram');
        if get(tbHist  ,'Value')
            set(histWin, 'Visible','on');
        else
            set(histWin, 'Visible','off');
        end
        set(tbHist, 'Enable','on','Tooltipstring','Show / hide histogram','String','','CData',histIcon);
        updateHistObjects;
        % Sort data for percentile function, called only once
        %         if isempty(dataPerc)
        %             sortData;
        %         end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function update2DHist(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        % Create and update 2D histogram
        if with2DHist && ishandle(matVis2DHist.figHandles.gui) && nargin && varargin{1} == tbHist && get(tbHist,'Value')  %  2Dhist exists already but is hidden: make visible
            set([matVis2DHist.figHandles.gui matVis2DHist.figHandles.imageWin], 'Visible','on');
        elseif with2DHist && ishandle(matVis2DHist.figHandles.gui) && nargin && varargin{1} == tbHist && ~get(tbHist,'Value')   % 2Dhist exists already and is visible: hide
            set([matVis2DHist.figHandles.gui matVis2DHist.figHandles.imageWin], 'Visible','off');
        elseif   ishandle(tbHist) && get(tbHist,'Value') % Update and/or create 2D hist
            % Create 2d hist from data and alpha map
            % 2D Histogram implementation from André Zeug
            dd = cat(2,currImVal{1}(:),currAlphaMapVal{1}(:));
            % if no weight is given, it is set to 1
            if size(dd,2)==2
                dd(:,3)=currAlphaMapVal{1}(:);
                %                 dd(:,3)=1;
            end
            dd(any(isnan(dd),2),:)=[];
            sz = size(dd);
            x_2DHist = linspace(minVal(1), maxVal(1),100);
            y_2DHist = linspace(minVal(nMat+1), maxVal(nMat+1),100);
            data2DHist = zeros(length(x_2DHist),length(y_2DHist));
            % fast method for lin spacing only
            dx = (x_2DHist(end)-x_2DHist(1))/(length(x_2DHist)-1);
            dy = (y_2DHist(end)-y_2DHist(1))/(length(y_2DHist)-1);
            l_x = length(x_2DHist);
            l_y = length(y_2DHist);
            for ii = 1:sz(1)
                xi = min(max(round((dd(ii,1)-x_2DHist(1))/dx)+1, 1), l_x);
                yi = min(max(round((dd(ii,2)-y_2DHist(1))/dy)+1, 1), l_y);
                data2DHist(xi,yi) = data2DHist(xi,yi) + dd(ii,3);
            end
            if ~with2DHist || ~ishandle(matVis2DHist.figHandles.gui) % 2D hist does not exist yet
                outStrct = exportOutputStrct;
                gp = get(gui,'Position');
                wp.gui = [gp(1) gp(2)-gp(4)-130-40]; % 130 is height of tooltip box
                wp.image = [gp(1)+gp(3)+20 gp(2)-gp(4)-40 gp(4) gp(4) ];
                %                 matVis2DHist =
                %                 matVis(data2DHist,'dimScale',[get(sldMax,'Min') get(sldMax,'Max');minVal(2) maxVal(2)],'dimNames',{'Data values';'Alpha values'},'matNames',{[varName{1} ': 2D Hist']},'startPar',{'aspRatio';0;'windowPositions';wp;'objVisibility';1;'windowVisibility';[1 0 0];'calledAs2DHist';{outStrct.fctHandles.loadNewSettings;double([min(dd(:,1)) max(dd(:,1))]);double([min(dd(:,2)) max(dd(:,2))])}});
                zv(1,1) = find(x_2DHist>=cmMinMax(1,1),1) ;
                zv(1,2) = find(x_2DHist>=cmMinMax(1,2),1)-find(x_2DHist>=cmMinMax(1,1),1)+1;
                zv(2,1) = find(y_2DHist>=cmMinMax(nMat+1,1),1) ;
                zv(2,2) = find(y_2DHist>=cmMinMax(nMat+1,2),1)-find(y_2DHist>=cmMinMax(nMat+1,1),1)+1;
                matVis2DHist = matVis(data2DHist,'dimScale',[minVal(1) maxVal(1);minVal(2) maxVal(2)],...
                                                 'dimNames',{'Data values';'Alpha values'},...
                                                 'matNames',{[varName{1} ': 2D Hist']},...
                                                 'startPar',{'zoomVal';zv;'aspRatio';0;'windowPositions';wp;'objVisibility';1;'windowVisibility';[1 0 0];'calledAs2DHist';{outStrct.fctHandles.loadNewSettings;double([min(dd(:,1)) max(dd(:,1))]);double([min(dd(:,2)) max(dd(:,2))])}});
                with2DHist = 1;
            else
                matVis2DHist.fctHandles.loadNewData(data2DHist);
            end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

% Function to sort data in order to provide percentage-based contrast or
% histrogram operations: NOT YET IMPLEMENTED
%     function sortData
%         tmp = sort(data{1}(:));
%         for ii=1:100
%             dataPerc(ii) = tmp(round(numel(tmp)*ii/100));
%         end
%         clear tmp
%     end

%Switch between configurations for gui and windows
    function switchConfig(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        set(btConfig, 'UserData', mod(get(btConfig, 'UserData')+1  ,2));
        switch get(btConfig, 'UserData')
            % Never used the "current config" option, thus disabled
            % (also made problems with the gui size when switching between tooltip and compact mode)
            %             case 0
            %                 setConfig(currConfig);
            %                 set(btConfig, 'String', 'CurrConf');
            case 0
                %                 currConfig = getCurrConfig;
                setConfig(customConfig);
                set(btConfig, 'String', 'CustConf');
            case 1
                setConfig(defaultConfig);
                set(btConfig, 'String', 'DefConf');
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Retrieve configuration from current state
    function currConfig = getCurrConfig
        if debugMatVis, debugMatVisFcn(1); end
        % Window Properties
        %Window Visibility
        currConfig.winVis.imageWin = get(tbWin(1), 'Value');       %Default: 1
        currConfig.winVis.zoomWin  = get(tbWin(2), 'Value');       %Default: 1
        currConfig.winVis.plotWin  = get(tbWin(3), 'Value');       %Default: 1
        %Window Position
        %For one data set
        currConfig.winPos.gui = get(gui, 'Position');
        %         currConfig.winPos.gui(3:4) = defaultConfig.winPos.gui(3:4);
        for ii=1:nMat
            currConfig.winPos.imageWin(ii,:) = get(imageWin(ii), 'Position');       %Default: [337, 545,  450, 450];
            currConfig.winPos.zoomWin(ii,:)  = get(zoomWin(ii), 'Position');       %Default: [800, 545,  450, 450];
            currConfig.winPos.plotWin(ii,:)  = get(plotWin(ii), 'Position');       %Default: [5,    10, 1250, 500];
        end
        
        % Image and Zoom Window Options
        %Aspect Ratio 1:1
        currConfig.aspectRatio = get(tbAspRatio, 'Value');           %Default: 1
        %Colorbar Display
        currConfig.colorbar = get(tbColorbar,'Value');              %Default: 0
        %Colormap
        currConfig.colormap = get(popLut, 'Value');              %Default: 1 (gray)
        %Gamma
        currConfig.gamma  = currGamma;
        %RGB Mode
        currConfig.RGB = get(tbSwitchRGB, 'Value');
        %Colormap Mode (Global, Local or Manual)
        currConfig.colormapMode = get(get(bg_colormap, 'SelectedObject'), 'String');   %Default: 'Global'
        %Lines Visibility
        currConfig.lineVis = get(tbShowObjects, 'Value');               %Default: 1
        %Menu Bar
        currConfig.menuBarVis = get(tbMenuBars, 'Value');            %Default: 0
        % Plot Options
        %Plot Dimensions
        currConfig.plotDim = plotDim;           %Default: [1 2]
        %Plot Zoom
        currConfig.plotZoom = get(tbPlotsXLim, 'Value');              %Default: 0
        %Scale Plot
        currConfig.plotScale = get(tbPlotsYLim, 'Value');       %Default: 0
        %Marker
        currConfig.marker = get(tbMarker, 'Value');                %Default: 0
        %Average
        currConfig.plotMean = get(btMean, 'UserData');              %Default: 0 (no averaging)
        % Tooltips
        currConfig.tooltips = get(tbTooltips, 'Value');   %Default: 1 (display tooltips)
        % Histogram update during playback
        currConfig.playHist = get(tb_playHist, 'Value');  % Default: 0 (don't update)
        % Histogram update during change of position indicators
        currConfig.moveHist = get(tb_moveHist, 'Value');  % Default: 0 (don't update)
        % Link figure size / position
        currConfig.linkFigSize = get(tbLinkWin, 'Value'); % Default: 1 (link figures)
        if debugMatVis, debugMatVisFcn(2); end
    end

%Set configuration
    function setConfig(config)
        if debugMatVis, debugMatVisFcn(1); end
        % Window Properties
        guiPos = get(gui,'Position');
        if (~config.tooltips && get(tbTooltips, 'Value')) || config.tooltips && ~get(tbTooltips, 'Value')
            shiftGui = 1;
        else
            shiftGui = 0;
        end
        set(gui, 'Position', [config.winPos.gui(1) config.winPos.gui(2)+shiftGui*(config.winPos.gui(4)-guiPos(4)) guiPos(3) guiPos(4)]);  %
        % Resize GUI if necessary
        if get(tbTooltips, 'Value') && ~config.tooltips
            set(tbTooltips, 'Value',0);
            toggleTooltipDisplay;
        elseif ~get(tbTooltips, 'Value') && config.tooltips
            set(tbTooltips, 'Value',1);
            toggleTooltipDisplay;
        end
        %Window Visibility
        set(tbWin(1), 'Value', config.winVis.imageWin);      %Default: 1
        set(tbWin(2), 'Value', config.winVis.zoomWin);       %Default: 1
        set(tbWin(3), 'Value', config.winVis.plotWin);       %Default: 1
        %Window Position
        %Menu Bar
        set(tbMenuBars, 'Value', config.menuBarVis);            %Default: 0
        toggleMenuBars;
        %For one data set
        if monSizeMin(1) <= min([config.winPos.imageWin(:,1); config.winPos.zoomWin(:,1); config.winPos.plotWin(:,1)]) && ...
          monSizeMin(2) <= min([config.winPos.imageWin(:,2); config.winPos.zoomWin(:,2); config.winPos.plotWin(:,2)]) && ...
          monSizeMax(1) >= max(sum([config.winPos.imageWin(:,[1 3]);config.winPos.zoomWin(:,[1 3]);config.winPos.plotWin(:,[1 3])],2)) && ...
          monSizeMax(2) >= max(sum([config.winPos.imageWin(:,[2 4]);config.winPos.zoomWin(:,[2 4]);config.winPos.plotWin(:,[2 4])],2))
          for ii=1:nMat
            set(imageWin(ii), 'Position', config.winPos.imageWin(ii,:));      %Default: [337, 545,  450, 450];
            set(zoomWin(ii),  'Position', config.winPos.zoomWin(ii,:));       %Default: [800, 545,  450, 450];
            set(plotWin(ii),  'Position', config.winPos.plotWin(ii,:));       %Default: [5,    10, 1250, 500];
          end
        end
        % Image and Zoom Window Options
        %Aspect Ratio 1:1
        set(tbAspRatio, 'Value', config.aspectRatio);            %Default: 1
        %Colorbar Display
        set(tbColorbar,'Value', config.colorbar);                 %Default: 0
        %Colormap
        set(popLut, 'Value', config.colormap);              %Default: 1 (gray)
        %Gamma
        currGamma = config.gamma(1);
        if nMat > numel(currGamma) ||withAlpha
            if withAlpha
                currGamma = currGamma*ones(1,2*nMat);
            else
                currGamma = currGamma(1)*ones(1,nMat);
            end
        end
        %RGB Mode
        if config.RGB && ~withAlpha
            set(cmImage, 'String', 'Channel');
            set(tbSwitchRGB, 'Value', 1);
            switchRGB;
        end
        %Colormap Mode (Global, Local or Manual)
        if oldMATLAB
            set(bg_colormap, 'SelectedObject',...
                strcmp(get(get(bg_colormap,'Children'),'String'),config.colormapMode)'*get(bg_colormap,'Children'));   %Default: 'Global'
        else
          if ~(config.RGB && withAlpha)
            set(bg_colormap, 'SelectedObject',bg_colormap.Children(strcmp(get(bg_colormap.Children,'String'),config.colormapMode)))
          end
        end
        %Lines Visibility
        set(tbShowObjects, 'Value', config.lineVis);               %Default: 1
        
        % Plot Options
        %Plot Dimensions
        plotDim = config.plotDim;           %Default: [1 2]
        %Plot Zoom
        set(tbPlotsXLim, 'Value', config.plotZoom);              %Default: 0
        %Scale Plot
        set(tbPlotsYLim, 'Value', config.plotScale);       %Default: 0
        %Marker
        set(tbMarker, 'Value', config.marker);                %Default: 0
        %Average
        set(btMean, 'UserData', config.plotMean) ;              %Default: 0 (no averaging)
        set(tbTooltips, 'Value', config.tooltips);         %Default: 1 (display tooltips)
        % Update of histogram
        set(tb_playHist, 'Value',config.playHist);
        set(tb_moveHist, 'Value',config.moveHist);
        % Lin figure size / position
        set(tbLinkWin, 'Value', config.linkFigSize);
        %         toggleMoveHist
        %Update
        set(cbPlots, 'Value', 0);
        set(cbPlots(plotDim), 'Value', 1);
        windowVisibility;
        updateImages;
        updateColormap;
        showColorbar;
        drawPlots;
        cbCallback;
        toggleShowObjects;
        updateObjects;
        if debugMatVis, debugMatVisFcn(2); end
    end

%Save configuration in matVis directory for future starts of matVis
    function saveConfig(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        compName = strtrim(getenv('Computername'));
        customConfig = getCurrConfig;
        [p,f,e] = fileparts(which('matVis.m'));    %#ok
        try
            save(fullfile(p,['matVisConfig_' compName '.mat']), 'customConfig');
        catch
            save(fullfile(p,['matVisConfig.mat']), 'customConfig');
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%Data export to workspace
    function exportData(varargin)
      if debugMatVis, debugMatVisFcn(1); end
      exportWin = figure('Units', 'Pixel', 'Position', [300 300 200 25*nDim+150], 'Name', 'Export Data',...
        'MenuBar', 'none', 'Resize', 'off', 'NumberTitle', 'off', 'HandleVisibility', 'off');
      uicontrol(exportWin, 'Style', 'Text', 'Position', [30 25*nDim+105 140 40], ...
        'String','Specify intervals of data to be exported!','FontWeight', 'bold',...
        'BackgroundColor', get(exportWin, 'Color'), 'HorizontalAlignment', 'center');
      uicontrol(exportWin, 'Style', 'Text', 'Position', [30 25*nDim+75 140 40], ...
        'String','Use '':'' for colon operator.',...
        'BackgroundColor', get(exportWin, 'Color'), 'HorizontalAlignment', 'center');
      for ii = 1:nDim
        uicontrol(exportWin, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(exportWin, 'Color'),...
          'Position', [30 25*nDim-(ii-2)*25+40 30 17], 'String', dimNames(ii), ...
          'HorizontalAlignment', 'left' );
        exportTxt(ii) = uicontrol(exportWin, 'Style', 'Edit', 'Units', 'Pixel', ...
          'Position', [90 25*nDim-(ii-2)*25+40 80 20], 'String', ['[' num2str(zoomVal(ii,1)) ':' num2str(sum(zoomVal(ii,:))-1) ']']);    %#ok
      end
      uicontrol(exportWin, 'Style', 'Text', 'Position', [30 55 50 17], 'String','Var Name',...
        'BackgroundColor', get(exportWin, 'Color'), 'HorizontalAlignment', 'left');
      exportName = uicontrol(exportWin, 'Style', 'Edit',  'String', 'matVisExport',...
        'Position', [90 55 80 20]);
      exportBt = uicontrol(exportWin, 'Style', 'Pushbutton', 'Units', 'Pixel', ...
        'Position', [30 25 140 20], 'String', {'Export to Workspace'},'Callback',@exportNow);
      % newMatVisBt:
      uicontrol(exportWin, 'Style', 'Pushbutton', 'Units', 'Pixel', ...
        'Position', [30 5 140 20], 'String', 'Open in new matVis','Callback',@newMatVis);
      if debugMatVis, debugMatVisFcn(2); end
      function exportNow(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        ind = []; indC = false(1,nDim);
        for iii = 1:nDim
          ind{iii} = eval(get(exportTxt(iii), 'String'));    %#ok
          indC(iii) = length(ind{iii})==1;
        end
        assignin('base', [get(exportName, 'String'),'_index'], ind);
        expProp = [];lEA = 1;
        if ~isempty(dimScale)
          expProp{lEA} = 'dimScale';
          for iii=1:length(ind)
            expProp{lEA+1}(iii,:) = dimScale(iii,1) + diff(dimScale(iii,:))/(dim(iii)-1)*(ind{iii}([1 end])-1);
          end
          expProp{lEA+1}(indC,:) = [];
          lEA = lEA+2;
        end
        if ~isempty(dimNames)
          expProp{lEA} = 'dimNames';
          expProp{lEA+1} = dimNames;
          expProp{lEA+1}(indC) = [];
          lEA = lEA+2;
        end
        if ~isempty(dimUnits)
          expProp{lEA} = 'dimUnits';
          expProp{lEA+1} = dimUnits;
          expProp{lEA+1}(indC) = [];
          lEA = lEA+2;
        end
        if ~isempty(allNames)
          expProp{lEA} = 'matNames';
          expProp{lEA+1} = {allNames};
          lEA = lEA+2;
        end
        assignin('base', [get(exportName, 'String'),'_props'], expProp);
        if nMat == 1
          d = data{1}(ind{:});
          assignin('base', get(exportName, 'String'), d);
        else
          for iii = 1:nMat
            d = data{iii}(ind{:});
            if isempty(varName{iii})
              assignin('base', [get(exportName, 'String') '_set' num2str(iii,'%02d')] , d);
            else
              if fromFile
                [pn vn] = fileparts(varName{iii});
              else
                vn = varName({iii});
              end
              assignin('base', [get(exportName, 'String') '_' vn] , d);
            end
          end
        end
        exportCount = exportCount + 1;
        clear d;
        delete(exportWin);
        exportWin = [];
        set(btExport, 'Value', 0);
        if debugMatVis, debugMatVisFcn(2); end
      end
      function newMatVis(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        ind = []; indC = false(1,nDim);
        for iii = 1:nDim
          ind{iii} = eval(get(exportTxt(iii), 'String'));    %#ok
          indC(iii) = length(ind{iii})==1;
        end
        for iii = 1:numel(data)
          inputArg{iii} = data{iii}(ind{:});
        end
        if withAlpha
          inputArg{end+1} = 'alphaMap';
          inputArg{end+1} = alphaMap{1}(ind{:});
        end
        if ~isempty(dimScale)
          inputArg{end+1} = 'dimScale';
          lIA = length(inputArg)+1;
          for iii=1:length(ind)
            inputArg{lIA}(iii,:) = dimScale(iii,1) + diff(dimScale(iii,:))/(dim(iii)-1)*(ind{iii}([1 end])-1);
          end
          inputArg{lIA}(indC,:) = [];
        end
        if ~isempty(dimNames)
          inputArg{end+1} = 'dimNames';
          inputArg{end+1} = dimNames;
          inputArg{end}(indC) = [];
        end
        if ~isempty(dimUnits)
          inputArg{end+1} = 'dimUnits';
          inputArg{end+1} = dimUnits;
          inputArg{end}(indC) = [];
        end
        if ~isempty(allNames)
          inputArg{end+1} = 'matNames';
          inputArg{end+1} = {allNames};
        end
        %         for iii = 1:numberOptionalArgs
        %           switch  optionalArgIdentifier{iii}
        %             case 'alphaMap'
        %               inputArg{end+1} = 'alphaMap';
        %               inputArg{end+1} = alphaMap{1}(ind{:});
        %             case 'dimNames'
        %               inputArg{end+1} = 'dimNames';
        %               inputArg{end+1} = dimNames;
        %               inputArg{end+1}(indC) = [];
        %             case 'dimUnits'
        %               inputArg{end+1} = 'dimUnits';
        %               inputArg{end+1} = dimUnits;
        %               inputArg{end+1}(indC) = [];
        %             case 'matNames'
        %               inputArg{end+1} = 'matNames';
        %               inputArg{end+1} = allNames;
        %               inputArg{end+1}(indC) = [];
        %             case 'startPar'
        %               inputArg{end+1} = 'startPar';
        %               inputArg{end+1} = startPar;
        %             case 'dimScale'
        %               inputArg{end+1} = 'dimScale';
        %               lIA = length(inputArg)+1;
        %               for n = 1:length(ind)
        %                 inputArg{lIA}(n,:) = dimScale(n,1) + diff(dimScale(n,:))/(dim(n)-1)*(ind{iiii}([1 end])-1);
        %               end
        %               inputArg{lIA}(indC,:) = [];
        %           end
        %         end
        matVis(inputArg{:});
        delete(exportWin);
        exportWin = [];
        if debugMatVis, debugMatVisFcn(2); end
      end
   end
%% Line-profile manager
%     nProfiles: actual number of profiles
%     The structure roiList contains information for all Rois:
%         profileList( ).number:        number of profile; displayed as text (profileText) in image/zoom
%                                   window. The number of a profile doesn't change even if other profiles are deleted from the
%                                   list. If a new profile is created, it gets the number max([profileList.number])+1.
%         profileList( ).name:          Name of profile. Also never changes unless changed by user. The default includes the number
%                                   of the profile for "normal" profiles.
%         profileList( ).rectangle:     Coordinates of rectangle surrounding profile
%                                   (left,down, right, up). If custom
%                                   dimension scales are used, the values
%                                   are specified in terms of the scales
%                                   (not the pixels). However, upon
%                                   export/import, these values are
%                                   translated to/from pixel values!
%         profileList( ).corners:   "Corners" of profile, i.e. "clicked" positions  
%         profileList( ).coordinates:       Coordinates of polygon defining profile (interpolated from the corners to correspond to distances of 1 pixel). If custom
%                                   dimension scales are used, the values
%                                   are specified in terms of the scales
%                                   (not the pixels). However, upon
%                                   export/import, these values are
%                                   translated to/from pixel values!
%         profileList( ).index.x/.y:    x/y index of all pixels of profile
%         profileList( ).settings.xySel xy dimension during time of profile selection
%         profileList( ).settings.position Position during profile selection
%         profileList( ).settings.zoomVal ZoomWindow settings during profile selection


%   During export, the profileList structure is exported or saved, during
%   import, this structure is loaded and the profiles are reconstructed from
%   it. Thus it is necessary that the chosen file contains a structure of
%   this kind.
%   When changing xySel (and not simply exchanging them), all profile
%   information will be lost and the profile Manager will be reset (resetprofiles).

%%  Build profile Gui during first call or make it visible/unvisible during
%   later calls or from WindowCloseRequestFcn
    function profileGui(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        %Make profile Gui and all related objects invisible
        if nargin == 3 || ~get(tbProfile,'Value')
            set(profileWin, 'Visible','off');
            set(btMean, 'Enable','on');
            set(cbPlots(xySel),'Enable', 'on');
            set(tbProfile,'Value',0)
            if nProfiles >0
                set([profileLine.im(:)', profileLine.zoom(:)', profileText.im(:)', profileText.zoom(:)'],  'Visible', 'off');
            end
            set([tb_newProfile tbProfileShift tbProfileRotate tbProfileScale profileBtReplace], 'Value', 0);
            set([imageWin zoomWin],'WindowButtonMotionFcn',@mouseMotion,...
                'WindowButtonDownFcn',@buttonDownCallback,...
                'WindowButtonUpFcn','');
            drawPlots;
            if debugMatVis, debugMatVisFcn(2); end
            return
        end
        %Make profile Gui and all related objects visible
        set(profileWin, 'Visible','on');
        set(btMean, 'Enable','off');
        if nProfiles > 0 && get(profileBtShowNames, 'Value')
            set([profileLine.im(:)', profileLine.zoom(:)', profileText.im(:)', profileText.zoom(:)'], 'Visible', 'on');
        elseif nProfiles > 0
            set([profileLine.im(:)', profileLine.zoom(:)', profileText.im(:)', profileText.zoom(:)'], 'Visible', 'off');
        end
        set(cbPlots(xySel),  'Value', 0 ,'Enable', 'off'); %necessary when call from WindowCloseRequestFcn
        plotSel(xySel) = 0;
        drawPlots;
        %   drawProfiles;
        profileStruct.trace.inwork = 0;
        %Create profile GuiWindow (only during first call)
        if isempty(profileWin)
            gp = get(gui,'Position');
            % Profile trace window
            profileTraceWin = figure('Position',[ gp(1)+140  gp(2)-452   560   420], 'Name', 'Profile view', 'MenuBar', 'none', 'NumberTitle', 'off');
            %profile Manager Window with text
            profileWin =  figure('Units', 'Pixel', 'Name', 'Profile Manager',...
                'MenuBar', 'none', 'NumberTitle', 'off','CloseRequestFcn', {@profileGui,0}, 'HandleVisibility', 'off',...
                'WindowStyle','normal', 'Position', [gp(1) gp(2)-362 130 330]);
            warning('off','MATLAB:HandleGraphicresizeRoiGuis:ObsoletedProperty:JavaFrame');
            profile_jf = get(profileWin,'JavaFrame');
            profile_jf.setFigureIcon(javax.swing.ImageIcon(im2java(uint8(icon_profile))));
            warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
            txt_profileList = uicontrol(profileWin, 'Style', 'Text', 'Position', [10 305 110 15], ...
                'String','List of profiles','FontWeight', 'bold',...
                'BackgroundColor', get(profileWin, 'Color'), 'HorizontalAlignment', 'center');
            %Listbox for profiles
            profileListbox = uicontrol(profileWin, 'Style', 'listbox', 'Units', 'Pixel', ...
                'BackgroundColor', get(profileWin, 'Color'),'Value',0,...
                'Position', [10 153 110 147], 'String', '','Max',1, 'Callback', @updateProfileSelection,...
                'Tag', 'List of profiles. Multiple selections are not (yet) possible. The Names of the profiles can be changed (see buttons below), the numbers in front of the names are generated by matVIs and can not be edited.');
            %Panel for display of profile properties: probably not
            %necessary
%             profilePropertiesPanel = uipanel(profileWin,'Units','Pixel','Position', [130,25,180,275],'Title', 'Selected profile',...
%                 'BackgroundColor',get(profileWin,'Color'),'FontWeight', 'bold','TitlePosition','centertop',...
%                 'Tag','Shows the currently selected profile.');
%             profileName = uicontrol('Parent',profilePropertiesPanel, 'Style', 'Text', 'Position', [5 210 170 30], ...
%                 'BackgroundColor', get(profileWin, 'Color'), 'HorizontalAlignment', 'center','FontWeight', 'bold',...
%                 'String', 'profile Name');
            %Axes for display of profile
            profileStruct.gui.profile.ax  = subplot(4,1,[1 3],'Parent',profileTraceWin);
            profileStruct.gui.traceim.ax  = subplot(4,1,4,'Parent',profileTraceWin);
            %set(profileStruct.gui.traceim.ax,'Color','k'); % required for alphaMap
            %Add Poly profile Button
            tb_newProfile = uicontrol(profileWin, 'Style', 'Togglebutton','Position', [10,125,32,22],...
                'CData', roiPolygonBt,'Callback',{@getNewProfile,'one'},...
                'ToolTipString', 'Add polygonal/free hand profile. Left click for single profile (''push button''), right click for multiple selections (''toggle button'') ',...
                'Tag', 'Add polygonal/free hand profile. Draw straight lines using short clicks and free hand selections by keeping the mouse button pressed. End selection by right or double click. Left click for single profile, right click for multiple selections');
            %Copy profile Button
            tb_copyProfile = uicontrol(profileWin, 'Style', 'Togglebutton','Position', [49,125,32,22],...
                'String', 'Copy','Callback',@copyProfile,'Enable','off',...
                'ToolTipString', 'Copy currently selected profile as new profile.',...
                'Tag', 'Copy currently selected profile as new profile.');
            
            %Shift profile Button
            tbProfileShift = uicontrol(profileWin, 'Style', 'Togglebutton','Position', [88,125,32,22], 'Enable', 'off',...
                'CData', roiShiftBt,'Callback',@shiftprofile,...
                'ToolTipString', ' Move selected profile with mouse or keyboard. Use control, shift or alt to move by 5, 10 or 20 px, respectively. ',...
                'Tag', 'Move selected profile(s) with mouse or keyboard. Use control, shift or alt to move by 5, 10 or 20 px, respectively.  Note that this is a toggle button.',...
                'Value',0);
           
            %Replace profile
            profileBtReplace = uicontrol(profileWin, 'Style', 'Togglebutton','Position', [10,100,32,22],...
                'String', 'Replace','Tag', 'Replace the selected profile. The position in the list and the name will be preserved. Note that this is a toggle button.',...
                'FontSize', 6, 'Enable', 'off','ToolTipString','Replace selected profile by new selection');
            
            %Rename profile Button
            profileBtRename = uicontrol(profileWin, 'Style', 'Pushbutton','Position', [49,100,32,22],...
                'String', 'Rename','Callback', @renameProfile,'ToolTipString', ' Rename selected profile ',...
                'FontSize', 6, 'Enable', 'off','ToolTipString','Rename Selected profile','Tag','Rename the selected profile.');
            
            %Show / hide profile names button
            profileBtShowNames = uicontrol(profileWin, 'Style', 'Togglebutton','Position', [88,100,32,22],...
                'String', 'ShowNames','Callback', @showProfileNames,'ToolTipString', ' Show or hide profile names in figure windows ',...
                'Tag', 'Show or hide the profile names in the figure windows.','FontSize', 6, 'Value', 1,'Enable','off');
            
            %Delete profile Button
            bt_deleteProfile = uicontrol(profileWin, 'Style', 'Pushbutton','Position', [10,75,32,22],...
                'CData',roiDeleteBt,'Callback', @deleteProfile,'ToolTipString', 'Delete selected profile(s).','Tag', 'Delete the selected profile(s).', 'Enable', 'off');
            %Export profile Button
            bt_profileExport = uicontrol(profileWin, 'Style', 'Pushbutton','Position', [49,75,32,22],...
                'CData', roiExportBt,'Callback', @exportProfiles,'Enable','on','ToolTipString', ' Export selected profile(s) ',...
                'Tag', 'Export the selected profile(s) to the Matlab workspace or to a mat-file.', 'Enable', 'off');
            %Import profile Button
            uicontrol(profileWin, 'Style', 'Pushbutton','Position', [88,75,32,22],...
                'CData', roiImportBt,'Callback', @importProfiles,'Enable','on','ToolTipString', ' Import profiles ',...
                'Tag', 'Import profiles.');
            if nDim>2
                %Popup for Dimension for profile data export
                n = dimNames;
                n(xySel) = [];
                profilePopDataExport = uicontrol(profileWin, 'Style', 'Popup','Position', [10,49,71  ,22],...
                    'Enable','on','String',n,...
                    'Tag', 'Select dimension for profile data export.');
                %Export profile Data Button
                bt_exportprofileData = uicontrol(profileWin, 'Style', 'Pushbutton','Position', [88,50,32,22],...
                    'CData', roiExportDataBt,'Callback', @exportprofileData,'Enable','on','ToolTipString', ' Export profile data along specified dimension ',...
                    'Tag', 'Export profile data along specified dimension into the workspace.', 'Enable', 'off');
            end
            profileStruct.show.width = 5;
            edt_lineWidt = uicontrol(profileWin, 'Style', 'Edit','Position', [88,25,32,22],...
                'String', num2str(profileStruct.show.width),'Callback', @changeProfileWidth,'Enable','off','ToolTipString', ' Line width in pixels ',...
                'Tag', 'Line width in pixels (only integers allowed)');
            pop_profileProjection= uicontrol(profileWin, 'Style', 'Popup','Position', [10,25,71  ,22],...
                    'Enable','on','String',{'Mean','Sum','Max','Min'},...
                    'Tag', 'Select projection method.','Callback',@imgTraceUpdate);
            uicontrol(profileWin, 'Style', 'Text','Position', [10,0,60,22],'String','Line length: ','BackgroundColor',get(profileWin,'color'), 'HorizontalAlignment','left');
            txt_profileLength = uicontrol(profileWin, 'Style', 'Text','Position', [70,0,50,22],'String','--','BackgroundColor',get(profileWin,'color'), 'HorizontalAlignment','left');
            %Fill profile Axes with current image
            set(profileWin, 'HandleVisibility', 'on','ResizeFcn', @resizeProfileWin,'WindowbuttonMotionFcn',@mouseMotion);
%             if customDimScale
%                 profileImage = imagesc(dimScale(xySel(2),:),dimScale(xySel(1),:),currIm{1});
%             else
%                 profileImage = imagesc(currIm{1});
%             end
%             axis image;
            set(profileStruct.gui.profile.ax,'FontSize',8);
            profileLine.profile = line(0,0,'Parent', profileStruct.gui.profile.ax,'Color','w');
            set(profileWin, 'HandleVisibility', 'off');
            %Set first entry of profile list (so it is not empty - will be
            %replaced by 1 when first profile is selected)
            profileList(1).number = 0;%stopHere
            updateColormap;
            set(profileWin, 'ResizeFcn',[]);adjustGuiSize(profileWin);set(profileWin,'ResizeFcn', @resizeProfileWin);
            % The following fields are set here temporarily and might have
            % to be implemented dynamically
            profileStruct.show.bw = true;  % should be false for 3d data!?
            profileStruct.trace.algor = 'pchip';
            profileStruct.gui.tracewidth = 1; % has to be implemented as edit field
            
            profileStruct.imgsz = size(currIm{1});
            profileStruct.show.chsz = 1;
            profileStruct.show.scale = 'abs';
            profileStruct.show.colors = [1 0 0]; % see SPECTRUM functions
            profileStruct.instance = 1; % not clear what this is used for
            profileStruct.gui.bgcolor = [.9 .9 .9];
        end
        function resizeProfileWin(varargin)
          if debugMatVis, debugMatVisFcn(1); end
            adjustGuiSize(profileWin,-1);
            newPos = get(profileWin, 'Position');
            if any(newPos(3:4) < [130 250])
                newPos = [newPos(1) newPos(2) 130 250];
                set(profileWin, 'Position', newPos);
            end
            set(profileListbox, 'Position', [10 153 newPos(3)-20 newPos(4)-183]);
            set(txt_profileList, 'Position', [10 newPos(4)-25 newPos(3)-20 15]);
%             set(profilePropertiesPanel, 'Position', [130 25 newPos(3)-140 newPos(4)-35]);
            set(profileName, 'Position', [5 newPos(4)-90 newPos(3)-150 30]);
%             set(profileStruct.gui.profile.ax, 'Position', [25 80 newPos(3)-175 newPos(4)-160]);
            adjustGuiSize(profileWin);
          if debugMatVis, debugMatVisFcn(2); end
        end
        function getNewProfile(varargin)
          if debugMatVis, debugMatVisFcn(1); end
            if nargin == 3 && strcmp(varargin{3},'one')
                if get(tb_newProfile, 'Value') %If previously selected button is not pressed
%                     set(profileStruct.gui.settrace.h,'Value',0);set(findobj('tag',sprintf('trace_%.0f',profileStruct.instance)),'Enable','off');
                    set([zoomWin imageWin],'Pointer','crosshair');
                    set([zoomWin imageWin],'WindowButtonDownFcn',@profileButtonDown);
                    set([zoomWin imageWin],'WindowButtonUpFcn',[])
                    set([zoomWin imageWin],'WindowButtonMotionFcn',[])
                    set([zoomWin imageWin],'WindowButtonUpFcn',[]);
                    if currProfile
                        set([profileLine.im profileLine.zoom], 'LineWidth',1,'Color','w');              %'Color', 'b');
                        set([profileText.im profileText.zoom], 'Color','w','FontWeight','normal')
                    end
                else % otherwise: unpress and quit
                    set(tb_newProfile, 'Value', 0);
                    set([zoomWin imageWin],'WindowButtonMotionFcn',@mouseMotion);
                    set([zoomWin imageWin],'WindowButtonDownFcn',@buttonDownCallback);
                    set([zoomWin imageWin],'WindowButtonUpFcn','');
                    set([zoomWin imageWin],'WindowButtonUpFcn',@zoomKeyReleaseFcn);
                    if currProfile
                        set(profileLine.im(:,currProfile,:), 'LineWidth',1,'Color','r');              %'Color', 'b');
                        set(profileLine.zoom(:,currProfile,:), 'LineWidth',1,'Color','r');            %'Color', 'b');
                        set([profileText.im(:,currProfile,:) profileText.zoom(:,currProfile,:)], 'Color','r','FontWeight','bold')
                    end
                end
            end
          if debugMatVis, debugMatVisFcn(2); end
        end
        function profileButtonDown(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            mousepos = get(myGca,'CurrentPoint');
            mousepos = mousepos(1,1:2);
%             currUnits = get(myGca, 'units');
%             set(myGca,'Units','Pixels');
%             axpos = get(myGca,'Position');
%             set(myGca,'Units',currUnits);
%             Xpos = get(myGca,'XLim');
%             Ypos = get(myGca,'YLim');
%             Xval = interp1([axpos(1) axpos(1)+axpos(3)],Xpos,mousepos(1),'linear');
%             Yval = interp1([axpos(2) axpos(2)+axpos(4)],fliplr(Ypos),mousepos(2),'linear');
%             mousepos = round([Xval Yval]);
            switch get(myGcf,'SelectionType')
                case 'normal'  % add points
                    set(myGcf,'WindowButtonMotionFcn',@imgTraceMove)
                    if ~profileStruct.trace.inwork
                        delete(findobj('tag',sprintf('traceplot_%.0f',profileStruct.instance)))
                        profileStruct.trace.points = [];
                        profileStruct.trace.points = mousepos;
                        profileStruct.trace.inwork = true;
                        profileStruct.trace.complete = false;
                        set(findobj('tag',sprintf('save_%.0f',profileStruct.instance)),'Enable','off');
%                         set(profileStruct.gui.resettrace.h,'Enable','on');
                    else
                        profileStruct.trace.points = [profileStruct.trace.points; mousepos];
                    end
                case 'alt'
                    if profileStruct.trace.inwork&&size(profileStruct.trace.points,1)>1
                        profileStruct.trace.points = profileStruct.trace.points(1:end-1,:);
                    end
                case 'open'
                    set(myGcf,'WindowButtonMotionFcn',[])
                    profileStruct.trace.inwork = false;
                    profileStruct.trace.complete = true;
                    set(findobj('tag',sprintf('save_%.0f',profileStruct.instance)),'Enable','on');
            end
            imgTraceUpdate([profileStruct.trace.points])
            if strcmp(get(myGcf,'SelectionType'),'open')
                addProfileToList([profileStruct.trace.points], profileStruct.trace.mat,'new');
                set(tb_newProfile, 'Value', 0);
                set([zoomWin imageWin],'WindowButtonMotionFcn',@mouseMotion);
                set([zoomWin imageWin],'WindowButtonDownFcn',@buttonDownCallback);
                set([zoomWin imageWin],'WindowButtonUpFcn','');
                set([zoomWin imageWin],'WindowButtonUpFcn',@zoomKeyReleaseFcn);
            end
            if debugMatVis, debugMatVisFcn(2); end
        end
        function addProfileToList(newPts,newMat,newOld)
          if debugMatVis, debugMatVisFcn(1); end
            nProfiles = nProfiles+1;
            numberProfile = nProfiles;
            profileList(numberProfile).number = numberProfile; %max([profileList.number])+1;
            profileList(numberProfile).name = sprintf('Prof %02d', nProfiles);
            profileList(numberProfile).points = newPts;
            profileList(numberProfile).smoothedPoints = imgTraceSmooth(newPts);
            profileList(numberProfile).mat = newMat;
            profileList(numberProfile).settings.xySel = xySel;
            profileList(numberProfile).settings.position = currPos;
            profileList(numberProfile).settings.zoomVal = zoomVal;
            for ii=1:nMat
                hold(imAx(ii),'on');
                profileLine.im(ii,numberProfile,1) = plot(profileList(numberProfile).smoothedPoints(:,1),profileList(numberProfile).smoothedPoints(:,2),'Color','w','tag',sprintf('traceplot_%.0f',profileStruct.instance),'Parent',imAx(ii),'LineWidth',2);
                profileLine.im(ii,numberProfile,2) = plot(profileList(numberProfile).mat(1,:,1),profileList(numberProfile).mat(1,:,2),'Color','w','tag',sprintf('traceplot_%.0f',profileStruct.instance),'Parent',imAx(ii),'LineStyle',':','LineWidth',2);
                profileLine.im(ii,numberProfile,3) = plot(profileList(numberProfile).mat(end,:,1),profileList(numberProfile).mat(end,:,2),'Color','w','tag',sprintf('traceplot_%.0f',profileStruct.instance),'Parent',imAx(ii),'LineStyle',':','LineWidth',2);
                profileText.im(ii,numberProfile) = text(max(profileList(numberProfile).smoothedPoints(:,1))+5,mean(profileList(numberProfile).smoothedPoints(:,2)),profileList(numberProfile).name,'Parent',imAx(ii),'FontWeight','bold');
                hold(zoomAx(ii),'on');
                profileLine.zoom(ii,numberProfile,1) = plot(profileList(numberProfile).smoothedPoints(:,1),profileList(numberProfile).smoothedPoints(:,2),'Color','w','tag',sprintf('traceplot_%.0f',profileStruct.instance),'Parent',zoomAx(ii),'LineWidth',2);
                profileLine.zoom(ii,numberProfile,2) = plot(profileList(numberProfile).mat(1,:,1),profileList(numberProfile).mat(1,:,2),'Color','w','tag',sprintf('traceplot_%.0f',profileStruct.instance),'Parent',zoomAx(ii),'LineStyle',':','LineWidth',2);
                profileLine.zoom(ii,numberProfile,3) = plot(profileList(numberProfile).mat(end,:,1),profileList(numberProfile).mat(end,:,2),'Color','w','tag',sprintf('traceplot_%.0f',profileStruct.instance),'Parent',zoomAx(ii),'LineStyle',':','LineWidth',2);
                profileText.zoom(ii,numberProfile) = text(max(profileList(numberProfile).smoothedPoints(:,1))+5,mean(profileList(numberProfile).smoothedPoints(:,2)),profileList(numberProfile).name,'Parent',zoomAx(ii),'FontWeight','bold');
            end
            profileStruct.trace.points = [];
            delete([profileStruct.trace.plot1 profileStruct.trace.plot2 profileStruct.trace.plot3]);
            if strcmp(newOld,'new')
                s = get(profileListbox, 'String');
                s{size(s,1)+1} = [num2str(size(get(profileListbox,'String'),1)+1  ,'%03d'),': ',profileList(numberProfile).name];
                if numberProfile == 1
                    set(profileListbox,'Value', numberProfile,'String', s);
                else
                    set(profileListbox,'String', s, 'Value', numberProfile);
                end
            end
            updateProfileSelection;
            set([bt_deleteProfile bt_profileExport edt_lineWidt],'Enable','on');
          if debugMatVis, debugMatVisFcn(2); end
        end
        function updateProfileSelection(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            currProfile = get(profileListbox, 'Value');
            %Update Profile lines and numbers in image/zoom windows
            set(profileLine.im, 'LineWidth',1,'Color','w');              %'Color', 'b');
            set(profileLine.zoom, 'LineWidth',1,'Color','w');            %'Color', 'b');
            set(profileLine.im(:,currProfile,:), 'Color','r');  % , 'LineWidth',3);  
            set(profileLine.zoom(:,currProfile,:), 'Color','r'); % , 'LineWidth',3); %
            set(profileText.im, 'FontWeight', 'normal','Color','w');     %'Color', 'b');
            set(profileText.zoom, 'FontWeight', 'normal','Color','w');   %'Color', 'b');
            set(profileText.im(:,currProfile), 'FontWeight', 'bold', 'Color','r');
            set(profileText.zoom(:,currProfile),'FontWeight', 'bold',  'Color','r');
            set(txt_profileLength, 'String', num2str(size(profileList(currProfile).smoothedPoints,1)));
            imgTraceUpdate;
            if debugMatVis, debugMatVisFcn(2); end
        end
        function imgTraceMove(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            mousepos = get(myGca,'CurrentPoint');
            mousepos = mousepos(1,1:2);
%             set(myGca,'Units','Pixels');
%             axpos = get(myGca,'Position');
%             set(myGca,'Units','Normalized');
%             Xpos = get(myGca,'XLim');
%             Ypos = get(myGca,'YLim');
%             Xval = interp1([axpos(1) axpos(1)+axpos(3)],Xpos,mousepos(1),'linear');
%             Yval = interp1([axpos(2) axpos(2)+axpos(4)],fliplr(Ypos),mousepos(2),'linear');
%             mousepos = round([Xval Yval])
            imgTraceUpdate([profileStruct.trace.points; mousepos])
            if debugMatVis, debugMatVisFcn(2); end
        end
        function changeProfileWidth(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            profileStruct.show.width = double(uint16(str2double(get(edt_lineWidt,'String'))));
            set(edt_lineWidt,'String', num2str(profileStruct.show.width));
            for iii = 1:numel(profileList)
                profileList(iii).smoothedPoints = imgTraceSmooth(profileList(iii).points);
                profileList(iii).mat = imgTraceMatrix(profileList(iii).smoothedPoints);
                for ii = 1:nMat
                    set([profileLine.im(ii,iii,2) profileLine.zoom(ii,iii,2)], 'XData',profileList(iii).mat(1,:,1),'Ydata',profileList(iii).mat(1,:,2));
                    set([profileLine.im(ii,iii,3) profileLine.zoom(ii,iii,3)], 'XData',profileList(iii).mat(end,:,1),'Ydata',profileList(iii).mat(end,:,2));
                end
            end
            imgTraceUpdate;
            if debugMatVis, debugMatVisFcn(2); end
        end
        function imgTraceUpdate(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            linpoints = [];
            if nargin == 1  % called with linpoints as input
                linpoints = varargin{1};
            elseif profileList(1).number > 0;  % called as function callback or without input
                linpoints = profileList(currProfile).points;
            end
            if ~isempty(linpoints) && numel(unique(linpoints(:,1)))>1
              tracepoints = imgTraceSmooth(linpoints);
              profileStruct.trace.mat = imgTraceMatrix(tracepoints);
              profileStruct.trace.img = [];
              if withAlpha
                profileStruct.trace.imgA = [];
              end
              for nt=1:profileStruct.show.chsz
                profileStruct.trace.img(:,:,nt) = interp2(double(currIm{1}(:,:,nt)),profileStruct.trace.mat(:,:,1),profileStruct.trace.mat(:,:,2));
                if withAlpha
                  profileStruct.trace.imgA(:,:,nt) = interp2(double(currAlphaMap{1}(:,:,nt)),profileStruct.trace.mat(:,:,1),profileStruct.trace.mat(:,:,2));
                end
              end
              if withAlpha
                profileStruct.trace.profile = squeeze(imgTraceProfile(profileStruct.trace.img,profileStruct.trace.imgA));
              else
                profileStruct.trace.profile = squeeze(imgTraceProfile(profileStruct.trace.img));
              end
              profileStruct.show.profile = profileStruct.trace.profile';
              %profileStruct.show.profile= profileStruct.show.profile(~isnan(profileStruct.show.profile))'; % added by Stephan check where the NaN values come from
              if regexp(profileStruct.show.scale,'norm')
                for np=1:size(profileStruct.trace.profile,2)
                  profileStruct.show.profile(:,np)=profileStruct.trace.profile(:,np)/max(profileStruct.trace.profile(:,np));
                end
              end
              profilesz = size(profileStruct.trace.img);
              if numel(profilesz) == 2 % inserted by sj, not clear whether correct, might be due to only single color channel
                profilesz(3) = 1;
              end
              if withAlpha
                profileStruct.trace.RGB = profileStruct.trace.img;
              else
                profileStruct.trace.RGB = reshape(profileStruct.trace.img,[prod(profilesz(1:end-1)) profilesz(end)]);
                for nt=1:profileStruct.show.chsz
                  profileStruct.trace.RGB(:,nt)=profileStruct.trace.RGB(:,nt)/max(profileStruct.trace.RGB(:,nt));
                end
                profileStruct.trace.RGB = profileStruct.trace.RGB*single(profileStruct.show.colors);
                profileStruct.trace.RGB = reshape(profileStruct.trace.RGB,[profilesz(1:end-1) 3]);
                profileStruct.trace.RGB(profileStruct.trace.RGB>1)=1;
                profileStruct.trace.RGB(profileStruct.trace.RGB<0)=0;
              end
              if isfield(profileStruct.trace, 'plot1') && ishandle(profileStruct.trace.plot1)
                set(profileStruct.trace.plot1,'XData',tracepoints(:,1),'YData',tracepoints(:,2))
                set(profileStruct.trace.plot2,'XData',profileStruct.trace.mat(1,:,1),'YData',profileStruct.trace.mat(1,:,2))
                set(profileStruct.trace.plot3,'XData',profileStruct.trace.mat(end,:,1),'YData',profileStruct.trace.mat(end,:,2))
                if profileStruct.show.bw
                  set(profileStruct.trace.profh{1},'YData',profileStruct.show.profile(:,1));
                else
                  for nt=1:profileStruct.show.chsz
                    set(profileStruct.trace.profh{nt},'YData',profileStruct.show.profile(:,nt),'Color',profileStruct.show.colors(nt,:)*.8);
                  end
                end
                if withAlpha
                  set(profileStruct.trace.imgh,'CData',profileStruct.trace.img,'AlphaData',profileStruct.trace.imgA)
                else
                  set(profileStruct.trace.imgh,'CData',profileStruct.trace.img);%profileStruct.trace.RGB)
                end
              else
                if profileStruct.trace.inwork  % create/update profile line only if a new line is currently selected (not when function is called from updateProfileSelection)
                  hold(myGca,'on');
                  profileStruct.trace.plot1 = plot(tracepoints(:,1),tracepoints(:,2),'Color','r','tag',sprintf('traceplot_%.0f',profileStruct.instance),'Parent',myGca,'LineWidth',2);
                  profileStruct.trace.plot2 = plot(profileStruct.trace.mat(1,:,1),profileStruct.trace.mat(1,:,2),'Color','r','tag',sprintf('traceplot_%.0f',profileStruct.instance),'Parent',myGca,'LineStyle',':','LineWidth',2);
                  profileStruct.trace.plot3 = plot(profileStruct.trace.mat(end,:,1),profileStruct.trace.mat(end,:,2),'Color','r','tag',sprintf('traceplot_%.0f',profileStruct.instance),'Parent',myGca,'LineStyle',':','LineWidth',2);
                end
                if profileStruct.show.bw
                  profileStruct.trace.profh{1} = plot(profileStruct.show.profile,'Parent',profileStruct.gui.profile.ax,'Color','k','tag',sprintf('traceplot_%.0f',profileStruct.instance));
                else
                  for nt=1:profileStruct.show.chsz
                    profileStruct.trace.profh{nt} = plot(profileStruct.show.profile(:,nt),'Parent',profileStruct.gui.profile.ax,'Color',profileStruct.show.colors(nt,:)*.8,'tag',sprintf('traceplot_%.0f',profileStruct.instance));
                  end
                end
                profileStruct.trace.imgh = imagesc(profileStruct.trace.img,'Parent',profileStruct.gui.traceim.ax,'tag',sprintf('traceplot_%.0f',profileStruct.instance));
                if withAlpha
                  set(profileStruct.trace.imgh,'AlphaData',profileStruct.trace.imgA);
                  set(profileStruct.gui.traceim.ax,'Color','k');
                end
              end
              set(profileStruct.gui.traceim.ax,'box','on', 'CLim', cmMinMax(1,:),...
                'XLim',[.5 size(tracepoints,1)-.5],'YLim',[.5 profileStruct.show.width+.5]);
              if profileStruct.show.width<6
                set(profileStruct.gui.traceim.ax,'YTick',1:profileStruct.show.width,...
                  'YTickLabel',(profileStruct.show.width-1)/2:-1:-1*(profileStruct.show.width-1)/2)
              else
                set(profileStruct.gui.traceim.ax,'YTick',[0 .5 1]*(profileStruct.show.width-1)+1,...
                  'YTickLabel',[(profileStruct.show.width-1)/2 0 -(profileStruct.show.width-1)/2])
              end
              set(profileStruct.gui.profile.ax,'box','on',...
                'XTickLabel',[],'Color',profileStruct.gui.bgcolor,...
                'XLim',[.5 size(tracepoints,1)-.5]);
              maxProfile = max(profileStruct.show.profile(:));  % Used to be nanmax
              minProfile = min(profileStruct.show.profile(:));   % Used to be nanmin
              if isfinite(maxProfile) && isfinite(minProfile) && maxProfile>minProfile
                set(profileStruct.gui.profile.ax,'YLim',[minProfile maxProfile]);
              end
            end
            if debugMatVis, debugMatVisFcn(2); end
        end
        function imgTraceRedraw(varargin)
            if debugMatVis > 1, debugMatVisFcn(1); end
            if isfield(profileStruct.trace,'points')&&~isempty(profileStruct.trace.points)
                imgTraceUpdate(profileStruct.trace.points);
            end
            if debugMatVis > 1, debugMatVisFcn(2); end
        end
        function smoothpoints = imgTraceSmooth(tracepoints)
            if debugMatVis > 1, debugMatVisFcn(1); end
            if numel(unique(tracepoints(:,1)))>1
                x = tracepoints(:,1);
                y = tracepoints(:,2);
                dist = [0; cumsum(sqrt((x(2:end)-x(1:end-1)).^2+(y(2:end)-y(1:end-1)).^2))];
                [~,ind,~] = unique(dist);
                dist = dist(ind);
                x = x(ind);
                y = y(ind);
                distt = (0:ceil(dist(end)))';
                xp = interp1(dist,x,distt,profileStruct.trace.algor);
                yp = interp1(dist,y,distt,profileStruct.trace.algor);
                dist = [0; cumsum(sqrt((xp(2:end)-xp(1:end-1)).^2+(yp(2:end)-yp(1:end-1)).^2))];
                distt = (0:ceil(dist(end)))';
                xp = interp1(dist,xp,distt,'linear','extrap');
                yp = interp1(dist,yp,distt,'linear','extrap');
                %                 % sj: Most of the time, the last yp value is nan, to avoid
                %                 % this:
                %                 nanIdx = isnan(yp);
                %                 xp(nanIdx) = [];
                %                 yp(nanIdx) = [];
                smoothpoints = [xp,yp];
            else
                smoothpoints = tracepoints;
            end
            if debugMatVis > 1, debugMatVisFcn(2); end
        end
        function perpendicularmat = imgTraceMatrix(tracepoints)
            if debugMatVis > 1, debugMatVisFcn(1); end
            xp = tracepoints(:,1);
            yp = tracepoints(:,2);
            % width of area to use for profile
%             profileStruct.show.width = get(profileStruct.gui.tracewidth.h,'Value');
            lw = profileStruct.show.width;
            w = (0:lw-1)-(lw-1)/2;
            % length of profile
            T = length(xp);
            t = 1:T;
            % calc slope
            t1 = max([t-1; ones(1,T)]);
            t2 = min([t+1; T*ones(1,T)]);
            p  = [xp(t) yp(t)];
            p1 = [xp(t1) yp(t1)];
            p2 = [xp(t2) yp(t2)];
            dp = p2-p1;
            clear t t1 t2 p p1 p2
            % norm slope
            dp = dp./repmat(sqrt(dp(:,1).^2+dp(:,2).^2),[1 2]);
            % perform matrix reshaping to allow rotation
            dp = dp';
            dp = dp(:)*w;
            dp = reshape(dp', lw,2,[]);
            dp = reshape(permute(dp,[1 3 2]),[],2);
            dp = dp*[0 1;-1 0]; % rotation vector
            % reshaping for x and y position
            xi = reshape(dp(:,1),lw,[])+repmat(xp',[lw 1]);
            yi = reshape(dp(:,2),lw,[])+repmat(yp',[lw 1]);
            perpendicularmat = cat(3,xi,yi);
            if debugMatVis > 1, debugMatVisFcn(2); end
        end
        
        function deleteProfile(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            profileSel = get(profileListbox, 'Value');
            profileList(profileSel) = [];
            s = get(profileListbox, 'String');
            s(profileSel) = [];
            for ii=1:size(s,1)
                s{ii}(1:3) = num2str(ii,'%03d');
            end
            set(profileListbox,'Value',1);
            set(profileListbox, 'String',s);
            delete(profileText.im(:,profileSel));
            delete(profileText.zoom(:,profileSel));
            delete(squeeze(profileLine.im(:,profileSel,:)));
            delete(squeeze(profileLine.zoom(:,profileSel,:)));
            profileText.im(:,profileSel) = [];
            profileText.zoom(:,profileSel) = [];
            profileLine.im(:,profileSel,:) = [];
            profileLine.zoom(:,profileSel,:) = [];
            nProfiles = nProfiles - size(profileSel,2);
            if nProfiles > 0
                set(profileListbox, 'Value',1);
                currProfile = 1;
                updateProfileSelection;
            else
                set(profileListbox, 'Value',0, 'String', '');
                set([profileBtRename tbProfileShift profileBtReplace bt_profileExport bt_deleteProfile bt_exportProfileData] , 'Enable', 'off');
                if ~isempty(bt_exportProfileData)
                    set([edt_lineWidt bt_exportProfileData], 'Enable','off')
                end
            end
            if debugMatVis, debugMatVisFcn(2); end
        end
        function exportProfiles(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            profileSel = get(profileListbox, 'Value');
            profileSel = 1:numel(profileList);  % This is only temporary since multiple selection in the listbox is currently not allowed. Thus all the profiles will be exported.
            if isempty(profileSel)
                warndlg('No ROIs selected. Select ROI first before exporting it/them.','No ROI available');
            else
                profileListExp = profileList;
                
                q = questdlg('Export profiles to workspace or save as .mat file?','Profile Export', 'To File', 'To Workspace', 'Cancel', 'To File');
                if strcmp(q,'To Workspace')
                    assignin('base', 'matVisProfileExport', profileListExp(profileSel));
                elseif strcmp(q, 'To File')
                    [f,p] = uiputfile('.mat','Choose folder and filename to save profiles!');
                    if f == 0
                      if debugMatVis, debugMatVisFcn(2); end
                      return
                    end
                    matVisProfileExport =  profileListExp(profileSel); %#ok
                    save([p,f],'matVisProfileExport');
                end
                clear matVisProfileExport profileListExp
                display('Profile export complete!');
            end
            if debugMatVis, debugMatVisFcn(2); end
        end
        function importProfiles(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            [f,p] = uigetfile('.mat','Choose file to load profiles!');
            if isequal(f,0)
              if debugMatVis, debugMatVisFcn(2); end
              return
            end
            try
                matVisProfileExport = cell2mat(struct2cell(load([p,f])));
            catch    %#ok
                errordlg('Chosen file does not contain profiles exported using matVis (''matVisProfileExport'')!');
                if debugMatVis, debugMatVisFcn(2); end
                return
            end
%             if customDimScale
%                 for ii=1:numel(matVisProfileExport)
%                     % Convert corners to scale
%                     for jj = 1:size(matVisProfileExport(ii).corners,2)
%                         matVisProfileExport(ii).corners(:,jj) = scaledLocation(matVisProfileExport(ii).corners(:,jj)');
%                     end
%                     % Convert rectangle to scale
%                     matVisProfileExport(ii).rectangle = [max(dimScale(xySel(2),1)+pxWidth(xySel(2))  ,-pxWidth(xySel(2))+min(matVisProfileExport(ii).corners(1,:))),...
%                         max(dimScale(xySel(1),1)+pxWidth(xySel(1))  ,-pxWidth(xySel(1))+min(matVisProfileExport(ii).corners(2,:))),...
%                         min(dimScale(xySel(2),2),pxWidth(xySel(2))+max(matVisProfileExport(ii).corners(1,:))),...
%                         min(dimScale(xySel(1),2),pxWidth(xySel(1))+max(matVisProfileExport(ii).corners(2,:)))];
%                 end
%             end
            if nProfiles  % nProfiles == 0 if there are no profiles existing
                q = questdlg('Append imported profiles or overwrite existing profiles?','Append or overwrite?','Overwrite','Append','Overwrite');
                if strcmp(q,'Append')
                    profileList(end+1:end+numel(matVisProfileExport)) = matVisProfileExport;
                else
                    profileList = matVisProfileExport;
                end
                delete([profileLine.im profileLine.zoom profileText.im profileText.zoom])
            else
                profileList = matVisProfileExport;
            end
            nProfiles = size(profileList,2);
            s = [];
            for ii=1:nProfiles
                s{ii} = profileList(ii).name;
            end
            s = [];
            set(profileListbox, 'String', '', 'Value', 1);
            for ii=1:nProfiles
                s{ii} = [num2str(ii,'%03d'),': ',profileList(ii).name];
            end
            for ii=1:nMat
                for pp = 1:numel(profileList)
                    hold(imAx(ii),'on');
                    profileLine.im(ii,pp,1) = plot(profileList(pp).smoothedPoints(:,1),profileList(pp).smoothedPoints(:,2),'Color','w','tag',sprintf('traceplot_%.0f',profileStruct.instance),'Parent',imAx(ii),'LineWidth',2);
                    profileLine.im(ii,pp,2) = plot(profileList(pp).mat(1,:,1),profileList(pp).mat(1,:,2),'Color','w','tag',sprintf('traceplot_%.0f',profileStruct.instance),'Parent',imAx(ii),'LineStyle',':','LineWidth',2);
                    profileLine.im(ii,pp,3) = plot(profileList(pp).mat(end,:,1),profileList(pp).mat(end,:,2),'Color','w','tag',sprintf('traceplot_%.0f',profileStruct.instance),'Parent',imAx(ii),'LineStyle',':','LineWidth',2);
                    profileText.im(ii,pp) = text(max(profileList(pp).smoothedPoints(:,1))+5,mean(profileList(pp).smoothedPoints(:,2)),profileList(pp).name,'Parent',imAx(ii),'FontWeight','bold');
                    hold(zoomAx(ii),'on');
                    profileLine.zoom(ii,pp,1) = plot(profileList(pp).smoothedPoints(:,1),profileList(pp).smoothedPoints(:,2),'Color','w','tag',sprintf('traceplot_%.0f',profileStruct.instance),'Parent',zoomAx(ii),'LineWidth',2);
                    profileLine.zoom(ii,pp,2) = plot(profileList(pp).mat(1,:,1),profileList(pp).mat(1,:,2),'Color','w','tag',sprintf('traceplot_%.0f',profileStruct.instance),'Parent',zoomAx(ii),'LineStyle',':','LineWidth',2);
                    profileLine.zoom(ii,pp,3) = plot(profileList(pp).mat(end,:,1),profileList(pp).mat(end,:,2),'Color','w','tag',sprintf('traceplot_%.0f',profileStruct.instance),'Parent',zoomAx(ii),'LineStyle',':','LineWidth',2);
                    profileText.zoom(ii,pp) = text(max(profileList(pp).smoothedPoints(:,1))+5,mean(profileList(pp).smoothedPoints(:,2)),profileList(pp).name,'Parent',zoomAx(ii),'FontWeight','bold');
                end
            end
            set(profileListbox, 'String',s,'Value', 1);
            set([profileBtRename tbProfileShift profileBtReplace bt_profileExport bt_deleteProfile bt_exportProfileData] , 'Enable', 'on');
            updateProfileSelection;
            if debugMatVis, debugMatVisFcn(2); end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function profile = imgTraceProfile(A, AA)
        if debugMatVis, debugMatVisFcn(1); end
        %profile = mean(A,1);
        val = get(pop_profileProjection,'Value');
        list = get(pop_profileProjection,'String');
        profileProj = list{val};
        sz1 = size(A);
        switch profileProj
            case 'Mean'
                if withAlpha
                    profile = sum(A.*AA,1, 'omitnan')./sum(AA.*~isnan(A),1, 'omitnan'); % weighted mean ratio
                else
                    profile = mean(A,1, 'omitnan'); % used to be nanmean
                end
            case 'Sum'
                if withAlpha
                    profile = sum(A.*AA,1, 'omitnan')./sum(AA.*~isnan(A),1, 'omitnan')*sz1(1); % weighted mean ratio % used to be nansum
                else
                    profile = sum(A,1, 'omitnan'); % used to be nansum
                end
            case 'Min'
                profile = min(A,[],1);
            case 'Max'
                if withAlpha
                    [~,AAInd] = max(AA, [], 1); %stopHere%squeeze() % used to be nanmax
                    AAInd  = sz1(1)*(0:sz1(2)-1) + AAInd(:)';%(1:prod(sz1(1:2)))' + prod(sz1(1:2)) * (AAInd(:)-1);
                    profile  = A(AAInd);
                else
                    profile = max(A,[],1); % used to be nanmax
                end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateProfileProperties
        if debugMatVis, debugMatVisFcn(1); end
        if debugMatVis, warning(sprintf('This function call is currently (more or less) a copy of ''imgTraceUpdate'',\nwhich is most likely not optimal. :-)')); end
        % This is currently (more or less) a copy of imgTraceUpdate, which most likely is
        % not optimal.
        linpoints = profileList(currProfile).points;
        if numel(unique(linpoints(:,1)))>1
            tracepoints = profileList(currProfile).smoothedPoints;
            profileStruct.trace.mat = profileList(currProfile).mat;
            profileStruct.trace.img = [];
            for nt=1:profileStruct.show.chsz
              profileStruct.trace.img(:,:,nt) = interp2(double(currIm{1}(:,:,nt)),profileStruct.trace.mat(:,:,1),profileStruct.trace.mat(:,:,2));
              if withAlpha
                profileStruct.trace.imgA(:,:,nt) = interp2(double(currAlphaMap{1}(:,:,nt)),profileStruct.trace.mat(:,:,1),profileStruct.trace.mat(:,:,2));
              end
            end
            if withAlpha
              profileStruct.trace.profile = squeeze(imgTraceProfile(profileStruct.trace.img,profileStruct.trace.imgA));
            else
              profileStruct.trace.profile = squeeze(imgTraceProfile(profileStruct.trace.img));
            end
            profileStruct.show.profile = profileStruct.trace.profile';
            %profileStruct.show.profile= profileStruct.show.profile(~isnan(profileStruct.show.profile))'; % added by Stephan check where the NaN values come from
            if regexp(profileStruct.show.scale,'norm')
                for np=1:size(profileStruct.trace.profile,2)
                    profileStruct.show.profile(:,np)=profileStruct.trace.profile(:,np)/max(profileStruct.trace.profile(:,np));
                end
            end
            profilesz = size(profileStruct.trace.img);
            if numel(profilesz) == 2 % inserted by sj, not clear whether correct, might be due to only single color channel
                profilesz(3) = 1;
            end
            profileStruct.trace.RGB = reshape(profileStruct.trace.img,[prod(profilesz(1:end-1)) profilesz(end)]);
            for nt=1:profileStruct.show.chsz
                profileStruct.trace.RGB(:,nt)=profileStruct.trace.RGB(:,nt)/max(profileStruct.trace.RGB(:,nt));
            end
            profileStruct.trace.RGB = profileStruct.trace.RGB*single(profileStruct.show.colors);
            profileStruct.trace.RGB = reshape(profileStruct.trace.RGB,[profilesz(1:end-1) 3]);
            profileStruct.trace.RGB(profileStruct.trace.RGB>1)=1;
            profileStruct.trace.RGB(profileStruct.trace.RGB<0)=0;
            if isfield(profileStruct.trace, 'plot1') && ishandle(profileStruct.trace.plot1)
                if profileStruct.show.bw
                    set(profileStruct.trace.profh{1},'YData',profileStruct.show.profile(:,1));
                else
                    for nt=1:profileStruct.show.chsz
                        set(profileStruct.trace.profh{nt},'YData',profileStruct.show.profile(:,nt),'Color',profileStruct.show.colors(nt,:)*.8);
                    end
                end
                set(profileStruct.trace.imgh,'CData',profileStruct.trace.RGB)
            else
                if profileStruct.show.bw
                    profileStruct.trace.profh{1} = plot(profileStruct.show.profile,'Parent',profileStruct.gui.profile.ax,'Color','k','tag',sprintf('traceplot_%.0f',profileStruct.instance));
                else
                    for nt=1:profileStruct.show.chsz
                        profileStruct.trace.profh{nt} = plot(profileStruct.show.profile(:,nt),'Parent',profileStruct.gui.profile.ax,'Color',profileStruct.show.colors(nt,:)*.8,'tag',sprintf('traceplot_%.0f',profileStruct.instance));
                    end
                end
                profileStruct.trace.imgh = imagesc(profileStruct.trace.img,'Parent',profileStruct.gui.traceim.ax,'tag',sprintf('traceplot_%.0f',profileStruct.instance));%profileStruct.trace.RGB
                if withAlpha
                  set(profileStruct.trace.imgh,'AlphaData',profileStruct.trace.imgA);
                  set(profileStruct.gui.traceim.ax,'Color','k');
                end
            end
            set(profileStruct.gui.traceim.ax,'box','on','YTickLabel',(profileStruct.show.width-1)/2:-1:-1*(profileStruct.show.width-1)/2,'XLim',[.5 length(tracepoints)-.5],'YLim',[.5 profileStruct.show.width+.5]);%,'Color',profileStruct.gui.bgcolor
            set(profileStruct.gui.profile.ax,'box','on','XTickLabel',[],'Color',profileStruct.gui.bgcolor,'XLim',[.5 length(tracepoints)-.5],'YLim',[0 max(profileStruct.show.profile(:))]);
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

%% Region of interest manager
%     nRois: actual number of "ROIs"
%     The structure roiList contains information for all Rois:
%         roiList( ).number:        number of Roi; displayed as text (roiText) in image/zoom
%                                   window. The number of a Roi doesn't change even if other Rois are deleted from the
%                                   list. Rois from the calculator receive zero as entry (as they don't appear in the
%                                   image windows. This is used to identify them as "virtual Rois" in the calculations
%                                   for the plots. If a new Roi is created, it gets the number max([roiList.number])+1.
%         roiList( ).name:          Name of Roi. Also never changes unless changed by user. The default includes the number
%                                   of the Roi for "normal" rois.
%         roiList( ).rectangle:     Coordinates of rectangle surrounding Roi
%                                   (left,down, right, up). If custom
%                                   dimension scales are used, the values
%                                   are specified in terms of the scales
%                                   (not the pixels). However, upon
%                                   export/import, these values are
%                                   translated to/from pixel values!
%         roiList( ).corners:       Coordinates of polygon defining Roi. If custom
%                                   dimension scales are used, the values
%                                   are specified in terms of the scales
%                                   (not the pixels). However, upon
%                                   export/import, these values are
%                                   translated to/from pixel values!
%         roiList( ).index.x/.y:    x/y index of all pixels of Roi
%         roiList( ).settings.xySel xy dimension during time of ROI selection
%         roiList( ).settings.position Position during selection of ROI
%         roiList( ).settings.zoomVal ZoomWindow settings during selection of ROI
%         roiList( ).virtualRoi:    Binary value defining whether Roi is
%                                   "virtual", ii.e. based on calculation
%                                   of two other Rois: not yet implemented
%                                   (if at all)
%         roiList( ).calc1:         Exists only if Roi is virtual Roi.
%                                   First Roi used for calculation
%         roiList( ).calc2:         Exists only if Roi is virtual Roi.
%                                   Second Roi used for calculation
%         roiList( ).op:            Exists only if Roi is virtual Roi.
%                                   Defines operation used to calculate Roi
%                                   (+,-,*,/)
%                                   Not yet implemented!!!

%   During export, the roiList structure is exported or saved, during
%   import, this structure is loaded and the rois are reconstructed from
%   it. Thus it is necessary that the chosen file contains a structure of
%   this kind.
%   When changing xySel (and not simply exchanging them), all Roi
%   information will be lost and the Roi Manager will be reset (resetRois).

%%  Build Roi Gui during first call or make it visible/unvisible during
%   later calls or from WindowCloseRequestFcn
    function roiGui(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        %Make Roi Gui and all related objects invisible
        if nargin == 3 || ~get(tbRoi,'Value')
            set(roiWin, 'Visible','off');
            set(btMean, 'Enable','on');
            set(cbPlots(xySel),'Enable', 'on');
            set(tbRoi,'Value',0)
            if nRois >0
                set([roiLine.im, roiLine.zoom, roiText.im, roiText.zoom],  'Visible', 'off');
            end
            set([tb_newRoi tbRoiShift tbRoiRotate tbRoiScale roiBtReplace], 'Value', 0);
            set([imageWin zoomWin],'WindowButtonMotionFcn',@mouseMotion,...
                'WindowButtonDownFcn',@buttonDownCallback,...
                'WindowButtonUpFcn','');
            drawPlots;
            if get(roiBtRename, 'Value')
              set(tempRoiPropWin, 'Visible','off');
            end
            if debugMatVis, debugMatVisFcn(2); end
            return
        end
        
        %Make Roi Gui and all related objects visible
        set(roiWin, 'Visible','on');
        set(btMean, 'Enable','off');
        if get(roiBtRename, 'UserData')
          set(tempRoiPropWin, 'Visible','on');
        end
        if nRois > 0 && get(roiBtShowNames, 'Value')
            set([roiLine.im, roiLine.zoom, roiText.im, roiText.zoom], 'Visible', 'on');
        elseif nRois > 0
            set([roiLine.im, roiLine.zoom, roiText.im, roiText.zoom], 'Visible', 'off');
        end
        set(cbPlots(xySel),  'Value', 0 ,'Enable', 'off'); %necessary when call from WindowCloseRequestFcn
        plotSel(xySel) = 0;
        drawPlots;
        drawRois;
        %Create ROI GuiWindow (only during first call)
        if isempty(roiWin)
            %Roi Manager Window with text
            gp = get(gui,'Position');
            roiWin =  figure;
            set(roiWin,'Units', 'Pixel', 'Name', 'ROI Manager',...
                'MenuBar', 'none', 'NumberTitle', 'off','CloseRequestFcn', {@roiGui,0}, 'HandleVisibility', 'off',...
                'WindowStyle','normal', 'Position', [gp(1) gp(2)-362 320 330]);
            warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
            roi_jf = get(roiWin,'JavaFrame');
            roi_jf.setFigureIcon(javax.swing.ImageIcon(im2java(uint8(icon_roi))));
            warning('on','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
            txt_roiList = uicontrol(roiWin, 'Style', 'Text', 'Position', [10 305 110 15], ...
                'String','List of ROIs','FontWeight', 'bold',...
                'BackgroundColor', get(roiWin, 'Color'), 'HorizontalAlignment', 'center');
            %Listbox for Rois
            roiListbox = uicontrol(roiWin, 'Style', 'listbox', 'Units', 'Pixel', ...
                'BackgroundColor', get(roiWin, 'Color'),'Value',0,...
                'Position', [10 153 110 147], 'String', '','Max',3, 'Callback', @listboxCallback,...
                'Tag', 'List of ROIs. Multiple selections are possible. The Names of the ROIs can be changed (see buttons below), the numbers in front of the names are generated by matVIs and can not be edited.');
            %Panel for display of Roi properties
            roiPropertiesPanel = uipanel(roiWin,'Units','Pixel','Position', [130,25,180,275],'Title', 'Selected ROI',...
                'BackgroundColor',get(roiWin,'Color'),'FontWeight', 'bold','TitlePosition','centertop',...
                'Tag','Shows the currently selected ROI and various properties (number of pixels, minimum, maximum and mean value). In case multiple ROIs are selected, only the first selected ROIs and its properties are shown here.');
            roiName = uicontrol('Parent',roiPropertiesPanel, 'Style', 'Text', 'Position', [5 210 170 30], ...
                'BackgroundColor', get(roiWin, 'Color'), 'HorizontalAlignment', 'center','FontWeight', 'bold',...
                'String', 'ROI Name');
            uicontrol(roiWin, 'Style', 'Text', 'Position', [152 75 40 15], ...
                'String','# Px:', 'BackgroundColor', get(roiWin, 'Color'), 'HorizontalAlignment', 'left');
            uicontrol(roiWin, 'Style', 'Text', 'Position', [152 60 40 15], ...
                'String','Min:', 'BackgroundColor', get(roiWin, 'Color'), 'HorizontalAlignment', 'left');
            uicontrol(roiWin, 'Style', 'Text', 'Position', [152 45 40 15], ...
                'String','Max:', 'BackgroundColor', get(roiWin, 'Color'), 'HorizontalAlignment', 'left');
            uicontrol(roiWin, 'Style', 'Text', 'Position', [152 30 55 15], ...
                'String','Mean:', 'BackgroundColor', get(roiWin, 'Color'), 'HorizontalAlignment', 'left');
            roiSize = uicontrol(roiWin, 'Style', 'Text', 'Position', [180 75 35 15], ...
                'BackgroundColor', get(roiWin, 'Color'), 'HorizontalAlignment', 'right');
            roiMin = uicontrol(roiWin, 'Style', 'Text', 'Position', [180 60 35 15], ...
                'BackgroundColor', get(roiWin, 'Color'), 'HorizontalAlignment', 'right');
            roiMax = uicontrol(roiWin, 'Style', 'Text', 'Position', [180 45 35 15], ...
                'BackgroundColor', get(roiWin, 'Color'), 'HorizontalAlignment', 'right');
            roiMean = uicontrol(roiWin, 'Style', 'Text', 'Position', [180 30 35 15], ...
                'BackgroundColor', get(roiWin, 'Color'), 'HorizontalAlignment', 'right');
            %roiFrame = uicontrol(roiWin,'Style','Frame',  'Position', [160,150,150,150],'Visible','off');
            %Axes for display of Roi
            roiAxes  = axes('Parent', roiPropertiesPanel, 'Unit','pixel','Position', [25 80 145 140]);
            
            %Add Rectangle Roi Button
            tb_newRoi(1) = uicontrol(roiWin, 'Style', 'Togglebutton','Position', [10,125,32,22],...
                'CData', roiRectBt,'Callback',{@getNewRoi,1,'one'}, ...
                'ToolTipString', ' Add rectangular ROI. Left click for single ROI (''push button''), right click for multiple selections (''toggle button'') ',...
                'Tag',  'Add rectangular ROI. Left click for single ROI (''push button''), right click for multiple selections (''toggle button'')');
            %Add Ellipse Roi Button
            tb_newRoi(2) = uicontrol(roiWin, 'Style', 'Togglebutton','Position', [49,125,32,22],...
                'CData', roiEllipseBt,'Callback',{@getNewRoi,2,'one'},...
                'ToolTipString', ' Add ellipsoid ROI. Left click for single ROI (''push button''), right click for multiple selections (''toggle button'')',...
                'Tag', 'Add ellipsoid ROI. Left click for single ROI (''push button''), right click for multiple selections (''toggle button'')');
            %Add Poly Roi Button
            tb_newRoi(3) = uicontrol(roiWin, 'Style', 'Togglebutton','Position', [88,125,32,22],...
                'CData', roiPolygonBt,'Callback',{@getNewRoi,3,'one'},...
                'ToolTipString', 'Add polygonal/free hand ROI. Left click for single ROI (''push button''), right click for multiple selections (''toggle button'') ',...
                'Tag', 'Add polygonal/free hand ROI. Draw straight lines using short clicks and free hand selections by keeping the mouse button pressed. End selection by right or double click. Left click for single ROI, right click for multiple selections');
            
            %Shift Roi Button
            tbRoiShift = uicontrol(roiWin, 'Style', 'Togglebutton','Position', [10,100,32,22], 'Enable', 'off',...
                'CData', roiShiftBt,'Callback',@shiftRoi,...
                'ToolTipString', ' Move selected Roi with mouse or keyboard. Use control, shift or alt to move by 5, 10 or 20 px, respectively. ',...
                'Tag', 'Move selected ROI(s) with mouse or keyboard. Use control, shift or alt to move by 5, 10 or 20 px, respectively.  Note that this is a toggle button.',...
                'Value',0);
            %Rotate Roi Button
            tbRoiRotate = uicontrol(roiWin, 'Style', 'Togglebutton','Position', [49,100,32,22], 'Enable', 'off',...
                'CData', roiRotateBt,'Callback',@rotateRoi,'ToolTipString', ' Rotate selected Roi ',...
                'Tag', 'Rotate selected ROI(s).  Note that this is a toggle button.');
            %Scale Roi Button
            tbRoiScale = uicontrol(roiWin, 'Style', 'Togglebutton','Position', [88,100,32,22], 'Enable', 'off',...
                'CData', roiScaleBt,'Callback',@scaleRoi,'ToolTipString', ' Scale selected Roi ','Tag', 'Scale the selected ROI(s). Note that this is a toggle button.');
            
            %Replace Roi
            roiBtReplace = uicontrol(roiWin, 'Style', 'Togglebutton','Position', [10,75,32,22],...
                'String', 'Replace','Tag', 'Replace the selected ROI. The position in the list and the name will be preserved. Note that this is a toggle button.',...
                'FontSize', 6, 'Enable', 'off','ToolTipString','Replace selected ROI by new selection');
            
            %Rename Roi Button
            roiBtRename = uicontrol(roiWin, 'Style', 'Togglebutton','Position', [49,75,32,22],...
                'String', 'ROI upd.','Callback', @updateRoiSettings,'FontSize', 6, 'Enable', 'off','UserData',0,...
                'ToolTipString','Update settings for selected ROI','Tag','Update settings for selected ROI');
            
            %Show / hide Roi names button
            roiBtShowNames = uicontrol(roiWin, 'Style', 'Togglebutton','Position', [88,75,32,22],...
                'String', 'ShowNames','Callback', @showRoiNames,'ToolTipString', ' Show or hide Roi names in figure windows ',...
                'Tag', 'Show or hide the ROI names in the figure windows.','FontSize', 6, 'Value', 1);
            
            %Delete Roi Button
            bt_deleteRoi = uicontrol(roiWin, 'Style', 'Pushbutton','Position', [10,50,32,22],...
                'CData',roiDeleteBt,'Callback', @deleteRoi,'ToolTipString', 'Delete selected ROI(s).','Tag', 'Delete the selected ROI(s).', 'Enable', 'off');
            %Export Roi Button
            bt_roiExport = uicontrol(roiWin, 'Style', 'Pushbutton','Position', [49,50,32,22],...
                'CData', roiExportBt,'Callback', @exportRois,'Enable','on','ToolTipString', ' Export selected Roi(s) ',...
                'Tag', 'Export the selected ROI(s) to the Matlab workspace or to a mat-file.', 'Enable', 'off');
            %Import Roi Button
            uicontrol(roiWin, 'Style', 'Pushbutton','Position', [88,50,32,22],...
                'CData', roiImportBt,'Callback', @importRois,'Enable','on','ToolTipString', ' Import Rois ',...
                'Tag', 'Import ROIs.');
            %Copy Roi Button
            bt_copyRoi = uicontrol(roiWin, 'Style', 'Pushbutton','Position', [10,24,32,22],...
                'CData', copyRoiBt,'Callback', @copyRoi,'Enable','on','ToolTipString', ' Copy Roi ',...
                'Tag', 'Copy ROI.', 'Enable','off');
            if nDim>2
                %Popup for Dimension for roi data export
                n = dimNames;
                n(xySel) = [];
                roiPopDataExport = uicontrol(roiWin, 'Style', 'Popup','Position', [49,24,32  ,22],...
                    'Enable','on','String',n,...
                    'Tag', 'Select dimension for ROI data export.','Fontsize',8);
                %Export Roi Data Button
                bt_exportRoiData = uicontrol(roiWin, 'Style', 'Pushbutton','Position', [88,25,32,22],...
                    'CData', roiExportDataBt,'Callback', @exportRoiData,'Enable','on','ToolTipString', ' Export Roi data along specified dimension ',...
                    'Tag', 'Export ROI data along specified dimension into the workspace.', 'Enable', 'off');
            end
            jump2ROIPos_cb = uicontrol(roiWin, 'Style', 'checkbox','Position', [10,2,250  ,22],...
                    'Value',1,'String','Jump to ROI selection position',...
                    'Tag', 'Jump to position in data set at which ROI was defined.');
            %             %Roi Calculator
            %              uicontrol(roiWin, 'Style', 'Text', 'Position', [90 28 110 15], ...
            %                     'String','ROI Calculator','FontWeight', 'bold',...
            %                     'BackgroundColor', get(roiWin, 'Color'), 'HorizontalAlignment', 'center');
            %             roiCalc1 = uicontrol(roiWin, 'Style', 'Popup','Position', [10,10,71  ,15],'String',{'Roi1';'Roi2';'Roi3'});
            %             roiCalcOp = uicontrol(roiWin, 'Style', 'Popup','Position', [88,10,30,15],'String',{'+';'-';'.*';'./'});
            %             roiCalc2 = uicontrol(roiWin, 'Style', 'Popup','Position', [130,10,71  ,15],'String',{'Roi1';'Roi2';'Roi3'});
            %             uicontrol(roiWin, 'Style', 'Text', 'Position', [201 6 8 15], ...
            %                 'String','=','FontWeight', 'bold',...
            %                 'BackgroundColor', get(roiWin, 'Color'), 'HorizontalAlignment', 'center');
            %             roiCalcRes = uicontrol(roiWin, 'Style', 'edit','Position', [213,4,60,20],'String','Roi1 + Roi2');
            %             %Add calculated term to list of Rois
            %             roiCalcAdd = uicontrol(roiWin, 'Style', 'PushButton','Position', [280,4,30,20],'String','Add');
            %Fill Roi Axes with current image
            set(roiWin, 'HandleVisibility', 'on','ResizeFcn', @resizeRoiWin,'WindowButtonDownFcn',@buttonDownRoiGui, 'WindowbuttonMotionFcn',@mouseMotion);
            if customDimScale
                roiImage = imagesc(dimScale(xySel(2),:),dimScale(xySel(1),:),currIm{1});
            else
                roiImage = imagesc(currIm{1});
            end
            if withAlpha
              set(roiImage, 'AlphaData', currAlphaMap{1},'AlphaDataMapping','scaled');
            end
            axis image;
            set(roiAxes,'FontSize',8,'Color','k');
            roiLine.roi = line(0,0,'Parent', roiAxes,'Color','w');
            for iii = 1:numel(data)
                roiCenterIndicator(1,iii) = text(0,0,'+','Parent',zoomAx(iii),'Color','r','FontWeight','b','FontSize',max([15,max(dim(xySel))/100]),'Visible','off','horizontalAlignment','center','FontName','Times');
                roiCenterIndicator(2,iii) = text(0,0,'+','Parent',imAx(iii),'Color','r','FontWeight','b','FontSize',max([15,max(dim(xySel))/100]),'Visible','off','horizontalAlignment','center','FontName','Times');
            end
            set(roiWin, 'HandleVisibility', 'off');
            %Set first entry of roi list (so it is not empty - will be
            %replaced by 1 when first Roi is selected)
            roiList(1).number = 0;
            updateColormap;
            set(roiWin, 'ResizeFcn',[]);adjustGuiSize(roiWin);set(roiWin,'ResizeFcn', @resizeRoiWin);
        end
        
        function resizeRoiWin(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            adjustGuiSize(roiWin,-1);
            newPos = get(roiWin, 'Position');
            if any(newPos(3:4) < [180 250])
                newPos = [newPos(1) newPos(2) 180 250];
                set(roiWin, 'Position', newPos);
            end
            set(roiListbox, 'Position', [10 153 110 newPos(4)-183]);
            set(txt_roiList, 'Position', [10 newPos(4)-25 110 15]);
            set(roiPropertiesPanel, 'Position', [130 25 newPos(3)-140 newPos(4)-35]);
            set(roiName, 'Position', [5 newPos(4)-90 newPos(3)-150 30]);
            set(roiAxes, 'Position', [25 80 newPos(3)-175 newPos(4)-160]);
            adjustGuiSize(roiWin);
            if debugMatVis, debugMatVisFcn(2); end
        end
        
        function buttonDownRoiGui(varargin)
			if debugMatVis, debugMatVisFcn(1); end
            if strcmp(get(roiWin,'SelectionType'),'alt')
                %Check for mouse position to enable fake 'button right-clicks'
                p1 = round(get(myGcf,'CurrentPoint'));
                for ii=1:3
                    btPosNewRoi(:,ii) = get(tb_newRoi(ii), 'Position'); %#ok
                    if (p1(1) > btPosNewRoi(1,ii) && ...
                            p1(1) < btPosNewRoi(1,ii) + btPosNewRoi(3,ii) && ...
                            p1(2) > btPosNewRoi(2,ii) && ...
                            p1(2) < btPosNewRoi(2,ii) + btPosNewRoi(4,ii))
                        if get(tb_newRoi(ii), 'Value')
                            set(tb_newRoi(ii), 'Value', 0);
                        else
                            set(tb_newRoi, 'Value', 0);
                            set(tb_newRoi(ii), 'Value', 1);
                            getNewRoi(0,0,ii);
                            if debugMatVis, debugMatVisFcn(2); end
                            return
                        end
                    end
                end
            elseif strcmp(get(roiWin,'SelectionType'),'extend')
                set(tempWin, 'HandleVisibility', 'on');
                if ~isempty(tempWin) && any(get(0,'Children') == tempWin)
                    set(tempWin, 'HandleVisibility', 'on');
                    close(tempWin);
                    tempWin = [];
                end
                tempWin = figure('Position',[gp(1) gp(2)-92 150 50],...
                    'Name', 'Enter ROI width','MenuBar', 'none', 'WindowStyle','normal', ...
                    'Resize', 'off', 'NumberTitle','off', 'HandleVisibility', 'off', 'CloseRequestFcn',@closeTempWin);
                uicontrol(tempWin, 'Style', 'Text', 'Position', [10 30 130 15], ...
                    'String','Specify ROI radius','FontWeight', 'bold',...
                    'BackgroundColor', get(tempWin, 'Color'), 'HorizontalAlignment', 'left');
                uicontrol(tempWin, 'Style', 'Edit', 'Position', [40 10 50 15], ...
                    'String',num2str(newRoiSize),'FontWeight', 'bold',...
                    'HorizontalAlignment', 'center', 'Callback', @updateRoiSize);
            end
            if debugMatVis, debugMatVisFcn(2); end
        end
        
        function updateRoiSize(varargin)
            newRoiSize = str2double(get(varargin{1},'string'));
        end
        
        %Choose new Roi (rectangular, ellipsoid or polygonal/free)
        function getNewRoi(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            roi = [];
            % Check if called with left mouse button ("push mode" = single ROI selection)
            if nargin == 4 && strcmp(varargin{4},'one')
                if ~get(tb_newRoi(varargin{3}), 'Value') %If previously selected button is pressed: unpress and quit
                    set(tb_newRoi, 'Value', 0);
                    set([zoomWin imageWin],'WindowButtonMotionFcn',@mouseMotion);
                    set([zoomWin imageWin],'WindowButtonDownFcn',@buttonDownCallback);
                    set([zoomWin imageWin],'WindowButtonUpFcn','');
                    if debugMatVis, debugMatVisFcn(2); end
                    return
                end
                set(tb_newRoi, 'Value', 0);
            end
            switch varargin{3}
                case 1                  %Select Rectangle
                    set([imageWin zoomWin],'WindowButtonDownFcn',@selectRectRoi);
                case 2                  %Select Ellipse
                    set([imageWin zoomWin],'WindowButtonDownFcn',@selectEllipseRoi);
                case 3                  %Select Polygon
                    tempLine = [];
                    roiSelLine = [];
                    set([imageWin zoomWin],'WindowButtonDownFcn',@selectPolyRoi);
            end
            function selectRectRoi(varargin)
                if debugMatVis, debugMatVisFcn(1); end
                roi = [];           %Clear old roi position
                p = get(myGca, 'CurrentPoint');
                p = p(1 ,[1 2]);
                if customDimScale
                    p = pixelLocation(p);
                end
                p1 = p(1 ,[1 2]);
                % If middle mouse button is used, create square ROI with
                % width 2*newRoiSize+1
                if strcmp(get(myGcf,'SelectionType'),'extend')
                    p = round(p1);
                    p1 = p-newRoiSize;
                    p2 = p+newRoiSize;
                else % Otherwise select rectangle
                    set(myGcf, 'HandleVisibility', 'on');
                    rbbox;
                    set(myGcf, 'HandleVisibility', 'off');
                    p = get(myGca, 'CurrentPoint');
                    p = p(1 ,[1 2]);
                    if customDimScale
                        p = pixelLocation(p);
                    else
                        p = round(p);
                    end
                    p2 = p(1 ,[1 2]);
                end
                roi(1,:) = [p1(1) p1(1) p2(1) p2(1) p1(1)];
                roi(2,:) = [p1(2) p2(2) p2(2) p1(2) p1(2)];
                nRois = nRois + 1;
                if customDimScale
                    roi(1,:) = round(roi(1,:)/pxWidth(xySel(1)))*pxWidth(xySel(1));
                    roi(2,:) = round(roi(2,:)/pxWidth(xySel(2)))*pxWidth(xySel(2));
                    for ii = 1:2
                        r = roi(ii,:);
                        if any(diff(r))
                            r(r==min(r)) = min(r)-pxWidth(xySel(ii))/2;
                            r(r==max(r)) = max(r)+pxWidth(xySel(ii))/2;
                            roi(ii,:) = r;
                        else % Single line or pixel
                            if ii==1;
                                roi(ii,:) = roi(ii,:)+pxWidth(ii)*[-0.5 -0.5 0.5 0.5 -0.5];
                            else
                                roi(ii,:) = roi(ii,:)+pxWidth(ii)*[-0.5 0.5 0.5 -0.5 -0.5];
                            end
                        end
                    end
                else
                    roi  = round(roi);
                    for ii = 1:2
                        r = roi(ii,:);
                        if any(diff(r))
                            r(r==min(r)) = min(r)-0.5;
                            r(r==max(r)) = max(r)+0.5;
                            roi(ii,:) = r;
                        else % Single line or pixel
                            if ii==1;
                                roi(ii,:) = roi(ii,:)+[-0.5 -0.5 0.5 0.5 -0.5];
                            else
                                roi(ii,:) = roi(ii,:)+[-0.5 0.5 0.5 -0.5 -0.5];
                            end
                        end
                    end
                end
                addNewRoi(roi,nRois,'new');
                if debugMatVis, debugMatVisFcn(2); end
            end
            function selectEllipseRoi(varargin)
                if debugMatVis, debugMatVisFcn(1); end
                roi = [];           %Clear old roi position
                p = get(myGca,'CurrentPoint');
                p = p(1  ,1:2);
                p1 = p;
                % If middle mouse button is used, create circular ROI with
                % diameter 2*newRoiSize+1
                if strcmp(get(myGcf,'SelectionType'),'extend')
                    if customDimScale
                        p = pixelLocation(p);
                    end
                    p1 = p;
                    p = round(p1);
                    p1 = p-newRoiSize;
                    p2 = p+newRoiSize;
                    x = min(p2(1),p1(1));
                    y = min(p2(2),p1(2));
                    w = abs(p2(1)-p1(1));
                    h = abs(p2(2)-p1(2));
                    t = 0:2*pi/(w+h):2*pi;
                    roi(1  ,:) = x+w/2+w/2*sin(t);
                    roi(2,:) = y+h/2+h/2*cos(t);
                    roi(:,end+1) = roi(:,1);
                    nRois = nRois + 1;
                    addNewRoi(roi,nRois,'new');
                else % Otherwise select ellipse
                    set(myGcf,'WindowButtonUpFcn',@finishEllipse);
                    set(myGcf, 'HandleVisibility', 'on');
                    tempEllipse = rectangle('Position',[p1(1) p1(2) .1 .1],'Curvature',[1 1],'EdgeColor','w');
                    tempRect = rectangle('Position',[p1(1) p1(2) .1 .1],'EdgeColor',[0.8 0.8 0.8], 'LineStyle',':');
                    set(myGcf, 'HandleVisibility', 'off');
                    set(myGcf,'WindowButtonMotionFcn',@updateRoiRect);
                end
                function updateRoiRect(varargin)
                    if debugMatVis, debugMatVisFcn(1); end
                    p = get(myGca,'CurrentPoint');
                    p = p(1  ,1:2);
%                     if customDimScale
%                         p = pixelLocation(p);
%                     end
                    x = min(p(1),p1(1));
                    y = min(p(2),p1(2));
                    w = abs(p(1)-p1(1));
                    h = abs(p(2)-p1(2));
                    if w == 0 || h == 0
                      if debugMatVis, debugMatVisFcn(2); end
                      return
                    end
                    set(tempEllipse, 'Position',[x y w h]);
                    set(tempRect, 'Position',[x y w h]);
                    if debugMatVis, debugMatVisFcn(2); end
                end
                function finishEllipse(varargin)
                    if debugMatVis, debugMatVisFcn(1); end
                    p2 = get(myGca, 'CurrentPoint');
                    p2 = p2(1 ,1:2);
                    if customDimScale
                        p2 = pixelLocation(p2);
                        p1 = pixelLocation(p1);
                    end
                    set(myGcf,'WindowButtonMotionFcn',@mouseMotion);
                    set(myGcf,'WindowButtonUpFcn','');
                    x = min(p2(1),p1(1));
                    y = min(p2(2),p1(2));
                    w = abs(p2(1)-p1(1));
                    h = abs(p2(2)-p1(2));
                    t = 0:2*pi/(w+h):2*pi;
                    roi(1  ,:) = x+w/2+w/2*sin(t);
                    roi(2,:) = y+h/2+h/2*cos(t);
                    roi(:,end+1) = roi(:,1);
                    delete(tempEllipse);
                    delete(tempRect);
                    nRois = nRois + 1;
                    addNewRoi(roi,nRois,'new');
                    if debugMatVis, debugMatVisFcn(2); end
                end
                if debugMatVis, debugMatVisFcn(2); end
            end
            function selectPolyRoi(varargin)
                if debugMatVis, debugMatVisFcn(1); end
                if ~any(myGcf==[imageWin zoomWin])
                  if debugMatVis, debugMatVisFcn(2); end
                  return
                end
                set(myGcf,'WindowButtonMotionFcn',@freeDraw);
                set(myGcf,'WindowButtonUpFcn',@stopFreeDraw);
                function freeDraw(varargin)
                    if debugMatVis >1, debugMatVisFcn(1); end
                    p = get(myGca, 'CurrentPoint');
                    p = p(1 ,1:2);
                    %fprintf('p(%f,%f)',p)
                    %if customDimScale
                    %    p = pixelLocation(p);
                    %    fprintf('\tp(%f,%f)',p)
                    %end
                    %fprintf('\n')
                    set(tempLine,'XData',[p(1  ,1),p(1  ,1)],'YData',[p(1  ,2),p(1  ,2)]);
                    roi(:,end+1) = p(1  ,1:2);
                    set(roiSelLine,'XData',roi(1  ,:),'YData',roi(2,:));
                    if debugMatVis >1, debugMatVisFcn(2); end
                end
                function stopFreeDraw(varargin)
                  if debugMatVis, debugMatVisFcn(1); end
                    % For some reason plotWin can appear as myGcf, even
                    % though varargin{1} == zoomWin/imageWin. Not clear how this can
                    % happen (stopFreeDraw can never become callback function
                    % for plotWin, since there's a check at the beginning
                    % of selectPolyRoi, the only place that stopFreeDraw is
                    % set. In fact, stopFreeDraw is never set as callback
                    % for plotWin, but updateLine is (resulting in errors
                    % when the mouse is moved over the plotWin during ROI
                    % selection.
                    if ~any(myGcf==[imageWin zoomWin])
                      if debugMatVis, debugMatVisFcn(2); end
                      return
                    end
                    if strcmp(get(myGcf,'SelectionType'),'normal')
                        set(myGcf,'WindowButtonMotionFcn',@updateLine);
                    end
                    if debugMatVis, debugMatVisFcn(2); end
                end
                %Select Polygon
                if strcmp(get(myGcf,'SelectionType'),'normal')
                    if isempty(tempLine)
                        p = get(myGca, 'CurrentPoint');
                        p = p(1 ,1:2);
                        %if customDimScale
                        %    p = pixelLocation(p);
                        %end
                        set(myGcf, 'HandleVisibility', 'on');
                        tempLine = line(p(1),p(2),'Color','w');
                        roiSelLine = line(p(1),p(2),'Color','w');
                        set(myGcf, 'HandleVisibility', 'off');
                    end
                    p = get(myGca, 'CurrentPoint');
                    p = p(1 ,1:2);
%                     if customDimScale
%                         p = pixelLocation(p);
%                     end
                    roi(:,end+1) = p(1  ,1:2);
                    set(roiSelLine,'XData',roi(1  ,:),'YData',roi(2,:));
                    %End Polygon Selection and extract ROI properties
                else
                    %Reset Window Mouse Function
                    set(myGcf,'WindowButtonMotionFcn',@mouseMotion);
                    if ~get(tb_newRoi(3), 'Value')
                        set(myGcf,'WindowButtonDownFcn',@buttonDownCallback);
                        set(myGcf,'WindowButtonUpFcn','');
                    end
                    %Complete polygon and round values
                    roi(:,end+1) = roi(:,1);
                    %Find pixels that were selected multiple time in a row
                    %and delete them
                    a = roi - circshift(roi,[0 1]);
                    b = intersect(find(a(1  ,:)==0),find(a(2,:) == 0));
                    b(1) = [];
                    roi(:,b) = [];
                    delete(roiSelLine);
                    roiSelLine = [];
                    delete(tempLine);
                    tempLine = [];
                    nRois = nRois + 1;
                    if customDimScale
                        for ii=1:size(roi,2)
                            roi(:,ii) = pixelLocation(roi(:,ii)');
                        end
                    end
                    addNewRoi(roi,nRois,'new');
                    roi = [];           %Clear old roi position
                end
                function updateLine(varargin)
                    if debugMatVis > 1, debugMatVisFcn(1); end
                    p = get(myGca, 'CurrentPoint');
                    p = p(1  ,1:2);
%                     if customDimScale
%                         p = pixelLocation(p);
%                     end
                    set(tempLine,'XData',[roi(1  ,end),p(1  ,1)],'YData',[roi(2,end),p(1  ,2)]);
                    if debugMatVis > 1, debugMatVisFcn(2); end
                end
                if debugMatVis, debugMatVisFcn(2); end
            end
            if debugMatVis, debugMatVisFcn(2); end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function addNewRoi(roi,numberRoi,newOld)
        if debugMatVis, debugMatVisFcn(1); end
        %Add selected Roi to list and create related objects
        % Check if get(roiBtReplace,'Value') == 1 and if so update numberRoi to current
        % selection
        if get(roiBtReplace,'Value')
            numberRoi = get(roiListbox, 'Value');
            nRois = nRois - 1;
        end
        
        [X, Y] = meshgrid(1:dim(xySel(2)), 1:dim(xySel(1)));
        % Matlab inpolygon function very slow. Use alternatives from file
        % exchange server:
        % tic
        % roiList(numberRoi).mask = inpolygon(Y,X,roi(2,:),roi(1  ,:));  % Slow Matlab function
        if exist('inpoly','file') == 3  % Check if mex-version (from Sebastien Paris) is available on path
            inregion = inpoly(cat(2,Y(:),X(:))',roi([2 1],:));
        else  % otherwise use "normal" inpoly.m function ( which I copied as a nested function into matVis (see below, including credits)
            inregion = inpoly_mv(cat(2,Y(:),X(:)),roi([2 1],:)');
        end
        % Didn't get this - supposedly faster - version to run (Matlab crashes). Maybe try again
        % later.
        % if exist('insidepoly','file') == 2  % Check if version of Bruno Luong is available on path
        %      inregion = insidepoly_sglengine(cat(2,Y(:),X(:))',roi([2 1],:));
        % else  % otherwise use "normal" inpoly.m function ( which I copied as a nested function into matVis (see below, including credits)
        %      inregion = inpoly_mv(cat(2,Y(:),X(:)),roi([2 1],:)');
        % end
        roiList(numberRoi).mask = false(dim(xySel(1)),dim(xySel(2)));
        roiList(numberRoi).mask(inregion) = 1;
        [maskXIdx, maskYIdx] = find(roiList(numberRoi).mask);
        roiList(numberRoi).centerOfGravity = [mean(maskXIdx), mean(maskYIdx)];
        %         toc
        if ~any(roiList(numberRoi).mask(:))
            roiList(numberRoi) = [];
            nRois = nRois - 1;
            warndlg('ROI completely outside of image. ROI will not be added / changed.');
            if strcmp(newOld, 'old')
                set(roiLine.im(:,numberRoi),'XData',roiList(numberRoi).corners(1  ,:),...
                    'YData',roiList(numberRoi).corners(2,:));
                set(roiLine.zoom(:,numberRoi),'XData',roiList(numberRoi).corners(1  ,:),...
                    'YData',roiList(numberRoi).corners(2,:));
                if customDimScale
                    set(roiText.im(:,numberRoi),'Position',...
                        [max(roiList(numberRoi).corners(1  ,:))+(dimScale(xySel(2),2)-dimScale(xySel(2),1))*0.01  ,mean(roiList(numberRoi).corners(2,:))] );
                    set(roiText.zoom(:,numberRoi),'Position',...
                        [max(roiList(numberRoi).corners(1  ,:))+(dimScale(xySel(2),2)-dimScale(xySel(2),1))*0.01  ,mean(roiList(numberRoi).corners(2,:))]);
                else
                    set(roiText.im(:,numberRoi),'Position',...
                        [max(roiList(numberRoi).corners(1  ,:))+dim(xySel(2))*0.01  ,mean(roiList(numberRoi).corners(2,:))] );
                    set(roiText.zoom(:,numberRoi),'Position',...
                        [max(roiList(numberRoi).corners(1  ,:))+dim(xySel(2))*0.01  ,mean(roiList(numberRoi).corners(2,:))]);
                
                end
            end
            if debugMatVis, debugMatVisFcn(2); end
            return
        end
        %Calculate index of each pixel in xySel plane
        [roiList(numberRoi).index.x roiList(numberRoi).index.y] = find(roiList(numberRoi).mask);
        if customDimScale
            roiList(numberRoi).corners = [];
            for ii = 1:size(roi,2)
                roiList(numberRoi).corners(:,ii) = scaledLocation(roi(:,ii)');
            end
        else
            roiList(numberRoi).corners = roi;
        end
        roiList(numberRoi).cornersPixelDim = roi;
        %Determine rectangle surrounding ROI pixels in polygon
        %with border of one pixel around rounded values (left,down, right, up)
        if customDimScale
            roiRegion = [max(dimScale(xySel(2),1)+pxWidth(xySel(2))  ,-pxWidth(xySel(2))+min(roiList(numberRoi).corners(1,:))),...
                max(dimScale(xySel(1),1)+pxWidth(xySel(1))  ,-pxWidth(xySel(1))+min(roiList(numberRoi).corners(2,:))),...
                min(dimScale(xySel(2),2),pxWidth(xySel(2))+max(roiList(numberRoi).corners(1,:))),...
                min(dimScale(xySel(1),2),pxWidth(xySel(1))+max(roiList(numberRoi).corners(2,:)))];
        else
            roiRegion = [round(max(1  ,-1+min(roiList(numberRoi).corners(1,:)))),...
                round(max(1  ,-1+min(roiList(numberRoi).corners(2,:)))),...
                round(min(dim(xySel(2)),1+max(roiList(numberRoi).corners(1  ,:)))),...
                round(min(dim(xySel(1)),1+max(roiList(numberRoi).corners(2,:))))];
        end
        roiList(numberRoi).rectangle = roiRegion;
        if numel(get(roiListbox,'Value'))==1 %~any(cell2mat(get([tbRoiShift tbRoiRotate tbRoiScale], 'Value')))
          roiList(numberRoi).settings.xySel = xySel;
          roiList(numberRoi).settings.position = currPos;
          % replace currPos of projection dim with position of max value
          % along projection dim
          if projMethod
            %updatePlotValues;
            %Get plot values for fs
            for m=1:nDim                                                    %Fill indices of all dimension with current-position-vectors of length of roi-size
              plotIndex{m} = currPos(m) * ones(size(roiList(numberRoi).index.x,1),1);
            end
            plotIndex{xySel(1)} = roiList(numberRoi).index.x;                       %Fill xySel dimension indices with roi indices
            plotIndex{xySel(2)} = roiList(numberRoi).index.y;
            plotIndex{projDim} = ones(size(roiList(numberRoi).index.x,1),1);        %Fill plot-dimension with ones (for first point)
            if get(tbSwitchRGB, 'Value') && (get(cmStretchRGBMean, 'Value') || get(cmStretchRGBMax, 'Value'))
              plotIndex{rgbDim} = ones(size(roiList(numberRoi).index.x,1),1);        %Fill plot-dimension with ones (for first point)
            end
            plotIndex = sub2ind(dim, plotIndex{:});                         %Determine linear index of roi pixels for first point
            plotIndex = repmat(plotIndex, [1 dim(projDim)]);                %Replicate linear index
            deltaIndex = prod(dim(1:projDim-1));                             %Difference value of indices along plotted dimension
            plotIndex = plotIndex + repmat(deltaIndex * (0:dim(projDim)-1),[size(roiList(numberRoi).index.x,1) 1]);   %Extend to all other points by adding deltaInd for each step
            if get(tbSwitchRGB, 'Value') && (get(cmStretchRGBMean, 'Value') || get(cmStretchRGBMax, 'Value'))
              plotIndex = repmat(plotIndex(:), [1 dim(rgbDim)]);                %Replicate linear index
              deltaIndex = prod(dim(1:rgbDim-1));                             %Difference value of indices along plotted dimension
              plotIndex = plotIndex + repmat(deltaIndex * (0:dim(rgbDim)-1),[size(roiList(numberRoi).index.x,1)*dim(projDim) 1]);   %Extend to all other points by adding deltaInd for each step
              plotIndex = reshape(plotIndex,[size(roiList(numberRoi).index.x,1) dim(projDim) dim(rgbDim)]);
              pV = mean(sum(data{1}(plotIndex),3, 'omitnan'),1, 'omitnan');
            else
              if withAlpha
                pV = mean(alphaMap{1}(plotIndex),    1, 'omitnan');
              else
                pV = mean(data{1}(plotIndex),1, 'omitnan');
              end
              
            end
            [~,pVMaxPos] = max(pV);
            roiList(numberRoi).settings.position(projDim) = pVMaxPos;
          end
          roiList(numberRoi).settings.zoomVal = zoomVal;
        end
        %Draw lines and add to list if Roi is new
        if strcmp(newOld, 'new')
            tb_newRoi_val = get(tb_newRoi, 'Value');
            if ~any([tb_newRoi_val{:}])
                set([imageWin zoomWin],'WindowButtonDownFcn',@buttonDownCallback);
            end
            set(tbRoiShift, 'Value', 0);
            if numberRoi == 1
                roiList(numberRoi).number = 1;
            elseif ~get(roiBtReplace,'Value')
                roiList(numberRoi).number = max([roiList.number])+1;
            end
            %Remove lines and text if Roi is replaced (instead of newly
            %defined)
            if get(roiBtReplace,'Value')
                delete(roiLine.im(:,numberRoi),roiText.im(:,numberRoi),roiLine.zoom(:,numberRoi),roiText.zoom(:,numberRoi));
            else
                roiList(numberRoi).name = ['Roi',num2str(roiList(numberRoi).number)];
            end
            %Draw final polygons and numbers in image and zoom
            %windows and delete temporary lines
            %             set(0, 'CurrentFigure', imageWin(1));
            %             set(imageWin(1), 'HandleVisibility', 'on');
            for ii = 1:nMat
                roiLine.im(ii,numberRoi) = line(roiList(numberRoi).corners(1,:),roiList(numberRoi).corners(2,:),'Color','w','Parent',imAx(ii));
                roiLine.zoom(ii,numberRoi) = line(roiList(numberRoi).corners(1,:),roiList(numberRoi).corners(2,:),'Color','w','Parent',zoomAx(ii));
                if customDimScale
                    roiText.im(ii,numberRoi) = text(max(roiList(numberRoi).corners(1,:))+(dimScale(xySel(2),2)-dimScale(xySel(2),1))*0.01  ,mean(roiList(numberRoi).corners(2,:)),roiList(numberRoi).name,'Color','w','HorizontalAlignment','left','Parent',imAx(ii));
                    roiText.zoom(ii,numberRoi) = text(max(roiList(numberRoi).corners(1,:))+(dimScale(xySel(2),2)-dimScale(xySel(2),1))*0.01  ,mean(roiList(numberRoi).corners(2,:)),roiList(numberRoi).name,'Color','w','HorizontalAlignment','left','Parent',zoomAx(ii));
                else
                    roiText.im(ii,numberRoi) = text(max(roiList(numberRoi).corners(1,:))+dim(xySel(2))*0.01  ,mean(roiList(numberRoi).corners(2,:)),roiList(numberRoi).name,'Color','w','HorizontalAlignment','left','Parent',imAx(ii));
                    roiText.zoom(ii,numberRoi) = text(max(roiList(numberRoi).corners(1,:))+dim(xySel(2))*0.01  ,mean(roiList(numberRoi).corners(2,:)),roiList(numberRoi).name,'Color','w','HorizontalAlignment','left','Parent',zoomAx(ii));
                end
            end
            %             set(imageWin(1), 'HandleVisibility', 'off');
            %             set(0, 'CurrentFigure', zoomWin(1));
            %             set(zoomWin(1), 'HandleVisibility', 'on');
            
            %             set(zoomWin(1), 'HandleVisibility', 'off');
            %Hide text if toggle button is off
            if ~get(roiBtShowNames, 'Value')
                set([roiText.im roiText.zoom], 'Visible', 'off');
            end
            if ~get(roiBtReplace,'Value')
                %Add to listbox if Roi is new
                s = get(roiListbox, 'String');
                s{size(s,1)+1} = [num2str(size(get(roiListbox,'String'),1)+1  ,'%03d'),': ',roiList(numberRoi).name];
                if numberRoi == 1
                    set(roiListbox,'Value', numberRoi,'String', s);
                else
                    set(roiListbox,'String', s, 'Value', numberRoi);
                end
            end
            updateRoiSelection(numberRoi);
            set([roiBtRename tbRoiShift tbRoiRotate tbRoiScale roiBtReplace bt_roiExport bt_deleteRoi] , 'Enable', 'on');
            %             if isempty(bt_exportRoiData)
            set(bt_exportRoiData, 'Enable','on')
            set(bt_copyRoi, 'Enable', 'on');
            %             end
            set(roiBtReplace,'Value',0)
        else
%             set(roiText.im(numberRoi),'Position',...
%                 [max(roiList(numberRoi).corners(1  ,:))+1  ,mean(roiList(numberRoi).corners(2,:))]);
%             set(roiText.zoom(numberRoi),'Position',...
%                 [max(roiList(numberRoi).corners(1  ,:))+1  ,mean(roiList(numberRoi).corners(2,:))]);
        end
        drawPlots;
        if debugMatVis, debugMatVisFcn(2); end
    end
    function copyRoi(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        copyRoiIdx = get(roiListbox, 'Value');
        for iii = 1:numel(copyRoiIdx)
            nRois = nRois + 1;
            addNewRoi(roiList(copyRoiIdx(iii)).cornersPixelDim,nRois,'new');
            roiList(nRois).settings = roiList(copyRoiIdx(iii)).settings;
        end
        set(roiListbox, 'Value',nRois-numel(copyRoiIdx)+1:nRois);
        updateRoiSelection;
        if debugMatVis, debugMatVisFcn(2); end
    end
    function deleteRoi(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        roiSel = get(roiListbox, 'Value');
        roiList(roiSel) = [];
        s = get(roiListbox, 'String');
        s(roiSel) = [];
        for ii=1:size(s,1)
            s{ii}(1:3) = num2str(ii,'%03d');
        end
        set(roiListbox,'Value',1);
        set(roiListbox, 'String',s);
        delete(roiText.im(:,roiSel));
        delete(roiText.zoom(:,roiSel));
        delete(roiLine.im(:,roiSel));
        delete(roiLine.zoom(:,roiSel));
        roiText.im(:,roiSel) = [];
        roiText.zoom(:,roiSel) = [];
        roiLine.im(:,roiSel) = [];
        roiLine.zoom(:,roiSel) = [];
        nRois = nRois - size(roiSel,2);
        if nRois > 0
            set(roiListbox, 'Value',min(roiSel(1),nRois));
            updateRoiSelection(min(roiSel(1),nRois));
        else
            set(roiListbox, 'Value',0, 'String', '');
            set([roiBtRename tbRoiShift tbRoiRotate tbRoiScale roiBtReplace bt_roiExport bt_deleteRoi bt_exportRoiData] , 'Enable', 'off');
            if ~isempty(bt_exportRoiData)
                set(bt_exportRoiData, 'Enable','off')
            end
            set(bt_copyRoi, 'Enable', 'off');
        end
        drawPlots;
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateRoiSettings(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if nargin == 3 || ~get(roiBtRename,'Value')
          set(tempRoiPropWin, 'Visible','off');
          set(roiBtRename,'Value',0)
          if debugMatVis, debugMatVisFcn(2); end
          return
        end
        currRoi = get(roiListbox, 'Value');
        ROIpos = roiList(currRoi).settings.position;
        nROIDim = length(ROIpos);
        ROIZoomVal = roiList(currRoi).settings.zoomVal;
        %ROIZoomVal(end+1:nROIDim,2) = dim(size(ROIZoomVal,1)+1:nROIDim);
        %ROIZoomVal(end+1:nROIDim,1) = 1; 
        for np = 1:min(nROIDim,nDim);
          ROIStr{np} = sprintf('Dim%d:  %4d %s %4d  %4d:%4d %s %4d:%4d',...
            np,ROIpos(np),char(hex2dec('2192')),currPos(np),...
            ROIZoomVal(np,1),sum(ROIZoomVal(np,:))-1,char(hex2dec('2192')),zoomVal(np,1),sum(zoomVal(np,:))-1);
        end
        gp = get(gui,'Position');
        if isempty(tempRoiPropWin)
          tempRoiPropWin = figure('Position', [gp(1), max(50,gp(2)-444-25*nROIDim), 320, 25*nROIDim + 50],...
            'MenuBar', 'none','NumberTitle', 'off', 'Name', 'ROI Properties!');%
          
          uicontrol(tempRoiPropWin,'Style', 'Text', 'Position', [12 25*nROIDim+26 55 18],...
            'String','ROI name:');
          etxt_newRoiName = uicontrol(tempRoiPropWin,'Style', 'Edit', 'Position', [80 25*nROIDim+25 150 20],...
            'String',roiList(currRoi).name,'HorizontalAlignment','left');
          uicontrol(tempRoiPropWin,'Style', 'Text', 'Position', [12 35 290 15*nROIDim+15],...
            'String','         Position        Zoom Settings','HorizontalAlignment','left','FontName','Courier');%
          for np = 1:min(nROIDim,nDim);
            if projMethod  && np == find(projDim == sort([xySel projDim])); myCol = [.5 .5 .5]; else myCol = [0 0 0]; end
            txt_RoiPropWinTxt(np) = uicontrol(tempRoiPropWin,'Style', 'Text', 'Position', [12 35+15*(nROIDim-np) 290 15],...
              'String',ROIStr{np},'HorizontalAlignment','left','FontName','Courier','ForegroundColor',myCol);%
          end
          btn1 = uicontrol(tempRoiPropWin,'Style', 'pushbutton', 'Position', [12 12 100 20],...
            'String','Close','Callback', {@updateRoiSettings,0});%ROIStr
          btn2 = uicontrol(tempRoiPropWin,'Style', 'pushbutton', 'Position', [200 12 100 20],...
            'String','Update','Callback', {@getNewRoiProps,1});%ROIStr
        else
          set(etxt_newRoiName,'String',roiList(currRoi).name)
          for np = 1:min(nROIDim,nDim);
            if projMethod  && np == find(projDim == sort([xySel projDim])); myCol = [.5 .5 .5]; else myCol = [0 0 0]; end
            set(txt_RoiPropWinTxt(np),'String',ROIStr{np},'ForegroundColor',myCol);%
          end
          set(tempRoiPropWin, 'Visible','on');
        end
        function getNewRoiProps(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            roiList(currRoi).name = get(etxt_newRoiName, 'String');
            roiList(currRoi).settings.xySel    = xySel;
            for np = 1:min(nROIDim,nDim);
              if ~(projMethod  && np == find(projDim == sort([xySel projDim])))
                roiList(currRoi).settings.position(np) = currPos(np);
              end
            end
            roiList(currRoi).settings.zoomVal  = zoomVal;
            set(roiText.im(:,currRoi), 'String', roiList(currRoi).name);
            set(roiText.zoom(:,currRoi), 'String', roiList(currRoi).name);
            s = get(roiListbox, 'String');
            s{currRoi} = [s{currRoi}(1:5) roiList(currRoi).name];
            set(roiListbox, 'String', s);
            updateRoiSelection
            %set(roiBtRename,'Value',0)
            %updateRoiSettings;
            if debugMatVis, debugMatVisFcn(2); end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateRoiSelection(numberRoi)
        if debugMatVis, debugMatVisFcn(1); end
        if nRois == 0
            if debugMatVis, debugMatVisFcn(2); end
            return
        end
        if nargin == 0
            numberRoi = get(roiListbox, 'Value');
        end

        if length(numberRoi) > 1
            set(roiBtReplace, 'Enable','off')
        else
            set(roiBtReplace, 'Enable','on')
        end
        if isempty(numberRoi)
            set(roiLine.im, 'LineWidth',1,'Color','w');
            set(roiLine.zoom, 'LineWidth',1,'Color','w');
            set(roiText.im, 'FontWeight', 'normal','Color','w');
            set(roiText.zoom, 'FontWeight', 'normal','Color','w');
            if debugMatVis, debugMatVisFcn(2); end
            return
        end
        %Set roi name above roiAxes
        rN = '';
        for ii=1:length(numberRoi)
            rN = [rN roiList(numberRoi(ii)).name '; '];  %#ok
        end
        rN(end-1:end) = [];
        set(roiName, 'String', rN);
        axis(roiAxes, [roiList(numberRoi(1)).rectangle(1)-0.5,...
            roiList(numberRoi(1)).rectangle(3)+0.5,...
            roiList(numberRoi(1)).rectangle(2)-0.5,...
            roiList(numberRoi(1)).rectangle(4)+0.5]);
        %Update Roi lines and numbers in image/zoom windows
        set(roiLine.roi, 'XData', roiList(numberRoi(1)).corners(1  ,:),...
            'YData', roiList(numberRoi(1)).corners(2,:));
        set(roiLine.im, 'LineWidth',1,'Color','w');              %'Color', 'b');
        set(roiLine.zoom, 'LineWidth',1,'Color','w');            %'Color', 'b');
        set(roiLine.im(:,numberRoi), 'LineWidth',3);   %'Color','r');
        set(roiLine.zoom(:,numberRoi), 'LineWidth',3); %'Color','r');
        set(roiText.im, 'FontWeight', 'normal','Color','w');     %'Color', 'b');
        set(roiText.zoom, 'FontWeight', 'normal','Color','w');   %'Color', 'b');
        set(roiText.im(:,numberRoi), 'FontWeight', 'bold'); %'Color','r');
        set(roiText.zoom(:,numberRoi),'FontWeight', 'bold'); % 'Color','r');
        if numel(numberRoi) == 1 
          if customDimScale
            set(roiCenterIndicator(:),'Visible','on','Position',[roiList(numberRoi).centerOfGravity([2 1]).*pxWidth(xySel([2 1])),0]);
          else
            set(roiCenterIndicator(:),'Visible','on','Position',[roiList(numberRoi).centerOfGravity([2 1]),0]);
          end
        else
            set(roiCenterIndicator(:),'Visible','off');
        end
        % Color for selected ROIs
        colorcodePlots = plotColors(numel(numberRoi));
        for ii = 1:numel(numberRoi)
            set([roiLine.im(:,numberRoi(ii)) roiLine.zoom(:,numberRoi(ii)) roiText.zoom(:,numberRoi(ii)) roiText.im(:,numberRoi(ii))], 'Color', colorcodePlots(ii,:));
        end
        if length(numberRoi) > 1
            set(roiBtRename, 'Enable', 'off');
        else
            set(roiBtRename, 'Enable', 'on');
        end
        % Set ROI properties that are indepedent of selected dimension
        set(roiSize, 'String', num2str(sum(roiList(numberRoi(1)).mask(:))));
        set(roiImage, 'CData',currIm{1}, 'AlphaData', 5/9*(0.8+roiList(numberRoi(1)).mask));
        set(roiAxes, 'Color','k');
        updateRoiProperties(0);
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateRoiProperties(plotUpdate)
        if debugMatVis, debugMatVisFcn(1); end
        numberRoi = get(roiListbox, 'Value');
        if numberRoi % Check in case function is called before a ROI has been created
            %Display in roiWin
            set(roiImage, 'CData',currIm{1});
            if withAlpha
              set(roiImage, 'AlphaData', currAlphaMap{1});
            end
            %Update Roi properties
            currRoi = currImVal{1}(roiList(numberRoi(1)).mask);
            if withAlpha
              currRoiAlpha = currAlphaMapVal{1}(roiList(numberRoi(1)).mask);
            end
            if isinteger(data{1})
                set(roiMin, 'String', num2str(min(currRoi(:))));     % Used to be nanmin
                set(roiMax, 'String', num2str(max(currRoi(:))));    % Used to be nanmax
                set(roiMean, 'String', num2str(mean(currRoi(:),'omitnan'),'%6.0f'));  % Used to be nanmean
            else
              if withAlpha
                set(roiMin, 'Position', [185 60 140 15], 'String',sprintf('%7.3f (without weight)',min(currRoi(:))), 'HorizontalAlignment', 'left');% Used to be nanmin
                set(roiMax, 'Position', [185 45 140 15], 'String', sprintf('%7.3f (without weight)',max(currRoi(:))), 'HorizontalAlignment', 'left');% Used to be nanmax
                set(roiMean, 'Position', [185 30 140 15], ...
                  'String', sprintf('%7.3f (alpha weighted)',sum(currRoi(:).*currRoiAlpha(:), 'omitnan')./sum(currRoiAlpha(:).*~isnan(currRoi(:)), 'omitnan') ), ... % Used to be nansum
                  'HorizontalAlignment', 'left');
              else
                set(roiMin, 'String', num2str(min(currRoi(:)),'%6.2f'));   % Used to be nanmin
                set(roiMax, 'String', num2str(max(currRoi(:)),'%6.2f'));   % Used to be nanmax
                set(roiMean, 'String', num2str(mean(currRoi(:), 'omitnan'),'%6.2f')); % Used to be nanmean
              end
            end
            if plotUpdate
                if isempty(subPlotHandles)
                    drawPlots;
                else
                    updatePlots;
                end
            end
            if get(roiBtRename,'Value'); updateRoiSettings; end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function listboxCallback(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        numberRoi = get(roiListbox, 'Value');
        if numel(numberRoi) == 1 && get(jump2ROIPos_cb,'Value')
          ROIPos = roiList(numberRoi).settings.position;
          matchDims = min([length(currPos) length(ROIPos)]);
          matchDims = find([ROIPos(1:matchDims)>dim(1:matchDims) 1],1,'first')-1; % checks if saved ROIpos is larger than data dimension and reduces MATCHDIMS, in case
          currPos(1:matchDims) = ROIPos(1:matchDims);
          currPos(currPos(1:matchDims)>dim(1:matchDims))=1;
          zoomVal(1:matchDims) = roiList(numberRoi).settings.zoomVal(1:matchDims); % checks if saved zoomVals are larger than data dimension and reduces MATCHDIMS, in case
          zoomVal(zoomVal((1:matchDims),2) > dim(1:matchDims)',1) = 1; 
          zoomVal(zoomVal((1:matchDims),2) > dim(1:matchDims)',2) = dim(zoomVal((1:matchDims),2)' > dim(1:matchDims));
          updateZoom(1,1,1:nDim,'ROIupdate');
          updateSelection(1:length(currPos));
        end
        drawPlots;
        
        updateRoiSelection(numberRoi);

        if debugMatVis, debugMatVisFcn(2); end
    end
    function drawRois
        if debugMatVis, debugMatVisFcn(1); end
        try
            delete(roiLine.zoom)
            delete(roiText.zoom)
        catch    %#ok
        end
        try
            delete(roiLine.im)
            delete(roiText.im)
        catch    %#ok
        end
        roiLine.zoom = [];
        roiText.zoom = [];
        roiLine.im = [];
        roiText.im = [];
        for ii=1:nRois
            for jj=1:nMat
                roiLine.im(jj,ii) = line(roiList(ii).corners(1  ,:),roiList(ii).corners(2,:),'Color','w','Parent',imAx(jj));
                roiText.im(jj,ii) = text(max(roiList(ii).corners(1  ,:))+1  ,mean(roiList(ii).corners(2,:)),roiList(ii).name,'Color','w','HorizontalAlignment','left','Parent',imAx(jj));
                roiLine.zoom(jj,ii) = line(roiList(ii).corners(1  ,:),roiList(ii).corners(2,:),'Color','w','Parent',zoomAx(jj));
                roiText.zoom(jj,ii) = text(max(roiList(ii).corners(1  ,:))+1  ,mean(roiList(ii).corners(2,:)),roiList(ii).name,'Color','w','HorizontalAlignment','left','Parent',zoomAx(jj));
            end
        end
        updateRoiSelection;
        if debugMatVis, debugMatVisFcn(2); end
    end
    function showRoiNames(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if get(roiBtShowNames, 'Value')
            set([roiText.im roiText.zoom], 'Visible', 'on');
        else
            set([roiText.im roiText.zoom], 'Visible', 'off');
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function shiftRoi(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if get(tbRoiShift, 'Value')
            set(imageWin,'WindowButtonDownFcn',@startShift);
            set(zoomWin,'WindowButtonDownFcn',@startShift);
            set(roiWin, 'KeyPressFcn', @shiftRoiWithKeyboard);
            set(tbRoiRotate, 'Value', 0);
            set(tbRoiScale, 'Value', 0);
        else
            set(imageWin,'WindowButtonDownFcn',@buttonDownCallback);
            set(zoomWin,'WindowButtonDownFcn',@buttonDownCallback);
            set(roiWin, 'KeyPressFcn', '');
        end
        function shiftRoiWithKeyboard(src,evnt)  %#ok
            if debugMatVis > 1, debugMatVisFcn(1); end
            switch evnt.Key
              case 'uparrow'
                shiftVector = [0 -1];
              case 'downarrow'
                shiftVector = [0 1];
              case 'leftarrow'
                shiftVector = [-1 0];
              case 'rightarrow'
                shiftVector = [1 0];
              otherwise
                if debugMatVis, debugMatVisFcn(2); end
                return
            end
            if ~isempty(evnt.Modifier)
                switch evnt.Modifier{1}
                    case 'control'
                        shiftVector = 5 * shiftVector;
                    case 'shift'
                        shiftVector = 10 * shiftVector;
                    case 'alt'
                        shiftVector = 20 * shiftVector;
                end
            end
            selRois = get(roiListbox,'Value');
            %             if any(roiList(selRois).corners(1,:)+shiftVector(1)<0) || ...
            %                any(roiList(selRois).corners(2,:)+shiftVector(2)<0) || ...
            %                any(roiList(selRois).corners(1,:)+shiftVector(1)>dim(xySel(1))+1) || ...
            %                any(roiList(selRois).corners(2,:)+shiftVector(2)>dim(xySel(2))+1)
            %                 return
            %             end
            set(roiWin, 'Name', 'BUSY!'); drawnow;
            for ii=1:size(selRois,2)
                roi(1,:) = roiList(selRois(ii)).corners(1,:)+shiftVector(1);
                roi(2,:) = roiList(selRois(ii)).corners(2,:)+shiftVector(2);
                addNewRoi(roi, selRois(ii),'old');
                roi = [];
            end
            updateRoiSelection(selRois);
            
            
            %             selRois = get(roiListbox,'Value');
            for ii=1:size(selRois,2)
                set(roiLine.im(:,selRois(ii)),'XData',roiList(selRois(ii)).corners(1,:),...
                    'YData',roiList(selRois(ii)).corners(2,:));
                set(roiLine.zoom(:,selRois(ii)),'XData',roiList(selRois(ii)).corners(1,:),...
                    'YData',roiList(selRois(ii)).corners(2,:));
                set(roiText.im(:,selRois(ii)),'Position',...
                    [max(roiList(selRois(ii)).corners(1  ,:))+1  ,mean(roiList(selRois(ii)).corners(2,:))]);
                set(roiText.zoom(:,selRois(ii)),'Position',...
                    [max(roiList(selRois(ii)).corners(1  ,:))+1  ,mean(roiList(selRois(ii)).corners(2,:))]);
            end
            set(roiWin, 'Name', 'ROI Manager');
            if debugMatVis > 1, debugMatVisFcn(2); end
        end
        function startShift(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            p = get(myGca, 'CurrentPoint');
            p = p(1  ,1:2);
            selRois = get(roiListbox,'Value');
            set(myGcf, 'WindowButtonMotionFcn',{@makeShift,p,selRois});
            set(myGcf, 'WindowButtonUpFcn',{@endShift,selRois,[0 0]});
            if debugMatVis, debugMatVisFcn(2); end
        end
        function makeShift(hO,h,p,selRois)  %#ok
            if debugMatVis, debugMatVisFcn(1); end
            p2 = get(myGca, 'CurrentPoint');
            p2 = p2(1  ,1:2);
            if customDimScale
                shiftVector = round((p2 - p)./pxWidth(xySel([2 1]))).*pxWidth(xySel([2 1]));
            else
                shiftVector = round(p2 - p);
            end
            for ii=1:size(selRois,2)
                set(roiLine.im(:,selRois(ii)),'XData',roiList(selRois(ii)).corners(1  ,:)+shiftVector(1),...
                    'YData',roiList(selRois(ii)).corners(2,:)+shiftVector(2));
                set(roiLine.zoom(:,selRois(ii)),'XData',roiList(selRois(ii)).corners(1  ,:)+shiftVector(1),...
                    'YData',roiList(selRois(ii)).corners(2,:)+shiftVector(2));
                set(roiText.im(:,selRois(ii)),'Position',...
                    [max(roiList(selRois(ii)).corners(1  ,:))+1  ,mean(roiList(selRois(ii)).corners(2,:))] + shiftVector);
                set(roiText.zoom(:,selRois(ii)),'Position',...
                    [max(roiList(selRois(ii)).corners(1  ,:))+1  ,mean(roiList(selRois(ii)).corners(2,:))] + shiftVector);
            end
            set(myGcf, 'WindowButtonUpFcn',{@endShift,selRois,shiftVector});
            if debugMatVis, debugMatVisFcn(2); end
        end
        function endShift(hO,h,selRois,shiftVector)  %#ok
            if debugMatVis, debugMatVisFcn(1); end
            set(myGcf, 'WindowButtonMotionFcn',@mouseMotion);
            set(myGcf, 'WindowButtonUpFcn','');
            for ii=1:size(selRois,2)
                roi(1,:) = roiList(selRois(ii)).corners(1,:)+shiftVector(1);
                roi(2,:) = roiList(selRois(ii)).corners(2,:)+shiftVector(2);
                if customDimScale
                    for jj=1:size(roi,2)
                        roi(:,jj) = pixelLocation(roi(:,jj)');
                    end
                end
                addNewRoi(roi, selRois(ii),'old');
                roi = [];
            end
            updateRoiSelection(selRois);
            set(tbRoiShift, 'Value', 0);
            if debugMatVis, debugMatVisFcn(2); end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function rotateRoi(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if get(tbRoiRotate, 'Value')
            set(imageWin,'WindowButtonDownFcn',@startRotation);
            set(zoomWin,'WindowButtonDownFcn',@startRotation);
            set(tbRoiShift, 'Value', 0);
            set(tbRoiScale, 'Value', 0);
        else
            set(imageWin,'WindowButtonDownFcn',@buttonDownCallback);
            set(zoomWin,'WindowButtonDownFcn',@buttonDownCallback);
        end
        function startRotation(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            p = get(myGca, 'CurrentPoint');
            p = p(1  ,1:2);
            selRois = get(roiListbox,'Value');
            set(myGcf, 'WindowButtonMotionFcn',{@makeRotation,p,selRois});
            set(myGcf, 'WindowButtonUpFcn',{@endRotation,selRois,p,p});
            if debugMatVis, debugMatVisFcn(2); end
        end
        function makeRotation(hO,h,p,selRois)  %#ok
            if debugMatVis, debugMatVisFcn(1); end
            p2 = get(myGca, 'CurrentPoint');
            p2 = p2(1  ,1:2);
            if numel(selRois)>1
              minCornerROIs = roiList(selRois(1)).corners(:  ,1);
              maxCornerROIs = minCornerROIs;
              for ii=1:size(selRois,2)
                minCornerROIs(1) = min([minCornerROIs(1) roiList(selRois(ii)).corners(1  ,1:end-1)]);
                minCornerROIs(2) = min([minCornerROIs(2) roiList(selRois(ii)).corners(2  ,1:end-1)]);
                maxCornerROIs(1) = max([maxCornerROIs(1) roiList(selRois(ii)).corners(1  ,1:end-1)]);
                maxCornerROIs(2) = max([maxCornerROIs(2) roiList(selRois(ii)).corners(2  ,1:end-1)]);
              end
              meanCornerROIs = (minCornerROIs + maxCornerROIs)/2;
              rotAngle = pi * sum(p-p2) / sum(maxCornerROIs - minCornerROIs);
              rotMat = [cos(rotAngle) sin(rotAngle);-sin(rotAngle) cos(rotAngle)];
              for ii=1:size(selRois,2)
                newCoord = [];
                newCoord(1,:) = roiList(selRois(ii)).corners(1,:) - meanCornerROIs(1);
                newCoord(2,:) = roiList(selRois(ii)).corners(2,:) - meanCornerROIs(2);
                newCoord = rotMat * newCoord ;
                newCoord(1,:) = newCoord(1,:)                     + meanCornerROIs(1);
                newCoord(2,:) = newCoord(2,:)                     + meanCornerROIs(2);
                set(roiLine.im(  :,selRois(ii)),'XData',newCoord(1  ,:),'YData',newCoord(2,:));
                set(roiLine.zoom(:,selRois(ii)),'XData',newCoord(1  ,:),'YData',newCoord(2,:));
              end
            else
              for ii=1:size(selRois,2)
                rotAngle = pi * sum(p-p2) / (roiList(selRois(ii)).rectangle(3)-roiList(selRois(ii)).rectangle(1)+roiList(selRois(ii)).rectangle(4)-roiList(selRois(ii)).rectangle(2));
                rotMat = [cos(rotAngle) sin(rotAngle);-sin(rotAngle) cos(rotAngle)];
                newCoord = [];
                newCoord(1,:) = roiList(selRois(ii)).corners(1,:) - mean(roiList(selRois(ii)).corners(1,1:end-1));
                newCoord(2,:) = roiList(selRois(ii)).corners(2,:) - mean(roiList(selRois(ii)).corners(2,1:end-1));
                newCoord = rotMat * newCoord ;
                newCoord(1,:) = newCoord(1,:)                     + mean(roiList(selRois(ii)).corners(1,1:end-1));
                newCoord(2,:) = newCoord(2,:)                     + mean(roiList(selRois(ii)).corners(2,1:end-1));
                set(roiLine.im(  :,selRois(ii)),'XData',newCoord(1  ,:),'YData',newCoord(2,:));
                set(roiLine.zoom(:,selRois(ii)),'XData',newCoord(1  ,:),'YData',newCoord(2,:));
              end
            end
            set(myGcf, 'WindowButtonUpFcn',{@endRotation,selRois,p,p2});
            if debugMatVis, debugMatVisFcn(2); end
        end
        function endRotation(hO,h,selRois,p,p2)  %#ok
            if debugMatVis, debugMatVisFcn(1); end
            set(myGcf, 'WindowButtonMotionFcn',@mouseMotion);
            set(myGcf, 'WindowButtonUpFcn','');
            if numel(selRois)>1
              minCornerROIs = roiList(selRois(1)).corners(:  ,1);
              maxCornerROIs = minCornerROIs;
              for ii=1:size(selRois,2)
                minCornerROIs(1) = min([minCornerROIs(1) roiList(selRois(ii)).corners(1  ,1:end-1)]);
                minCornerROIs(2) = min([minCornerROIs(2) roiList(selRois(ii)).corners(2  ,1:end-1)]);
                maxCornerROIs(1) = max([maxCornerROIs(1) roiList(selRois(ii)).corners(1  ,1:end-1)]);
                maxCornerROIs(2) = max([maxCornerROIs(2) roiList(selRois(ii)).corners(2  ,1:end-1)]);
              end
              meanCornerROIs = (minCornerROIs + maxCornerROIs)/2;
              rotAngle = pi * sum(p-p2) / sum(maxCornerROIs - minCornerROIs);
              rotMat = [cos(rotAngle) sin(rotAngle);-sin(rotAngle) cos(rotAngle)];
              for ii=1:size(selRois,2)
                newCoord = [];
                newCoord(1,:) = roiList(selRois(ii)).corners(1,:) - meanCornerROIs(1);
                newCoord(2,:) = roiList(selRois(ii)).corners(2,:) - meanCornerROIs(2);
                newCoord = rotMat * newCoord ;
                newCoord(1,:) = newCoord(1,:)                     + meanCornerROIs(1);
                newCoord(2,:) = newCoord(2,:)                     + meanCornerROIs(2);
                if customDimScale
                    for jj=1:size(newCoord,2)
                        newCoord(:,jj) = pixelLocation(newCoord(:,jj)');
                    end
                end
                addNewRoi(newCoord, selRois(ii),'old');
              end
            else
              for ii=1:size(selRois,2)
                rotAngle = pi * sum(p-p2) / (roiList(selRois(ii)).rectangle(3)-roiList(selRois(ii)).rectangle(1)+roiList(selRois(ii)).rectangle(4)-roiList(selRois(ii)).rectangle(2));
                rotMat = [cos(rotAngle) sin(rotAngle);-sin(rotAngle) cos(rotAngle)];
                newCoord = [];
                newCoord(1,:) = roiList(selRois(ii)).corners(1,:) - mean(roiList(selRois(ii)).corners(1,1:end-1));
                newCoord(2,:) = roiList(selRois(ii)).corners(2,:) - mean(roiList(selRois(ii)).corners(2,1:end-1));
                newCoord = rotMat * newCoord ;
                newCoord(1,:) = newCoord(1,:)                     + mean(roiList(selRois(ii)).corners(1,1:end-1));
                newCoord(2,:) = newCoord(2,:)                     + mean(roiList(selRois(ii)).corners(2,1:end-1));
                if customDimScale
                    for jj=1:size(newCoord,2)
                        newCoord(:,jj) = pixelLocation(newCoord(:,jj)');
                    end
                end
                addNewRoi(newCoord, selRois(ii),'old');
              end
            end
            updateRoiSelection(selRois);
            set(tbRoiRotate, 'Value', 0);
            if debugMatVis, debugMatVisFcn(2); end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function scaleRoi(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if get(tbRoiScale, 'Value')
            set(imageWin,'WindowButtonDownFcn',@startScaling);
            set(zoomWin,'WindowButtonDownFcn',@startScaling);
            set(tbRoiRotate, 'Value', 0);
            set(tbRoiShift, 'Value', 0);
        else
            set(imageWin,'WindowButtonDownFcn',@buttonDownCallback);
            set(zoomWin,'WindowButtonDownFcn',@buttonDownCallback);
        end
        function startScaling(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            p = get(myGca, 'CurrentPoint');
            p = p(1  ,1:2);
%             if customDimScale
%                 p = pixelLocation(p);
%             end
            selRois = get(roiListbox,'Value');
            set(myGcf, 'WindowButtonMotionFcn',{@makeScaling,p,selRois});
            set(myGcf, 'WindowButtonUpFcn',{@endScaling,selRois,[0 0]});
            if debugMatVis, debugMatVisFcn(2); end
        end
        function makeScaling(hO,h, p, selRois)  %#ok
            if debugMatVis, debugMatVisFcn(1); end
            p2 = get(myGca, 'CurrentPoint');
            p2 = p2(1  ,1:2);
%             if customDimScale
%                 p2 = pixelLocation(p2);
%             end
            for ii=1:size(selRois,2)
                scalingFct(1) = 1 + (p2(1)-p(1)) / (roiList(selRois(ii)).rectangle(3)-roiList(selRois(ii)).rectangle(1));
                scalingFct(2) = 1 + (p2(2)-p(2)) / (roiList(selRois(ii)).rectangle(4)-roiList(selRois(ii)).rectangle(2));
                newCoord(1,:) = scalingFct(1) * (roiList(selRois(ii)).corners(1,:) - mean(roiList(selRois(ii)).corners(1,1:end-1)));
                newCoord(2,:) = scalingFct(2) * (roiList(selRois(ii)).corners(2,:) - mean(roiList(selRois(ii)).corners(2,1:end-1)));
                newCoord(1,:) = newCoord(1,:) + mean(roiList(selRois(ii)).corners(1,1:end-1));
                newCoord(2,:) = newCoord(2,:) + mean(roiList(selRois(ii)).corners(2,1:end-1));
                set(roiLine.im(:,  selRois(ii)),'XData',newCoord(1  ,:),'YData',newCoord(2,:));
                set(roiLine.zoom(:,selRois(ii)),'XData',newCoord(1  ,:),'YData',newCoord(2,:));
                newCoord = [];
            end
            set(myGcf, 'WindowButtonUpFcn',{@endScaling,selRois,p,p2});
            if debugMatVis, debugMatVisFcn(2); end
        end
        function endScaling(hO,h,selRois,p,p2)  %#ok
            if debugMatVis, debugMatVisFcn(1); end
            set(myGcf, 'WindowButtonMotionFcn',@mouseMotion);
            set(myGcf, 'WindowButtonUpFcn','');
            for ii=1:size(selRois,2)
                scalingFct(1) = 1 + (p2(1)-p(1)) / (roiList(selRois(ii)).rectangle(3)-roiList(selRois(ii)).rectangle(1));
                scalingFct(2) = 1 + (p2(2)-p(2)) / (roiList(selRois(ii)).rectangle(4)-roiList(selRois(ii)).rectangle(2));
                newCoord(1,:) = scalingFct(1) * (roiList(selRois(ii)).corners(1,:) - mean(roiList(selRois(ii)).corners(1,1:end-1)));
                newCoord(2,:) = scalingFct(2) * (roiList(selRois(ii)).corners(2,:) - mean(roiList(selRois(ii)).corners(2,1:end-1)));
                newCoord(1,:) = newCoord(1,:) + mean(roiList(selRois(ii)).corners(1,1:end-1));
                newCoord(2,:) = newCoord(2,:) + mean(roiList(selRois(ii)).corners(2,1:end-1));
                if customDimScale
                    for jj=1:size(newCoord,2)
                        newCoord(:,jj) = pixelLocation(newCoord(:,jj)');
                    end
                end
                addNewRoi(newCoord, selRois(ii),'old');
                newCoord = [];
            end
            updateRoiSelection(selRois);
            set(tbRoiScale, 'Value', 0);
            if debugMatVis, debugMatVisFcn(2); end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function exportRois(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        roiSel = get(roiListbox, 'Value');
        if isempty(roiSel)
            warndlg('No ROIs selected. Select ROI first before exporting it/them.','No ROI available');
        else
            roiListExp = roiList;
            % Convert corners to pixel values
            if customDimScale
                for ii=1:numel(roiSel)
                    for jj = 1:size(roiListExp(roiSel(ii)).corners,2)
                        roiListExp(roiSel(ii)).corners(:,jj) = pixelLocation(roiListExp(roiSel(ii)).corners(:,jj)');
                    end
                    roiListExp(roiSel(ii)).rectangle = [round(max(1  ,-1+min(roiListExp(roiSel(ii)).corners(1,:)))),...
                        round(max(1  ,-1+min(roiListExp(roiSel(ii)).corners(2,:)))),...
                        round(min(dim(xySel(2)),1+max(roiListExp(roiSel(ii)).corners(1  ,:)))),...
                        round(min(dim(xySel(1)),1+max(roiListExp(roiSel(ii)).corners(2,:))))];
                end
            end
            % Remove corners outside image (is this necessary?)
            for ii=1:numel(roiSel)
                roiListExp(roiSel(ii)).corners(1,roiListExp(roiSel(ii)).corners(1,:) > dim(xySel(2))) = dim(xySel(2));
                roiListExp(roiSel(ii)).corners(2,roiListExp(roiSel(ii)).corners(2,:) > dim(xySel(1))) = dim(xySel(1));
                roiListExp(roiSel(ii)).corners(1,roiListExp(roiSel(ii)).corners(1,:) < 0) = 0;
                roiListExp(roiSel(ii)).corners(2,roiListExp(roiSel(ii)).corners(2,:) < 0) = 0;
            end
            
            q = questdlg('Export Rois to workspace or save as .mat file?','ROI Export', 'To File', 'To Workspace', 'Cancel', 'To File');
            if strcmp(q,'To Workspace')
                assignin('base', 'matVisRoiExport', roiListExp(roiSel));
            elseif strcmp(q, 'To File')
                [f,p] = uiputfile('.mat','Choose folder and filename to save rois!');
                if f == 0
                  if debugMatVis, debugMatVisFcn(2); end
                  return
                end
                matVisRoiExport =  roiListExp(roiSel); %#ok
                save([p,f],'matVisRoiExport');
            end
            clear matVisRoiExport roiListExp
            display('Roi export complete!');
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function importRois(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        [f,p] = uigetfile('.mat','Choose file to load rois!');
        if isequal(f,0)
          if debugMatVis, debugMatVisFcn(2); end
          return
        end
        try
            matVisRoiExport = cell2mat(struct2cell(load([p,f])));
        catch    %#ok
            errordlg('Chosen file does not contain rois exported using matVis (''matVisRoiExport'')!');
            if debugMatVis, debugMatVisFcn(2); end
            return
        end
        if customDimScale
            for ii=1:numel(matVisRoiExport)
                % Convert corners to scale
                for jj = 1:size(matVisRoiExport(ii).corners,2)
                    matVisRoiExport(ii).corners(:,jj) = scaledLocation(matVisRoiExport(ii).corners(:,jj)');
                end
                % Convert rectangle to scale
                matVisRoiExport(ii).rectangle = [max(dimScale(xySel(2),1)+pxWidth(xySel(2))  ,-pxWidth(xySel(2))+min(matVisRoiExport(ii).corners(1,:))),...
                    max(dimScale(xySel(1),1)+pxWidth(xySel(1))  ,-pxWidth(xySel(1))+min(matVisRoiExport(ii).corners(2,:))),...
                    min(dimScale(xySel(2),2),pxWidth(xySel(2))+max(matVisRoiExport(ii).corners(1,:))),...
                    min(dimScale(xySel(1),2),pxWidth(xySel(1))+max(matVisRoiExport(ii).corners(2,:)))];
            end
        end
        if nRois  % nRois == 0 if there are no rois existing
            q = questdlg('Append imported ROIs or overwrite existing ROIs?','Append or overwrite?','Overwrite','Append','Overwrite');
            if strcmp(q,'Append')
                roiList(end+1:end+numel(matVisRoiExport)) = matVisRoiExport;
            else
                roiList = matVisRoiExport;
            end
        else
            roiList = matVisRoiExport;
        end
        nRois = size(roiList,2);
        s = [];
        for ii=1:nRois
            s{ii} = roiList(ii).name;
        end
        s = [];
        set(roiListbox, 'String', '', 'Value', 1);
        for ii=1:nRois
            s{ii} = [num2str(ii,'%03d'),': ',roiList(ii).name];
        end
        set(roiListbox, 'String',s,'Value', 1);
        set([roiBtRename tbRoiShift tbRoiRotate tbRoiScale roiBtReplace bt_roiExport bt_deleteRoi bt_exportRoiData] , 'Enable', 'on');
        drawPlots;
        drawRois;
        if debugMatVis, debugMatVisFcn(2); end
    end
    function exportRoiData(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        selRois = get(roiListbox, 'Value');
        if isempty(selRois)
            warndlg('No ROIs selected. Select ROI first before exporting ROI data.','No ROI available');
        else
            s = get(roiPopDataExport,'String');
            exportDim = find(strcmp(dimNames, s{get(roiPopDataExport,'Value')}));
            for ii = 1:nDim
                imIndex{ii} = currPos(ii);  %#ok
            end
            deltaIndex = prod(dim(1:exportDim-1));
            for ii = 1:nMat
                for kk=1:size(selRois,2);
                    dataIndex = imIndex;
                    for m=1:nDim                                                    %Fill indices of all dimension with current-position-vectors of length of roi-size
                        dataIndex{m} = dataIndex{m} * ones(size(roiList(selRois(kk)).index.x,1),1);
                    end
                    dataIndex{xySel(1)} = roiList(selRois(kk)).index.x;                       %Fill xySel dimension indices with roi indices
                    dataIndex{xySel(2)} = roiList(selRois(kk)).index.y;
                    dataIndex{exportDim} = ones(size(roiList(selRois(kk)).index.x,1),1);        %Fill plot-dimension with ones (for first point)
                    dataIndex = sub2ind(dim, dataIndex{:});                         %Determine linear index of roi pixels for first point
                    dataIndex = repmat(dataIndex, [1 dim(exportDim)]);                %Replicate linear index
                    dataIndex = dataIndex + repmat(deltaIndex * (0:dim(exportDim)-1),[size(roiList(selRois(kk)).index.x,1) 1]);   %Extend to all other points by adding deltaInd for each step
                    exportValues(kk,:) = mean(data{ii}(dataIndex),1,'omitnan');   %#ok
                end
                assignin('base',['roiData_set' num2str(ii,'%02d') '_',dimNames{exportDim}], exportValues);
            end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function resetRois(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        delete(roiLine.zoom);
        delete(roiLine.im);
        delete(roiLine.roi(ishandle(roiLine.roi)));
        delete(roiText.zoom);
        delete(roiText.im);
        currRoi = [];                            %Current Roi
        roiMin = [];                             %Min of current Roi
        roiMax = [];                             %Max of current Roi
        roiMean = [];                            %Mean of current Roi
        roiSize = [];                            %Number of pixels of Roi
        roiList = [];                            %List of Rois
        roiLine.zoom = [];                            %Handle to lines indicating Rois
        roiLine.im = [];                            %Handle to lines indicating Rois
        roiLine.roi = [];                            %Handle to lines indicating Rois
        nRois = 0;                               %Number of Rois
        roiText.zoom = [];                             %Handle to text (numbers) for Rois
        roiText.im = [];                            %Handle to text (numbers) for Rois
        roiList(1).number = 0;
        %drawRois;
        set(roiListbox, 'String','','Value',0);
        set(roiMin, 'String','');
        set(roiMax, 'String','');
        set(roiMean, 'String','');
        set(roiSize, 'String','');
        if ~isempty(roiWin)
            set(0, 'CurrentFigure', roiWin);
            set(roiWin, 'HandleVisibility', 'on');
            if customDimScale
                roiImage = imagesc(dimScale(xySel(2),:),dimScale(xySel(1),:),currIm{1},'AlphaData',1);
            else
                roiImage = imagesc(currIm{1});
            end
            roiLine.roi = line(0,0,'Parent', roiAxes,'Color','w');
            axis image;
            set(roiWin, 'HandleVisibility', 'off');
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function [cn,on] = inpoly_mv(p,node,edge,TOL)
        if debugMatVis > 1, debugMatVisFcn(1); end
        % Function of Darren Engwirda from Matlab File exchange server. See below for credit.
        %  INPOLY: Point-in-polygon testing.
        %
        % Determine whether a series of points lie within the bounds of a polygon
        % in the 2D plane. General non-convex, multiply-connected polygonal
        % regions can be handled.
        %
        % SHORT SYNTAX:
        %
        %   in = inpoly(p,node);
        %
        %   p   : The points to be tested as an Nx2 array [x1 y1; x2 y2; etc].
        %   node: The vertices of the polygon as an Mx2 array [X1 Y1; X2 Y2; etc].
        %         The standard syntax assumes that the vertices are specified in
        %         consecutive order.
        %
        %   in  : An Nx1 logical array with IN(i) = TRUE if P(i,:) lies within the
        %         region.
        %
        % LONG SYNTAX:
        %
        %  [in,on] = inpoly(p,node,edge);
        %
        %  edge: An Mx2 array of polygon edges, specified as connections between
        %        the vertices in NODE: [n1 n2; n3 n4; etc]. The vertices in NODE
        %        do not need to be specified in connsecutive order when using the
        %        extended syntax.
        %
        %  on  : An Nx1 logical array with ON(i) = TRUE if P(i,:) lies on a
        %        polygon edge. (A tolerance is used to deal with numerical
        %        precision, so that points within a distance of
        %        eps^0.8*norm(node(:),inf) from a polygon edge are considered "on"
        %        the edge.
        %
        % EXAMPLE:
        %
        %   polydemo;       % Will run a few examples
        %
        % See also INPOLYGON
        
        % The algorithm is based on the crossing number test, which counts the
        % number of times a line that extends from each point past the right-most
        % region of the polygon intersects with a polygon edge. Points with odd
        % counts are inside. A simple implementation of this method requires each
        % wall intersection be checked for each point, resulting in an O(N*M)
        % operation count.
        %
        % This implementation does better in 2 ways:
        %
        %   1. The test points are sorted by y-value and a binary search is used to
        %      find the first point in the list that has a chance of intersecting
        %      with a given wall. The sorted list is also used to determine when we
        %      have reached the last point in the list that has a chance of
        %      intersection. This means that in general only a small portion of
        %      points are checked for each wall, rather than the whole set.
        %
        %   2. The intersection test is simplified by first checking against the
        %      bounding box for a given wall segment. Checking against the bbox is
        %      an inexpensive alternative to the full intersection test and allows
        %      us to take a number of shortcuts, minimising the number of times the
        %      full test needs to be done.
        %
        %   Darren Engwirda: 2005-2007
        %   Email          : d_engwirda@hotmail.com
        %   Last updated   : 23/11/2007 with MATLAB 7.0
        %
        % Problems or suggestions? Email me.
        
        %% ERROR CHECKING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if nargin<4
            TOL = 1.0e-12;
            if nargin<3
                edge = [];
                if nargin<2
                    error('Insufficient inputs');
                end
            end
        end
        nnode = size(node,1);
        if isempty(edge)                                                           % Build edge if not passed
            edge = [(1:nnode-1)' (2:nnode)'; nnode 1];
        end
        if size(p,2)~=2
            error('P must be an Nx2 array.');
        end
        if size(node,2)~=2
            error('NODE must be an Mx2 array.');
        end
        if size(edge,2)~=2
            error('EDGE must be an Mx2 array.');
        end
        if max(edge(:))>nnode || any(edge(:)<1)
            error('Invalid EDGE.');
        end
        
        %% PRE-PROCESSING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        n  = size(p,1);
        nc = size(edge,1);
        
        % Choose the direction with the biggest range as the "y-coordinate" for the
        % test. This should ensure that the sorting is done along the best
        % direction for long and skinny problems wrt either the x or y axes.
        dxy = max(p,[],1)-min(p,[],1);
        if dxy(1)>dxy(2)
            % Flip co-ords if x range is bigger
            p = p(:,[2,1]);
            node = node(:,[2,1]);
        end
        tol = TOL*min(dxy);
        
        % Sort test points by y-value
        [y,ii] = sort(p(:,2));
        x = p(ii,1);
        
        %% MAIN LOOP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cn = false(n,1);     % Because we're dealing with mod(cn,2) we don't have
        % to actually increment the crossing number, we can
        % just flip a logical at each intersection (faster!)
        on = cn;
        for kk = 1:nc         % Loop through edges
            
            % Nodes in current edge
            n1 = edge(kk,1);
            n2 = edge(kk,2);
            
            % Endpoints - sorted so that [x1,y1] & [x2,y2] has y1<=y2
            %           - also get xmin = min(x1,x2), xmax = max(x1,x2)
            y1 = node(n1,2);
            y2 = node(n2,2);
            if y1<y2
                x1 = node(n1,1);
                x2 = node(n2,1);
            else
                yt = y1;
                y1 = y2;
                y2 = yt;
                x1 = node(n2,1);
                x2 = node(n1,1);
            end
            if x1>x2
                xmin = x2;
                xmax = x1;
            else
                xmin = x1;
                xmax = x2;
            end
            
            % Binary search to find first point with y<=y1 for current edge
            if y(1)>=y1
                start = 1;
            elseif y(n)<y1
                start = n+1;
            else
                lower = 1;
                upper = n;
                for jj = 1:n
                    start = round(0.5*(lower+upper));
                    if y(start)<y1
                        lower = start;
                    elseif y(start-1)<y1
                        break;
                    else
                        upper = start;
                    end
                end
            end
            
            % Loop through points
            for jj = start:n
                % Check the bounding-box for the edge before doing the intersection
                % test. Take shortcuts wherever possible!
                
                Y = y(jj);   % Do the array look-up once & make a temp scalar
                if Y<=y2
                    X = x(jj);   % Do the array look-up once & make a temp scalar
                    if X>=xmin
                        if X<=xmax
                            
                            % Check if we're "on" the edge
                            on(jj) = on(jj) || (abs((y2-Y)*(x1-X)-(y1-Y)*(x2-X))<tol);
                            
                            % Do the actual intersection test
                            if (Y<y2) && ((y2-y1)*(X-x1)<(Y-y1)*(x2-x1))
                                cn(jj) = ~cn(jj);
                            end
                            
                        end
                    elseif Y<y2   % Deal with points exactly at vertices
                        % Has to cross edge
                        cn(jj) = ~cn(jj);
                    end
                else
                    % Due to the sorting, no points with >y
                    % value need to be checked
                    break
                end
            end
            
        end
        
        % Re-index to undo the sorting
        cn(ii) = cn|on;
        on(ii) = on;
        if debugMatVis > 1, debugMatVisFcn(2); end
    end      % inpoly()
% Sets name of main gui to 'Busy' during long calculations (projections
% and RGB image calculations)
    function busy(in)
        if debugMatVis > 1, debugMatVisFcn(1); end
        % in == 1: matVis busy
        % in == 0: matVis not busy
        if in
            isBusy = isBusy + 1;
            set(gui, 'Name', 'Busy ...');
            drawnow;
        else
            isBusy = isBusy - 1;
            % if busy reaches 0 reset matVis name
            if isBusy < 1
                isBusy = 0;  % to avoid accidental negative values
                set(gui, 'Name', ['matVis: ', allNames]);
                %                 drawnow;
            end
        end
        if debugMatVis > 1, debugMatVisFcn(2); end
    end
    function out = scaleMinMax(in,scaleMin, scaleMax, minIn, maxIn, trunc)
      % scaleMin, scaleMax: final scale limits 
      % minIn, maxIn:       source scale limits
      % trunc:              truncation switch, false: no truncation, 1: truncation to scale limits, -1: outlayers are replaced by NaN
      if debugMatVis > 1, debugMatVisFcn(1); end
      if nargin == 1
        scaleMin = 0;
        scaleMax = 1;
      end
      if nargin < 4 || isempty(minIn)
        minIn = min(in(:));
        maxIn = max(in(:));
      end
      
      out = (in-minIn) / (maxIn-minIn) * (scaleMax - scaleMin) + scaleMin;
      
      if nargin == 6 && trunc
        switch trunc
          case -1
            out(out<scaleMin | out>scaleMax) = NaN;
          case 1
            out(out<scaleMin) = scaleMin;
            out(out>scaleMax) = scaleMax;
        end
      end
      if debugMatVis > 1, debugMatVisFcn(2); end
    end
    function myGcf = myGcf
        if debugMatVis > 2, debugMatVisFcn(1); end
        set(0, 'ShowHiddenHandles', 'on');
        myGcf = gcf;
        set(0, 'ShowHiddenHandles', 'off');
        if debugMatVis > 2, debugMatVisFcn(2); end
    end
    function myGca = myGca
        if debugMatVis > 2, debugMatVisFcn(1); end
        set(0, 'ShowHiddenHandles', 'on');
        myGca = gca;
        set(0, 'ShowHiddenHandles', 'off');
        if debugMatVis > 2, debugMatVisFcn(2); end
    end
    function closeTempWin(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        set(tempWin, 'HandleVisibility', 'on');
        delete(tempWin);
        tempWin = [];
        if debugMatVis, debugMatVisFcn(2); end
    end
    function openMatVisGuide(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        [w flag] = urlread('http://www.colors-and-contrasts.com/');
        if ~flag
            msgbox('Internet connection or webhost not available. Connect to the internet or try later.');
        else
            system('start http://www.colors-and-contrasts.com/Documents/matVisGuide.pdf');
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateMatVis(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        m = msgbox('Looking for updates ...');
        % [w flag] = urlread('http://www.colors-and-contrasts.com/');
        % if ~flag
        %     set([bt_matVisGuide bt_updateMatVis], 'Enable','off');
        % end
        [newMatVis flag] = urlread('http://www.colors-and-contrasts.com/Documents/matVis.m');
        if ~flag
            delete(m);
            msgbox('Internet connection or webhost not available. Connect to the internet or try later.');
        else
            w = findstr(newMatVis,'versionNumber');
            newVersionNumber = str2double(newMatVis(w(1)+15:w(1)+20));
            delete(m);
            if newVersionNumber > versionNumber
                questString = {'Save new version and rename old version','Save new version and overwrite old version','Cancel'};
                a = questdlg({'New version of matVis available.' ['Current version: ' num2str(versionNumber,'%6.3f')] ['New version: ' num2str(newVersionNumber)] 'Update is recommended!'},...
                    'matVis update',questString{1},questString{2},questString{3},questString{1});
                currMatVis = which('matVis');
                mvDir = fileparts(currMatVis);
                switch a
                    case questString{1}
                        movefile(currMatVis, [mvDir filesep 'matVis_v' num2str(floor(versionNumber)) '-' num2str(1000*(versionNumber-floor(versionNumber)),'%03.0f') '.m'],'f');
                        urlwrite('http://www.colors-and-contrasts.com/Documents/matVis.m',[mvDir filesep 'matVis.m']);
                    case questString{2}
                        urlwrite('http://www.colors-and-contrasts.com/Documents/matVis.m',[mvDir filesep 'matVis.m']);
                    case questString{3}
                      if debugMatVis, debugMatVisFcn(2); end
                      return;
                end
                % Opening data in updated matVis does not work for some reasons
                % (works though using breakpoints).
                %             a = questdlg('Open data in new matVis?','','Yes','No','Yes');
                %             if strcmp(a,'Yes')
                %                 rehash;
                %                 pause(0.1);
                %                 matVis(data{1});
                %                 closeGui;
                %             end
                a = questdlg('View log file?' , 'Log file','Yes','No','Yes');
                if strcmp(a,'Yes')
                    system('start http://www.colors-and-contrasts.com/Documents/matVis_logFile.pdf')
                end
            else
                a = questdlg(['Your matVis is up to date. Current version: ' num2str(versionNumber,'%6.3f')],'','View log-file','Close','Close');
                if strcmp(a,'View log-file')
                    system('start http://www.colors-and-contrasts.com/Documents/matVis_logFile.pdf')
                end
            end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
% Load new data into matVis. Handle to this function can be returned
% during initial calling of matVis, and can then be used to exchange
% the loaded data set(s).
    function loadNewData(d,mode,replaceIdx,varargin)
        if debugMatVis, debugMatVisFcn(1); end
        % If single data set is given as matrix instead of cell, convert
        % into cell first for consistency.
        if ~iscell(d)
            d = {d};
        end
        % Check if number of data sets, number of dimensions and length of
        % dimensions is consistent with currently loaded data. If not, give
%         % error message and return.
        if numel(d) ~= nMat
            errordlg('New data cannot be loaded into existing matVis GUI: Number of data sets not correct!');
        end
       
        %         if any(size(data{1}) ~= size(d{1}))
        %             errordlg('New data cannot be loaded into existing matVis GUI: Size of data set not correct!');
        %         end
        % Exchange data and clear input variable
        if nargin == 1
            mode = 'replaceData';
        end
        switch mode
            case 'replaceData'
                if numel(size(data{1})) ~= numel(size(d{1}))
                    errordlg('New data cannot be loaded into existing matVis GUI: Wrong number of dimensions!');
                end
                data =  d;
            case 'replaceFrame'
                if size(data{1},1) ~= size(d{1},1) || size(data{1},2) ~= size(d{1},2)
                    errordlg('New data cannot be loaded into existing matVis GUI: Wrong number of dimensions!');
                end
                data{1}(:,:,replaceIdx) =  d{1};
            case 'append'
                if numel(dim) == 2
                    errordlg('Cannot append frame to 2D data set!');  % Possible solution: restart matVis automatically. This would break the possibility to use the function handles, as the new matVis wouldn't be linked to the old handles.
                end
                if size(data{1},1) ~= size(d{1},1) || size(data{1},2) ~= size(d{1},2)
                    errordlg('New data cannot be loaded into existing matVis GUI: Wrong number of dimensions!');
                end
                data{1}(:,:,end+1:end+size(d{1},3)) =  d{1};
        end
        clear d;
        newDim = size(data{1});
        if any(newDim ~= dim)
            dim = newDim;
            set(sld_down, 'Value', 1);
            set(etxt_down, 'String', 1);
            for ii = 1:numel(dim)
                if get(sld(ii), 'Value') > dim(ii)
                    set(sld(ii), 'Value', dim(ii));
                    set(etxt(ii), 'String', dim(ii));
                end
                set(sld_up(ii), 'Value', dim(ii));
                set(etxt_up, 'String', dim(ii));
                set([sld(ii) sld_down(ii) sld_up(ii)], 'Max', dim(ii));
                set(dimSize(ii), 'String', [' / ', num2str(dim(ii))]);
                zoomVal(ii,:) = [1 dim(ii)];
                currPos(ii) = get(sld(ii), 'Value');
                plotXLim(ii,:) = [1 dim(ii)];
            end
            zoomValXY = [1,1,dim(xySel(2)),dim(xySel(1))];
            drawPlots;
            updateImages;
            updateZoom;
        else
            % Update matVis
            updateSelection;
            updatePlots;
        end
        if nargin > 1
            loadNewSettings(varargin{1:end});
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function loadNewSettings(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if isMatVis2DHist
            % Check if function was called from matVis' own updateColormap. In
            % this case it is likely (?) that it was called as an update of the
            % 2D histogram
            fc = dbstack;
            if any(strfind(fc(2).name,'updateColormap')) || any(strfind(fc(2).name,'updateAlphaContrast'))
                calledFrom2DHist = 1;
            end
        end
        %Read optional arguments
        for ii = 1:2:nargin-1
            optionalArgIdentifier = varargin{ii};
            val = varargin{ii+1};
            switch  optionalArgIdentifier
                case 'alphaMap'
                    alphaMap{1} = squeeze(val);  %#ok
                    updateImages;
                case 'updateParameter' % Parameters that should be updated according to properties of new data set
                    if ~iscell(val)
                        val = {val};
                    end
                    for jj=1:numel(val)
                        switch val{jj}
                            case 'sldMinMax'
                                for kk=1:nMat
                                    [maxVal(kk) maxValInd(kk)] = max(data{kk}(:));         %Maximum Data Value
                                    [minVal(kk) minValInd(kk)] = min(data{kk}(:));         %Minimum Data Value
                                    cmMinMax(kk,:) = [minVal(kk) maxVal(kk)]'; %Colormap Limits
                                end
                        end
                    end
                    set(sldLimMin, 'String',num2str(minVal(1)));
                    set(sldLimMax, 'String',num2str(maxVal(1)));
                    updateSldLim;
                    set(sldMin, 'Value', get(sldMin, 'Min'));
                    set(sldMax, 'Value', get(sldMax, 'Max'));
                    if ~rgbCount
                        updateGuiHistVal;
                    end
                    updateColormap;
                case 'sldMinMax' % Give fixed values to be used for limits and values of colormap (sliders)
                    minVal = val(:,1);
                    maxVal = val(:,2);
                    set(sldLimMin, 'String',num2str(minVal(1)));
                    set(sldLimMax, 'String',num2str(maxVal(1)));
                    updateSldLim;
                    set(sldMin, 'Value', get(sldMin, 'Min'));
                    set(sldMax, 'Value', get(sldMax, 'Max'));
                    if ~rgbCount
                        updateGuiHistVal;
                    end
                    updateColormap;
                case 'cmMinMax'
                    cmMinMax = val;
                    set(sldMin, 'Value', cmMinMax(1,1));
                    set(sldMax, 'Value', cmMinMax(1,2));
                    if ~rgbCount
                        updateGuiHistVal;
                    end
                    updateColormap;
                    if with2DHist
                        updateAlphaContrast;
                    end
                    calledFrom2DHist = 0;
                case 'zoomVal'
                    zoomVal = val;
                    zoomValXY([1,3]) = zoomVal(xySel(2),:);
                    zoomValXY([2,4]) = zoomVal(xySel(1),:);
                    updateZoom;
                    calledFrom2DHist = 0;
            end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function updateOutputStrct
        if debugMatVis, debugMatVisFcn(1); end
        % Figure handles
        out.figHandles.gui = gui;
        out.figHandles.imageWin = imageWin;
        out.figHandles.zoomWin = zoomWin;
        out.figHandles.plotWin = plotWin;
        % Function handles
        out.fctHandles.loadNewData = @loadNewData;
        out.fctHandles.loadNewSettings = @loadNewSettings;
        out.fctHandles.exportOutputStrct = @exportOutputStrct;
        out.fctHandles.closeGui = @closeGui;
        % Settings
        out.settings.zoom = zoomVal;
        out.settings.zoomXY = zoomValXY;
        % Data (file or variable) name
        out.dataName = varName;
        % ROIs
        if ~isempty(roiWin)
            out.figHandles.roiWin = roiWin;
            out.roiManager.roiList = roiList;
            out.roiManager.roiSelection = get(roiListbox, 'Value');
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
    function newOutputStrct = exportOutputStrct
        if debugMatVis, debugMatVisFcn(1); end
        updateOutputStrct;
        newOutputStrct = out;
        if debugMatVis, debugMatVisFcn(2); end
    end
%% Movie recording functions
    function vidGenerator(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        if ~get(tbRecord, 'Value')
            set(tbPlotsYLim, 'Value',0);
            vidGenerator_showHide(0,0,'hide');
            set(tbLinkWin, 'Enable','on');
            if get(tbLinkWin, 'Value')
                set(imageWin, 'ResizeFcn', {@resizeWin, 'image'});
                set(zoomWin, 'ResizeFcn', {@resizeWin, 'zoom'});
                set(plotWin, 'ResizeFcn', {@resizeWin, 'plot'});
            end
        elseif ~isfield(movdata,'gui');
            % Some general paramters
            gp = get(gui, 'Position');
            pos0=[gp(1)+gp(3)+20 gp(2)]; % Position of vidGenerator windows x0 y0
            res = 3; % internal parameter print resolution for 'hardcopy'
            set([imageWin, zoomWin, plotWin], 'ResizeFcn','');
            set(tbLinkWin, 'Enable','off');
            set(tbPlotsYLim, 'Value',1);
            vidGenerator_initialize;
        else
            set([imageWin, zoomWin, plotWin], 'ResizeFcn','');
            set(tbLinkWin, 'Enable','off');
            set(tbPlotsYLim, 'Value',1);
            vidGenerator_showHide(0,0,'show');
        end
            
        
        function vidGenerator_initialize
            if debugMatVis, debugMatVisFcn(1); end
            % Test for available movdata.codecs
            profiles=VideoWriter.getProfiles();
            movdata.codecs=[];
            for iii=1:length(profiles);
                movdata.codecs{iii,1}=profiles(iii).Name;
            end
            
            % Test for ffmpeg
            [status,result] = system('ffmpeg');
            if ~strcmp(result(1:14),'ffmpeg version');
                hffmpegwarn = questdlg({'Could not find FFMPEG in the system path.';'Codecs will be limited to the ones which come with Matlab.';'Press Ok to visit Source.'},'FFMPEG Warning','Ok','Skip Information','Ok');
                switch hffmpegwarn
                    case 'Ok'
                        [status result] = system('start http://www.ffmpeg.org');
                end
            else
                if ~any(strcmp(movdata.codecs,'MPEG-4'));
                    movdata.codecs{end+1,1}='FFMPEG > MPEG-4';
                end
                movdata.codecs{end+1,1}='FFMPEG > Theora';
            end
            
            % Test for export_fig
            % Stephan: Not clear what's it for, since export_fig is not
            % used !?
%             exitstartup = false;
%             if ~exist('export_fig','file');
%                 hexport_figwarn = msgbox({'export_fig.m is required for VidGenerator.';'Press Ok to visit Source.'},'export_fig Warning','warn');
%                 uiwait(hexport_figwarn);
%                 [status result] = system('start http://www.mathworks.se/matlabcentral/fileexchange/23629-exportfig');
%                 exitstartup = true;
%             end
            
            % Set movdata.FilterSpecs
            movdata.FilterSpec = [];
            if any(strcmp(movdata.codecs,'Archival'));movdata.FilterSpec{end+1,1}='*.mj2';movdata.FilterSpec{end,2}='Archival (*.mj2)';end
            if any(strcmp(movdata.codecs,'Motion JPEG 2000'));movdata.FilterSpec{end+1,1}='*.mj2';movdata.FilterSpec{end,2}='Motion JPEG 2000 (*.mj2)';end
            if any(strcmp(movdata.codecs,'Uncompressed AVI'));movdata.FilterSpec{end+1,1}='*.avi';movdata.FilterSpec{end,2}='Uncompressed AVI (*.avi)';end
            if any(strcmp(movdata.codecs,'Motion JPEG AVI'));movdata.FilterSpec{end+1,1}='*.avi';movdata.FilterSpec{end,2}='Motion JPEG AVI (*.avi)';end
            if any(strcmp(movdata.codecs,'MPEG-4'));movdata.FilterSpec{end+1,1}='*.mp4';movdata.FilterSpec{end,2}='MPEG-4 (*.mp4)';end
            if any(strcmp(movdata.codecs,'FFMPEG > MPEG-4'));movdata.FilterSpec{end+1,1}='*.mp4';movdata.FilterSpec{end,2}='MPEG-4 (*.mp4)';end
            if any(strcmp(movdata.codecs,'FFMPEG > Theora'));movdata.FilterSpec{end+1,1}='*.ogg';movdata.FilterSpec{end,2}='Theora (*.ogg)';end
            
            % Set File Name from presettings
            if ~isfield(movdata,'set.movname'); movdata.set.movname = 'movie_file_name';end
            if ~isfield(movdata,'set.movdir'); movdata.set.movdir = pwd;end
            if ~isfield(movdata,'set.movcodec'); movdata.set.movcodec = 'Uncompressed AVI';end
            movdata.set.fileext = [];
            for iii=1:length(movdata.FilterSpec);
                if regexp(movdata.FilterSpec{iii,2},movdata.set.movcodec,'ONCE')==1;
                    movdata.set.fileext = movdata.FilterSpec{iii,1}(2:end);
                    break
                end
            end
            movdata.set.movstring=[movdata.set.movdir '\' movdata.set.movname movdata.set.fileext];
            
            % Set used codec from presettings
            if ~isfield(movdata,'set.codecval');
                for iii=1:length(movdata.codecs);
                    if ~isempty(regexp(movdata.codecs{iii},movdata.set.movcodec,'ONCE'));
                        movdata.set.codecval = iii;
                        break
                    end
                end
            end
            if strcmp(movdata.set.movcodec,'Uncompressed AVI'); movdata.set.movquality = 100; end
            
            % set Quality/compression from presettings
            if strcmp(movdata.set.movcodec,'Motion JPEG 2000') || strcmp(movdata.set.movcodec,'Archival');
                movdata.set.qualstring = 'Compress:';
                movdata.set.qualval = movdata.set.movcompress;
            else
                movdata.set.qualstring = 'Quality [%]:';
                movdata.set.qualval = movdata.set.movquality;
            end
            vidGenerator_generateGUI;
            if debugMatVis, debugMatVisFcn(2); end
        end
        function vidGenerator_generateGUI(varargin)
          if debugMatVis, debugMatVisFcn(1); end
            if isfield(movdata,'gui') && ~isempty(movdata.gui.hvidgen)
                delete(movdata.gui.hvidgen);
                delete(movdata.prev.hprev);
                movdata.gui.hvidgen = [];
            end
            % Find visible windows
            oldMovParts = movdata.parts;
            movdata.parts = [];
            for iii = 1:nMat
                if get(tbWin(1), 'Value')
                    movdata.parts(end+1) = imageWin(iii);
                end
                if get(tbWin(2), 'Value')
                    movdata.parts(end+1) = zoomWin(iii);
                end
                if get(tbWin(3), 'Value')
                    movdata.parts(end+1) = plotWin(iii);
                end
            end
            if strcmp(get(histWin,'Visible'),'on')
                movdata.parts(end+1) = histWin;
            end
            if with2DHist
                if strcmp(get(matVis2DHist.figHandles.zoomWin,'Visible'),'on')
                    movdata.parts(end+1) = matVis2DHist.figHandles.zoomWin;
                end
                if strcmp(get(matVis2DHist.figHandles.imageWin,'Visible'),'on')
                    movdata.parts(end+1) = matVis2DHist.figHandles.imageWin;
                end
                if strcmp(get(matVis2DHist.figHandles.plotWin,'Visible'),'on')
                    movdata.parts(end+1) = matVis2DHist.figHandles.plotWin;
                end
            end
            if ~isequal(oldMovParts, movdata.parts)
                movdata.set.parts = [];
                for iii=1:length(movdata.parts);
                    movdata.set.parts{iii}.handle = movdata.parts(iii); % include input figure handle in presetting structure
                    movdata.set.parts{iii}.name = get(movdata.parts(iii),'name');
                    movdata.set.parts{iii}.status = 1;
                    set(movdata.set.parts{iii}.handle, 'Color',movdata.set.movPartsbg);
                    for kkk=1:numel(movdata.set.parts)
                        al = findobj(movdata.set.parts{kkk}.handle ,'Type','Axes');
                        for jjj=1:numel(al)
                            set(al(jjj),'XColor',movdata.set.movLabelColor,'YColor',movdata.set.movLabelColor);
                            set(get(al(jjj),'XLabel'),'Color',movdata.set.movLabelColor);
                            set(get(al(jjj),'YLabel'),'Color',movdata.set.movLabelColor);
                        end
                    end
                end
                nRows = ceil(length(movdata.parts)/2);
                if length(movdata.parts)==1
                    nCols = 1;
                else
                    nCols = 2;
                end
                ct = 0;
                defaultTitleHeight = 0.9;
                for kkk = 1:nRows
                    for jjj=1:nCols
                        ct = ct +1;
                        movdata.set.parts{ct}.pos = [(jjj-1)/nCols defaultTitleHeight*(1-kkk/nRows) 1/nCols defaultTitleHeight/nRows];
                    end
                end
                if mod(length(movdata.parts),2) && nRows > 1
                    movdata.set.parts(ct) = [];
                    movdata.set.parts{ct-1}.pos = [0 0 1 defaultTitleHeight/nRows];
                end
            end
            %% Evaluate existing movie parts
            for iii=1:length(movdata.parts);
                set(movdata.parts(iii),'Units','pixels');
                movdata.set.parts{iii}.origpos = get(movdata.parts(iii),'Position');
            end
           %% Open GUI
           nparts=length(movdata.parts);
           movdata.gui.hvidgen=figure('Units','pixels','Color',[.94 .94 .94],'MenuBar','none','HandleVisibility', 'off',...
                'Name','Video Generator','NumberTitle','off','Position',[pos0(1) pos0(2) 495 259.3+nparts*15.6],...
                'Resize','off','tag','movdata.gui.hvidgen', 'CloseRequestFcn', {@vidGenerator_showHide,'hide'});  %pos0(2)-(259.3+nparts*15.6)
            % Stephan: Not sure about the function of exitstartup
            %                 if exitstartup;
            %                     vidGenerator_hide;
            %                 else
            % Movie Input Panel
            movdata.gui.hinput=uipanel('Parent',movdata.gui.hvidgen,'Units','pixels','Title','Movie Input',...
                'Position',[11 141.4 476 114.4+nparts*15.6],'tag','movdata.gui.hinput');
            row1=85+nparts*15.6; % rowposition for panel 1
            row2=65+nparts*15.6-3;
            row3=45+nparts*15.6-3;
            % Movie Size
            movdata.gui.hinput1 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','HorizontalAlignment','left',...
                'Position',[6 row1 141 13],'String','Movie Frame Size (in pixels):',...
                'Style','text','tag','movdata.gui.hinput1');
            movdata.gui.hinput2 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','BackgroundColor',[1 1 1],...
                'Callback',@vidGenerator_preview,'HorizontalAlignment','right','Position',[151 row1-2.5 50 15.6],...
                'String',movdata.set.movsize(1),'Style','edit','tag','movdata.gui.hinput2');
            movdata.gui.hinput3 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','HorizontalAlignment','left',...
                'Position',[201 row1 10 13],'String','x','Style','text','tag','movdata.gui.hinput3');
            movdata.gui.hinput4 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','BackgroundColor',[1 1 1],...
                'Callback',@vidGenerator_preview,'HorizontalAlignment','right','Position',[211 row1-2.5 50 15.6],...
                'String',movdata.set.movsize(2),'Style','edit','tag','movdata.gui.hinput4');
            % Parts background color
            movdata.gui.hinput12 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','HorizontalAlignment','right',...
                'Position',[310 row2-1.3+5 60 14.3],'String','Parts color:','Style','text','tag','movdata.gui.hinput5');
            movdata.gui.hinput13 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','BackgroundColor',movdata.set.movPartsbg,...
                'Callback',{@vidGenerator_bgcolor,'axes'},'Position',[456-80 row2-2.5+7 14 14],'String',[],'Style','pushbutton','tag','movdata.gui.hinput6');
            % Figure background color
            movdata.gui.hinput5 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','HorizontalAlignment','right',...
                'Position',[310 row1-1.3 60 14.3],'String','Figure color:','Style','text','tag','movdata.gui.hinput5');
            movdata.gui.hinput6 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','BackgroundColor',movdata.set.movbg,...
                'Callback',{@vidGenerator_bgcolor,'figure'},'Position',[456-80 row1-2.5+2 14 14],'String',[],'Style','pushbutton','tag','movdata.gui.hinput6');
            % Title color
            movdata.gui.hinput14 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','HorizontalAlignment','right',...
                'Position',[390 row1-1.3 60 14.3],'String','Title color:','Style','text','tag','movdata.gui.hinput5');
            movdata.gui.hinput15 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','BackgroundColor',movdata.set.movTitleColor,...
                'Callback',{@vidGenerator_bgcolor,'title'},'Position',[456 row1-2.5+2 14 14],'String',[],'Style','pushbutton','tag','movdata.gui.hinput6');
            % Label color
            movdata.gui.hinput16 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','HorizontalAlignment','right',...
                'Position',[390 row2-1.3+5 60 14.3],'String','Label color:','Style','text','tag','movdata.gui.hinput5');
            movdata.gui.hinput17 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','BackgroundColor',movdata.set.movLabelColor,...
                'Callback',{@vidGenerator_bgcolor,'label'},'Position',[456 row2-2.5+7 14 14],'String',[],'Style','pushbutton','tag','movdata.gui.hinput6');
            % Update preview
            movdata.gui.hinput7 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','Callback',@vidGenerator_generateGUI,...
                'Position',[370 row3-2.5+2 100 20],'String','Update preview','Style','pushbutton','tag','movdata.gui.hinput7','FontSize',8);
            % Scaling
            movdata.gui.hinput8 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','Position',[6 mean([row2,row3])+3 50 13],...
                'String','Scaling:','Style','text','tag','movdata.gui.hinput8','HorizontalAlignment','left');
            movdata.gui.hinput9 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','BackgroundColor',[1 1 1],...
                'Callback',@vidGenerator_preview,'HorizontalAlignment','right','Position',[50 mean([row2,row3]) 60 20.8],...
                'String',{'100 %';'75 %';'50 %'},'Style','popup','tag','movdata.gui.hinput9','value',movdata.set.prevscale);
            % Units
            movdata.gui.hinput10 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','HorizontalAlignment','left',...
                'Position',[115 mean([row2,row3])+3 122 13],'String','Units for Movie Parts:','Style','text','tag','movdata.gui.hinput10');
            movdata.gui.hinput11 = uicontrol('Parent',movdata.gui.hinput,'Units','pixels','BackgroundColor',[1 1 1],...
                'Callback',@vidGenerator_preview,'Position',[225 mean([row2,row3]) 80 20.8],'String',{'normalized'},...
                'Style','popup','tag','movdata.gui.hinput11','value',movdata.set.partunits);
            
            % Panel - movie parts
            movdata.gui.hparts = uipanel('Parent',movdata.gui.hinput,'Units','pixels','Title','Available Movie Parts',...
                'Position',[6 7.5 462.5 39+nparts*15.6],'tag','movdata.gui.hparts');
            % input: title
            coltitle=7.8+nparts*15.6;
            if movdata.set.titlestatus==true;
                enablestat = 'on';
            else
                enablestat = 'off';
            end
            movdata.gui.hparts1 = uicontrol('Parent',movdata.gui.hparts,'Units','pixels','Callback',@vidGenerator_enable,...
                'Position',[9 coltitle-1.3 47 13],'String',[],...
                'Style','checkbox','Value',movdata.set.titlestatus,'tag','movdata.gui.hparts1');
            movdata.gui.hparts2 = uicontrol('Parent',movdata.gui.hparts,'BackgroundColor',[1 1 1],'Callback',@vidGenerator_preview,...
                'Units','pixels','HorizontalAlignment','left','Position',[26 coltitle-2.5 130 15.6],...
                'String',movdata.set.title,'Style','edit','tag','movdata.gui.hparts2','enable',enablestat);
            movdata.gui.hparts3 = uicontrol('Parent',movdata.gui.hparts,'Units','pixels','enable',enablestat,...
                'HorizontalAlignment','left','Position',[166 coltitle 93 13],...
                'String','Position in Frame:','Style','text','tag','movdata.gui.hparts3');
            movdata.gui.hparts4 = uicontrol('Parent',movdata.gui.hparts,'Units','pixels','BackgroundColor',[1 1 1],...
                'Callback',@vidGenerator_preview,'HorizontalAlignment','right','Position',[253.5 coltitle-2.3 50 15.6],...
                'String',movdata.set.titlepos(1),'Style','edit','tag','movdata.gui.hparts4','enable',enablestat);
            movdata.gui.hparts5 = uicontrol('Parent',movdata.gui.hparts,'Units','pixels','BackgroundColor',[1 1 1],...
                'Callback',@vidGenerator_preview,'HorizontalAlignment','right','Position',[303.5 coltitle-2.3 50 15.6],...
                'String',movdata.set.titlepos(2),'Style','edit','tag','movdata.gui.hparts5','enable',enablestat);
            movdata.gui.hparts6 = uicontrol('Parent',movdata.gui.hparts,'Units','pixels',...
                'HorizontalAlignment','left','Position',[358.5 coltitle 80 13],...
                'String','FontSize (pt):','Style','text','tag','movdata.gui.hparts6','enable',enablestat);
            movdata.gui.hparts7 = uicontrol('Parent',movdata.gui.hparts,'Units','pixels','BackgroundColor',[1 1 1],...
                'Callback',@vidGenerator_preview,'HorizontalAlignment','right','Position',[428.5 coltitle-2.3 25 15.6],...
                'String',movdata.set.titlesize,'Style','edit','tag','movdata.gui.hparts7','enable',enablestat);
            % Loop over image inputs
            for iii=1:nparts;
                col=coltitle-15.6*iii;
                if movdata.set.parts{iii}.status==true;
                    enablestat = 'on';
                else
                    enablestat = 'off';
                end
                movdata.gui.hpartsim1{iii} = uicontrol('Parent',movdata.gui.hparts,'Units','pixels','Callback',@vidGenerator_enable,...
                    'Position',[9 col-1.3 47 13],'String',[],...
                    'Style','checkbox','Value',movdata.set.parts{iii}.status,'tag',['movdata.gui.hpartsim1{' num2str(iii) '}']);
                movdata.gui.hpartsim2{iii} = uicontrol('Parent',movdata.gui.hparts,'enable',enablestat,...
                    'Units','pixels','HorizontalAlignment','left','Position',[26 col 130 13],...
                    'String',movdata.set.parts{iii}.name,'Style','text','tag',['movdata.gui.hpartsim2{' num2str(iii) '}']);
                movdata.gui.hpartsim3{iii} = uicontrol('Parent',movdata.gui.hparts,'enable',enablestat,'Units','pixels',...
                    'HorizontalAlignment','left','Position',[166 col 93 13],...
                    'String','Position in Frame:','Style','text','tag',['movdata.gui.hpartsim3{' num2str(iii) '}']);
                movdata.gui.hpartsim4{iii} = uicontrol('Parent',movdata.gui.hparts,'enable',enablestat,'Units','pixels','BackgroundColor',[1 1 1],...
                    'Callback',@vidGenerator_preview,'HorizontalAlignment','right','Position',[253.5 col-2.5 50 15.6],...
                    'String',movdata.set.parts{iii}.pos(1),'Style','edit','tag',['movdata.gui.hpartsim4{' num2str(iii) '}']);
                movdata.gui.hpartsim5{iii} = uicontrol('Parent',movdata.gui.hparts,'enable',enablestat,'Units','pixels','BackgroundColor',[1 1 1],...
                    'Callback',@vidGenerator_preview,'HorizontalAlignment','right','Position',[303.5 col-2.5 50 15.6],...
                    'String',movdata.set.parts{iii}.pos(2),'Style','edit','tag',['movdata.gui.hpartsim5{' num2str(iii) '}']);
                movdata.gui.hpartsim6{iii} = uicontrol('Parent',movdata.gui.hparts,'enable',enablestat,'Units','pixels','BackgroundColor',[1 1 1],...
                    'Callback',@vidGenerator_preview,'HorizontalAlignment','right','Position',[353.5 col-2.5 50 15.6],...
                    'String',movdata.set.parts{iii}.pos(3),'Style','edit','tag',['movdata.gui.hpartsim6{' num2str(iii) '}']);
                movdata.gui.hpartsim7{iii} = uicontrol('Parent',movdata.gui.hparts,'enable',enablestat,'Units','pixels','BackgroundColor',[1 1 1],...
                    'Callback',@vidGenerator_preview,'HorizontalAlignment','right','Position',[403.5 col-2.5 50 15.6],...
                    'String',movdata.set.parts{iii}.pos(4),'Style','edit','tag',['movdata.gui.hpartsim7{' num2str(iii) '}']);
            end
            % Movie Output Panel
            movdata.gui.houtput=uipanel('Parent',movdata.gui.hvidgen,'Units','pixels','Title','Movie Output',...
                'Position',[11 37.4 476 104],'tag','movdata.gui.houtput');
            movdata.gui.houtput1 = uicontrol('Parent',movdata.gui.houtput,'Units','pixels','BackgroundColor',[1 1 1],...
                'Callback',@vidGenerator_updatefile,'HorizontalAlignment','left','Position',[9 65 395 15.6],...
                'String',[movdata.set.movdir '\' movdata.set.movname movdata.set.fileext],'Style','edit','tag','movdata.gui.houtput1');
            movdata.gui.houtput2 = uicontrol('Parent',movdata.gui.houtput,'Units','pixels','Callback',@vidGenerator_savefile,...
                'Position',[405 64.3 60 18],'String','Save as','Style','pushbutton','tag','movdata.gui.houtput2');
            movdata.gui.houtput3 = uicontrol('Parent',movdata.gui.houtput,'Units','pixels','HorizontalAlignment','left',...
                'Position',[10 42 38 13],'String','Codec:','Style','text','tag','movdata.gui.houtput3');
            movdata.gui.houtput4 = uicontrol('Parent',movdata.gui.houtput,'Units','pixels','BackgroundColor',[1 1 1],...
                'Callback',@vidGenerator_choosecodec,'Position',[49 39 202 20],'String',movdata.codecs,...
                'Style','popup','tag','movdata.gui.houtput4','value',movdata.set.codecval);
            movdata.gui.houtput5 = uicontrol('Parent',movdata.gui.houtput,'Units','pixels','HorizontalAlignment','left',...
                'Position',[261 42 77 13],'String','FrameRate [f/s]:','Style','text','tag','movdata.gui.houtput5');
            movdata.gui.houtput6 = uicontrol('Parent',movdata.gui.houtput,'Units','pixels','BackgroundColor',[1 1 1],...
                'Callback',@vidGenerator_FCQ,'HorizontalAlignment','right','Position',[342 38 30 22],...
                'String',movdata.set.movFrameRate,'Style','edit','tag','movdata.gui.houtput6');
            movdata.gui.houtput7 = uicontrol('Parent',movdata.gui.houtput,'Units','pixels','HorizontalAlignment','left',...
                'Position',[377 42 62 13],'String',movdata.set.qualstring,'Style','text','tag','movdata.gui.houtput7');
            movdata.gui.houtput8 = uicontrol('Parent',movdata.gui.houtput,'Units','pixels','BackgroundColor',[1 1 1],...
                'Callback',@vidGenerator_FCQ,'HorizontalAlignment','right','Position',[434 38 30 22],...
                'String',movdata.set.qualval,'Style','edit','tag','movdata.gui.houtput8');
            if strcmp(movdata.set.movcodec,'Uncompressed AVI');set([movdata.gui.houtput7 movdata.gui.houtput8],'String',num2str(100),'Enable','off');end
            movdata.gui.houtput9 = uicontrol('Parent',movdata.gui.houtput,'Units','pixels','Enable','off',...
                'HorizontalAlignment','left','Position',[10 15 70 13],'String','ffmpeg:',...
                'Style','text','tag','movdata.gui.houtput9');
            movdata.gui.houtput10 = uicontrol('Parent',movdata.gui.houtput,'Units','pixels','BackgroundColor',[1 1 1],...
                'Callback',[],'Enable','off','HorizontalAlignment','left',...
                'Position',[50 12.5 338 15.6],'String',[],'Style','edit','tag','movdata.gui.houtput10');
            movdata.gui.houtput11 = uicontrol('Parent',movdata.gui.houtput,'Units','pixels','Callback',[],...
                'Enable','off','Position',[395 11.5 70 18],'String','More Options',...
                'Style','pushbutton','tag','movdata.gui.houtput11');
            movdata.gui.hvidgen1 = uicontrol('Parent',movdata.gui.hvidgen,'Units','pixels','Callback',@vidGenerator_record,...
                'CData',icon_Record,'Position',[10 10 22 22],'String',[],'Style','togglebutton','tag','movdata.gui.hvidgen1','value',movdata.rec);
%             movdata.gui.hvidgen2 = uicontrol('Parent',movdata.gui.hvidgen,'Units','pixels','Callback',@vidGenerator_showHide,...
%                 'Position',[60 10 50 22],'String','Hide','Style','pushbutton','tag','movdata.gui.hvidgen2');
            movdata.gui.hvidgen3 = uicontrol('Parent',movdata.gui.hvidgen,'Units','pixels','HorizontalAlignment','left',...
                'Position',[35 12 500 15],'String','Adjust settings, then press button to activate recording.',...
                'Style','text','tag','movdata.gui.hvidgen3');
            % %                 end
            %             else
            %                 vidGenerator_hide
            adjustGuiSize(movdata.gui.hvidgen);
            vidGenerator_preview; 
            if debugMatVis, debugMatVisFcn(2); end
        end
        %% Draw Preview
        function vidGenerator_preview(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            movdata.set.partunits = get(movdata.gui.hinput11,'value');
            prevunits = get(movdata.gui.hinput11,'String');
            prevunits = prevunits{movdata.set.partunits};
            movdata.set.movsize(1) = str2num(get(movdata.gui.hinput2,'String'));
            movdata.set.movsize(2) = str2num(get(movdata.gui.hinput4,'String'));
            movdata.set.prevscale = get(movdata.gui.hinput9,'Value');
            switch movdata.set.prevscale
                case 1
                    prevscale=1;
                case 2
                    prevscale=.75;
                case 3
                    prevscale=.5;
            end
            movdata.set.movbg = get(movdata.gui.hinput6,'BackgroundColor');
            scalestring = get(movdata.gui.hinput9,'String');
            scalestring = scalestring{get(movdata.gui.hinput9,'Value')};
            prevname = [movdata.set.movname ' - ' num2str(movdata.set.movsize(1)) ' x ' num2str(movdata.set.movsize(2)) ' - ' scalestring];
            prevpos = [pos0(1)+512 pos0(2) movdata.set.movsize(1)*prevscale movdata.set.movsize(2)*prevscale];   %pos0(2)-movdata.set.movsize(2)*prevscale
            
%             if isfield(movdata.prev,hprev)
%                 pos=get(movdata.prev.hprev,'Position');
%                 prevpos(1)=pos(1);
%                 prevpos(2)=pos(2)+pos(4)-movdata.set.movsize(2)*prevscale;
%                 close(movdata.prev.hprev);
%             end
            % create preview window
            if ~(isfield(movdata,'prev') && ishandle(movdata.prev.hax))
                movdata.prev.hprev = figure('Position',prevpos,'Color',movdata.set.movbg,'MenuBar','none','HandleVisibility', 'off',...
                    'Units','pixels','NumberTitle','off','Resize','off','Name',['Movie preview. File name: ' prevname],'tag','movdata.prev.hprev', 'CloseRequestFcn', {@vidGenerator_showHide,'hide'});
                movdata.prev.hax = axes('Position',[0 0 1 1],'Color',movdata.set.movbg,'XTick',[],'YTick',[],'Parent',movdata.prev.hprev);
                axis(movdata.prev.hax,'off');
            end
            
            %         % create title line in preview
            %         if get(movdata.gui.hparts1,'Value')==1;
            %           movdata.set.title = get(movdata.gui.hparts2,'String');
            %           movdata.set.titlepos = [str2num(get(movdata.gui.hparts4,'String')) str2num(get(movdata.gui.hparts5,'String'))];
            %           movdata.set.titlesize = str2num(get(movdata.gui.hparts7,'String'));
            %           prevtitlepos = [0 movdata.set.titlepos(2)-(movdata.set.titlesize*1.5/movdata.set.movsize(2))/2 1 (movdata.set.titlesize*1.5/movdata.set.movsize(2))];
            %           movdata.prev.prevtitleaxes = axes('Parent',movdata.prev.hprev,'Position',prevtitlepos,'Color',movdata.set.movbg);
            %           movdata.prev.prevtitletext = text('Position',[movdata.set.titlepos(1) .5],'Parent',movdata.prev.prevtitleaxes,'Color',1-movdata.set.movbg,'HorizontalAlignment','center','String',movdata.set.title,'FontSize',movdata.set.titlesize*prevscale);
            %           axis off;
            %         end
            
            % Create Title Data
            if movdata.set.titlestatus == true;
                movdata.titleframe.handle = figure('Position',[1 1 movdata.set.movsize(1) movdata.set.titlesize*1.5],...
                    'MenuBar','none','Units','pixels','NumberTitle','off','Color',movdata.set.movbg,...
                    'Resize','off','Name','Frame for Title','tag','movdata.frame','Visible','off');
                movdata.titleframe.titleaxes = axes('Parent',movdata.titleframe.handle,'Position',[0 0 1 1],'Color',movdata.set.movbg);
                movdata.titleframe.titletext = text('Parent',movdata.titleframe.titleaxes,'Position',[movdata.set.titlepos(1) .5],'Color', movdata.set.movTitleColor,'HorizontalAlignment','center','String',movdata.set.title,'FontSize',movdata.set.titlesize,'FontWeight','bold');
                axis off;
                drawnow;
                if movdata.fastCaptureMode == true;
                    movdata.titleframe = vidGenerator_myHardcopyFigureSet(movdata.titleframe, res);
                    movdata.titleframe.CData = vidGenerator_myHardcopy(movdata.titleframe, res);
                    vidGenerator_myHardcopyFigureReset(movdata.titleframe)
                else
                    movdata.titleframe.CData = export_fig(movdata.titleframe.handle,'-nocrop','-q100','-RGB');%,'-painters'
                    if ndims(movdata.titleframe.CData)==2;
                        movdata.titleframe.CData = repmat(movdata.titleframe.CData,[1 1 3]);
                    end
                end
                close(movdata.titleframe.handle);
            end
            
            %         % create parts in preview (and adjust original)
            %         for nprevparts=1:nparts;
            %           if get(movdata.gui.hpartsim1{nprevparts},'Value')==1;
            %             movdata.set.parts{nprevparts}.pos = [str2num(get(movdata.gui.hpartsim4{nprevparts},'String'))...
            %               str2num(get(movdata.gui.hpartsim5{nprevparts},'String'))...
            %               str2num(get(movdata.gui.hpartsim6{nprevparts},'String'))...
            %               str2num(get(movdata.gui.hpartsim7{nprevparts},'String'))];
            %             if get(movdata.gui.hinput11,'Value') == 1;
            %               part{nprevparts}.abs = movdata.set.parts{nprevparts}.pos.*[movdata.set.movsize(1) movdata.set.movsize(2) movdata.set.movsize(1) movdata.set.movsize(2)];
            %             else
            %               part{nprevparts}.abs = movdata.set.parts{nprevparts}.pos;
            %             end
            %             % reset figure Position
            %             movdata.set.parts{nprevparts}.origpos = get(movdata.set.parts{nprevparts}.handle,'Position');
            %             movdata.set.parts{nprevparts}.origpos(3:4) = part{nprevparts}.abs(3:4);
            %             set(movdata.set.parts{nprevparts}.handle,'Units','pixels','Position',movdata.set.parts{nprevparts}.origpos);
            %             figure(movdata.set.parts{nprevparts}.handle);
            %             if movdata.fastCaptureMode == true;
            %               movdata.set.parts{nprevparts} = vidGenerator_myHardcopyFigureSet(movdata.set.parts{nprevparts}, res);
            %               movdata.set.parts{nprevparts}.CData = vidGenerator_myHardcopy(movdata.set.parts{nprevparts}, res);
            %               vidGenerator_myHardcopyFigureReset(movdata.set.parts{nprevparts})
            %             else
            %               % get image
            %               movdata.set.parts{nprevparts}.CData = export_fig(movdata.set.parts{nprevparts}.handle,'-nocrop','-q100','-RGB');
            %             end
            %             % draw preview
            %             movdata.prev.hprevpart{nprevparts} = axes(...
            %               'Units',prevunits,'Position',movdata.set.parts{nprevparts}.pos,...
            %               'XColor',1-movdata.set.movbg,'YColor',1-movdata.set.movbg,...
            %               'XTick',[],'YTick',[],'Box','on','Layer','Top',...
            %               'Parent',movdata.prev.hprev,'Color',movdata.set.movbg);
            %             image('Parent',movdata.prev.hprevpart{nprevparts},'CData',movdata.set.parts{nprevparts}.CData,'CDataMapping','scaled');
            %             if ndims(movdata.set.parts{nprevparts}.CData)==2;
            %               colormap(grey);
            %             end
            %             axis(movdata.prev.hprevpart{nprevparts},'ij');
            %             axis(movdata.prev.hprevpart{nprevparts},'tight');
            %           end
            %         end
            %         clear nprevparts
            
            % Create Parts Data
            for np = 1:length(movdata.set.parts)
                movdata.set.parts{np}.pos = [str2num(get(movdata.gui.hpartsim4{np},'String'))...
                    str2num(get(movdata.gui.hpartsim5{np},'String'))...
                    str2num(get(movdata.gui.hpartsim6{np},'String'))...
                    str2num(get(movdata.gui.hpartsim7{np},'String'))];
                if get(movdata.gui.hinput11,'Value') == 1;
                    part{np}.abs = movdata.set.parts{np}.pos.*[movdata.set.movsize(1) movdata.set.movsize(2) movdata.set.movsize(1) movdata.set.movsize(2)];
                else
                    part{np}.abs = movdata.set.parts{np}.pos;
                end
                % reset figure Position
                movdata.set.parts{np}.origpos = get(movdata.set.parts{np}.handle,'Position');
                prevFigPos = get(movdata.prev.hprev,'Position');
                movdata.set.parts{np}.origpos(1:2) = [prevFigPos(1)+prevFigPos(3)*movdata.set.parts{np}.pos(1) prevFigPos(2)-prevFigPos(4)-50++prevFigPos(4)*movdata.set.parts{np}.pos(2)];
                movdata.set.parts{np}.origpos(3:4) = part{np}.abs(3:4);
                % Font scaling as set in the settings of the OS has to be
                % taken into account, as the getframe commad returns the
                % scaled image (e.g. 1.25x larger than screen pixels dimensions), while get(h, 'Position') returns the image
                % dimensions as if there was no scaling.
                set(movdata.set.parts{np}.handle,'Units','pixels','Position',screenSizeScaling*movdata.set.parts{np}.origpos);
                figure(movdata.set.parts{np}.handle);
                %Get CData
                if movdata.set.parts{np}.status == true;
                    % two figure capture options are implemented
                    % fast: directly uses 'hardcopy' and presets the figure
                    % properies before scanning
                    % slow: uses export_fig
                    if movdata.fastCaptureMode == true;
                        movdata.set.parts{np} = vidGenerator_myHardcopyFigureSet(movdata.set.parts{np}, res);
                        movdata.set.parts{np}.CData = vidGenerator_myHardcopy(movdata.set.parts{np},res);
                    else
                        movdata.set.parts{np}.CData = export_fig(movdata.set.parts{np}.handle,'-nocrop','-q100','-RGB');%,'-painters'
                        if ndims(movdata.set.parts{np}.CData)==2;
                            movdata.set.parts{np}.CData = repmat(movdata.set.parts{np}.CData,[1 1 3]);
                        end
                    end
                end
            end
            
            % Create frame matrix
            movdata.prev.CData = uint8(repmat(reshape(255*movdata.set.movbg,[1 1 3]),[movdata.set.movsize(2) movdata.set.movsize(1)]));
            % fill in titleframe
            if movdata.set.titlestatus == true;
                movdata.set.title = get(movdata.gui.hparts2,'String');
                movdata.set.titlepos = [str2num(get(movdata.gui.hparts4,'String')) str2num(get(movdata.gui.hparts5,'String'))];
                movdata.set.titlesize = str2num(get(movdata.gui.hparts7,'String'));
                movdata.titleframe.posy(1) = max([1,round(movdata.set.movsize(2)-(movdata.set.titlepos(2)*movdata.set.movsize(2)+size(movdata.titleframe.CData,1)))]);
                movdata.titleframe.posy(2) = round(movdata.titleframe.posy(1)+size(movdata.titleframe.CData,1)-1);
                movdata.prev.CData(movdata.titleframe.posy(1):movdata.titleframe.posy(2),:,:) = movdata.titleframe.CData;
            end
            % fill in movie parts
            for np = 1:length(movdata.set.parts)
                if movdata.set.parts{np}.status == true;
                    movdata.set.parts{np}.posx(1) = round(movdata.set.parts{np}.pos(1)*movdata.set.movsize(1)+1);
                    movdata.set.parts{np}.posx(2) = round(movdata.set.parts{np}.posx(1) + size(movdata.set.parts{np}.CData,2)-1);
                    movdata.set.parts{np}.posy(1) = round(movdata.set.movsize(2)-(movdata.set.parts{np}.pos(2)*movdata.set.movsize(2)+size(movdata.set.parts{np}.CData,1))+1);
                    movdata.set.parts{np}.posy(2) = round(movdata.set.parts{np}.posy(1) + size(movdata.set.parts{np}.CData,1)-1);
                    movdata.prev.CData(movdata.set.parts{np}.posy(1):movdata.set.parts{np}.posy(2),movdata.set.parts{np}.posx(1):movdata.set.parts{np}.posx(2),:) = movdata.set.parts{np}.CData;
                end
            end
            cla(movdata.prev.hax);
            movdata.prev.himage = image(movdata.prev.CData,'Parent',movdata.prev.hax);
            axis(movdata.prev.hax,'off');
            figure(movdata.gui.hvidgen);
			if debugMatVis, debugMatVisFcn(2); end
        end
        
        %% Enable/Disable MovieParts
        function vidGenerator_enable(varargin)
          if debugMatVis, debugMatVisFcn(1); end
            nparts=length(movdata.parts);
            movdata.set.titlestatus = get(movdata.gui.hparts1,'Value');
            if get(movdata.gui.hparts1,'Value')==1;
                set(movdata.gui.hparts2,'Enable','on');
                set(movdata.gui.hparts3,'Enable','on');
                set(movdata.gui.hparts4,'Enable','on');
                set(movdata.gui.hparts5,'Enable','on');
                set(movdata.gui.hparts6,'Enable','on');
                set(movdata.gui.hparts7,'Enable','on');
            else
                set(movdata.gui.hparts2,'Enable','off');
                set(movdata.gui.hparts3,'Enable','off');
                set(movdata.gui.hparts4,'Enable','off');
                set(movdata.gui.hparts5,'Enable','off');
                set(movdata.gui.hparts6,'Enable','off');
                set(movdata.gui.hparts7,'Enable','off');
            end
            for nen=1:nparts;
                movdata.set.parts{nen}.status = get(movdata.gui.hpartsim1{nen},'Value');
                if get(movdata.gui.hpartsim1{nen},'Value')==1;
                    set(movdata.gui.hpartsim2{nen},'Enable','on');
                    set(movdata.gui.hpartsim3{nen},'Enable','on');
                    set(movdata.gui.hpartsim4{nen},'Enable','on');
                    set(movdata.gui.hpartsim5{nen},'Enable','on');
                    set(movdata.gui.hpartsim6{nen},'Enable','on');
                    set(movdata.gui.hpartsim7{nen},'Enable','on');
                else
                    set(movdata.gui.hpartsim2{nen},'Enable','off');
                    set(movdata.gui.hpartsim3{nen},'Enable','off');
                    set(movdata.gui.hpartsim4{nen},'Enable','off');
                    set(movdata.gui.hpartsim5{nen},'Enable','off');
                    set(movdata.gui.hpartsim6{nen},'Enable','off');
                    set(movdata.gui.hpartsim7{nen},'Enable','off');
                end
            end
            vidGenerator_preview
            clear nen;
            figure(movdata.gui.hvidgen);
			if debugMatVis, debugMatVisFcn(2); end
        end
        
        %% Choose BackgroundColor
        function vidGenerator_bgcolor(varargin)
          if debugMatVis, debugMatVisFcn(1); end
            switch varargin{3}
                case 'figure'
                    movdata.set.movbg = uisetcolor(get(movdata.gui.hinput6,'BackgroundColor'),'Choose Color for Movie Background');
                    set(movdata.gui.hinput6,'BackgroundColor',movdata.set.movbg);
                case 'axes'
                    movdata.set.movPartsbg = uisetcolor(get(movdata.gui.hinput13,'BackgroundColor'),'Choose Color for Movie Parts');
                    set(movdata.gui.hinput13,'BackgroundColor',movdata.set.movPartsbg);
                    for iii=1:numel(movdata.set.parts)
                        set(movdata.set.parts{iii}.handle,'Color',movdata.set.movPartsbg);
                    end
                case 'title'
                    movdata.set.movTitleColor = uisetcolor(get(movdata.gui.hinput15,'BackgroundColor'),'Choose Color for Movie Title');
                    set(movdata.gui.hinput15,'BackgroundColor',movdata.set.movTitleColor);
                case 'label'
                    movdata.set.movLabelColor = uisetcolor(get(movdata.gui.hinput17,'BackgroundColor'),'Choose Color for Axes Labels');
                    set(movdata.gui.hinput17,'BackgroundColor',movdata.set.movLabelColor);
                    for iii=1:numel(movdata.set.parts)
                        al = findobj(movdata.set.parts{iii}.handle ,'Type','Axes');
                        for jjj=1:numel(al)
                            set(al(jjj),'XColor',movdata.set.movLabelColor,'YColor',movdata.set.movLabelColor);
                            set(get(al(jjj),'XLabel'),'Color',movdata.set.movLabelColor);
                            set(get(al(jjj),'YLabel'),'Color',movdata.set.movLabelColor);
                        end
                    end
            end
            vidGenerator_preview
            figure(movdata.gui.hvidgen);
            if debugMatVis, debugMatVisFcn(2); end
        end
        
        %% Choose File to write
        function vidGenerator_savefile(varargin)
          if debugMatVis, debugMatVisFcn(1); end
            [movname,movdir,movtype]=uiputfile(movdata.FilterSpec,'Choose File for Movie Output',get(movdata.gui.houtput1,'String'));
            movdata.set.movname = movname(1:end-4);
            movdata.set.movdir = movdir(1:end-1);
            movdata.set.fileext = movdata.FilterSpec{movtype,1}(2:end);
            movdata.set.movstring = [movdir movname];
            set(movdata.gui.houtput1,'String',movdata.set.movstring);
            scalestring = get(movdata.gui.hinput9,'String');
            scalestring = scalestring{get(movdata.gui.hinput9,'Value')};
            set(movdata.prev.hprev,'name',[movdata.set.movname ' - (' num2str(movdata.set.movsize(1)) 'x' num2str(movdata.set.movsize(2)) ') - PREVIEW (' scalestring ')']);
            
            movdata.set.movcodec = movdata.FilterSpec{movtype,2}(1:end-8);
            codecpopup = get(movdata.gui.houtput4,'String');
            for iii = 1:length(codecpopup);
                if strcmp(codecpopup{iii},movdata.set.movcodec) || strcmp(codecpopup{iii},['FFMPEG > ' movdata.set.movcodec]);
                    setpopupvalue = iii;
                    break
                end
            end
            set(movdata.gui.houtput4,'Value',setpopupvalue);
            switch movdata.set.movcodec
                case 'Motion JPEG 2000'
                    set(movdata.gui.houtput7,'String','Compress:','Enable','on');
                    set(movdata.gui.houtput8,'String',num2str(10),'Enable','on');
                case 'Archival'
                    set(movdata.gui.houtput7,'String','Compress:','Enable','off');
                    set(movdata.gui.houtput8,'String',num2str(10),'Enable','off');
                case 'Uncompressed AVI'
                    set(movdata.gui.houtput7,'String','Quality [%]:','Enable','off');
                    set(movdata.gui.houtput8,'String',num2str(100),'Enable','off');
                otherwise
                    set(movdata.gui.houtput7,'String','Quality [%]:','Enable','on');
                    set(movdata.gui.houtput8,'String',num2str(75),'Enable','on');
            end
            clear codecpopup setpopupvalue movname movdir movtype;
            figure(movdata.gui.hvidgen);
			if debugMatVis, debugMatVisFcn(2); end
        end
        
        %% Update File String
        function vidGenerator_updatefile(varargin)
          if debugMatVis, debugMatVisFcn(1); end
            movdata.set.movstring = get(movdata.gui.houtput1,'String');
            dirpos = regexp(movdata.set.movstring,'\');
            pointpos = regexp(movdata.set.movstring, '\.');
            movdata.set.movname = movdata.set.movstring(dirpos(end)+1:pointpos(end)-1);
            scalestring = get(movdata.gui.hinput9,'String');
            scalestring = scalestring{get(movdata.gui.hinput9,'Value')};
            set(movdata.prev.hprev,'name',[movdata.set.movname ' - (' num2str(movdata.set.movsize(1)) 'x' num2str(movdata.set.movsize(2)) ') - PREVIEW (' scalestring ')']);
            if isempty(dirpos); dirpos = 0; end
            if ~isdir(movdata.set.movstring(1:dirpos(end)));
                movdata.set.movdir = pwd;
                movdata.set.movstring = [pwd '\' movdata.set.movstring(dirpos(end)+1:end)];
                set(movdata.gui.houtput1,'String',movdata.set.movstring);
            end
            pointpos = regexp(movdata.set.movstring, '\.');
            if isempty(pointpos); pointpos = length(movdata.set.movstring); end
            if length(movdata.set.movstring)-3~=pointpos(end);
                movdata.set.fileext = '.avi';
                movdata.set.movstring = [movdata.set.movstring(1:pointpos(end)-1) '.avi'];
                set(movdata.gui.houtput1,'String',movdata.set.movstring);
                pointpos = regexp(movdata.set.movstring, '\.');
            end
            movdata.set.fileext = movdata.set.movstring(pointpos:end);
            movdata.set.movcodec = [];
            for iii = 1:length(movdata.FilterSpec);
                if regexp(movdata.FilterSpec{iii,1},movdata.set.fileext)==2;
                    movdata.set.movcodec = movdata.FilterSpec{iii,2}(1:end-8);
                    break
                end
            end
            if isempty(movdata.set.movcodec);
                movdata.set.movcodec = 'Uncompressed AVI';
                movdata.set.fileext = '.avi';
                movdata.set.movstring = strrep(movdata.set.movstring,movdata.set.fileext,'.avi');
                set(movdata.gui.houtput1,'String',movdata.set.movstring);
            end
            codecpopup = get(movdata.gui.houtput4,'String');
            for iii = 1:length(codecpopup);
                if strcmp(codecpopup{iii},movdata.set.movcodec) || strcmp(codecpopup{iii},['FFMPEG > ' movdata.set.movcodec]);
                    setpopupvalue = iii;
                    break
                end
            end
            set(movdata.gui.houtput4,'Value',setpopupvalue);
            switch movdata.set.movcodec
                case 'Motion JPEG 2000'
                    set(movdata.gui.houtput7,'String','Compress:','Enable','on');
                    set(movdata.gui.houtput8,'String',num2str(10),'Enable','on');
                case 'Archival'
                    set(movdata.gui.houtput7,'String','Compress:','Enable','off');
                    set(movdata.gui.houtput8,'String',num2str(10),'Enable','off');
                case 'Uncompressed AVI'
                    set(movdata.gui.houtput7,'String','Quality [%]:','Enable','off');
                    set(movdata.gui.houtput8,'String',num2str(100),'Enable','off');
                otherwise
                    set(movdata.gui.houtput7,'String','Quality [%]:','Enable','on');
                    set(movdata.gui.houtput8,'String',num2str(75),'Enable','on');
            end
            figure(movdata.gui.hvidgen);
			if debugMatVis, debugMatVisFcn(2); end
        end
        
        %% Choose Codec
        function vidGenerator_choosecodec(varargin)
          if debugMatVis, debugMatVisFcn(1); end
            popupvalue = get(movdata.gui.houtput4,'Value');
            codecpopup = get(movdata.gui.houtput4,'String');
            movdata.set.movcodec = codecpopup{popupvalue};
            for iii = 1:length(movdata.FilterSpec);
                if regexp(movdata.FilterSpec{iii,2},strrep(movdata.set.movcodec,'FFMPEG > ',''))==1;
                    movdata.set.fileext = movdata.FilterSpec{iii,1}(2:end);
                end
            end
            movdata.set.movstring = get(movdata.gui.houtput1,'String');
            pointpos = regexp(movdata.set.movstring, '\.');
            if length(movdata.set.movstring)-3==pointpos;
                movdata.set.movstring = strrep(movdata.set.movstring,movdata.set.movstring(pointpos:end),movdata.set.fileext);
            else
                movdata.set.movstring = [movdata.set.movstring movdata.set.fileext];
            end
            set(movdata.gui.houtput1,'String',movdata.set.movstring);
            switch movdata.set.movcodec
                case 'Motion JPEG 2000'
                    set(movdata.gui.houtput7,'String','Compress:','Enable','on');
                    set(movdata.gui.houtput8,'String',num2str(10),'Enable','on');
                case 'Archival'
                    set(movdata.gui.houtput7,'String','Compress:','Enable','off');
                    set(movdata.gui.houtput8,'String',num2str(10),'Enable','off');
                case 'Uncompressed AVI'
                    set(movdata.gui.houtput7,'String','Quality [%]:','Enable','off');
                    set(movdata.gui.houtput8,'String',num2str(100),'Enable','off');
                otherwise
                    set(movdata.gui.houtput7,'String','Quality [%]:','Enable','on');
                    set(movdata.gui.houtput8,'String',num2str(75),'Enable','on');
            end
            clear popupvalue codecpopup setcodec movstring pointpos
            figure(movdata.gui.hvidgen);
			if debugMatVis, debugMatVisFcn(2); end
        end
        
        %% Set FrameRate/CompressionRate/Quality
        function vidGenerator_FCQ(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            movdata.set.movFrameRate = str2num(get(movdata.gui.houtput6,'String'));
            switch movdata.set.movcodec
                case 'Motion JPEG 2000'
                    movdata.set.movcompress = str2num(get(movdata.gui.houtput8,'String'));
                case 'Motion JPEG AVI'
                    movdata.set.movquality = str2num(get(movdata.gui.houtput8,'String'));
                case 'MPEG-4'
            end
			if debugMatVis, debugMatVisFcn(2); end
        end
        
        %% Activate Record
        function vidGenerator_record(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            movdata.rec = get(movdata.gui.hvidgen1,'Value');
            if movdata.rec
                vidWriter('initiate');
                set(movdata.gui.hvidgen1,'CData',icon_Stop);
                set(btPlayAll , 'CData',0.9*icon_RecordAll);
                set(btPlayZoom, 'CData',0.9*icon_RecordZoom);
            else
                vidWriter('terminate');
                set(movdata.gui.hvidgen1,'CData',icon_Record);
                set(btPlayAll , 'CData',arrowAll);
                set(btPlayZoom , 'CData',arrowZoom);
            end
			if debugMatVis, debugMatVisFcn(2); end
        end
        
        %% Close VidGenerator
        function vidGenerator_showHide(varargin)
            if debugMatVis, debugMatVisFcn(1); end
            switch varargin{3}
                case 'hide'
                    set([movdata.gui.hvidgen movdata.prev.hprev], 'Visible','off');
                    set(tbRecord, 'Value',0);
                    if withAlpha
                        c = [0 0 0];
                    else
                        c = get(0,'DefaultFigureColor');
                    end
                    for iii=1:length(movdata.parts);
                        set(movdata.set.parts{iii}.handle, 'Color',c);
                        al = findobj(movdata.set.parts{iii}.handle ,'Type','Axes');
                        for jjj=1:numel(al)
                            set(al(jjj),'XColor',movdata.set.movLabelColor,'YColor',movdata.set.movLabelColor);
                            set(get(al(jjj),'XLabel'),'Color','k');
                            set(get(al(jjj),'YLabel'),'Color','k');
                        end
                    end
                case 'show'
                    set([movdata.gui.hvidgen movdata.prev.hprev], 'Visible','on');
                    set(tbRecord, 'Value',1);
                    for iii=1:length(movdata.parts);
                        set(movdata.set.parts{iii}.handle, 'Color',movdata.set.movPartsbg);
                        al = findobj(movdata.set.parts{iii}.handle ,'Type','Axes');
                        for jjj=1:numel(al)
                            set(al(jjj),'XColor',movdata.set.movLabelColor,'YColor',movdata.set.movLabelColor);
                            set(get(al(jjj),'XLabel'),'Color',movdata.set.movLabelColor);
                            set(get(al(jjj),'YLabel'),'Color',movdata.set.movLabelColor);
                        end
                    end
            end
			if debugMatVis, debugMatVisFcn(2); end
        end
    end

    function vidWriter(writetype)
        if debugMatVis, debugMatVisFcn(1); end
        res = 3;
        switch writetype
            case 'initiate'
                % disable options dialog
                movdata.h = findobj(movdata.gui.hvidgen,'Enable','on');
                set(movdata.h,'Enable','off');
                set([movdata.gui.hvidgen1 movdata.gui.hvidgen3],'Enable','on');
                set(movdata.gui.hvidgen3,'Position',[35 8 500 27], 'String', {'Use matVis controls to scroll through dimensions while recording.';'Each change of a position selection will produce a new movie frame. Click to end recording.'});
                movprofile = movdata.set.movcodec;
                switch movdata.set.movcodec
                    case 'FFMPEG > MPEG-4'
                        movprofile = 'Uncompressed AVI';
                        movdata.set.movstring = strrep(movdata.set.movstring,'.mp4','.avi');
                    case 'FFMPEG > Theora'
                        movprofile = 'Uncompressed AVI';
                        movdata.set.movstring = strrep(movdata.set.movstring,'.ogg','.avi');
                end
%                 profile clear;
%                 profile on;
                prevunits = get(movdata.gui.hinput11,'String');
                prevunits = prevunits{movdata.set.partunits};
                try
                    movdata.mov = VideoWriter(movdata.set.movstring,movprofile); % initiate movie object
                catch erms
                    warndlg(['Moviefile could not be generated: ]' erms.message '. Select different file/folder name.']);
                    if debugMatVis, debugMatVisFcn(2); end
                    return
                end
                movdata.mov.FrameRate = movdata.set.movFrameRate; % set framerate
                switch movprofile
                    case 'Motion JPEG AVI'
                        movdata.mov.Quality = movdata.set.movquality; % set quality
                    case 'MPEG-4'
                        movdata.mov.Quality = movdata.set.movquality; % set quality
                    case 'Motion JPEG 2000'
                        movdata.mov.CompressionRatio = movdata.set.movcompress; % set compression
                end
                
                % Create Title Data
                if movdata.set.titlestatus == true;
                    movdata.titleframe.handle = figure('Position',[1 1 movdata.set.movsize(1) movdata.set.titlesize*1.5],...
                        'Color',movdata.set.movbg,'MenuBar','none','Units','pixels','NumberTitle','off','Color',movdata.set.movbg,...
                        'Resize','off','Name','Frame for Title','tag','movdata.frame','Visible','off');
                    movdata.titleframe.titleaxes = axes('Parent',movdata.titleframe.handle,'Position',[0 0 1 1],'Color',movdata.set.movbg);
                    movdata.titleframe.titletext = text('Parent',movdata.titleframe.titleaxes,'Position',[movdata.set.titlepos(1) .5],'Color',movdata.set.movTitleColor,'HorizontalAlignment','center','String',movdata.set.title,'FontSize',movdata.set.titlesize,'FontWeight','bold');
                    axis off;
                    drawnow;
                    if movdata.fastCaptureMode == true;
                        movdata.titleframe = vidGenerator_myHardcopyFigureSet(movdata.titleframe, res);
                        movdata.titleframe.CData = vidGenerator_myHardcopy(movdata.titleframe, res);
                        vidGenerator_myHardcopyFigureReset(movdata.titleframe)
                    else
                        movdata.titleframe.CData = export_fig(movdata.titleframe.handle,'-nocrop','-q100','-RGB');%,'-painters'
                        if ndims(movdata.titleframe.CData)==2;
                            movdata.titleframe.CData = repmat(movdata.titleframe.CData,[1 1 3]);
                        end
                    end
                    close (movdata.titleframe.handle);
                end
                
                % Create Parts Data
                for np = 1:length(movdata.set.parts)
                    %Get CData
                    if movdata.set.parts{np}.status == true;
                        % two figure capture options are implemented
                        % fast: directly uses 'hardcopy' and presets the figure
                        % properies before scanning
                        % slow: uses export_fig
                        if movdata.fastCaptureMode == true;
                            movdata.set.parts{np} = vidGenerator_myHardcopyFigureSet(movdata.set.parts{np}, res);
                            movdata.set.parts{np}.CData = vidGenerator_myHardcopy(movdata.set.parts{np},res);
                        else
                            movdata.set.parts{np}.CData = export_fig(movdata.set.parts{np}.handle,'-nocrop','-q100','-RGB');%,'-painters'
                            if ndims(movdata.set.parts{np}.CData)==2;
                                movdata.set.parts{np}.CData = repmat(movdata.set.parts{np}.CData,[1 1 3]);
                            end
                        end
                    end
                end
                
                % Create frame matrix
                movdata.frame.CData = uint8(repmat(reshape(255*movdata.set.movbg,[1 1 3]),[movdata.set.movsize(2) movdata.set.movsize(1)]));
                % fill in titleframe
                if movdata.set.titlestatus == true;
                    movdata.titleframe.posy(1) = round(movdata.set.movsize(2)-(movdata.set.titlepos(2)*movdata.set.movsize(2)+size(movdata.titleframe.CData,1)));
                    movdata.titleframe.posy(2) = round(movdata.titleframe.posy(1)+size(movdata.titleframe.CData,1)-1);
                    movdata.frame.CData(movdata.titleframe.posy(1):movdata.titleframe.posy(2),:,:) = movdata.titleframe.CData;
                end
                % fill in movie parts
                for np = 1:length(movdata.set.parts)
                    if movdata.set.parts{np}.status == true;
                        movdata.set.parts{np}.posx(1) = round(movdata.set.parts{np}.pos(1)*movdata.set.movsize(1)+1);
                        movdata.set.parts{np}.posx(2) = round(movdata.set.parts{np}.posx(1) + size(movdata.set.parts{np}.CData,2)-1);
                        movdata.set.parts{np}.posy(1) = round(movdata.set.movsize(2)-(movdata.set.parts{np}.pos(2)*movdata.set.movsize(2)+size(movdata.set.parts{np}.CData,1))+1);
                        movdata.set.parts{np}.posy(2) = round(movdata.set.parts{np}.posy(1) + size(movdata.set.parts{np}.CData,1)-1);
                        movdata.frame.CData(movdata.set.parts{np}.posy(1):movdata.set.parts{np}.posy(2),movdata.set.parts{np}.posx(1):movdata.set.parts{np}.posx(2),:) = movdata.set.parts{np}.CData;
                    end
                end
                
                % open VideoWriter Object
                open(movdata.mov);
                % write first frame to Object
                writeVideo(movdata.mov,movdata.frame.CData);
                
            case 'append'
                % update frame parts from updated parts
                for np = 1:length(movdata.set.parts)
                    if movdata.set.parts{np}.status == true;
                        if movdata.fastCaptureMode == true;
                            movdata.set.parts{np}.CData = vidGenerator_myHardcopy(movdata.set.parts{np},res);
                        else
                            movdata.set.parts{np}.CData = export_fig(movdata.set.parts{np}.handle,'-nocrop','-q100','-RGB');%,'-painters'
                            if ndims(movdata.set.parts{np}.CData)==2;
                                movdata.set.parts{np}.CData = repmat(movdata.set.parts{np}.CData,[1 1 3]);
                            end
                        end
                        % update frame Matrix
                        movdata.frame.CData(movdata.set.parts{np}.posy(1):movdata.set.parts{np}.posy(2),movdata.set.parts{np}.posx(1):movdata.set.parts{np}.posx(2),:) = movdata.set.parts{np}.CData;
                    end
                end
                writeVideo(movdata.mov,movdata.frame.CData);
                
            case 'terminate'
                set(movdata.gui.hvidgen3, 'String', 'Saving movie file. Please wait ...');
                drawnow;
                close(movdata.mov);
                switch movdata.set.movcodec
                    case 'FFMPEG > MPEG-4'
                        [status,result] = system(['ffmpeg -i ' movdata.set.movstring ' -y ' strrep(movdata.set.movstring,'.avi','.mp4')]);
                        delete(movdata.set.movstring);
                    case 'FFMPEG > Theora'
                        [status,result] = system(['ffmpeg -i ' movdata.set.movstring ' -vcodec libtheora -y ' strrep(movdata.set.movstring,'.avi','.ogg')]);
                        delete(movdata.set.movstring);
                end
                if movdata.fastCaptureMode == true;
                    for np = 1:length(movdata.set.parts)
                        vidGenerator_myHardcopyFigureReset(movdata.set.parts{np})
                    end
                end
                % disable options dialog
                set(movdata.h,'Enable','on');
                set(movdata.gui.hvidgen3, 'Position',[35 12 500 15],'String', 'Adjust settings, then press button to activate recording.');
%                 profile off;
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

    function hh = vidGenerator_myHardcopyFigureSet(hh,rres)
        if debugMatVis, debugMatVisFcn(1); end
        % MATLAB "feature": axes limits can change when printing
        hh.Hlims = findall(hh.handle, 'Type', 'axes');
        hh.Xlims = make_cell(get(hh.Hlims, 'XLimMode'));
        hh.Ylims = make_cell(get(hh.Hlims, 'YLimMode'));
        hh.Zlims = make_cell(get(hh.Hlims, 'ZLimMode'));
        hh.XTicks = make_cell(get(hh.Hlims, 'XTickMode'));
        hh.YTicks = make_cell(get(hh.Hlims, 'YTickMode'));
        hh.ZTicks = make_cell(get(hh.Hlims, 'ZTickMode'));
        
        % Set all axes limit modes to manual, so the limits can't change
        set(hh.Hlims, 'XLimMode', 'manual', 'YLimMode', 'manual', 'ZLimMode', 'manual');
        set(hh.Hlims, 'XTickMode', 'manual', 'YTickMode', 'manual', 'ZTickMode', 'manual');
        
        % Set to print exactly what is there
        hh.old_IH_mode = get(hh.handle, 'InvertHardcopy');
        set(hh.handle, 'InvertHardcopy', 'off');
        
        % get figure size
        old_U_mode = get(hh.handle, 'Units');
        set(hh.handle, 'Units', 'pixels');
        px = get(hh.handle, 'Position');
        hh.px = [px([4 3])*rres 3];
        
        set(hh.handle, 'Units', old_U_mode);
        % set PaperPositionMode
        hh.old_PPM_mode = get(hh.handle, 'PaperPositionMode');
        set(hh.handle, 'PaperPositionMode', 'auto');
        % Helper function
        function A = make_cell(A)
            if ~iscell(A)
                A = {A};
            end
        end
        if debugMatVis, debugMatVisFcn(2); end
    end

    % extended hardcopy function
    function A = vidGenerator_myHardcopy(hh,rres)
        if debugMatVis, debugMatVisFcn(1); end
        % '-painters' is the preferd renderer but can only capture visible figures (figures do not have to be updated by drawnow!!!)
        % as workaround we use the slower '-opengl' renderer for invisible figures
        if strcmp(get(hh.handle,'Visible'),'on'); % this is not time consuming (about 10^-5s)
            A = getframe(hh.handle);%hardcopy(hh.handle,'temp.tif','-dpainters');%,-dopengl',sprintf('-r%d',96*rres)'
            A = A.cdata;
        else
            % figure capture
            A = getframe(hh.handle);%hardcopy(hh.handle,'temp.tif','-dopengl',sprintf('-r%d',96*rres));%,'-dpainters'
            A = A.cdata;
            % cuts dimension
            if ~isequal(size(A), hh.px)
                % Correct the output size
                %                 A = A(1:min(end,hh.px(1)),1:min(end,hh.px(2)),:);  %
                %                 Correction doesn't appear to (alway) work
                % Corrected version (Stephan Junek)
                if size(A,1) < hh.px(1)
                    A(end+1:round(hh.px(1)),:,:) = 0;
                elseif size(A,1) > hh.px(1)
                    A = A(1:round(hh.px(1)),:,:);
                end
                if size(A,2) < hh.px(2)
                    A(:,end+1:round(hh.px(2)),:) = 0;
                elseif size(A,2) > hh.px(2)
                    A = A(:,1:round(hh.px(2)),:);
                end
            end
            % does resizing (requires image processing toolbox :-(
            A = imresize(A, 1/rres, 'bilinear');
        end
        if debugMatVis, debugMatVisFcn(2); end
    end
%% resets figure properties
    function vidGenerator_myHardcopyFigureReset(hh)
        if debugMatVis, debugMatVisFcn(1); end
        % Reset the axes limit modes
        for a = 1:numel(hh.Hlims)
            set(hh.Hlims(a), 'XLimMode', hh.Xlims{a}, 'YLimMode', hh.Ylims{a}, 'ZLimMode', hh.Zlims{a});
            set(hh.Hlims(a), 'XTickMode', hh.XTicks{a}, 'YTickMode', hh.YTicks{a}, 'ZTickMode', hh.ZTicks{a});
        end
        set(hh.handle, 'PaperPositionMode', hh.old_PPM_mode);
        set(hh.handle, 'InvertHardcopy', hh.old_IH_mode);
        if debugMatVis, debugMatVisFcn(2); end
    end
    function y = nansum(x,dim)
        if debugMatVis>1, debugMatVisFcn(1); end
        x(isnan(x)) = 0;
        if nargin == 1
            y = sum(x); 
        else
            y = sum(x,dim);
        end
        if debugMatVis>1, debugMatVisFcn(2); end
    end
    function m = nanmean(x,dim)
        if debugMatVis>1, debugMatVisFcn(1); end
        % Find NaNs and set them to zero
        nans = isnan(x);
        x(nans) = 0;
        if nargin == 1 % let sum deal with figuring out which dimension to use
            % Count up non-NaNs.
            n = sum(~nans);
            n(n==0) = NaN; % prevent divideByZero warnings
            % Sum up non-NaNs, and divide by the number of non-NaNs.
            m = sum(x) ./ n;
        else
            % Count up non-NaNs.
            n = sum(~nans,dim);
            n(n==0) = NaN; % prevent divideByZero warnings
            % Sum up non-NaNs, and divide by the number of non-NaNs.
            m = sum(x,dim) ./ n;
        end
        if debugMatVis>1, debugMatVisFcn(2); end
    end
%% Close all windows
    function closeGui(varargin)
        if debugMatVis, debugMatVisFcn(1); end
        % Set flag to stop execution of calculation of global histogram
        shuttingDown = 1;
        %         try
        %             delete(timerHist);
        %         catch %#ok
        %         end
        isPlaying = 0;  %Stop play to avoid error message
        if with2DHist && ishandle(matVis2DHist.figHandles.gui)
            matVis2DHist.fctHandles.closeGui();
        end
        if isfield(movdata,'gui')
            try delete(movdata.gui.hvidgen);catch; end %#ok
            try delete(movdata.prev.hprev); catch; end %#ok
        end
        try
            delete(exportWin); catch    %#ok
        end
        try
            delete(profileWin); catch    %#ok
        end
        try
            delete(profileTraceWin); catch    %#ok
        end
        try
            delete(tempWin); catch    %#ok
        end
        try
            delete(imageWin); catch    %#ok
        end
        try
            delete(zoomWin); catch    %#ok
        end
        try
            delete(plotWin); catch    %#ok
        end
        try
            delete(roiWin); catch    %#ok
        end
        try
            delete(tifParFig); catch    %#ok
        end
        try
            delete(histWin); catch    %#ok
        end
        try
            delete(tempRoiWin); catch    %#ok
        end
        try
            delete(tempRoiPropWin); catch    %#ok
        end
        
        delete(gui);
        if debugMatVis, debugMatVisFcn(2); end
    end
    function  debugMatVisFcn(inOut)
        inOutStr = {sprintf('%sStart ',char(9484)),sprintf('%sEnd   ',char(9492))'};
        ST = dbstack;
        %ST.name
        fctLevel = length(ST)-1;
        fctName = strrep(strrep(ST(2).name,'matVis/',''),'/','_');          % name of executed function
        if inOut == 1
          if fctLevel == 1 && isstruct(fctCount)
            fnames = fieldnames(fctCount);
            for ii = 1:length(fnames)
              fctCount.(fnames{ii}).count = 0;
            end
          end
          if strcmp(fctName,'updateRoiSelection') && inOut == 1      % to stop at specific function name
            %stopHere
          end
          if isfield(fctCount,fctName)                                      % counts for executed function
            fctCount.(fctName).count = fctCount.(fctName).count+1;          % counts within callback
            fctCount.(fctName).countTot = fctCount.(fctName).countTot+1;    % total calls since matVis start
          else
            fctCount.(fctName).count = 1;
            fctCount.(fctName).countTot = 1;
            fctCount.(fctName).callFct = [];
            fctCount.(fctName).countCallFct = 0;
          end
          % stores calling subfunction name
          if fctLevel > 1
            callFctName = strrep(strrep(ST(3).name,'matVis/',''),'/','_');  
            if isfield(fctCount.(fctName).callFct,callFctName)
              fctCount.(fctName).callFct.(callFctName).count    = fctCount.(fctName).callFct.(callFctName).count + 1;
              %fctCount.(fctName).callFct.(callFctName).countTot = fctCount.(fctName).callFct.(callFctName).countTot + 1;
            else
              fctCount.(fctName).callFct.(callFctName).count = 1;
              %fctCount.(fctName).callFct.(callFctName).countTot = 1;
              fctCount.(fctName).countCallFct = fctCount.(fctName).countCallFct + 1;
            end
          end
          fctNameHLp = regexp(fctName,'_');
          if isempty(fctNameHLp), fctNameHL=fctName; else fctNameHL=fctName(fctNameHLp+1:end); end
          fctCount.(fctName).time  = tic;
          timeStr = [];
          fctNameStr = sprintf('<a href="matlab:matlab.desktop.editor.openAndGoToFunction(which(''matVis.m''),''%s'');">%s</a>',fctNameHL, fctName);%fctName
          if fctLevel > 1
            fctNameStr = sprintf('%s/<a href="matlab:matlab.desktop.editor.openAndGoToLine(which(''matVis.m''),%d);">line%d</a>',fctNameStr,ST(3).line,ST(3).line);
          end
        else
          timeStr = sprintf('(execution time: %.3f s)',toc(fctCount.(fctName).time));
          fctNameStr = sprintf('%s',fctName);%fctName
        end
        fprintf('%s%s %d/%d:%s %s\n',repmat(sprintf('|\t'), 1, fctLevel-1), inOutStr{inOut},fctCount.(fctName).count,fctCount.(fctName).countTot, fctNameStr, timeStr);
        if fctLevel == 1 && inOut == 2
          fprintf('\n')
          assignin('base', 'matVisFunctionCalls', fctCount);  % saves fctCount to matlab workspace when complete function call for later analysis
        end
    end
    if debugMatVis, debugMatVisFcn(2); end
end

%% To do:
% - alphaMap: Open subset in new matVis does not export ‚alphaMap‘
% - alphaMap: ROIdata are not updated while 'playing'
% - alphaMap: Plot window could show (weighted mean) of ratio data AND (scaled?) alpha data
% - Find max in zoom region?
% - Extend more functions for several data sets (global hist, ...)
% - Apply filter to alpha map? Checkbox or toggle button?
% - 2D histogram feature:
%    - limits from slider limits of contrast sliders? in this case, the
%    dimScale property has to be update-able
%    - extend from data/alpha-case to data1/data2 case. not clear
%      yet how best to implement
%    - enable global histogram display AND 2D histogram feature (i.e don't
%    use the same button depending on input)
% - "mask drawing" module? (like cleanStack-function)
% - refresh button? re-load variable of same name into matVis (without
% making any changes to the settings)
% - add new features to custom config: directions of x,y-axis, filter
% setting, linked window size/position


%% Known bugs:
% - Handling of figure icons not 'clean' (e.g. one of them might appear as
%   icon of Matlab figure container)
% - Produces complete Matlab crash when data with alpha maps are loaded and
%   a video is created with the video generator
% - Clicking twice in a row on toggle button of the colomap button
%   group will mess up the 'SelectedObject' property of the button group,
%   resulting in an error in the switch expression (updateColormap at 5549)

%% For further bug reports or feature requests send an e-mail to
%  Stephan Junek (sjunek@gwdg.de)