function [png] = PNGentropy(dat0, sampint, nsymb, name)

dat = dat0(:,1:sampint:end,:); % downsample
clear dat0
warning('off', 'MATLAB:DELETE:FileNotFound');

png = nan(size(dat,1),1);
n = size(dat,2)*size(dat,3);
num = log(n)/log(nsymb); % same normalization as LZ76 (entropy rate)
%num = log(n)/log(nsymb + 1); % ask Pedro if we should do it this way?
normfact = num/n;
dpi =  1000;
dpm = dpi * 39.3701; % Convert inches to meters

% Store filenames for cleanup
filenames = cell(size(dat,1), 1);

% remove byte size of header by saving a blank image
try 
    file_info = dir('./recyclemepls/PNG_REF_IMAGE.png');
    refval = file_info.bytes;
catch
    refimg = uint8(zeros(10,10));
    filename = './recyclemepls/PNG_REF_IMAGE.png';
    imwrite(refimg, filename, 'png', 'XResolution', dpm, 'YResolution', dpm, 'ResolutionUnit', 'meter');
    file_info = dir('./recyclemepls/PNG_REF_IMAGE.png');
    refval = file_info.bytes;
end

parfor ich = 1:size(dat,1)
    X = squeeze(dat(ich,:,:));
    tmp = discretize(X,quantile(X(:),0:1/nsymb:1));
    img = uint8((tmp - 1) / (nsymb - 1) * 255);

    % Generate unique filename for each iteration
    % @Dimitrios update the file path for your computer 
    filename = sprintf('./recyclemepls/%s_CH%i_deleteme.png', name, ich);
    filenames{ich} = filename;

    % Write with DPI metadata
    imwrite(img', filename, 'png', 'XResolution', dpm, 'YResolution', dpm, 'ResolutionUnit', 'meter');
    % Get file size
    file_info = dir(filenames{ich});
    png(ich) = (file_info.bytes - refval) * normfact;
    fprintf('.')
end

% Clean up files after parfor loop
recycle('off'); % Turn off recycle for cleanup
for i = 1:length(filenames)
    if exist(filenames{i}, 'file')
        delete(filenames{i});
    end
end
recycle('on'); % Turn recycle back on
fprintf('\n')
close all
end
