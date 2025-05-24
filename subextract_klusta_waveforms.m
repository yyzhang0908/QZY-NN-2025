function subextract_klusta_waveforms(ds, wp, trg, clus, nbefore, nafter, ygap)
wave_color = 'c';

tmin = 0;
tmax = max(ds.SpkTs);
if isnumeric(trg) && length(trg) == 2
    tmin = max(min(trg), tmin);
    tmax = min(max(trg), tmax);
end

if isempty(clus)
    clus = ds.ClusterCh(:, 1);
else
    clus(~ismember(clus, ds.ClusterCh(:, 1))) = [];
end

[fpath, fname] = fileparts(ds.LFPPath);
op = fullfile(wp, fname);
if ~isfolder(op); mkdir(op); end

wave_ts = -nbefore:1:nafter;
xlim = [min(wave_ts), max(wave_ts)];

wn = length(wave_ts);
nch = size(ds.ChGrp, 1);
nclu = length(clus);

[start_ts, ~, start_ts_valid, stop_ts_valid] = ...
    get_chunk_info(fullfile(fpath, 'process.prm'), tmax);

opt_str = sprintf('%.0fs_%.0fs_L%d_R%d', tmin, tmax, nbefore, nafter);

for cc = 1:nclu
    clu = clus(cc);
    spkt = ds.SpkTs(ds.SpkCluster == clu);
    spkt(spkt < min(trg) | spkt > max(trg)) = [];
    if ~any(spkt); continue; end
    spk_ts = round(spkt * ds.SpkFs);
    spk_chunk = arrayfun(@(x)(find(start_ts_valid <= x & stop_ts_valid >= x, 1)), spk_ts);

    chunk_ls = unique(spk_chunk);
    chunk_ls(chunk_ls == 0) = [];
    if isempty(chunk_ls); continue; end
    nchunk = length(chunk_ls);
    waves = cell(nchunk, 1);
    for kk = 1:nchunk
        da = readNPY(fullfile(fpath, '.spikedetekt\group_all\filtered', ...
            sprintf('chunk_%d.npy', start_ts_valid(chunk_ls(kk))-1)));
        ts = spk_ts(spk_chunk == chunk_ls(kk));
        nspk = numel(ts);
        wts = ts - start_ts(chunk_ls(kk)) + 2 + wave_ts;
        waves{kk} = reshape(da(wts(:), :), nspk, wn, nch);
    end
    waves = cat(1, waves{:});

    gstr = ds.ChGrp{cat(1, ds.ChGrp{:, 1}) == ds.ClusterCh(clu, 2), 2};
    gch = cat(1, ds.ChGrp{strcmp(ds.ChGrp(:, 2), gstr), 1});
    ngrp = length(gch);
    figure; set(gcf, 'Position', [100, 100, 300, 600]); hold on;
    for ii = 1:ngrp
        plot(xlim, [1, 1] * ygap * ii, 'Color', [.5, .5, .5], 'LineWidth', 1);
        plot(wave_ts, waves(:, :, gch(ii)) + ii * ygap, ...
            'Color', wave_color, 'LineWidth', 0.3);
        text(min(xlim) - range(xlim) * 0.1, ygap * ii, ...
            num2str(gch(ii)));
    end
    set(gca, 'Visible', 'off', 'XLim', xlim, 'YLim', [0, ngrp + 1] * ygap);
    text(mean(wave_ts), (ngrp + 1) * ygap, sprintf('%s\n%s\nCluster #%d', ...
        fname, opt_str, clu), 'Interpreter', 'none', 'FontSize', 10, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'baseline');
    fp = fullfile(op, sprintf('%s_cluster#%.3d', opt_str, clu));
    print(gcf, fp, '-djpeg');
    print(gcf, fp, '-dpdf', '-r300', '-vector');
    saveas(gcf, fp, 'fig');
end
end

function [ibeg, iend, ibego, iendo] = get_chunk_info(prmfile, tmax)
fid = fopen(prmfile, 'r');
prm_str = fread(fid, [1, inf], '*char');
fclose(fid);

fs = str2double(regexp(prm_str, ...
    '(?<=sample_rate = )\S+(?=,)', 'match', 'once'));
chunk_sec = str2double(regexp(prm_str, ...
    '(?<=chunk_size_seconds = )\S+(?=,)', 'match', 'once'));
overlap_sec = str2double(regexp(prm_str, ...
    '(?<=chunk_overlap_seconds = )\S+(?=,)', 'match', 'once'));
chunkn = round(chunk_sec * fs);
overlapn = round(overlap_sec * fs);

max_ts = round(tmax * fs);
nchunk = ceil(max_ts / (chunkn - overlapn));

[ibeg, iend, ibego, iendo] = deal(nan(nchunk, 1));
ibeg(1) = 1;
iend(1) = 1 + chunkn;
ibego(1) = ibeg(1);
iendo(1) = iend(1) - floor(overlapn / 2);
kk = 1;
while iend(kk) - overlapn + chunkn < max_ts
    kk = kk + 1;
    ibeg(kk) = iend(kk - 1) - overlapn;
    iend(kk) = ibeg(kk) + chunkn;
    ibego(kk) = iendo(kk - 1);
    iendo(kk) = iend(kk) - floor(overlapn / 2);
end
ibeg(kk + 1) = iend(kk) - overlapn;
iend(kk + 1) = max_ts;
ibego(kk + 1) = iendo(kk);
iendo(kk + 1) = iend(kk + 1);
end

function data = readNPY(filename)
% Function to read NPY files into matlab.
% *** Only reads a subset of all possible NPY files, specifically N-D arrays of certain data types.
% See https://github.com/kwikteam/npy-matlab/blob/master/tests/npy.ipynb for
% more.
%

[shape, dataType, fortranOrder, littleEndian, totalHeaderLength, ~] = readNPYheader(filename);

if littleEndian
    fid = fopen(filename, 'r', 'l');
else
    fid = fopen(filename, 'r', 'b');
end

try

    [~] = fread(fid, totalHeaderLength, 'uint8');

    % read the data
    data = fread(fid, prod(shape), [dataType '=>' dataType]);

    if length(shape)>1 && ~fortranOrder
        data = reshape(data, shape(end:-1:1));
        data = permute(data, length(shape):-1:1);
    elseif length(shape)>1
        data = reshape(data, shape);
    end

    fclose(fid);

catch me
    fclose(fid);
    rethrow(me);
end
end

function [arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename)
% function [arrayShape, dataType, fortranOrder, littleEndian, ...
%       totalHeaderLength, npyVersion] = readNPYheader(filename)
%
% parse the header of a .npy file and return all the info contained
% therein.
%
% Based on spec at http://docs.scipy.org/doc/numpy-dev/neps/npy-format.html

fid = fopen(filename);

% verify that the file exists
if (fid == -1)
    if ~isempty(dir(filename))
        error('Permission denied: %s', filename);
    else
        error('File not found: %s', filename);
    end
end

try
    
    dtypesMatlab = {'uint8','uint16','uint32','uint64','int8','int16','int32','int64','single','double', 'logical'};
    dtypesNPY = {'u1', 'u2', 'u4', 'u8', 'i1', 'i2', 'i4', 'i8', 'f4', 'f8', 'b1'};
    
    
    magicString = fread(fid, [1 6], 'uint8=>uint8');
    
    if ~all(magicString == [147,78,85,77,80,89])
        error('readNPY:NotNUMPYFile', 'Error: This file does not appear to be NUMPY format based on the header.');
    end
    
    majorVersion = fread(fid, [1 1], 'uint8=>uint8');
    minorVersion = fread(fid, [1 1], 'uint8=>uint8');
    
    npyVersion = [majorVersion minorVersion];
    
    headerLength = fread(fid, [1 1], 'uint16=>uint16');
    
    totalHeaderLength = 10+headerLength;
    
    arrayFormat = fread(fid, [1 headerLength], 'char=>char');
    
    % to interpret the array format info, we make some fairly strict
    % assumptions about its format...
    
    r = regexp(arrayFormat, '''descr''\s*:\s*''(.*?)''', 'tokens');
    dtNPY = r{1}{1};    
    
    littleEndian = ~strcmp(dtNPY(1), '>');
    
    dataType = dtypesMatlab{strcmp(dtNPY(2:3), dtypesNPY)};
        
    r = regexp(arrayFormat, '''fortran_order''\s*:\s*(\w+)', 'tokens');
    fortranOrder = strcmp(r{1}{1}, 'True');
    
    r = regexp(arrayFormat, '''shape''\s*:\s*\((.*?)\)', 'tokens');
    shapeStr = r{1}{1}; 
    arrayShape = str2num(shapeStr(shapeStr~='L'));

    
    fclose(fid);
    
catch me
    fclose(fid);
    rethrow(me);
end
end
