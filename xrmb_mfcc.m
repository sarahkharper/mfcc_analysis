function [outdata] = xrmb_mfcc(filename, window, overlap, cfnum, delt);
%Function designed to get mfccs from specific time points in an audio file.
%These time points presumably correspond to some meaningful articulatory or 
%acoustic landmark in a token of interest.
%
%INPUT:
%FILENAME = name of .csv file containing list of paired time points (time
%in ms) and .wav file names to use in MFCC calculation and extraction
%WINDOW = Length of window to be used in MFCC calculation (default to 25 ms)
%OVERLAP = Number of overlapping samples between adjacent windows (default to 5
%ms)
%CFNUM = # of mfccs to be calculated
%DELTA = Distance from window containing articulatory landmark at which to
%compare change

%Read data file containing timepoints at which mfccs should be extracted
data = readtable(filename);
cfstr = num2str(cfnum);
shift = num2str(window - overlap);

for i = 1:length(data.SUBJ),
        i
        %Get timepoints at which formants will be extracted
        maxC = data.MAXC_TIME_ms(i);

        %Find existing .wav file (and cd to directory) or (if none found)
        %create .wav file from .mat file
        fn = sprintf('%s_%s', char(data.SUBJ(i)), char(data.TASK(i)));
        wavname = strcat(fn, '_scaled.wav');
        if exist(wavname, 'file'),
            dir = which(wavname);
            dir = dir(1:end-numel(wavname));
            cd(dir);
        else
            matname = strcat(fn, '.mat');
            matdat = load(matname);
            fns = fieldnames(matdat);
            matdat = matdat.(fns{3});
            dir = which(matname);
            dir = dir(1:end-numel(matname));
            cd(dir);
            audio_raw = matdat(1).SIGNAL;
            sr = matdat(1).SRATE;
            %scale raw audio data from .mat file and write it to a .wav file
            audio_scale = 2*(audio_raw-min(audio_raw))/(max(audio_raw)-min(audio_raw))-1;
            audiowrite(wavname, audio_scale,sr);
        end

        %load wav file
        [wav,fs] = audioread(wavname);

        %convert window length in ms to samples
        wm = 1000/window; %calculate portion of a second represented by window ms
        window_samp = round(fs/wm);
        %convert overlap in ms to samples
        om = 1000/overlap; %calculate portion of a second represented by overlap ms
        overlap_samp = round(fs/om);

    %     %calculate number of windows that fit in audio file
    %     wintot = round((length(wav) - overlap_samp)/(window_samp - overlap_samp));

        %Find existing MFCC file (and cd to directory) or (if none found)
        %create MFCC file from .wav file
        fn = sprintf('%s_%s', char(data.SUBJ(i)), char(data.TASK(i)));
        mfccname = sprintf('%s_%smfcc_%sshift.csv', fn, cfstr, shift);
        cd('/Users/sarahharper/Dropbox/Research/Dissertation/MFCC Analysis');
        if exist(mfccname, 'file'),
            deltname = sprintf('%s_%smfcc_%sshift.csv', fn, 
            dir = which(mfccname);
            dir = dir(1:end-numel(mfccname));
            cd(dir);
            coeffs = csvread(mfccname);  %read in mfcc file as matrix
            delta = csvread(
        else
            % calculate mfccs over entire audio file (default # of mfccs = 13)
            %cd('/Users/sarahharper/Dropbox/Research/Dissertation/MFCC Analysis');
            [coeffs, delta, deltaDelta, ~] = mfcc(wav, fs, 'WindowLength', window_samp, 'OverlapLength', overlap_samp, 'NumCoeffs', cfnum);
            writematrix(coeffs, mfccname); %save mfccs as .csv file
            fprintf('created %s\n', mfccname);
        end

        %find number of the window that corresponds to the time of the articulatory
        %landmark
        audC = maxC; %+ 15;
        windata = round((audC - (window/2))/(window-overlap)) + 1;
        if windata > size(coeffs,1),
            windata = size(coeffs,1);
        end

        %extract mfcc information at window corresponding to time of
        %articulatory landmark and window at indicated distance(s)
        MFCC_maxc = coeffs(windata,:);
        MFCC_delta = delta(windata,:);
        MFCC_deltaDelta = deltaDelta(windata,:);
        for j = 1:size(coeffs,2),
            jstr = num2str(j-1);
            dlt = num2str(d);
            mxnm = sprintf('MFCC%s', jstr);
            dlnm = sprintf('MFCC%s_delta', jstr);
            ddnm = sprintf('MFCC%s_deltaDelta', jstr);
            %add mfcc @ MAXC to table
            outdata(i).(mxnm) = MFCC_maxc(j);
            %add delta
            outdata(i).(dlnm) = MFCC_delta(j);
            outdata(i).(ddnm) = MFCC_deltaDelta(j);
            for h = 1:length(delt),
                d = delt(h);
                MFCC_prev = coeffs(windata - d,:);
                vrnm = sprintf('MFCC%s_prev%s', jstr, dlt);
                %add mfcc @ distanced window
                outdata(i).(vrnm) = MFCC_prev(j);
        end
    end
end

%Save time point MFCCs to data file
out = struct2table(outdata);
tfinal = [data, out];
[~, name, ~] = fileparts(filename);
outname = sprintf("%s_%sMFCC_%sSHIFT_new.csv", name, cfstr, shift);
writetable(tfinal,outname);
end

