%% AnalyzeAudioScript.m
% By Dariusz Mrugala
% This script is designed to analyze WAV (or other uncompressed audio) for
% microphone positioning. The user selects two or more audio files to
% analyze. These audiofiles are recorded simultaneously and contain an
% identical number of samples. By convention, microphone 1 or audiofile 1
% should correspond to a L hand side microphone (either at collarbone or
% ear). Microphone 2 (audiofile 2) should correspond to a R hand side
% microphone (also at either collarbone or ear). Finally microphone 3
% corresponds to a microphone at the wrist.

%phase replacement is optimized for Fs = 16000 and N = 512
%% Cleanup 
clc, clear% -except x_vertex_orig
%% Global Variables - you can modify these
%note number of files selected (uigetfile; line 62) must match number of microphones in
%microphone_locs or this code will crash
%stmpath = "C:\Users\Dariusz_admin\Documents\Podcasting192\Audio Files\"; %change this to your audio file path for convenience
stmpath = "C:\Users\dariusz.DESKTOP-JP65O25\Dropbox\PC (2)\Documents\focus rite audio files 11-22-22";
Fs_target = 48000; %downsample to this frequency, must be factor of original sampling rate
l_buff = 512; %l_buff is N or buffer/FFT size. Window time is l_buff/Fs_target ie 2048/48000 = 0.032 seconds (per slice)
Enable_welch = false;
HPF_cutoff = 200; %currently unused
Voice_activity_detection_freq = 187.5; %Hz; sets threshold for spectral entropy
Noise_cancellation = false; %enable for noise cancellation code, not implemented
Noise_cancel_ch = [3 4]; % only used if Noise_cancellation is true
%speaker_loc = [0 -0.005 0.006]; % speaker coordinate in meter (old trials)
speaker_loc = [0.13 0.04 -0.21]; % ie relative coordinates of speakers mouth
distance_interval = .3; %meters
max_distance = 3; %meters
speech_window_time = .5; %seconds
%compute additional sound sources based on speaker location
%ie to a precision to a tenth of a meter, create an acoustic map, should be
% length x length by number of 
non_speaker_locs_x = [speaker_loc(1)-max_distance:distance_interval:speaker_loc(1)-distance_interval speaker_loc(1)+distance_interval:distance_interval:speaker_loc(1)+max_distance];
non_speaker_locs_y = [speaker_loc(2)-max_distance:distance_interval:speaker_loc(2)-distance_interval speaker_loc(2)+distance_interval:distance_interval:speaker_loc(2)+max_distance];
non_speaker_locs_z = [speaker_loc(3)-max_distance:distance_interval:speaker_loc(3)-distance_interval speaker_loc(3)+distance_interval:distance_interval:speaker_loc(3)+max_distance];
[locs_x, locs_y, locs_z] =  meshgrid(non_speaker_locs_x, non_speaker_locs_y, non_speaker_locs_z);
locs_x_reshape = reshape(locs_x, numel(locs_x),[]);
locs_y_reshape = reshape(locs_y, numel(locs_y),[]);
locs_z_reshape = reshape(locs_z, numel(locs_z),[]);
non_speaker_locs_total = [locs_x_reshape, locs_y_reshape, locs_z_reshape];
% initial trial
%microphones_locs = [.161 .471 .230; .357 .448 .227; .503 .448 .231; .173 .456 .130]; % old trials, microphone coordinates <x,y,z> 1 to 4 in meter
% December 2 trial
microphones_locs = [0.237 -0.1587765 -0.085; 0.008 -0.1587765 -0.09; 0.125 0.18334 0; 0.115	0 -0.32];
loc_diff = microphones_locs - speaker_loc;
c = 346; % speed of sound in meters/second 
buff_plot_num = 4; %set to multiple of 4% # of plots to see if buffered FFTs look ok, don't set this too high
% every channel is also plotted
%% Microphone & Source positioning
% Verification
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6359636
time = sqrt(sum((microphones_locs(:,:) - speaker_loc).^2,2))/c;
time_pairs = nchoosek(time,2); 
time_delay = time_pairs(:,1)-time_pairs(:,2); %delay from known sound source (mouth)
% TDOA for non-speaker for grid
% Adding TDOA for non-speaker
time_non_speaker = zeros([size(microphones_locs,1) size(non_speaker_locs_total, 1)]);
for i = 1:size(microphones_locs,1)
    time_non_speaker(i,:) = sqrt(sum((microphones_locs(i,:)-non_speaker_locs_total).^2,2))/c;
end
microphone_pairs = nchoosek([1:size(microphones_locs,1)], 2); %ie [1 2; 1 3; ... ; 1 n; 2 3 ; ... ; (n-1) n]
time_non_speaker_one = time_non_speaker(microphone_pairs(:, 1),:);
time_non_speaker_two = time_non_speaker(microphone_pairs(:, 2),:);
time_delay_non_speaker = time_non_speaker_one - time_non_speaker_two;

% code below calculates the tau sample points. Search the delta(kappa)
% array below 
% confirmed that this code below is equivalent
% sqrt(sum((mic_euc_pair_el1-mic_euc_pair_el2).^2,2))/c*Fs_target
% to code that follows
microphone_pairs = nchoosek([1:size(microphones_locs,1)], 2); %ie [1 2; 1 3; ... ; 1 n; 2 3 ; ... ; (n-1) n]
mic_euc = microphones_locs(microphone_pairs, :); %get pairwise matrix in 3D-space
mic_euc_pair_el1 = microphones_locs(microphone_pairs(:,1), :); %extract first element of pair, 
mic_euc_pair_el2 = microphones_locs(microphone_pairs(:,2), :); %extract second element of pair
mic_euc_pairs = mic_euc_pair_el1 - mic_euc_pair_el2; %subtract pair coordinates
mic_euc_pairs_distance = sqrt(sum(mic_euc_pairs.^2,2)); %Here calculate euclidean distance between n pairs 
tau_m = mic_euc_pairs_distance/c*Fs_target; % distance/speed of sound * sampling rate
%tau_m_3pairs = [tau_m-1 tau_m tau_m+1]; %use this code below for parabola interpolation calculation 
% n*(n-1)/2 pairs
% 
% function prototype
%fi = @(varargin)varargin{length(varargin)-varargin{1}};


fprintf("This script analyzes 2 to n microphones.\n");
fprintf("Only WAV files are accepted\n\n");

% Change and uncomment folder path here below if necessary

%% Data Input
cd(stmpath);

[file,stmpath,indx] = uigetfile('*.wav', 'MultiSelect', 'on'); %tic
% Code below is a fix for analysis when only one audio file is selected.
% In this case file will be of type char instead which breaks subsequent
% code which expects file to be of type cell.
if class(file) == 'char' %ensures conversion to cell if user selects just one audio file
    file = mat2cell(file,1,length(file));
end
% if max(size(file)) < 2
%     msg = ['Two or more audio files must be selected'];
%     ME = MException('MATLAB:mycode:userinput',msg);    
%     throw(ME);
% end
%tic 
%[file_path,name,ext] = fileparts(string(strcat(stmpath, char(file(1,:))))); %toc
% code speedup, read 1st file with audioread, then read the rest in a loop
% after pre-allocating an array.
[y, Fs] = audioread(strcat(stmpath, char(file(1))));
y_append = zeros(max(size(y)), max(size(file))-1,'double'); %for n audio files, pre-allocate n-1 arrays
Fs_append = zeros(1,3);
for i = 2:max(size(file)) %find number of files selected
    [y_append(:,i-1), Fs_append(i-1)] = audioread(strcat(stmpath, char(file(i)))); tic
end
y = [y y_append];
Fs = [Fs Fs_append];
clear y_append Fs_append;
%y_denoise = wdenoise(y);
Fs_orig = Fs;
Fs = Fs(1);
%T = 1/Fs;
L = length(y); 
%% Low pass filter - basic decimation
% Possibly implement custom FIR vs IIR filter later
if Fs > Fs_target + 1 %first check to make sure we sampled high enough
%     if mod(log2(Fs) - log2(Fs_target),1) == 0 % next check original sampling rate a multiple of 48k
%         dec_factor = log2(Fs) - log2(Fs_target); tic
%         y_out = zeros(ceil(max(size(y))/2^dec_factor), max(size(file)),'double'); %preallocate for speedup
%         for i = 1:max(size(file))
%             y_out(:,i) = decimate(y(:,i),2^dec_factor); %maybe replace with something custom later
%         end 
%         toc
%         y_orig = y;
%         y = y_out;
%         %clear y_out;
%         Fs = Fs_target;
    if mod(Fs, Fs_target) == 0 % support any factor of original sampling rate 
        dec_factor = Fs/Fs_target;
        y_out = zeros(ceil(max(size(y))/dec_factor), max(size(file)),'double'); %preallocate for speedup
        for i = 1:max(size(file))
            %y(:,i) = bandpass(y(:,i), [100 90000], Fs, 'ImpulseResponse','fir');
            y_out(:,i) = decimate(y(:,i),dec_factor); %maybe replace with something custom later
        end 
        %toc
        y_orig = y;
        y = y_out;
        %clear y_out;
        Fs = Fs_target;
    else
        fprintf("This frequency is unsupported\n\n");
    end
else
    fprintf("This frequency is too low for downsampling\n\n");
end
%% Placeholder for downsampling code




%% FFT full signal
%y_fft = fft(y(:,:), n);  
hann_window = hann(length(y), 'periodic'); %window first, apparently periodic useful for spectral analysis
y_window = hann_window.*y(:,:);
n = 2^nextpow2(size(y,1)); %used for padding signal with n - size(y,1) zeroes
y_fft = fft(y_window, n); %windowed signal padded
%% Buffered FFT - l_buff point buffer
% https://appliedacousticschalmers.github.io/scaling-of-the-dft/AES2020_eBrief/
% Window correction, section 4 equation (11)
% Still need to scale signals appropriately, accounting for windowing
hann_window_buff = hann(l_buff,'symmetric');
y_buffer = zeros(l_buff, ceil(max(size(y))/l_buff), max(size(file)));
for i = 1:max(size(file))
    y_buffer(:,:,i) = buffer(y(:,i), l_buff);
end
y_buff_win = y_buffer.*hann_window_buff;
y_buff_fft = fft(y_buff_win, l_buff); %/sum(hann_window_buff); % Implement other correction factors later
Fx_buffer = 1:Fs/l_buff:Fs/2+1; % (0:Fs/l_buff:Fs/2) 

%y_buffer_conj = conj(y_buffer);
%% STFT 
% STFT with Hann window of k points where 2^n = k & n is integer. Symmetric
% window, 50% overlap ie <= k/2.
%% Voice Activity Detection - Section 2
% https://appliedacousticschalmers.github.io/scaling-of-the-dft/AES2020_eBrief/#sec:scaling_broadband
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6359636
% Section 2.1 Spectral entropy
%Todo:1) fix indexing. attach 0 or NaN to omitted indices
%2) Assign a naive noise threshold solution, validate empirically and modify later
%   Since a pure sinusoid is being recorded in a relatively quiet room,
%   every frame should be a noise frame
%3) Calculate an SNR for assigned voice frames, that reclassifies voice
%frames as noise if SNR is <= 12, for example.
buff_power = abs(y_buff_fft(:,:,:)).^2/sum(hann_window_buff)^2; %Section 4, Equation 12 of link above
w_min = uint32(Voice_activity_detection_freq*l_buff/Fs+1);
w = w_min:1:l_buff/2;
buff_power_w = buff_power(w,:,:);
buff_power_w_sum(:,:) = sum(buff_power_w,1);
prob_pow = zeros(size(buff_power_w));
for i = 1:size(buff_power_w, 1)
    prob_pow(i, :, :) = squeeze(buff_power_w(i,:,:))./buff_power_w_sum;
end
H_t = squeeze((-1)*sum(prob_pow.*log2(prob_pow))); %Shannon's Entropy
t = 2:size(H_t, 1); %discard the first frame
H_t_fixed = max(H_t(t, :), H_t(t-1, :)); %Equation 3
bad_index = cell(1,size(H_t_fixed,2));
% Threshhold for noise frame set to < median - 1.5 absolute deviations
for i = 1:size(H_t_fixed, 2)
    bad_index{i} = uint32(find(H_t_fixed(:,i) < (median(H_t_fixed(:,i)) - 2*mad(H_t_fixed(:,i),1))));
end
%Empirically classify a frame as a noise frame
% Section 2.2 SNR Verification

% Implement later, not necessary

%% Compute TDOA - buffered signal
% https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6359636% 
% y_fft_conj = conj(y_fft); %For full length signal
y_buff_fft_conj = conj(y_buff_fft); %buffered signal, split into l_buff sizes
pairs_fft = nchoosek(1:max(size(file)),2); %We have n*(n-1)/2 unique microphones pairs
% In simplest case, for two microphones, Gx1x2f =% y_fft(:,1).*y_fft_conj(:,2);
% For all pairs, G_x1_u_x2_v_f = y_fft(:,pairs_fft(:,1)).*y_fft_conj(:,pairs_fft(:,2)); 
G_xu_xv_f = y_buff_fft(:,:, pairs_fft(:,1)).*y_buff_fft_conj(:,:, pairs_fft(:,2)); %for buffered signal
mag_G_xu_xv_f = abs(G_xu_xv_f); % Magnitude of signal
%change phase of Gxix2f
phase_G_xu_xv_f = angle(G_xu_xv_f);
k = uint32(l_buff/4+1:l_buff/2); %phase replacement, equation (8) of paper above
phase_G_xu_xv_f_pos(1:min(k)-1,:,:) = phase_G_xu_xv_f(1:min(k)-1, :, :);
phase_G_xu_xv_f_pos(k,:,:) = 2*phase_G_xu_xv_f(floor(k/2),:,:); % replace phase from [0 to pi)
%The FT of a real signal has the following relevant property: 
% angle X(jw) = (-1)*angle X(-jw)
%So correct the phase from 0 to Fs/2-Fs/N (as performed above)
%Then reconstruct the phase as [0, Fs/2-Fs/N]U[Fs/2]U(-Fs/2, -Fs/N] 
positive_pi = zeros(size(phase_G_xu_xv_f(1,:,:))); %which of these to use
positive_pi(1,:,:) = pi;
negative_pi = (-1)*positive_pi;
phase_G_xu_xv_f = [phase_G_xu_xv_f_pos; 
    zeros(size(phase_G_xu_xv_f(1,:,:))); 
    (-1)*phase_G_xu_xv_f_pos(l_buff/2:-1:2,:,:)]; %here indexing is decremented avoiding index 1 (DC)
G_xu_xv_f_fixed = mag_G_xu_xv_f.*exp(1j*phase_G_xu_xv_f);
R_yu_yv_tau = ifft(G_xu_xv_f_fixed./mag_G_xu_xv_f, 'symmetric'); %fixed Discrete pulse function delta(k) per paper above
delta = zeros(3, size(R_yu_yv_tau,2), size(R_yu_yv_tau,3));
kappa_guess = delta(1,:,:);  
for i = 1:max(size(tau_m)) % 0 to tau_m => 1 to ceil(tau_m + 1) in code
    [delta(2,:,i), kappa_guess(1,:,i)] = max(R_yu_yv_tau(1:uint32(tau_m(i)),:,i)); %ceiling of IDFT elements 0 to tau_m
end
kappa = permute([kappa_guess-1; kappa_guess; kappa_guess+1], [2 3 1]);
delta = permute(delta, [2 3 1]); 
for idx = 1:size(tau_m, 1)
    for kdx = 1:size(R_yu_yv_tau, 2)
        delta(kdx, idx, 3) = R_yu_yv_tau(kappa(kdx, idx, 3), kdx, idx);
        if kappa(kdx, idx, 1) ~= 0            
            delta(kdx, idx, 1) = R_yu_yv_tau(kappa(kdx, idx, 1), kdx, idx);
        else
            delta(kdx, idx, 1) = NaN;
        end
        %delta2(kdx, idx) = R_yu_yv_tau(kappa(kdx, idx, 2), kdx, idx);
    end
end
%a = tau_m_3pairs; %+1; %indexing in matlab starts at 1, housekeeping
%y_vertex = zeros(1,6);
x_vertex = zeros(size(R_yu_yv_tau,2), max(size(tau_m)));
x_vertex_testing = x_vertex;
sign_a = x_vertex; %tic %30 seconds for l_buff 2048, Fs_target 48000, 4 microphones, ~30 seconds
for i = 1:size(R_yu_yv_tau,2) %this outer for loop is timed at only .0107 seconds
    %b(:,:,i) = R_yu_yv_tau(ceil(a(i,:)),:,i)';
    %Will update code below with Vandermonde matrix or something different, fit and polyfit are
    %very slow
    for kdx = 1:max(size(tau_m)) %
        if kappa(i, kdx, 1) ~= 0
            %curve_fit = fit(squeeze(kappa(i, kdx, :)), squeeze(delta(i,kdx,:)), 'poly2'); % parabolic curve interpolation
            curve_fit = polyfit(squeeze(kappa(i, kdx, :)), squeeze(delta(i,kdx,:)), 2);
            %y_vertex(i,kdx) = curve_fit.p3 - curve_fit.p2^2/(4*curve_fit.p1); %c - b^2/4a
            %add code checking sign of a aka curve_fit.p1, if a is negative
            %(minima,unbounded), use kappa_guess instead. Decide what variable to use for
            %Xc. x_vertex, kappa_guess, etc.
            %sign_a(i, kdx) = sign(curve_fit.p1);
            sign_a(i, kdx) = sign(curve_fit(1));
            %x_vertex(i, kdx) = -curve_fit.p2/(2*curve_fit.p1); %the value Xc that maximizes
            x_vertex(i, kdx) = -curve_fit(2)/(2*curve_fit(1)); %the value Xc that maximizes
        else
            x_vertex(i, kdx) = NaN;
            x_vertex_testing(i, kdx) = delta(i, kdx, 2); %retain original guess
            %x_vertex(i, kdx) = delta(i, kdx, 2);
        end
    end
end %x_vertex is Xc in paper
toc
for i = 1:max(size(tau_m))
    fprintf("Microphone pair %d %d\nMean: %.6f\nMedian: %.6f\nStandard Deviation: %.6f\nMean Absolute deviation: %.6f\nMedian Absolute deviation: %.6f\n\n", microphone_pairs(i,1), microphone_pairs(i,2), mean(x_vertex(:,i),'omitnan'),median(x_vertex(:,i),'omitnan'), std(x_vertex(:,i),'omitnan'), mad(x_vertex(:,i),0), mad(x_vertex(:,i),1))
end
%% Direction Estimation (Source localization or Patient Identification)
%https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6359636
% Let first element of DV be the time delay for patient
DV = [time_delay time_delay_non_speaker]'; %Eq 9 generalized, (calculated above), theoretical
EV = x_vertex/Fs; %EV in seconds; x_vertex = distance/speed of sound * Fs

EV_DV_similarity = zeros(size(EV,1), size(DV,1));

for i = 1:size(DV,1)
    EV_DV_similarity(:,i) = sqrt(sum((EV(:,:)-DV(i,:)).^2,2));
end 
[EV_DV_similarity_val, EV_DV_similarity_index] = min(EV_DV_similarity,[],2);
%sum
speech_window_samples = uint32(speech_window_time*Fs/l_buff);
max_idx = floor(size(EV_DV_similarity_index,1)/speech_window_samples)-1;
idx_end = mod(size(EV_DV_similarity_index,1), speech_window_samples) + max_idx*speech_window_samples;
counter_v = 0; 
for i = 1:max_idx
    if sum(EV_DV_similarity_index(1+(i-1)*speech_window_samples:i*speech_window_samples) == 1) > (3*speech_window_samples/5)
        counter_v = counter_v + 1;
    end
end
percentage_spent_speaking = 100*counter_v / max_idx;
fprintf("\nSpeech from pt detected %f%% of time\n", percentage_spent_speaking);
%DV_rm = repmat(DV', [1 max(size(EV))])';
%geo_dis = sum((EV-DV).^2,2); %geometric distance?
geo_dis_log = log(geo_dis);
% ACMP after this
azimuth_angle_start = -pi;
azimuth_angle_finish = pi;
elevation_angle_start = -pi/2;
elevation_angle_finish = pi/2;
angle_increment = pi/36;
azimuth_angle = azimuth_angle_start:angle_increment:azimuth_angle_finish;
elevation_angle = elevation_angle_start:angle_increment:elevation_angle_finish;
gamma_angle = azimuth_angle;
two_d_matrix = azimuth_angle'*elevation_angle;
%[azimuth_conv, elevation_conv, r_conv] = cart2sph(microphones_locs(:,1), microphones_locs(:,2), microphones_locs(:,3));
%[azimuth_pt, elevation_pt, r_pt] = cart2sph(speaker_loc(1), speaker_loc(2), speaker_loc(3));
[azi_diff_pt, ele_diff_pt, r_diff_pt] = cart2sph(loc_diff(:,1), loc_diff(:,2), loc_diff(:,3)); %
% Acoustic Map, Citation 8 from primary paper
% DOI: 10.1109/HSCMA.2008.4538690     
% SRP-PHAT aka Global Coherence Field (GCF)
% Citations 9 and 10 from secondary paper speak more about GCF
% DOI:
% DOI:
% (respectively)
%% Continue power code
pow_y = abs(y_fft).^2; %get power by squaring the magnitude elementwise
% For windowed power signals, divide by sum of squares of window elements
% https://appliedacousticschalmers.github.io/scaling-of-the-dft/AES2020_eBrief/#sec:scaling_broadband
pow_y = pow_y/(hann_window'*hann_window);
% hann_window'*hann_window which is < length(y)
num_unique_pts = length(pow_y)/2+1; % fft is 2-sided; make one-sided index with 
% 0
pow_y = pow_y(1:num_unique_pts,:); %convert fft to 1-sided with index
pow_y = pow_y*2; %power correction factor for one-sided fft, technically don't multiply by two at 0, fix this later
pow_y(2:end-1,:) = pow_y(2:end-1,:)*2; %POWER correction factor for hann window (2)
%amplitude is now identical but energy is not (1.63)
Pxx1 = pow_y/Fs; %divide power by sampling rate to obtain PSD
Fx1 = (0:num_unique_pts-1)*Fs/n;
%phase_y = unwrap(angle(fftshift(y_fft))); %angle is also two sided
%phase_y = phase_y(length(phase)/2,:);
%[y,Fs] = audioread(filename);

%% FFT based on windowing 

%% Welch Periodogram
%https://www.mathworks.com/matlabcentral/answers/33653-psd-estimation-fft-vs-welch
if Enable_welch == true
    [Pxx2, Fx2] = pwelch(y, hann_window, [], n, Fs);
end
%% Plotting spectogram, PSD and Welch 
%Basic signal integrity information
%yticks(-150:2:-40);
% Spectrogram
hann_spect = hann(Fs/20); %toc
if Enable_welch == true
    for i = 1:length(file)
        figure
        subplot(2,1,1);
        plot(Fx1/1000, 10*log10(Pxx1(:,i)),'r:');
        xlabel('Frequency (kHz)')
        ylabel('Signal Power (dB) ')
        str = sprintf('Figure %d: Power Spectral Density with Hann Window', i);
        title(str);
        subplot(2,1,2);
        plot(Fx2/1000, 10*log10(Pxx2(:,i)),'r:');
        xlabel('Frequency (kHz)')
        ylabel('Signal Power (dB) ')
        str = sprintf('Figure %d: PWelch with Hann Window', i);
        title(str);
    end
else
    for i = 1:length(file)
        figure
        plot(Fx1/1000, 10*log10(Pxx1(:,i)),'r:');
        xlabel('Frequency (kHz)')
        ylabel('Signal Power (dB)')
        str = sprintf('Figure %d: Power Spectral Density with Hann Window', i);
        title(str);
        figure
        semilogx(Fx1, 10*log10(Pxx1(:,i)),'r:');
        xlabel('Frequency (kHz)')
        ylabel('Signal Power (dB)')
        str = sprintf('Figure %d: Power Spectral Density with Hann Window', i);
        title(str);
    end
end
% for i = 1:length(file)
%     figure
%     plot(Fx1, phase_y(:,i),'r:')
%     str = sprintf('Figure %d: Phase', i);
%     title(str);
% end
%% Sanity check for buffer code
% finish rest tomorrow
% plot channels different files
rand_plots = randsample(size(y_buff_fft,2), 4*ceil(buff_plot_num/4)); %in case user selects non-multiple of 4
for h = 1:max(size(file))
    for i = 1:ceil(buff_plot_num/4)
        %         figure
        %         plot(Fx_buffer/1000, 10*log10(y_buff_fft( 1:length(Fx_buffer), rand_plots(i), 1)),'r:');
        %         xlabel('Frequency (kHz)')
        %         ylabel('Signal Power (dB)')
        %         str = sprintf('Figure %d: Magnitude of Buffered Signal %d with Hann Window', i, rand_plots(i));
        %         title(str);
        figure
        subplot(2,2,1);
        semilogx(Fx_buffer, 10*log10(imag(y_buff_fft( 1:length(Fx_buffer), rand_plots(4*i-3), h))),'r:');
        xlabel('Frequency (kHz)')
        ylabel('Signal Power (dB)')
        str = sprintf('Figure %d: Imag part of Buffered Signal %d on Channel %d with Hann Window', i, rand_plots(4*i-3),h);
        title(str);
        subplot(2,2,2);
        semilogx(Fx_buffer, 10*log10(imag(y_buff_fft( 1:length(Fx_buffer), rand_plots(4*i-2), h))),'r:');
        xlabel('Frequency (kHz)')
        ylabel('Signal Power (dB)')
        str = sprintf('Figure %d: Imag part of Buffered Signal %d on Channel %d with Hann Window', i, rand_plots(4*i-2),h);
        title(str);
        subplot(2,2,3);
        semilogx(Fx_buffer, 10*log10(y_buff_fft( 1:length(Fx_buffer), rand_plots(4*i-1), h)),'r:');
        xlabel('Frequency (kHz)')
        ylabel('Signal Power (dB)')
        str = sprintf('Figure %d: Real part of Buffered Signal %d on Channel %d with Hann Window', i, rand_plots(4*i-1),h);
        title(str);
        subplot(2,2,4);
        semilogx(Fx_buffer, 10*log10(y_buff_fft( 1:length(Fx_buffer), rand_plots(4*i), h)),'r:');
        xlabel('Frequency (kHz)')
        ylabel('Signal Power (dB)')
        str = sprintf('Figure %d: Real part of Buffered Signal %d on Channel %d with Hann Window', i, rand_plots(4*i),h);
        title(str);
    end
end
%% Noise Cancellation


