function [wav_high , wav_low] = filterhjb(wav , f_cutoff , df , av_yes);

% function [wav_tot_high , wav_tot_low] = filterhjb(wav , f_cutoff , df , av_yes);
%
% This function divides a wave signal in a high frequency part and a low frequency part
% separated by a cutoff frequency.
% Input variables:
%   - wav       : wave signal, series in columns, must be even number of points;
%   - f_cutoff  : cut-off frequency;
%   - df        : frequency step;
%   - av_yes    : include average value yes = 1, no = 0 (default is no);
%
% Author: Henk Jan Bakkenes

if nargin <4
    av_yes = 0;
end

ii_cutoff = round(f_cutoff/df); % calculates integer for cut off frequency

A_wav = fft(wav);

A_high = zeros(size(wav));
A_low  = zeros(size(wav));

size(wav);
length(wav);
A_high(ii_cutoff+2:length(wav)/2,:) = A_wav(ii_cutoff+2:length(wav)/2,:);

if av_yes
    A_low(1,:) = 0.5*A_wav(1,:);
    A_low(2:ii_cutoff+1,:) = A_wav(2:ii_cutoff+1,:);
else
    A_low(2:ii_cutoff+1,:) = A_wav(2:ii_cutoff+1,:);
end

wav_high = 2*real(ifft(A_high)); % 2 compensates for the missing 'mirror' fourier coefficients
wav_low = 2*real(ifft(A_low));