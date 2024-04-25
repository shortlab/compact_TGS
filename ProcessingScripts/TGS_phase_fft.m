function [fft] = TGS_phase_fft(SAW_only,psd_out,truncate_fraction)
%   Generates a filtered fft of the SAW data. It is a revised version of
%   TGS_phase_fft that doesn't read the raw files. Instead, it is made to
%   work in the context of TGSPhaseAnalysis.
%   SAW_only: This should be the total_signal w/ an erfc fit subtracted
%   from it so it's macroscopically flat.
%   psd: if psd=1, save out power spectrum, else save out fft magnitude
%   truncate_fraction: What percentage of the signal should be analyzed (i.e,
%   first 70%)

%Boolean options for plotting and processing
derivative=1;
copy=0; % set true if you want to double the signal, mirroring it to artificially increase the frequency resolution. This is a quirk of the matlab way of doing power spectral densities
plotfft=0;
saveout=0;
cut_tails=1000; % Fraction of data removed from beginning and end of power spectral density spectra

% Shorten 'SAW_only' to 'SAW_only_truncated' so the FFT samples better data
% and normalise it to the maximum value
newlength=ceil(numel(SAW_only(:,1))*truncate_fraction);
SAW_only_truncated=zeros(newlength,2);
SAW_only_truncated(:,1)=SAW_only(1:newlength,1);
SAW_only_truncated(:,2)=SAW_only(1:newlength,2)/max(SAW_only(1:newlength,2));

if copy
    SAW_only_truncated_mirrored=zeros(newlength*2-1,2);
    for i=1:length(SAW_only_truncated(:,1))
        SAW_only_truncated_mirrored(i,1)=SAW_only_truncated(i,1);
        SAW_only_truncated_mirrored(i,2)=SAW_only_truncated(i,2);
    end
    for i=length(SAW_only_truncated(:,1))+2:length(SAW_only_truncated_mirrored(:,1))+1
        SAW_only_truncated_mirrored(i-1,2)=SAW_only_truncated((2*length(SAW_only_truncated(:,1))+1)-i,2);
        SAW_only_truncated_mirrored(i-1,1)=SAW_only_truncated(i-length(SAW_only_truncated(:,1)),1)+SAW_only_truncated(length(SAW_only_truncated(:,1)),1);
    end
    clear SAW_only_truncated
    SAW_only_truncated=SAW_only_truncated_mirrored;
end

% Time step info necessary for differentiation and flat padding
tstep=SAW_only_truncated(end,1)-SAW_only_truncated(end-1,1);

% If option selected, take transform of derivative of recorded signal
% to filter out DC even more than just the background subtraction. This
% derivative is then normalized by the maximum.
SAW_only_derivative=diff(SAW_only_truncated(:,2))/tstep;
SAW_only_derivative=SAW_only_derivative/max(SAW_only_derivative);
if derivative
    SAW_only_truncated=[SAW_only_truncated(1:length(SAW_only_derivative),1) SAW_only_derivative];
end

% Find the stuff we need to take the spectral profile
num_points=length(SAW_only_truncated(:,1));
sampling_rate=num_points/(SAW_only_truncated(end,1)-SAW_only_truncated(1,1));
padding=18; %magnitude of zero padding to increase resolution in power spectrum
padsize=2^padding-num_points-2; %more padding = smoother transform

%Only pad on the positive end
pad_val=0;
pad=zeros(padsize,1);
pad(1:end)=pad_val;
padded_end_time=SAW_only_truncated(end,1):tstep:SAW_only_truncated(end,1)+(padsize-1)*tstep; %this makes sure that the padded data has a sensible end point in time, consistent with the regular data.

SAW_only_padded=[SAW_only_truncated(:,1) SAW_only_truncated(:,2);padded_end_time' pad];

number_of_fft_points=length(SAW_only_padded(:,2));
%Find the Power Spectral density

%Use a hamming window and a Welchs method. Hamming does the best of the
%ones I've tried and Welch does slightly better than the normal
%periodogram.
[power_spectral_density,freq]=periodogram(SAW_only_padded(:,2),rectwin(number_of_fft_points),number_of_fft_points,sampling_rate); %periodogram method

psd_length_adjust=ceil(numel(power_spectral_density)*(1/5))-6*cut_tails; %trim off the end of the fft
power_spectral_density(1:cut_tails)=0;
power_spectral_density(psd_length_adjust:end)=0;
power_spectral_density=power_spectral_density/(max(power_spectral_density));

%Don't save out DC spike in FFT/PSD
if psd_out
    output_amplitude=power_spectral_density(1:end);
else
    output_amplitude=sqrt(power_spectral_density(1:end));
end

fft=[freq(1:end-1) sgolayfilt(output_amplitude(1:length(freq(1:end-1))),5,201)];

if saveout
    dlmwrite('dat_spec.txt',out);
end

if plotfft
    figure()
    axes('Position',[0 0 1 1],'xtick',[],'ytick',[],'box','on','handlevisibility','off','LineWidth',5)
    plot(freq(1:end-1),sgolayfilt(output_amplitude(1:length(freq(1:end-1)))/max(output_amplitude),5,201),'k','LineWidth',2);
    hold on
    xlim([0 1.0e9]);
    set(gca,...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',24,...
        'FontName','Times',...
        'LineWidth',3)
    ylabel({'Intensity [a.u.]'},...
        'FontUnits','points',...
        'FontSize',24,...
        'FontName','Times')
    xlabel({'Frequency [Hz]'},...
        'FontUnits','points',...
        'FontSize',24,...
        'FontName','Times')
end
end
