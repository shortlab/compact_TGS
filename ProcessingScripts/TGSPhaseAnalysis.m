%   Function to process TGS signals, fitting the TGS response equation - eq.
%   12 in Dennett et al. (2018) [https://doi.org/10.1063/1.5026429] - using
%   a Levenberg-Marquardt nonlinear least squares method.
%   Created by C.A. Dennett.
%   Modified by B.R. Dacus, A.P.C. Wylie, K. Zoubkova and S. Engebretson.


function [frequency_final,frequency_error,SAW_speed,thermal_diffusivity,thermal_diffusivity_err,acoustic_damping_constant, acoustic_damping_error, A, A_err, displacement_reflectance_ratio, displacement_reflectance_ratio_err, B, B_err, acoustic_phase, acoustic_phase_err, C, C_err] = TGSPhaseAnalysis(pos_file,neg_file,grating,start_point,two_SAW_frequencies,baselineBool,POSbaselineStr,NEGbaselineStr, delimiter, verbose, header_length, plot_things, mono_heterodyne)
%   Data is saved in two files, positive (with one heterodyne phase of pi/2) and
%   negative (with another of -pi/2), must provide both files
%
%   ------Inputs------
%   pos_file:                           positive phase TGS data file
%   neg_file:                           negative phase TGS data file
%   grating:                            calibrated grating spacing in um
%   start_point:                        provide integer between 1 and 4 to pick the null-point start
%                                       from which fit will begin
%   two_SAW_frequencies:                a boolean value, default value is 0 if only one acoustic mode
%                                       is present in the measurement. If two are present, provide 1 for both
%                                       frequencies and speeds to be output.
%   baselineBool:                       a boolean value telling the script if baseline subtraction for your
%                                       data is desired. 1 for 'yes do baseline subtraction', 0 for 'don't'.
%   POSbaselineStr:                     file for the positive baseline subtraction. Note that
%                                       this must always be supplied, even if no subtraction is
%                                       desired - supply a dummy file in this case.
%   NEGbaselineStr:                     file for the negative baseline subtraction. Again a file must be supplied.
%   delimiter:                          the delimiter of the scope output
%   verbose:                            whether to print outputs
%   header_length:                      the number of lines in the header
%                                       of the scope output
%   plot_things:                        whether to plot outputs
%   mono_heterodyne:                    whether the setup is single
%                                       heterodyne
%
%   -----Outputs-----
%   frequency_final:                    SAW frequency [Hz] from signal FFT
%   frequency_error:                    error on SAW frequency [Hz] (1 σ)
%   SAW_speed:                          SAW speed [ms^-1]
%   thermal_diffusivity:                thermal diffusivity [m^2s^-1]
%   thermal_diffusivity_err:            thermal diffusivity err [m^2s^-1] taken using 95%
%                                       confidence interval. The rest of the errors in the outputs use this
%                                       confidence interval-based method.
%   acoustic_damping_constant:          acoustic damping parameter [s]
%   acoustic_damping_error:             acoustic damping parameter error [s]
%   A:                                  diffracted response amplitude [Wm^-2]
%   A_err:                              diffracted response amplitude error [Wm^-2]
%   displacement_reflectance_ratio      ratio of diffraction contributions from thermoreflectance changes
%                                       to those of surface displacement [s^0.5]. This ratio is also 
%                                       referred to as beta occasionally in these scripts
%   displacement_reflectance_ratio_err: error on above ratio [s^0.5].
%   B:                                  signal sinusoid contribution amplitude [Wm^-2].
%   B_err:                              error on signal sinusoid contribution amplitude [Wm^-2].
%   acoustic_phase                      phase of signal sinusoin contribution [].
%   acoustic_phase_err                  error on phase of signal sinusoin contribution [].
%   C                                   fitting constant [Wm^-2].
%   C_err                               error on fitting constant [Wm^-2].

%%%%%%%%%%%%%%%%%%%
%USER SETTINGS
%These control lots of plot outputs and some fitting settings for signals
%that have more than one acoustic mode.
%%%%%%%%%%%%%%%%%%%
%Settings for various plotting and output options to be set by boolean arguments
find_max=0; %gives the option to manually locate the beginning time t_0
plot_everything=0; %plots a bunch of extra graphs at certain steps to check things
plot_trace=plot_things; %plots pos - neg file data
psd_out=1; % generates power spectrum density plots rather than FFT magnitudes
plot_psd=plot_things; %adds bonus lines in lorentzian_peak_fit plotting
plot_final=plot_things; % Generates the final plots with the fft inserts

%updated tstep parsing to better handle variable sampling rates that can happen if the scope is left in 4-channel mode rather than 2-channel mode.
%in general, this should never happen as 20 gigasamples per second is how the scope should be configured, but accidents happen and this renders
%those data as salvageable.
%timestep found by parsing first few time data points
% pos_aux = readvars(pos_file,"NumHeaderLines",header_length); %[[[Maybe roll this in when the first total_signal is created to minimise the number of times that the file is called]]]
% scope_timebase = pos_aux(2) - pos_aux(1);
no_pre_calc=0; %use this if you want to hard-code in your initial guesses for diffusivity and a few other parameters
amp_grat=0; % say 1 if you want to do ONLY the reflectance exponential fit for the data and nothing else
steel=0; % this makes the script call find_freq_steel, which fits a sin term directly to the SAW component, rather than determining the frequency via the peak of the FFT (old and not often used)
printFitsToFiles = 0; %outputs both thermal and full fit data to output file

%FFT settings for TGS_phase_fft.m
truncate_fraction=0.9;      % The fraction of 'SAW_only' that is used to calculate the FFT, can adjust this if SAWs only make up a very small portion of the signal

%if on, record the amplitude ratio of A/B in final fit at the end of tau
%vector (so it is four elements instead of 3)
record_amp_ratio=0; %[[[THIS MAY GET CHANGED FOR TAU REWORK]]]

spot_size=140e-6; %measured pump spot size of 140um - MIT, used to find the time for a SAW to leave the pump spot

derivative=0; %plots the derivative of the total signal

%%%%%%%%%%%%%%%%%%%
%Tau selection block. Allows you to select the final model using a variety
%of time constant schemes
%%%%%%%%%%%%%%%%%%%
%if fitting tau need to select one of three options for start point
tau_final_fit=1;
start_constant=1;
start_lorentz=0;
start_walkoff=0;
%if not fitting tau, need to select one of two constant schemes
fixed_tau=0;
fixed_lorentz=0;
fixed_walkoff=0;
%check schemes to make sure that the code will always run. default to
%fitting with constant values
if tau_final_fit==0 && fixed_tau==0
    tau_final_fit=1;
elseif tau_final_fit==1 && fixed_tau==1
    tau_final_fit=0;
end
if tau_final_fit==1 && start_constant==0 && start_lorentz==0 && start_walkoff==0
    start_constant=1;
end
if fixed_tau==1 && fixed_lorentz==0 && fixed_walkoff==0
    fixed_walkoff=1;
end
%%%%%%%%%%%%%%%%%%%
%%
%Kristyna addition to handle baseline arguments
% if nargin < 6
%     baselineBool = 0;
%     POSbaselineStr = NaN;
%     NEGbaselineStr = NaN;
% end
% 
% if nargin == 6 && baselineBool == 0
%     POSbaselineStr = NaN;
%     NEGbaselineStr = NaN;
% elseif nargin == 6 && baselineBool == 1 % can be either error or warning and no subtraction
%     error('For baselineBool = 1 it is necessary to provide positive and negative baseline signals.')
%     %     warning('baselineBool = 1, but positive and negative baseline signals not provided -> continuing without baseline subtraction')
%     %     baselineBool = 0
%     %     POSbaselineStr = NaN;
%     %     NEGbaselineStr = NaN;
% end
% 
% if nargin == 7 && baselineBool == 0
%     warning('baselineBool = 0 -> continuing without baseline subtraction')
%     POSbaselineStr = NaN;
%     NEGbaselineStr = NaN;
% elseif nargin == 7 && baselineBool == 1
%     error('For baselineBool = 1 it is necessary to provide BOTH positive and negative baseline signals.')
% end
% 
% if baselineBool == 1 && (isempty(POSbaselineStr) || isempty(NEGbaselineStr))
%     error('baselineBool = 1 but positive or negative baseline signal not provided')
% end
%End Kristyna addition to handle baseline arguments
if baselineBool == 0
    POSbaselineStr = NaN;
    NEGbaselineStr = NaN;
end

q=2*pi/(grating*10^(-6));

%How far from guessed values for diffusivity and beta do you vary in the
%end, d for diffusivity b for beta. Diffusivity is most important
percent_range_d=0.45;
percent_range_b=2.25;
if tau_final_fit && ~start_constant
    if start_walkoff
        percent_range_t=2.5;
    elseif start_lorentz
        percent_range_t=0.9;
    end
end

%Output tau will have multiple values, initialize variable.
%tau(1) is the lorentzian decay time
%tau(2) is the walk-off time
%tau(3) is the final value either used in or optimized by fit
acoustic_damping_constant=zeros(1,3);

% if nargin<5
%     two_SAW_frequencies=0;
% end

%Difference in file write format based on newer or older acquisition. hdr_len should be 16 for the Ge dataset

%read in data files for this procedure
pos=dlmread(pos_file,delimiter,header_length,0);
if mono_heterodyne == 0
    neg=dlmread(neg_file,delimiter,header_length,0);
else
    neg=pos;
    neg(:, 2) = 0;
end
if baselineBool
    posbas=dlmread(POSbaselineStr,delimiter,header_length,0);
    if mono_heterodyne == 0
        negbas=dlmread(NEGbaselineStr,delimiter,header_length,0);
    else
        negbas=posbas;
        negbas(:, 2) = 0;
    end
end
scope_timebase = pos(2, 1) - pos(1, 1);

%sometimes written data is off by one time step at the end, chop that off if they do not match
if length(pos(:,1))>length(neg(:,1))
    pos=pos(1:length(neg(:,1)),:);
    if baselineBool
        posbas=posbas(1:length(neg(:,1)),:);
    end
elseif length(neg(:,1))>length(pos(:,1))
    neg=neg(1:length(pos(:,1)),:);
    if baselineBool
        negbas=negbas(1:length(pos(:,1)),:);
    end
end

if baselineBool
    pos(:,2)=pos(:,2)-posbas(:,2);
    neg(:,2)=neg(:,2)-negbas(:,2);
end

%normalize each set of data to the zero level before the pump impulse
pos(:,2)=pos(:,2)-mean(pos(1:50,2));
neg(:,2)=neg(:,2)-mean(neg(1:50,2));

%%%%%Time indexing block%%%%%%%%
[pump_time_index,end_time]=findTimeIndex(pos,neg);
if verbose
    display(['Time Index is: ',num2str(pump_time_index)])
    display(['End Time is: ',num2str(end_time)])
end
time_naught=neg(pump_time_index,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end_index=floor(end_time/scope_timebase)-36;

%re-normalize data to end signal decayed state if grating decays entirely during collection window, if not, do not re-normalize
if grating<8
    base_index=floor(length(pos(:,2))/5);
    pre_signal_average=mean((pos(1:50,2)-neg(1:50,2)));
else
    pre_signal_average=0;
end

if plot_trace
    figure()
    plot(neg(:,1)*10^9,(pos(:,2)-neg(:,2)-pre_signal_average)*10^3,'-','Color',[0 0 0],'LineWidth',1.25)
    hold on
    xlim([0 (end_time/2)*10^9])
    set(gca,...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',16,...
        'FontName','Times',...
        'LineWidth',1)
    ylabel({'Amplitude [mV]'},...
        'FontUnits','points',...
        'FontSize',20,...
        'FontName','Times')
    xlabel({'Time [ns]'},...
        'FontUnits','points',...
        'FontSize',20,...
        'FontName','Times')
    txt3 = {['Time Index = ',num2str(pump_time_index)],['End Time = ',num2str(end_time)]};
    saveas(gcf,"TGS_Trace.png")
end

if start_point==0
    thermal_diffusivity=0;
    thermal_diffusivity_err=0;
    acoustic_damping_error = 0;
    A= 0;
    A_err= 0;
    displacement_reflectance_ratio= 0;
    displacement_reflectance_ratio_err= 0;
    B= 0;
    B_err= 0;
    acoustic_phase= 0;
    acoustic_phase_err= 0;
    C= 0;
    C_err= 0;
else

    if amp_grat==0
        total_signal=[pos(pump_time_index:end_index,1)-time_naught pos(pump_time_index:end_index,2)-neg(pump_time_index:end_index,2)-pre_signal_average];
        if derivative
            der_len=length(total_signal(:,1))-1;
            fixed_derivative=zeros(1,der_len);
            for jj=1:der_len
                fixed_derivative(jj)=(total_signal(jj+1,2)-total_signal(jj,2))/scope_timebase;
            end
            figure()
            plot(total_signal(1:der_len,1),fixed_derivative,'k-')
            title('This is the derivative of fixed short')
        end

        %if you don't want to automatically find the peak t_0 of the profile, setting find_max to true above will allow
        %you to select a region on the plot within with to search. Useful if there are initial transients.
        if find_max || plot_everything
            figure()
            plot(total_signal(:,1),total_signal(:,2),'k-')
            xlim([0 1.5*10^-7])
            title('this is fixed short');
            if find_max
                hold on
                [x_cord,~]=ginput(2);
                neg_x_cord=x_cord(1);
                pos_x_cord=x_cord(2);
                pos_x_ind=floor(pos_x_cord/scope_timebase);
                neg_x_ind=floor(neg_x_cord/scope_timebase);
                [~,max_index]=max(total_signal(neg_x_ind:pos_x_ind,2));
                time_max=total_signal(max_index+neg_x_ind,1);
                close(gcf)
            end
        end

        %Otherwise, find t_0 from the profile directly
        if ~find_max
            time_offset_index=20; %was 20 before
            [~,max_index]=max(total_signal(time_offset_index:end,2));
            max_index=max_index+time_offset_index-1;
            time_max=total_signal(max_index,1);

            start_time_phase=find_start_phase(total_signal(max_index:end,1),total_signal(max_index:end,2),start_point,grating);
            start_index_master=round(start_time_phase/scope_timebase)+1;
            start_time_master=total_signal(start_index_master,1);

            %Fitting parameters for initial naive fit
            LB=[0 0];
            UB=[1 5*10^-4]; % Increased to account for silver-alloys, which are REALLY HIGH thermal diffusivity
            ST=[.05 5*10^-6];

            OPS=fitoptions('Method','NonLinearLeastSquares','Lower',LB,'Upper',UB,'Start',ST);
            TYPE=fittype('A.*erfc(q*sqrt(k*(x+time_max)))','options',OPS,'problem',{'q','time_max'},'coefficients',{'A','k'});
            [f0,gof]=fit(total_signal(:,1),total_signal(:,2),TYPE,'problem',{q,time_max});

            thermal_diffusivity=f0.k;
            con_int_error=confint(f0,0.95);
            %factor of 2 makes the 1 sigma confidence interval come out
            thermal_diffusivity_err=[thermal_diffusivity-con_int_error(1,2) con_int_error(2,2)-thermal_diffusivity]/2;

            if plot_everything
                figure()
                plot(total_signal(:,1),total_signal(:,2),total_signal(:,1),f0(total_signal(:,1)))
                hold on
                title('First naive fit')
            end

            SAW_only=[total_signal(:,1) total_signal(:,2)-f0(total_signal(:,1))];
            [fft] = TGS_phase_fft(SAW_only,psd_out,truncate_fraction);
            [frequency_final,frequency_error,tau_lorentz]=lorentzian_peak_fit(fft,two_SAW_frequencies,plot_psd);
            acoustic_damping_constant(1)=tau_lorentz(1);

            SAW_speed=frequency_final*grating*10^(-6);
            walk_off=spot_size/(SAW_speed(1)*2);
            acoustic_damping_constant(2)=walk_off;
            peak_freq=frequency_final(1);

            %We'll call the parameter beta the ratio of the amplitudes of the
            %displacement versus temperature grating. beta is a small number in metals.
            %Note beta here eventually becomes the
            %displacement_reflectance_ratio output of the script

            for jj=1:10

                beta=q*sqrt(thermal_diffusivity/pi)*(q^2*thermal_diffusivity+1/(2*time_max))^(-1);

                start_time=start_time_master;
                start_index=start_index_master;

                %Conduct initial parameter estimation without using an sin(x) contribution to the fit

                LB1=[0 0];
                UB1=[1 5*10^-4];   % Increased to account for silver-alloys, which are REALLY HIGH thermal diffusivity
                ST1=[.05 10^-5];

                OPS1=fitoptions('Method','NonLinearLeastSquares','Lower',LB1,'Upper',UB1,'Start',ST1);
                TYPE1=fittype('A.*(erfc(q*sqrt(k*(x+start_time)))-beta*exp(-q^2*k*(x+start_time))./sqrt((x+start_time)))','options',OPS1,'problem',{'q','beta','start_time'},'coefficients',{'A','k'});
                [f1,gof]=fit(total_signal(start_index:end,1),total_signal(start_index:end,2),TYPE1,'problem',{q,beta,start_time});

                thermal_diffusivity=f1.k;
                con_int_error=confint(f1,0.95);
                %factor of 2 makes the 1 sigma confidence interval come out
                thermal_diffusivity_err=[thermal_diffusivity-con_int_error(1,2) con_int_error(2,2)-thermal_diffusivity]/2;

                if plot_everything
                    figure()
                    plot(total_signal(:,1),total_signal(:,2))
                    hold on
                    title(strcat('Fit number ',num2str(jj+1),' - fixed beta'))
                end
            end

            %If you've elected not to pre-compute, provide hard initial guesses for diffusivity and beta based on
            %the bulk value of Ge diffusivity. Set ranges for final fit.
            if no_pre_calc
                thermal_diffusivity=0.3636*10^-4;
                beta=2e-5;
                low_bound=[1e-5 0];
                up_bound=[1e-3 1e-4];
            else
                low_bound=[thermal_diffusivity*(1-percent_range_d) beta*(1-percent_range_b)];
                up_bound=[thermal_diffusivity*(1+percent_range_d) beta*(1+percent_range_b)];
                if percent_range_d>1
                    low_bound(1)=0;
                end
                if percent_range_b>1
                    low_bound(2)=0;
                end
            end

            start_time2=start_time_master;
            start_index2=start_index_master;

            if tau_final_fit
                if start_constant
                    low_t=5e-9;
                    up_t=7e-7;
                    start_tau=1e-9;
                elseif start_lorentz
                    low_t=acoustic_damping_constant(1)*(1-percent_range_t);
                    up_t=acoustic_damping_constant(1)*(1+percent_range_t);
                    if percent_range_t>1
                        low_t(1)=0;
                    end
                    start_tau=acoustic_damping_constant(1);
                    LB2=[0 low_bound(1) low_bound(2) 0 -2*pi low_t -5e-3];
                    UB2=[1 up_bound(1) up_bound(2) 10 2*pi up_t 5e-3];
                    ST2=[.05 thermal_diffusivity beta 0.05 0 acoustic_damping_constant(1) 0];
                elseif start_walkoff
                    low_t=acoustic_damping_constant(2)*(1-percent_range_t);
                    up_t=acoustic_damping_constant(2)*(1+percent_range_t);
                    if percent_range_t>1
                        low_t=0;
                    end
                    start_tau=acoustic_damping_constant(2);
                end

                LB2=[0 low_bound(1) low_bound(2) 0 -2*pi low_t -5e-3];
                UB2=[1 up_bound(1) up_bound(2) 10 2*pi up_t 5e-3];
                ST2=[.05 thermal_diffusivity beta 0.05 0 start_tau 0];

                OPS2=fitoptions('Method','NonLinearLeastSquares','Lower',LB2,'Upper',UB2,'Start',ST2);
                TYPE2=fittype('A.*(erfc(q*sqrt(k*(x+start_time)))-beta*exp(-q^2*k*(x+start_time))./sqrt((x+start_time)))+B.*sin(2*pi*(peak_freq)*(x+start_time)+p)*exp(-(x+start_time)/t)+D','options',OPS2,'problem',{'q','start_time','peak_freq'},'coefficients',{'A','k','beta','B','p','t','D'});
                [f2,gof]=fit(total_signal(start_index2:end,1),total_signal(start_index2:end,2),TYPE2,'problem',{q,start_time2,peak_freq});

            elseif fixed_tau
                if fixed_lorentz
                    pp_tau=acoustic_damping_constant(1);
                elseif fixed_walkoff
                    pp_tau=acoustic_damping_constant(2);
                end
                LB2=[0 low_bound(1) low_bound(2) 0 -2*pi -5e-3];
                UB2=[1 up_bound(1) up_bound(2) 10 2*pi 5e-3];
                ST2=[.05 thermal_diffusivity beta 0.05 0 0];

                OPS2=fitoptions('Method','NonLinearLeastSquares','Lower',LB2,'Upper',UB2,'Start',ST2);
                TYPE2=fittype('A.*(erfc(q*sqrt(k*(x+start_time)))-beta*exp(-q^2*k*(x+start_time))./sqrt((x+start_time)))+B.*sin(2*pi*(peak_freq)*(x+start_time)+p)*exp(-(x+start_time)/t)+D','options',OPS2,'problem',{'q','start_time','peak_freq','t'},'coefficients',{'A','k','beta','B','p','D'});
                [f2,gof]=fit(total_signal(start_index2:end,1),total_signal(start_index2:end,2),TYPE2,'problem',{q,start_time2,peak_freq,pp_tau});
            end

            if verbose
                display(f2)
            end

            thermal_diffusivity=f2.k;
            A= f2.A;
            displacement_reflectance_ratio= f2.beta;
            B= f2.B;
            acoustic_phase= f2.p;
            C= f2.D;

            con_int_error=confint(f2,0.95);
            %factor of 2 makes the 1 sigma confidence interval come out
            thermal_diffusivity_err=(thermal_diffusivity-con_int_error(1,2))/2;
            acoustic_damping_error = (f2.t - con_int_error(1,6))/2;
            A_err= (f2.A - con_int_error(1,1))/2;
            displacement_reflectance_ratio_err= (f2.beta - con_int_error(1,3))/2;
            B_err= (f2.B - con_int_error(1,4))/2;
            acoustic_phase_err= (f2.p - con_int_error(1,5))/2;
            C_err= (f2.D - con_int_error(1,7))/2;

            %final fit (on constant provided) version of tau
            acoustic_damping_constant(3)=f2.t;

            if record_amp_ratio
                gamma=f2.A/f2.B;
                acoustic_damping_constant=[acoustic_damping_constant gamma];
            end

            %first checks that diffusivity has not pegged to fit bounds, second checks that beta
            %has not pegged. If either do, it is a bad fit
            bad_alpha=isnan(thermal_diffusivity_err(1));
            bad_beta=isnan(con_int_error(1,3));
            if bad_alpha && bad_beta
                display(strcat('Bad fit for: ',pos_file,'~re: tau (likely)'))
            elseif bad_alpha && ~bad_beta
                display(strcat('Bad fit for: ',pos_file,'~re: alpha'))
            elseif ~bad_alpha && bad_beta
                display(strcat('Bad fit for: ',pos_file,'~re: beta'))
            end

            if plot_final
                %Plotting factor for generation of traces in Figure 6
                amp_factor=1;
                %%%%%%%%%%%%%%%
                %Block to reconstruct the best-fit model without the sinusoidal
                %contribution, for comparison
                if tau_final_fit
                    f_remove_sine=cfit(TYPE2,f2.A,f2.k,f2.beta,f2.B,f2.p,0,f2.D,q,start_time2,peak_freq);
                elseif fixed_tau
                    f_remove_sine=cfit(TYPE2,f2.A,f2.k,f2.beta,f2.B,f2.p,f2.D,q,start_time2,peak_freq,0);
                end
                %%%%%%%%%%%%%%%
                
                %This is where the final fit composite figure is made

                figure()
                plot((neg(:,1)-time_naught)*10^9,(pos(:,2)-neg(:,2)-pre_signal_average)*10^3/amp_factor,'-','Color','#464646','LineWidth',2.5,'DisplayName','Raw TGS trace')
                hold on
%               %%%%%%plot vertical line at start time
                plot([total_signal(start_index2,1) total_signal(start_index2,1)]*10^9,ylim,'--','Color','#71C371','LineWidth',2.5,'DisplayName','Start index')
                hold on
                plot([neg(pump_time_index,1)-time_naught neg(pump_time_index,1)-time_naught]*10^9,ylim,'--','Color','#7DD3DA','LineWidth',2.5,'DisplayName','Time index')
                hold on
                plot(total_signal(start_index2:end,1)*10^9,(f2(total_signal(start_index2:end,1)))*10^3/amp_factor,'--','Color','#C43838','LineWidth',2.5,'DisplayName','Full functional fit')
                hold on
                plot(total_signal(start_index2:end,1)*10^9,(f_remove_sine(total_signal(start_index2:end,1)))*10^3/amp_factor,'-','Color','#3E59DC','LineWidth',2.5,'DisplayName','Thermal fit')
                if printFitsToFiles
                    printFile1 = 'DataOfRaw.txt';
                    fid1 = fopen( printFile1, 'w' );
                    fprintf(fid1, 'time[ns] Raw_signal[mV]\n');
                    for i=1:length(neg(:,1))
                        fprintf(fid1, '%6f %6f\n',(neg(i,1)-time_naught)*10^9, (pos(i,2)-neg(i,2)-pre_signal_average)*10^3/amp_factor);
                    end
                    fclose(fid1);
                    
                    %this second print only goes part of the way through the file and then throws an error but it will do!
                    printFile2 = 'DataOfFits.txt';
                    fid2 = fopen( printFile2, 'w' );
                    thermalFunc = (f_remove_sine(total_signal(start_index2:end,1)))*10^3/amp_factor;
                    fprintf(fid2, 'fit_time[ns] Functional_fit[mV] Thermal_fit[mV]\n');
                    for i=start_index2:length(total_signal(:,1))
                        fprintf(fid2, '%6f %6f %6f\n', total_signal(i,1)*10^9, (f2(total_signal(i,1)))*10^3/amp_factor, thermalFunc(i+1-start_index2,1));
                    end
                    fclose(fid2);
                end

                hold on
                xlim([-5 end_time*10^9])
                set(gcf,'Position',[0 0 1520 880])
                hold on
                set(gca,...
                    'FontUnits','points',...
                    'FontWeight','normal',...
                    'FontSize',30,...
                    'FontName','Times',...
                    'LineWidth',1.5)
                ylabel({'Signal amplitude [mV]'},...
                    'FontUnits','points',...
                    'FontSize',40,...
                    'FontName','Times')
                xlabel({'Time [ns]'},...
                    'FontUnits','points',...
                    'FontSize',40,...
                    'FontName','Times')
                legend('Location','southeast')

                % Display the already-made FFT as an inset image
                if ~steel
                    axes('pos',[.48 .48 .42 .42])
                    imshow('TGS_FFT.png')
                    saveas(gcf,strcat(pos_file,"_TGS_Final_Fit.png"))
                end
            end
        end

    %This will only execute if amp grating is set to 0, then it just fits
    %with the exponential and no SAWs.
    else
        total_signal=[pos(pump_time_index:end_index,1)-time_naught pos(pump_time_index:end_index,2)-neg(pump_time_index:end_index,2)-pre_signal_average];
        LBamp=[0 0 0];
        UBamp=[1 5*10^-3 1];
        STamp=[.1 5*10^-5 .1];
        start_time=pump_time_index*scope_timebase;

        OPSamp=fitoptions('Method','NonLinearLeastSquares','Lower',LBamp,'Upper',UBamp,'Start',STamp);
        TYPEamp=fittype('(A./sqrt((x+start_time)))*exp(-q^2*k*(x+start_time))+D','options',OPSamp,'problem',{'q','start_time'},'coefficients',{'A','k','D'});
        [famp,gof]=fit(total_signal(pump_time_index:end,1),total_signal(pump_time_index:end,2),TYPEamp,'problem',{q,start_time});

        thermal_diffusivity=famp.k;
        amp_factor=1;

        if plot_final
            figure()
            plot((neg(:,1)-time_naught)*10^9,(pos(:,2)-neg(:,2)-pre_signal_average)*10^3/amp_factor,'k-','LineWidth',4,'DisplayName','Raw TGS trace')
            hold on
            plot(total_signal(pump_time_index:end,1)*10^9,(famp(total_signal(pump_time_index:end,1)))*10^3/amp_factor,'r--','LineWidth',4,'DisplayName','Full functional fit')
            hold on
            xline(total_signal(pump_time_index,1)*10^9,'-b','LineWidth',4,'DisplayName','time index');
            xlim([-5 end_time*10^9])
            set(gcf,'Position',[0 0 960 540])
            annotation('textbox',[0.02 0.01 0.5 0.03],'String',strcar('Time Index = ',num2str(pump_time_index)),'FontSize',25,'FontName','Times','FontWeight','bold','LineStyle','none')
            hold on
            set(gca,...
                'FontUnits','points',...
                'FontWeight','normal',...
                'FontSize',30,...
                'FontName','Times',...
                'LineWidth',1.5)
            ylabel({'Signal amplitude [mV]'},...
                'FontUnits','points',...
                'FontSize',40,...
                'FontName','Times')
            xlabel({'Time [ns]'},...
                'FontUnits','points',...
                'FontSize',40,...
                'FontName','Times')
            legend('Location','northeast')
            saveas(gcf,strcat(pos_file,"_TGS_Final_Fit.png"))
        end      
        if verbose
                display(famp)
        end
    end
end



