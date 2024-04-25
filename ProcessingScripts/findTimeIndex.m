function [pump_time_index,end_time] = findTimeIndex(pos,neg)
%%%%
% Finds the pump_time_index based on the index at which the 2nd derivative of
% fixed_short is a maximum. This is set to work for 20ns and 50ns a
% division, but it will give wrong values for a different subdivision. This
% is a revised version of findTimeIndex that doesn't read the .txt, but
% instead takes the signal matrixes to save opening a bunch of extra files.

plot_ti=0;   % Show the peaks of raw trace and time index point
plot_derv=0; % Show second derivative and peak fits

%This is the offset applied to the location of the second derivative to
%find the pump_index. Not sure why this is necessary.
off_20ns=23;
off_50ns=19;

%%%%%%%%% Make the total signal
total_signal=[pos(:,1) pos(:,2)-neg(:,2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%% Methods for finding time_index %%%%%%%%%%%%%%%%
%%%%%%%%% Find Time Index from Trace MAX
[max_trace,max_time]=max(total_signal(1:600,2));


%%%%%%%%% Find Time Index from Second Derivative
first = gradient(total_signal(:,2)) ;   % derivative of total_signal
second_derv = gradient(first);         % second derivative of total_signal
[~,derv_index] = max(second_derv(1:max_time));
minprom=5*max(second_derv(1:50));
[pos_val,pos_loc]=findpeaks(second_derv(1:max_time),'MinPeakProminence',minprom);

%%%%%%%%% Make sure the peak it is using is the first after the rise of the
%%%%%%%%% raw trace
if length(pos_loc)>0
    if pos_loc(1)<derv_index
        derv_index=pos_loc(1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Assign pump_time_index and end time values %%%%%%%%%%%%%%%%
%%%%%%%%% Determine indices for different methods and assign end times
time_len=length(pos(:,1))/1000;
if time_len < 5
    end_time=2e-7; %for 20ns base on scope
    derv_index = derv_index-off_20ns;
else
    end_time=5e-7; %for 50ns base on scope
    derv_index = derv_index-off_50ns;
end
pump_time_index = derv_index;

%%%%%%%%% Plot
if plot_ti
    figure()
    plot(total_signal(:,2),'k','LineWidth',3)
    xlim([0 max_time+50])
    ylim([-0.005 max_trace+0.005])
    hold on
    plot([derv_index derv_index],[-0.005 max_trace+0.005],'r--','LineWidth',3)
    hold on
    title('Raw Trace with Calculated Time Index');
    ylabel({'Signal Amplitude [mV]'})
    xlabel({'Time Index'})

    txt3 = {['Time Index = ',num2str(pump_time_index)],['End Time = ',num2str(end_time)]};
    text(10,0.9*max_trace,txt3,'FontSize',20)

    set(gca,...
        'FontUnits','points',...
        'FontWeight','normal',...
        'FontSize',24,...
        'FontName','Times',...
        'LineWidth',3)
    ylabel({'Raw Trace'},...
        'FontUnits','points',...
        'FontSize',24,...
        'FontName','Times')
    xlabel({'Time Index'},...
        'FontUnits','points',...
        'FontSize',24,...
        'FontName','Times')

    if plot_derv
        yyaxis right
        hold on
        plot(second_derv*10^4,'b','LineWidth',3)
        hold on
        plot(pos_loc,pos_val*10^4,'o','LineWidth',3)
        ylabel({'2^{nd} Derivative'},...
            'FontUnits','points',...
            'FontSize',24,...
            'FontName','Helvetica')
    end
end
end
