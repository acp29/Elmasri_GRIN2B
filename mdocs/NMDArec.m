function out = NMDArec(filepath, offset, makeplots, figID)

    % Offset should be in volts (e.g. 0.01)

    % Set default
    if nargin < 3
        makeplots = false;
    end
    
    % Specify settings (part 1)
    duration_stim_artifact = 0.003;    % default is 0.003 seconds
    R_s_final = 3e+06;                 % default is 3e+06 ohms (0 for 100% compensation)
    risetime_limits = [0.2, 0.5, 0.8]; % default is [0.2, 0.5, 0.8]
    set_nmda_peak_mean = 90;          % default is 90
    frac_decay_time = 0.5;             % default is 0.5
    nmda_fwhm_scale_factor = 30;       % default is 30 
    nmda_delay = -0.000;                % default is 0 seconds
    
    % Preallocate output structure
    out = {};
    out.filepath = filepath;
    out.decay_fraction = 0.5;
    pathref = pwd;
    
    % Get file and path information
    notes = '';
    holding = [];
    S = ephysIO(filepath);
    try
        temp = split(S.notes{10},':');
        sprintf(char(join(S.notes,'\n')));
        holding = eval(temp{end});
        if holding ~= 0.03
            return
        end
    catch
        holding = 0.03;
    end
    time = S.array(:,1);
    iwaves = S.array(:,2);
    sampInt = S.xdiff;

    % Correct holding potential for LJP (offset)
    holding = holding - offset;
    ref = [0.02];
    sz = 1;
    
    % Specify settings (part 2)
    % Set baseline cursor positions 
    b1 = dsearchn(time,0.095);  % default is 0.095 seconds (baseline start)
    b2 = dsearchn(time,0.099);  % default is 0.099 seconds (baseline end)
    % Set peak cursors for mixed EPSC
    p1 = dsearchn(time,0.1002+duration_stim_artifact-0.004); % default is 0.1034 seconds (peak start; mixed)
    p3 = dsearchn(time,0.105);  % default is 0.105 seconds (peak start; nmda)
    p4 = dsearchn(time,0.160);  % default is 0.160 seconds (peak end; nmda)
    
    % measure and record pre-stimulus baseline
    raw_current_baseline = mean(iwaves(b1:b2,:))';
    out.raw_current_baseline_pA = raw_current_baseline' * 1e12;
    out.raw_current_peaks_pA = [];
    out.raw_current_peaks_reversal_potential_mV = [];
    
    % subtract pre-stimulus baseline
    iwaves = iwaves - (ones(sz(1),1) * raw_current_baseline');
 
    % Interpolate stimulus artifact
    for i=1:sz
      iwaves(2500:2500+duration_stim_artifact/sampInt,i) = interp1q(time([2500;2500+duration_stim_artifact/sampInt]),iwaves([2500;2500+duration_stim_artifact/sampInt],i),time([2500:2500+duration_stim_artifact/sampInt]'));
    end
    raw_iwaves = iwaves;
      
    % Perform offline series resistance compensation
    [R_s, C_m, R_c] = wcp (time, iwaves(:,i), -0.005, 0.01, 0.02, sampInt);
    fraction = max(1 - (R_s_final / R_s), 0);
    Erev = -0.003; % determined from untransfected neurons
    % Method 2
    if fraction > 0
        iwaves(:,i) = rscomp (iwaves(:,i), sampInt, R_s_final, R_s, C_m, holding(i), Erev);
    end
    out.raw_series_resistance_Mohm = R_s * 1e-06;
    out.raw_cell_resistance_Mohm = R_c * 1e-06;
    out.raw_cell_capacitance_pF = C_m  * 1e+12;
    
    out.ampa_peak_ampl_nS = [];
    out.ampa_peak_current_pA = [];
    out.ampa_risetime_ms = [];
    out.ampa_decaytime_ms = [];
    out.ampa_fwhm_ms = [];
    out.ampa_decay_fit = [];
    out.ampa_charge_transfer_pC = [];
    
    % Carry out conversions
    gnmda = iwaves ./ (ones(sz(1),1) * holding - Erev);
    raw_gwaves = raw_iwaves ./ (ones(sz(1),1) * holding - Erev);
    
    % Perform NMDA-only conductance peak measurements 
    [nmda_conductance_peaks, nmda_peak_idx, nmda_conductance_smoothed] = measure_peaks (gnmda, p3, p4, set_nmda_peak_mean,'positive');
    
    % Calculate NMDA-EPSC peak
    out.nmda_peak_ampl_nS = nmda_conductance_peaks(end) * 1e+09;
    out.nmda_peak_ampl_pA = nmda_conductance_peaks(end) * holding(end) * 1e+12;
    nmda_peak_time = time(nmda_peak_idx(end));

    % Calculate NMDA-EPSC risetime
    [nmda_risetime, nmda_rise_times] = measure_risetime (time, nmda_conductance_smoothed(:,end), p3, nmda_peak_idx(end), risetime_limits * nmda_conductance_peaks(end));
    out.nmda_risetime_ms = nmda_risetime * 1e+03;
    
    % Calculate NMDA-EPSC half-decay time
    [nmda_decaytime, nmday_decay_time] = measure_decaytime(time, nmda_conductance_smoothed(:,end), nmda_peak_idx(end), nmda_conductance_peaks(end), frac_decay_time);
    out.nmda_decaytime_ms = nmda_decaytime * 1e+03;

    % Calculate NMDA-EPSC full-width half-maximum
    nmda_fwhm = nmda_peak_time + nmda_decaytime - nmda_rise_times(2);
    out.nmda_fwhm_ms = nmda_fwhm * 1e+03;
     
    out.ampa_nmda_ratio = [];
    
    % Measure decay of NMDA-EPSC
    % time of peak
    nmda_fit_start_time = nmda_peak_time + nmda_delay;
    nmda_fit_end_time = min(nmda_peak_time + nmda_fwhm * nmda_fwhm_scale_factor, 1.1);  
    [curve, out.nmda_decay_fit, nmda_indiv_curves, nmda_time] = fit_decay (time, gnmda(:,end),nmda_fit_start_time,nmda_fit_end_time);
    
    % Measure charge transfer of NMDA-EPSC
    f2 = dsearchn(time, nmda_fit_end_time);
    out.nmda_charge_transfer_pC = abs(1e12 * trapz(time(p3:f2),iwaves(p3:f2,1)));
       
    % Plot NMDA-EPSC conductance
    if makeplots
        area(time(p3:f2),gnmda(p3:f2,end),0,'EdgeColor','none','FaceColor',[0.90 0.90 1.0]);
        hold on;
        nh_rs = plot(time,gnmda(:,end),'k'); xlim([0.05,0.6]); ylim('auto');
        nh = plot(time,raw_gwaves(:,1),'b-');
        %plot(time,nmda_conductance_smoothed(:,end),'Color',[0.75,0.75,0.75]);
        plot(nmda_peak_time,nmda_conductance_peaks(end),'co','MarkerSize',8,'MarkerFaceColor','c') % peak point
        plot(nmda_rise_times(1),min(risetime_limits) * nmda_conductance_peaks(end),'go','MarkerSize',8,'MarkerFaceColor','g') % lower rise point
        plot(nmda_rise_times(3),max(risetime_limits) * nmda_conductance_peaks(end),'go','MarkerSize',8,'MarkerFaceColor','g') % upper rise point
        plot(nmda_time,curve,'r-','LineWidth',3) % decay fit
        for i = 1:out.nmda_decay_fit.n
            plot(nmda_time,nmda_indiv_curves,'r--') % decay fit
        end
        plot(nmda_rise_times(2),frac_decay_time*nmda_conductance_peaks(end),'yo','MarkerSize',8,'MarkerFaceColor','y') % half-rise point
        plot(nmday_decay_time,frac_decay_time*nmda_conductance_peaks(end),'yo','MarkerSize',8,'MarkerFaceColor','y') % half-decay point
        xlabel('Time (s)');ylabel('Conductance (S)');
        hleg = [nh(1) nh_rs(1)]; legend(hleg,'Rs uncompensated','Rs compensated')
        hold off      
    end
    
    out.raw_current_slopes_pA_ms = [];
    out.raw_current_slopes_reversal_potential_mV = [];
    
    % Save figure
    if nargin > 3
        h1 = gcf;
        sgtitle(['Recording ' num2str(figID)]);
        fig_name = cat(2,'./nmda_img/',num2str(figID),'.fig');
        savefig(h1,fig_name,'compact');
        png_name = cat(2,'./nmda_img/',num2str(figID),'.png');
        print(h1,png_name,'-dpng','-r100','-opengl');
    end
    

end


function [peaks, idx, Yf] = measure_peaks (Y, p1, p2, set_peak_mean, sign_peaks)
    
    % Peak measurement function 
    
    % p1 is the index of the first cursor position
    % p2 is the index of the second cursor position
    % set_peak_mean sets the number of points to average over the peak
     
    % Apply boxcar filter (see example at
    % https://uk.mathworks.com/help/matlab/ref/filter.html)
    b = (1 / set_peak_mean) * ones(1, set_peak_mean);
    a = 1;
    Yf = filter(b, a, Y); % apply boxcar filter
    
    % Center the moving average
    Yf = circshift(Yf,fix(-set_peak_mean/2));
    
    switch sign_peaks
        case {'positive'}
            
            peaks = max (Yf (p1:p2, :)); 
            
        case {'negative'}
            
            peaks = min (Yf (p1:p2, :)); 
            
        case {'both'}

            % Calculate extreme positive and negative peaks between the cursors
            p_neg = min (Yf (p1:p2, :)); 
            p_pos = max (Yf (p1:p2, :));
            P = [p_neg; p_pos];
    
            % Which absolute value of the positive and negative peaks is bigger
            absmax = max(abs(P)); 
            logic = (abs(P) == absmax);
    
            % Accomodate for 0 or NaN
            ridx = or(all(logic),~any(logic));
            logic(1,ridx)=1;
            logic(2,ridx)=0;
    
            % Return column peak values
            peaks = P(logic);
    
    end
    
    % Return row number of peak values
    n = numel(peaks);
    i = nan(n,1);
    for i = 1:n
        if ~isnan(peaks(i)) && (peaks(i) ~= 0) && ~isinf(peaks(i))
          tmp = find(Yf(p1:p2,i) == peaks(i));  
          idx(i,1) = p1 - 1 + tmp(1);
        end
    end
    
end

function [risetime, times] = measure_risetime (t, y, p1, p2, limits)
    
    % Rise-time measurement function 
    times = zeros(1,3);
    min_lim = find(y(p1:p2,1) > min(limits));
    mid_lim = find(y(p1:p2,1) > limits(2));
    max_lim = find(y(p1:p2,1) > max(limits));
    times(1) = interp1([y(p1+min_lim(1)-2);y(p1-1+min_lim(1))],[t(p1+min_lim(1)-2);t(p1-1+min_lim(1))],min(limits),'linear','extrap');
    times(2) = interp1([y(p1+mid_lim(1)-2);y(p1-1+mid_lim(1))],[t(p1+mid_lim(1)-2);t(p1-1+mid_lim(1))],limits(2),'linear','extrap');
    times(3) = interp1([y(p1+max_lim(1)-2);y(p1-1+max_lim(1))],[t(p1+max_lim(1)-2);t(p1-1+max_lim(1))],max(limits),'linear','extrap');
    if times(1) > times(3)
        times = fliplr(times); % do this to prevent an error later but results will look weird; adjust cursor positions to fix
    end   
    risetime = times(3) - times(1);
end

function [decaytime, decay_time] = measure_decaytime (t, y, peak_idx, peak_ampl, frac_decay_time)
    
    % Decay-time measurement function 
    decay_idx = find(y(peak_idx:end,1) < frac_decay_time * peak_ampl);
    decay_time = interp1([y(peak_idx+decay_idx(1)-2);y(peak_idx-1+decay_idx(1))],[t(peak_idx+decay_idx(1)-2);t(peak_idx-1+decay_idx(1))],frac_decay_time * peak_ampl,'linear','extrap');
    decaytime = decay_time - t(peak_idx(1));
    
end

function [fit, stats, indiv_curves, time] = fit_decay (t, y, fit_start, fit_end)
    
    % Set fit cursors
    f1 = dsearchn(t, fit_start);
    f2 = dsearchn(t, fit_end);
    time = t(f1:f2);
  
    % Scale the data
    t = t * 1e+03; % time in ms
    y = y * 1e+09; % current in nS
    Tn = 100;
    
    % Decay measurement function using Chebyshev algorithm
    % 2 exponential components
    try
        [fit, stats] = chebexp(t(f1:f2), y(f1:f2), 2, Tn); 
        if ~any(isreal(stats.tau))
            error('invoke catch statement if there is a harmonic component')
        end     
        if any(stats.a < 0)
            error('invoke catch statement if amplitude of any component is negative')
        end
        if any(stats.tau < 0)
            error('invoke catch statement if time constant of any component is negative')
        end
        if any((stats.a/sum(stats.a)) < 0.05)
            error('invoke catch statement if amplitude of any component is less than 5%') 
        end
        stats.n = 2;
    catch
        try
            [fit, stats] = chebexp(t(f1:f2), y(f1:f2), 1);  
            stats.n = 1;
        catch
            fit = zeros(size(t(f1:f2),1),1);
            stats.a = NaN;
            stats.tau = NaN;
            stats.offset = NaN;
            stats.n = 1;
            stats.wtau = NaN;
            stats.perc = NaN;
        end
    end
    indiv_curves = zeros(size(t(f1:f2),1),stats.n);
    for i = 1:stats.n
        indiv_curves(:,i) = stats.a(i) * exp(-(t(f1:f2)-t(f1))/stats.tau(i)) + stats.offset;
    end
    stats.wtau = sum(stats.a / sum(stats.a) .* stats.tau); % weighted tau
    stats.perc = 100 * stats.a / sum(stats.a);             % amplitude as percantage
    
    % Rescale output 
    fit = fit * 1e-09;
    indiv_curves = indiv_curves * 1e-09;
    
end
