function out = combiRec(filepath, offset, makeplots, figID)

    % Offset should be in volts (e.g. 0.01)

    % Set default
    if nargin < 3
        makeplots = false;
    end
    
    % Specify settings (part 1)
    duration_stim_artifact = 0.003;   % default is 0.003 seconds
    R_s_final = 3e+06;                 % default is 3e+06 ohms (0 for 100% compensation)
    risetime_limits = [0.2, 0.5, 0.8]; % default is [0.2, 0.5, 0.8]
    set_mixed_peak_mean = 30;          % default is 30
    set_nmda_peak_mean  = 300;          % default is 90
    frac_decay_time = 0.5;             % default is 0.5
    ampa_fwhm_scale_factor = 17;       % default is 10
    nmda_fwhm_scale_factor = 5;       % default is 30 
    ampa_delay = -0.000;                % default is 0 seconds
    nmda_delay = -0.002;                % default is 0 seconds
    
    % Preallocate output structure
    out = {};
    out.filepath = filepath;
    
    % Get file and path information
    pathref = pwd;
    [dirpath,fname,ext] = fileparts(filepath);
    filename = cat(2,fname,ext);                % e.g. "1.phy"
    dirstruct = split(dirpath,filesep);
    protocol = dirstruct{end}(1:end-4);         % e.g. "dual_mixed_eEPSC"
    rootdir = join(dirstruct(1:end-1),filesep); % e.g. "<path>/pair_000/"
    rootdir = rootdir{1};
    
    % Load data from channel
    chdir(rootdir)
    count = 0;
    iwaves = [];
    notes = '';
    holding = [];
    while true
        wavename = cat(2,protocol,'_00',num2str(count));
        if isfolder(wavename)
            chdir(wavename);
            try
                S = ephysIO(filename);
            catch
                chdir(pathref);
                break
            end
            if count == 0
                time = S.array(:,1);
                sampInt = S.xdiff;
            end
            iwaves{count+1} = S.array(:,2);
            if strcmpi(filepath(end-3:end),'.phy')
                try 
                    temp = split(S.notes{10},':');
                catch
                    sprintf(char(join(S.notes,'\n')));
                    chdir(pathref);
                    break               
                end
                holding(count+1) = eval(temp{end});
            elseif strcmpi(filepath(end-2:end),'.h5')
                temp = [-0.09, -0.07, -0.05, -0.03, -0.01, 0.01, 0.03];
                holding(count+1) = temp(count+1);
            end
            count = count+1;
            chdir('..')
        else
            break
        end
    end
    iwaves = cell2mat(iwaves);
    sz = size(iwaves);
    
    % waves sorted in ascending order of holding potential
    [holding, I] = sort(holding);
    iwaves = iwaves(:,I);
    
    % do not use incomplete recordings
    if count < 7
        chdir(pathref);
        return
    end
    if count > 7
        count = 7;
        iwaves(:,8:end)=[];
        holding(8:end)=[];
        sz(2) = count;
    end

    % Correct holding potential for LJP (offset)
    holding = holding - offset;
    ref = [-0.1,-0.08,-0.06,-0.04,-0.02,0.00,0.02];
    if sum(diff([holding;ref])) > 0
        chdir(pathref);
        return
    end
    
    % Specify settings (part 2)
    % Set baseline cursor positions 
    b1 = dsearchn(time,0.095);  % default is 0.095 seconds (baseline start)
    b2 = dsearchn(time,0.099);  % default is 0.099 seconds (baseline end)
    % Set peak cursors for mixed EPSC
    p1 = dsearchn(time,0.1002+duration_stim_artifact); % default is 0.1034 seconds (peak start; mixed)
    p2 = dsearchn(time,0.150);  % default is 0.160 seconds (peak end; mixed)
    p3 = dsearchn(time,0.103);  % default is 0.105 seconds (peak start; nmda)
    p4 = dsearchn(time,0.160);  % default is 0.160 seconds (peak end; nmda)
    
    % measure and record pre-stimulus baseline
    raw_current_baseline = mean(iwaves(b1:b2,:))';
    
    % subtract pre-stimulus baseline
    iwaves = iwaves - (ones(sz(1),1) * raw_current_baseline');
 
    % Interpolate stimulus artifact
    for i=1:count
      iwaves(2500:2500+duration_stim_artifact/sampInt,i) = interp1q(time([2500;2500+duration_stim_artifact/sampInt]),iwaves([2500;2500+duration_stim_artifact/sampInt],i),time([2500:2500+duration_stim_artifact/sampInt]'));
    end
    raw_iwaves = iwaves;
    
    % Perform raw (uncompensated and mixed) current peak measurements
    [raw_current_peaks, raw_peak_idx] = measure_peaks (iwaves, p1, p2, set_mixed_peak_mean,'both');

    % Estimate reversal potential by linear regression 
    % (through points at -20, 0 and +2) and find the root
    p = polyfit(holding(end-2:end)',1e+012*raw_current_peaks(end-2:end),1);
    Erev = roots(p);
    out.Erev = Erev;
    
    %if makeplots 
    %    subplot(2,2,2);
    %    plot(holding,raw_current_peaks,'ko-','MarkerFaceColor','k');
    %    title('Mixed I-V curve');xlabel('Voltage (V)');ylabel('Current (A)');
    %    h = gca; hold on; plot([0 0],h.YLim,'k-'); plot(h.XLim,[0 0],'k-'); 
    %    try
    %        plot(Erev,0,'r+','MarkerSize',10,'LineWidth',2)
    %    catch
    %    end
    %    hold off; set(gca,'YLim',h.YLim);set(gca,'XLim',h.XLim);
    %end 
    
    % measure and record pre-stimulus baseline
    mixed_current_baseline = mean(iwaves(b1:b2,:))';
    
    % subtract pre-stimulus baseline
    iwaves = iwaves - (ones(sz(1),1) * mixed_current_baseline');
    
    % Perform mixed current peak measurements 
    mixed_current_peaks = measure_peaks (iwaves, p1, p2, set_mixed_peak_mean,'both');
    
    % Carry out conversions
    raw_gwaves = raw_iwaves ./ (ones(sz(1),1) * holding - Erev);
    gwaves = iwaves ./ (ones(sz(1),1) * holding - Erev);
    gwaves(:,6) = NaN; % very noisy conductance estimates at/close to reversal potential
    gnmda  = gwaves - gwaves(:,1) * ones(1,sz(2));
    inmda = gnmda .* (ones(sz(1),1) * holding - Erev);
    
    % Perform mixed conductance peak measurements 
    [mixed_conductance_peaks, mixed_peak_idx, mixed_conductance_smoothed] = measure_peaks (gwaves, p1, p2, set_mixed_peak_mean,'positive');
    
    if numel(mixed_peak_idx) ~= numel(mixed_conductance_peaks)
        chdir(pathref);
        return
    end
    
    % Perform NMDA-only conductance slope measurements 
    [nmda_conductance_peaks, nmda_peak_idx, nmda_conductance_smoothed] = measure_peaks (gnmda, p3, p4, set_nmda_peak_mean,'positive');
    [nmda_risetime, nmda_rise_times] = measure_risetime (time, gnmda(:,end), p1, nmda_peak_idx(end), risetime_limits * nmda_conductance_peaks(end));

    % Measure and plot rising slope of nmda currents for more robust IV
    % customise risetime points to tweak cursors for slop measurements
    % ampa_rise_times = [0.1038 0.107 0.1048];
    % nmda_rise_times = [0.1038 0.107 0.1052];    
    %slope_start = dsearchn(time,linspace(ampa_rise_times(1),ampa_rise_times(1),7)');
    %slope_end = dsearchn(time,linspace(ampa_rise_times(3),nmda_rise_times(3),7)');
    slope_start = dsearchn(time,ones(7,1)*nmda_rise_times(1));
    slope_end = dsearchn(time,ones(7,1)*nmda_rise_times(3));
    linear_fits = zeros(7,2);
    for i = 1:7
        linear_fits(i,:) = polyfit(time(slope_start(i):slope_end(i)),inmda(slope_start(i):slope_end(i),i),1);
    end
    out.nmda_current_slopes_pA_ms = linear_fits(:,1)' * 1e+9;
    linear_fits(end-1,:) = 0;
    
    % Plot current slope as tangents on EPSC rise
    subplot(1,2,1);
    plot(time,inmda,'-','Color',[0.5,0.5,0.5]); xlim([0.08,0.2]); ylim('auto');title('NMDA currents')
    hold on
    for i = 1:7
       plot(time(slope_start(i):slope_end(i)),polyval(linear_fits(i,:),time(slope_start(i):slope_end(i))),'-r','LineWidth',1);
    end
    hold off
    subplot(1,2,2);
    plot(holding,linear_fits(:,1)','ko-','MarkerFaceColor','k');

    chdir(pathref);

    % Save figure
    if nargin > 3
        h1 = gcf;
        sgtitle(['Recording ' num2str(figID)]);
        %fig_name = cat(2,'./IV_img/',num2str(figID),'.fig');
        %savefig(h1,fig_name,'compact');
        png_name = cat(2,'./IV_img/',num2str(figID),'.png');
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