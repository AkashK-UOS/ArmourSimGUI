function armourSimGUI()
% type 'armourSimGUI' into the command window to open the GUI, configure
% settings and then hit run simulation box to start simulation. 

    %% Define Material Library 
    materials = { ...
        struct('name','Boron Carbide (B4C)',             'E',450e9, 'nu',0.17, 'rho',2520, 'role','Strike Face'), ...
        struct('name','Silicon Carbide (SiC)',           'E',410e9, 'nu',0.14, 'rho',3100, 'role','Intermediate'), ...
        struct('name','AR500 Steel',                     'E',210e9, 'nu',0.30, 'rho',7850, 'role','Backing'), ...
        struct('name','Alumina',                         'E',380e9, 'nu',0.22, 'rho',3900, 'role','Strike Face/Intermediate'), ...
        struct('name','Titanium Diboride (TiB2)',        'E',520e9, 'nu',0.20, 'rho',4500, 'role','Strike Face'), ...
        struct('name','Silicon Nitride (Si3N4)',         'E',310e9, 'nu',0.27, 'rho',3200, 'role','Intermediate'), ...
        struct('name','Zirconia Toughened Alumina (ZTA)','E',380e9,'nu',0.23,  'rho',3800, 'role','Intermediate'), ...
        struct('name','Tungsten Carbide (WC-Co)',        'E',600e9, 'nu',0.29, 'rho',15000,'role','Strike Face'), ...
        struct('name','Alumina/SiC Composite',           'E',400e9, 'nu',0.20, 'rho',3500, 'role','Intermediate'), ...
        struct('name','Kevlar Composite',                'E',70e9,  'nu',0.36, 'rho',1440, 'role','Backing'), ...
        struct('name','Dyneema Composite',               'E',50e9,  'nu',0.30, 'rho',970,  'role','Backing') ...
    };

    %% Define Projectile Library
    projectiles = { ...
        struct('name','7.62 mm NATO', 'E',210e9, 'nu',0.30, 'rho',7850, 'vel',854,  'radius',0.00381), ...
        struct('name','5.56 mm NATO', 'E',200e9, 'nu',0.30, 'rho',7850, 'vel',945,  'radius',0.00278), ...
        struct('name','.50 BMG',      'E',250e9, 'nu',0.30, 'rho',7850, 'vel',928,  'radius',0.00635) ...
    };
    projNames = cellfun(@(p)p.name, projectiles, 'UniformOutput', false);

    %% Create GUI Figure
    figure('Name','Armour Material & Simulation Parameter Selection', ...
               'Position',[100 100 800 700]);

    %% Material selection controls
    matNames = cellfun(@(m)[m.name ' (' m.role ')'], materials, 'UniformOutput', false);
    uicontrol('Style','text','Position',[20 650 160 20], 'String','Layer 1 (Strike Face) [mm]:','FontWeight','bold');
    popup1 = uicontrol('Style','popupmenu','Position',[190 650 300 20],'String',matNames);
    edit1  = uicontrol('Style','edit','Position',[500 650 60 20],'String','0');

    uicontrol('Style','text','Position',[20 620 160 20], 'String','Layer 2 (Intermediate) [mm]:','FontWeight','bold');
    popup2 = uicontrol('Style','popupmenu','Position',[190 620 300 20],'String',matNames);
    edit2  = uicontrol('Style','edit','Position',[500 620 60 20],'String','0');

    uicontrol('Style','text','Position',[20 590 160 20], 'String','Layer 3 (Backing) [mm]:','FontWeight','bold');
    popup3 = uicontrol('Style','popupmenu','Position',[190 590 300 20],'String',matNames);
    edit3  = uicontrol('Style','edit','Position',[500 590 60 20],'String','0');

    %% Projectile selection
    uicontrol('Style','text','Position',[20 560 160 20], 'String','Projectile:','FontWeight','bold');
    popupProj = uicontrol('Style','popupmenu','Position',[190 560 300 20],'String',projNames,'Value',1);

    %% Simulation Parameter Controls
    uicontrol('Style','text','Position',[20 530 300 20], 'String','Simulation Parameters:','FontWeight','bold');
    uicontrol('Style','text','Position',[20 500 150 20], 'String','Dispersion Factor (k):');
    edit_dispersion = uicontrol('Style','edit','Position',[180 500 100 20],'String','0.5');
    uicontrol('Style','text','Position',[20 470 150 20], 'String','Damping Coefficient (alpha):');
    edit_damping    = uicontrol('Style','edit','Position',[180 470 100 20],'String','2');
    uicontrol('Style','text','Position',[20 440 150 20], 'String','Time Pad:');
    edit_timepad    = uicontrol('Style','edit','Position',[180 440 100 20],'String','12');
    uicontrol('Style','text','Position',[20 410 150 20], 'String','Number of Points:');
    edit_npoints    = uicontrol('Style','edit','Position',[180 410 100 20],'String','2000');
    uicontrol('Style','text','Position',[20 380 150 20], 'String','Min Amplitude:');
    edit_minAmp     = uicontrol('Style','edit','Position',[180 380 100 20],'String','1e-14');
    uicontrol('Style','text','Position',[320 500 100 20],'String','Verbose:','FontWeight','bold');
    check_verbose   = uicontrol('Style','checkbox','Position',[420 500 20 20],'Value',0);

    %% Adaptive WaveEvents
    uicontrol('Style','text','Position',[20 350 160 20], 'String','Adaptive WaveEvents:','FontWeight','bold');
    check_adaptive  = uicontrol('Style','checkbox','Position',[190 350 20 20],'Value',0);
    uicontrol('Style','text','Position',[20 320 220 20], 'String','Fixed Cap (if unchecked):');
    edit_fixedCap   = uicontrol('Style','edit','Position',[240 320 80 20],'String','500000');

    %% Damage Controls
    uicontrol('Style','text','Position',[20 280 150 20], 'String','Damage Value:','FontWeight','bold');
    edit_damageVal  = uicontrol('Style','edit','Position',[180 280 100 20],'String','0.8');
    uicontrol('Style','text','Position',[300 280 150 20], 'String','Damage Layer:','FontWeight','bold');
    popup_damageLayer = uicontrol('Style','popupmenu','Position',[450 280 100 20],...
                                  'String',{'Layer 1','Layer 2','Layer 3'});

    %% Sensor positions (x,y) for 4 virtual sensors 
    defaultSensors = [ ...
        -0.10, -0.10; ...
         0.20, -0.05; ...
         0.15,  0.15; ...
        -0.05,  0.25  ...
    ];
    sensorEdits = gobjects(4,2);
    yStart      = 240;   
    ySpacing    = 30;    
    for sensorIdx = 1:4
        yPos = yStart - (sensorIdx-1)*ySpacing;
        uicontrol('Style','text','Position',[20  yPos 100 20], ...
                  'String',sprintf('Sensor %d (x,y):',sensorIdx));
        sensorEdits(sensorIdx,1) = uicontrol('Style','edit','Position',[130 yPos 60 20], ...
                  'String',num2str(defaultSensors(sensorIdx,1)));
        sensorEdits(sensorIdx,2) = uicontrol('Style','edit','Position',[200 yPos 60 20], ...
                  'String',num2str(defaultSensors(sensorIdx,2)));
    end

    %% Batch Analysis
    uicontrol('Style','text','Position',[320 350 150 20], 'String','Batch Analysis:','FontWeight','bold');
    handles.check_batch     = uicontrol('Style','checkbox','Position',[480 350 20 20],'Value',0);
    uicontrol('Style','text','Position',[580 500 120 20],'String','Param to Vary:','FontWeight','bold');
    handles.paramPopup      = uicontrol('Style','popupmenu','Position',[580 480 140 20],...
                          'String',{'Damping','Dispersion','Time Pad','N Points','Min Amp','WaveEventsMax'},...
                          'Value',1);
    uicontrol('Style','text','Position',[580 450 70 20], 'String','Start:');
    handles.edit_rangeStart = uicontrol('Style','edit','Position',[650 450 70 20],'String','0.5');
    uicontrol('Style','text','Position',[580 420 70 20], 'String','End:');
    handles.edit_rangeEnd   = uicontrol('Style','edit','Position',[650 420 70 20],'String','5');
    uicontrol('Style','text','Position',[580 390 70 20], 'String','#Steps:');
    handles.edit_rangeSteps = uicontrol('Style','edit','Position',[650 390 70 20],'String','5');
    guidata(gcf, handles);

    %% Run Simulation Button
    uicontrol('Style','pushbutton','Position',[300 20 150 40], ...
              'String','Run Simulation','Callback',@runSimulation);

    %% runSimulation callback
    function runSimulation(~,~)
        simStartTotal = tic;

        % Read materials & thicknesses
        sel1 = get(popup1,'Value'); sel2 = get(popup2,'Value'); sel3 = get(popup3,'Value');
        t1   = str2double(get(edit1,'String'));
        t2   = str2double(get(edit2,'String'));
        t3   = str2double(get(edit3,'String'));
        
        % Scale thickness from mm to m
        t1 = t1 * 1e-3;
        t2 = t2 * 1e-3;
        t3 = t3 * 1e-3;

        % Read simulation parameters
        k_factor   = str2double(get(edit_dispersion,'String'));
        alpha_damp = str2double(get(edit_damping,'String'));
        time_pad   = str2double(get(edit_timepad,'String'));
        N_points   = str2double(get(edit_npoints,'String'));
        minAmp     = str2double(get(edit_minAmp,'String'));
        verbose    = get(check_verbose,'Value');

        % Adaptive vs fixed waveEvents
        if get(check_adaptive,'Value')
            fixedCap = [];
        else
            fixedCap = str2double(get(edit_fixedCap,'String'));
        end

        % Damage parameters
        dVal        = str2double(get(edit_damageVal,'String'));
        damageLayer = get(popup_damageLayer,'Value');
        doBatch     = handles.check_batch.Value;

        % Read selected projectile
        selProj = get(popupProj,'Value');
        config.projectile = projectiles{selProj};

        % Read sensor positions
        config.sensorPositions = zeros(4,2);
        for i = 1:4
            x = str2double(get(sensorEdits(i,1),'String'));
            y = str2double(get(sensorEdits(i,2),'String'));
            config.sensorPositions(i,:) = [x, y];
        end

        % Build config struct
        config.E_layers_orig      = [materials{sel1}.E, materials{sel2}.E, materials{sel3}.E];
        config.nu_layers          = [materials{sel1}.nu, materials{sel2}.nu, materials{sel3}.nu];
        config.rho_layers_orig    = [materials{sel1}.rho, materials{sel2}.rho, materials{sel3}.rho];
        config.thick              = [t1, t2, t3];
        config.numLayers          = 3;
        config.k_factor           = k_factor;
        config.alpha_damp         = alpha_damp;
        config.time_pad           = time_pad;
        config.N_points           = N_points;
        if isempty(fixedCap)
            config.waveEventsMax = Inf;
        else
            config.waveEventsMax = fixedCap;
        end
        config.minAmp             = minAmp;
        config.verbose            = verbose;
        config.adaptiveWaveEvents = get(check_adaptive,'Value');
        config.damageVal          = dVal;
        config.damageLayer        = damageLayer;

        % Custom display of configuration 
        fprintf('Selected Configuration:\n');

        % projectile name
        fprintf('  %-18s : %s\n', 'projectile', config.projectile.name);
    
        % sensorPositions matrix
        fprintf('  %-18s :\n', 'sensorPositions');
        for r = 1:size(config.sensorPositions,1)
            fprintf('    [%8.4f  %8.4f]\n', config.sensorPositions(r,:));
        end

        % the rest
        fprintf('  %-18s : [%s]\n', 'E_layers_orig',      num2str(config.E_layers_orig,    ' %.4e'));
        fprintf('  %-18s : [%s]\n', 'nu_layers',          num2str(config.nu_layers,        ' %.4f'));
        fprintf('  %-18s : [%s]\n', 'rho_layers_orig',    num2str(config.rho_layers_orig,  ' %.0f'));
        fprintf('  %-18s : [%s]\n', 'thick',              num2str(config.thick,            ' %.4f'));
        fprintf('  %-18s : %d\n',   'numLayers',          config.numLayers);
        fprintf('  %-18s : %.4f\n', 'k_factor',           config.k_factor);
        fprintf('  %-18s : %.4f\n', 'alpha_damp',         config.alpha_damp);
        fprintf('  %-18s : %.4f\n', 'time_pad',           config.time_pad);
        fprintf('  %-18s : %d\n',   'N_points',           config.N_points);
        fprintf('  %-18s : %g\n',   'waveEventsMax',      config.waveEventsMax);
        fprintf('  %-18s : %.4e\n', 'minAmp',             config.minAmp);
        fprintf('  %-18s : %d\n',   'verbose',            config.verbose);
        fprintf('  %-18s : %d\n',   'adaptiveWaveEvents', config.adaptiveWaveEvents);
        fprintf('  %-18s : %.4f\n', 'damageVal',          config.damageVal);
        fprintf('  %-18s : %d\n',   'damageLayer',        config.damageLayer);

        % Run simulation or batch
        if ~doBatch
        % Single‐run: discard any outputs
        dis(config);
        else
        % Batch analysis: collect regression‐DI curves and visualize sensitivity
        % Retrieve GUI batch controls
        paramList    = get(handles.paramPopup, 'String');
        selParam     = handles.paramPopup.Value;
        startVal     = str2double(get(handles.edit_rangeStart, 'String'));
        endVal       = str2double(get(handles.edit_rangeEnd,   'String'));
        nSteps       = str2double(get(handles.edit_rangeSteps, 'String'));
    
        fprintf('\n=== BATCH MODE ENABLED ===\n');
        fprintf('Parameter to vary : %s\n', paramList{selParam});
        fprintf('Range             : [%.4g → %.4g]\n', startVal, endVal);
        fprintf('Number of steps   : %d\n\n', nSteps);

        % Sweep the selected parameter and store DI vs damageVal
        damageLevels = 0:0.1:1;  
        all_DI       = zeros(nSteps, numel(damageLevels));
        
        % Choose linear or logarithmic spacing depending on parameter
        if selParam == 5  % Min Amp is the 5th entry 
            % Logarithmic sweep for amplitude thresholds
            sweepVals = logspace(log10(startVal), log10(endVal), nSteps);
        else
            % Linear sweep for all other parameters
            sweepVals = linspace(startVal, endVal, nSteps);
        end

        for i = 1:nSteps
            v = sweepVals(i);
            switch selParam
                case 1, config.alpha_damp    = v;
                case 2, config.k_factor      = v;
                case 3, config.time_pad      = v;
                case 4, config.N_points      = round(v);
                case 5, config.minAmp        = v;
                case 6, config.waveEventsMax = v;
            end
            % capture the regression DI curve returned by dis()
            fprintf('Step %2d/%d: %s = %.4g →', i, nSteps, paramList{selParam}, v);
            DI_fit       = dis(config);    
            all_DI(i, :) = DI_fit(:)';
            fprintf('done\n');
        end
    
        % Plot overlaid DI curves vs damageVal
        figure('Color','w','Position',[200 200 600 400]);
        hold on;        
        % choose a colormap matching the number of sweep steps
        cmap = parula(nSteps);        
        % plot each DI curve colored by its parameter value
        for i = 1:nSteps
            plot(damageLevels, all_DI(i, :), 'Color', cmap(i, :), 'LineWidth', 1.2);
        end        
        xlabel('True Damage Level (damageVal)');
        ylabel('Composite Damage Index');
        title(sprintf('Regression-DI vs Damage for varying %s', paramList{selParam}));        
        % apply the custom colormap
        colormap(cmap);
        % force the colour axis to span actual parameter range
        clim([startVal, endVal]);
        % relabel the colourbar to show real parameter values
        cb = colorbar;
        cb.Label.String = paramList{selParam};
        % place ticks at each sweep value
        cb.Ticks = sweepVals;
        % format the tick labels as strings
        cb.TickLabels = arrayfun(@(v) sprintf('%.2g', v), sweepVals, 'UniformOutput', false);
        grid on;
    
        % Compute and plot RMSE of DI curves vs. the ideal 45-degree line
        errors = zeros(nSteps,1);
        for i = 1:nSteps
            DI_fit    = all_DI(i, :);                     
            errors(i) = sqrt(mean((DI_fit - damageLevels).^2));
        end
        
        figure('Color','w','Position',[200 200 600 400]);
        plot(sweepVals, errors, 'b-o', 'LineWidth', 1.5);
        xlabel(paramList{selParam});
        ylabel('RMSE(DI, ideal)');
        title(sprintf('Sensitivity of DI to %s (RMSE metric)', paramList{selParam}));
        grid on;
        end

    fprintf('Total Simulation Time = %.3f seconds.\n', toc(simStartTotal));
    end
end

%% Signal Processing
function DI_fit = dis(config)

    %% Generate Healthy and Damaged Signals
    tStart = tic;
    [t, healthySignal, c_layers] = simulateMasterSignal(0, config);
    timeHealthy = toc(tStart);
    fprintf('Time for healthy master signal generation = %.3f seconds.\n', timeHealthy);

    tStart = tic;
    [~, damagedSignal, ~] = simulateMasterSignal(config.damageVal, config);
    timeDamaged = toc(tStart);
    fprintf('Time for damaged master signal generation = %.3f seconds.\n', timeDamaged);

    %% Plot Master Signals
    figure('Color','w','Position',[100 100 600 800]);
    subplot(2,1,1);
        plot(t*1e6, healthySignal/1e6,'b','LineWidth',2);
        xlabel('Time (\mus)'); ylabel('Force (MN)');
        title('Healthy Master Signal'); grid on; xlim([0,45]);
    subplot(2,1,2);
        plot(t*1e6, damagedSignal/1e6,'r','LineWidth',2);
        xlabel('Time (\mus)'); ylabel('Force (MN)');
        title(sprintf('Damaged Master Signal (Damage = %.2f)', config.damageVal));
        grid on; xlim([0,45]);

    %% Virtual sensors & TDOA (use healthy master signal for triangulation)
    tStart = tic;
    numSensors      = 4;
    sensorPositions = config.sensorPositions;
    initialGuessLoc = [0, 0];             % starting point for triangulation
    c_surface       = c_layers(1);
    sensorSignals   = zeros(numSensors, length(t));

    % interpolate from the healthy master signal
    for iS = 1:numSensors
        d       = norm(sensorPositions(iS,:) - initialGuessLoc);
        tDelay  = d / c_surface;
        sensorSignals(iS,:) = interp1( t - tDelay, healthySignal, t, 'linear', 0 );
    end
    timeSensors = toc(tStart);
    fprintf('Time for virtual sensor processing = %.3f seconds.\n', timeSensors);

    % plot sensor waveforms
    figure('Color','w','Position',[100 100 600 800]);
    for iS = 1:numSensors
        subplot(numSensors,1,iS);
        plot(t*1e6, sensorSignals(iS,:)/1e6, 'LineWidth',1.5);
        xlabel('Time (\mus)');
        ylabel(sprintf('Sensor %d (MN)',iS));
        grid on; xlim([0,45]);
    end
    sgtitle('Virtual Sensor Signals (Healthy Case)');

    % TDOA Triangulation
    TOA = zeros(numSensors,1);
    for iS = 1:numSensors
        [~,idx] = max(abs(sensorSignals(iS,:)));
        TOA(iS) = t(idx);
    end
    t_ref = TOA(1);
    dT    = TOA - t_ref;

    % iterative least-squares method
    s_est = initialGuessLoc;  % start at the initial guess
    for it = 1:100
        H = zeros(numSensors-1,2);
        b = zeros(numSensors-1,1);
        d_ref = norm(s_est - sensorPositions(1,:));
        for iS = 2:numSensors
            d_i       = norm(s_est - sensorPositions(iS,:));
            H(iS-1,:) = (s_est - sensorPositions(iS,:))/d_i ...
                        - (s_est - sensorPositions(1,:))/d_ref;
            b(iS-1)   = c_surface * dT(iS) - (d_i - d_ref);
        end
        delta = H \ b;
        s_est = s_est + delta';
        if norm(delta) < 1e-6
            break;
        end
    end
    fprintf('Triangulated Strike Location: x=%.3f m, y=%.3f m\n', s_est(1), s_est(2));

    % plot sensor locations, initial guess, and triangulated strike
    figure('Color','w');
    hold on;
      plot(sensorPositions(:,1), sensorPositions(:,2), 'bo', ...
           'MarkerSize',6, 'LineWidth',1.5);
      plot(initialGuessLoc(1), initialGuessLoc(2), 'g*', ...
           'MarkerSize',8, 'LineWidth',1.5);
      plot(s_est(1), s_est(2), 'rx', ...
           'MarkerSize',8, 'LineWidth',1.5);
    hold off;
    legend('Sensors','Initial Guess','Triangulated Strike','Location','best');
    xlabel('X (m)');
    ylabel('Y (m)');
    title('TDOA Triangulation');
    grid on;


%% FFT Analysis 
tStart = tic;
Fs   = 1/(t(2)-t(1));
L    = length(t);
f    = Fs*(0:(L/2))/L;           % in Hz
f_MHz = f/1e6;                   % convert to MHz

FFT_h = fft(healthySignal);
P2h   = abs(FFT_h/L);
P1h   = P2h(1:L/2+1);
P1h(2:end-1) = 2*P1h(2:end-1);

FFT_d = fft(damagedSignal);
P2d   = abs(FFT_d/L);
P1d   = P2d(1:L/2+1);
P1d(2:end-1) = 2*P1d(2:end-1);

fprintf('Time for FFT computation = %.3f seconds.\n', toc(tStart));

figure('Color','w','Position',[100 100 800 600]);
% Full-range healthy
subplot(2,2,1);
    plot(f_MHz, P1h/1e6, 'b', 'LineWidth', 2);
    xlim([0, 3]);                
    xlabel('Frequency (MHz)');
    ylabel('|FFT(Healthy)| (MN)');
    title('Healthy FFT (Full Range)');
    grid on;
% Full-range damaged
subplot(2,2,2);
    plot(f_MHz, P1d/1e6, 'r', 'LineWidth', 2);
    xlim([0, 3]);
    xlabel('Frequency (MHz)');
    ylabel('|FFT(Damaged)| (MN)');
    title('Damaged FFT (Full Range)');
    grid on;
% Zoomed healthy
subplot(2,2,3);
    plot(f_MHz, P1h/1e6, 'b', 'LineWidth', 2);
    xlim([0, 1]);                
    xlabel('Frequency (MHz)');
    ylabel('|FFT(Healthy)| (MN)');
    title('Healthy FFT (Zoomed)');
    grid on;
% Zoomed damaged
subplot(2,2,4);
    plot(f_MHz, P1d/1e6, 'r', 'LineWidth', 2);
    xlim([0, 1]);
    xlabel('Frequency (MHz)');
    ylabel('|FFT(Damaged)| (MN)');
    title('Damaged FFT (Zoomed)');
    grid on;

%% PSD Analysis (0–20 MHz & 0–5 MHz) 
tStart = tic;
[PSD_h, f_psd] = pwelch(healthySignal, hamming(256), 128, 1024, Fs);
[PSD_d, ~    ] = pwelch(damagedSignal,  hamming(256), 128, 1024, Fs);
f_psd_MHz     = f_psd / 1e6;
fprintf('Time for PSD computation = %.3f seconds.\n', toc(tStart));

figure('Color','w','Position',[100 100 800 800]);

% Healthy PSD, Full Range 
subplot(2,2,1);
plot(f_psd_MHz, 10*log10(PSD_h), 'b', 'LineWidth', 1.5);
xlim([0, 20]);
xticks(0:5:20);
xlabel('Frequency (MHz)');
ylabel('Power (dB/Hz)');
title('Healthy PSD (0–20 MHz)');
grid on;

% Damaged PSD, Full Range 
subplot(2,2,2);
plot(f_psd_MHz, 10*log10(PSD_d), 'r', 'LineWidth', 1.5);
xlim([0, 20]);
xticks(0:5:20);
xlabel('Frequency (MHz)');
ylabel('Power (dB/Hz)');
title(sprintf('Damaged PSD (0–20 MHz, Damage=%.2f)', config.damageVal));
grid on;

% Healthy PSD, Zoomed
subplot(2,2,3);
plot(f_psd_MHz, 10*log10(PSD_h), 'b', 'LineWidth', 1.5);
xlim([0, 5]);
xticks(0:1:5);
xlabel('Frequency (MHz)');
ylabel('Power (dB/Hz)');
title('Healthy PSD (0–5 MHz)');
grid on;

% Damaged PSD, Zoomed  
subplot(2,2,4);
plot(f_psd_MHz, 10*log10(PSD_d), 'r', 'LineWidth', 1.5);
xlim([0, 5]);
xticks(0:1:5);
xlabel('Frequency (MHz)');
ylabel('Power (dB/Hz)');
title('Damaged PSD (0–5 MHz)');
grid on;

sgtitle('PSD Comparison');  

%% CWT Scalograms (0–20 MHz) with Mega‐Newton Scaling and Difference
% Compute and plot healthy, damaged, and difference scalograms with MN units

% 1) Compute CWT coefficients 
freqLimits = [0 20e6];
[Hh, freqs] = cwt(healthySignal, Fs, 'amor', 'FrequencyLimits', freqLimits);
[Hd, ~]      = cwt(damagedSignal,  Fs, 'amor', 'FrequencyLimits', freqLimits);

% 2) Prepare time axis in μs and cone of influence
t_us   = t * 1e6;                   % μs
coi_us = sqrt(2)./freqs * 1e6;      % μs

% 3) Frequency ticks (MHz)
fmax = 20e6;
ytHz = 0:5e6:fmax;                  % every 5 MHz

% Healthy CWT Scalogram
figure('Color','w','Position',[100 100 600 400]);
imagesc(t_us, freqs/1e6, abs(Hh)/1e6);  % divide by 1e6 is MN, freqs in MHz
axis xy; xlim([0, t_us(end)]); ylim([0, fmax/1e6]);
set(gca, 'YTick', ytHz/1e6, 'YTickLabel', ytHz/1e6);
xlabel('Time (\mus)');
ylabel('Frequency (MHz)');
title('Healthy CWT Scalogram');

commonCLim = clim;                     % capture for shared color scaling

cb = colorbar;
cb.Label.String  = '|CWT Coefficient| (MN)';
cb.TickDirection = 'out';

hold on
  plot(coi_us,           freqs/1e6, 'w--','LineWidth',1.5);
  plot(t_us(end)-coi_us, freqs/1e6, 'w--','LineWidth',1.5);
hold off

% Damaged CWT Scalogram
figure('Color','w','Position',[100 100 600 400]);
imagesc(t_us, freqs/1e6, abs(Hd)/1e6);
axis xy; xlim([0, t_us(end)]); ylim([0, fmax/1e6]);
set(gca, 'YTick', ytHz/1e6, 'YTickLabel', ytHz/1e6);
xlabel('Time (\mus)');
ylabel('Frequency (MHz)');
title(sprintf('Damaged CWT Scalogram (Damage = %.2f)', config.damageVal));

clim(commonCLim);

cb = colorbar;
cb.Label.String  = '|CWT Coefficient| (MN)';
cb.TickDirection = 'out';

hold on
  plot(coi_us,           freqs/1e6, 'w--','LineWidth',1.5);
  plot(t_us(end)-coi_us, freqs/1e6, 'w--','LineWidth',1.5);
hold off

% Difference Scalogram (Damaged – Healthy)
D = abs(Hd) - abs(Hh);

figure('Color','w','Position',[100 100 600 400]);
imagesc(t_us, freqs/1e6, D/1e6);
axis xy; xlim([0, t_us(end)]); ylim([0, fmax/1e6]);
set(gca, 'YTick', ytHz/1e6, 'YTickLabel', ytHz/1e6);
xlabel('Time (\mus)');
ylabel('Frequency (MHz)');
title('CWT Difference (Damaged – Healthy)');

diffMax = max(abs(D(:)))/1e6;           % symmetric limits in MN
clim([-diffMax diffMax]);

colormap(redbluecmap);

cb = colorbar;
cb.Label.String  = 'Δ|CWT Coefficient| (MN)';
cb.TickDirection = 'out';

hold on
  plot(coi_us,           freqs/1e6, 'k--','LineWidth',1.5);
  plot(t_us(end)-coi_us, freqs/1e6, 'k--','LineWidth',1.5);
hold off
     

%% Automatic Damage Index Calibration (with Composite DI, Prints, and Curvature-Based Threshold Regions)
% Computes and normalizes multiple damage metrics, fuses them into a single
% Damage Index (equal weight and regression weighted), prints all to CMD
% Shades “early warning”, “caution” and “danger” regions based on curvature detection.

tStart       = tic;
damageLevels = 0:0.1:1;
nLevels      = numel(damageLevels);

% Preallocate raw metrics
spectDiff   = zeros(1,nLevels);
timeShift     = zeros(1,nLevels);
ampRatio      = zeros(1,nLevels);
magRatio   = zeros(1,nLevels);
centroidShift = zeros(1,nLevels);

% Baseline healthy signal
[t, baseSig, ~] = simulateMasterSignal(0, config);
L  = length(baseSig);
dt = t(2) - t(1);
Fs = 1/dt;

% Spectral baseline
FFT_b = fft(baseSig);
P2b   = abs(FFT_b / L);
P1b   = P2b(1:floor(L/2)+1);
P1b(2:end-1) = 2*P1b(2:end-1);

% Time-domain envelope & first-echo baseline
env_h = abs(hilbert(baseSig));
[peaks_h, locs_h] = findpeaks(env_h, t, 'NPeaks',1,'SortStr','descend');
t0_h = locs_h(1);
A0_h = peaks_h(1);

% Envelope mag baseline
E_h = trapz(t, env_h);

% Spectral centroid baseline
freqs  = (0:floor(L/2)) * (Fs / L);
cent_h = freqs * P1b(:) / sum(P1b(:));

% Loop over damage levels
for i = 1:nLevels
    dVal = damageLevels(i);
    [t, curSig, ~] = simulateMasterSignal(dVal, config);

    % 1) Spectral DF
    FFT_c = fft(curSig);
    P2c   = abs(FFT_c / L);
    P1c   = P2c(1:floor(L/2)+1);
    P1c(2:end-1) = 2*P1c(2:end-1);
    spectDiff(i) = norm(P1b - P1c) / norm(P1b);

    % 2) ToF shift
    env_d = abs(hilbert(curSig));
    [peaks_d, locs_d] = findpeaks(env_d, t, 'NPeaks',1,'SortStr','descend');
    timeShift(i) = locs_d(1) - t0_h;

    % 3) Amplitude ratio
    ampRatio(i) = peaks_d(1) / A0_h;

    % 4) Envelope mag ratio
    E_d = trapz(t, env_d);
    magRatio(i) = E_d / E_h;

    % 5) Spectral centroid shift
    cent_d = freqs * P1c(:) / sum(P1c(:));
    centroidShift(i) = cent_d - cent_h;

    % Print raw metrics
    fprintf('dVal=%0.2f | SpectDF=%0.3f | AmpRa=%0.3f | ToF=%0.2fµs | MagRa=%0.3f | CenSh=%0.2fHz\n', ...
        dVal, spectDiff(i), ampRatio(i), 1e6*timeShift(i), magRatio(i), centroidShift(i));
end
fprintf('Calibration time: %.3f s\n', toc(tStart));

% Normalize each metric to [0,1]
normDF = (spectDiff - min(spectDiff)) ./ (max(spectDiff) - min(spectDiff));
normTS = (abs(timeShift)   - min(abs(timeShift)))   ./ (max(abs(timeShift))   - min(abs(timeShift)));
normAR = ((1 - ampRatio)   - min(1 - ampRatio))     ./ (max(1 - ampRatio)     - min(1 - ampRatio));
normED = (abs(magRatio - 1) - min(abs(magRatio - 1))) ./ (max(abs(magRatio - 1)) - min(abs(magRatio - 1)));
normCS = (abs(centroidShift)   - min(abs(centroidShift)))   ./ (max(abs(centroidShift))   - min(abs(centroidShift)));

% Composite Damage Index: Equal weight average
DI_equal = (normDF + normTS + normAR + normED + normCS) / 5;

% Composite Damage Index: Regression weighted
dl     = damageLevels(:);                        % 11×1 vector
M      = [normDF(:), normTS(:), normAR(:), normED(:), normCS(:)];  % 11×5
w      = M \ dl;                                 % learned weights
DI_fit = M * w;                                  % raw fitted DI
DI_fit = (DI_fit - min(DI_fit)) / (max(DI_fit) - min(DI_fit));  % re-normalize

% Print learned weights and R squared
mu_dl = mean(dl);
R2    = 1 - sum((DI_fit - dl).^2) / sum((dl - mu_dl).^2);
fprintf('Learned weights [SpectDF, |ToF|, 1-AmpRa, EnDev, |CenSh|]:\n');
disp(w');
fprintf('R^2 of fitted Damage Index = %.3f\n', R2);

% Curvature‐based threshold detection 
x    = damageLevels(:);    % 11×1 vector (0:0.1:1)
y    = DI_fit(:);          % 11×1 regression DI

% Compute first and second derivatives
dx   = diff(x);            
dy   = diff(y);            
d1   = dy ./ dx;           % first derivative (10×1)
d2   = diff(d1) ./ dx(2:end);  % second derivative (9×1)

% Compute curvature 
kappa = abs(d2) ./ (1 + d1(2:end).^2).^(3/2);  % 9×1

    % Find peaks in curvature
    [pks, locs_p] = findpeaks(kappa);

    if numel(locs_p) >= 2
        % pick two largest curvature peaks
        [~, sortIdx] = sort(pks, 'descend');
        top2 = sort(locs_p(sortIdx(1:2)));
    elseif isscalar(locs_p)
        % exactly one knee found: replicate its index and midpoint
        midIdx = ceil(length(kappa)/2);
        top2   = [locs_p; midIdx];
    else
        % no clear knees: fall back to fixed thresholds (0.40 & 0.70)
        t1idx = find(x >= 0.40,  1) - 1;
        t2idx = find(x >= 0.70,  1) - 1;
        top2  = [t1idx; t2idx];
    end

% map back to damageLevels 
th1 = x(top2(1) + 1);  
th2 = x(top2(2) + 1);

% Print threshold values
fprintf('Early‐warning threshold (damageVal) = %.2f\n', th1);
fprintf('Danger threshold = %.2f\n', th2);

% Plot with shaded threshold regions
figure('Color','w'); hold on;
yl = [0 1];

% Shade Early-Warning = green
patch([0 th1 th1 0],        [yl(1) yl(1) yl(2) yl(2)], 'g', 'FaceAlpha',0.1,'EdgeColor','none');
% Shade Caution = yellow
patch([th1 th2 th2 th1],    [yl(1) yl(1) yl(2) yl(2)], 'y', 'FaceAlpha',0.1,'EdgeColor','none');
% Shade Danger = red
patch([th2 1 1 th2],        [yl(1) yl(1) yl(2) yl(2)], 'r', 'FaceAlpha',0.1,'EdgeColor','none');

% Plot curves
plot(damageLevels, damageLevels, 'k--','LineWidth',1.5,'DisplayName','Ideal');
plot(damageLevels, DI_equal,    'b-o','LineWidth',1.5,'DisplayName','Equal-weight DI');
plot(damageLevels, DI_fit,      'r-s','LineWidth',1.5,'DisplayName','Regression-DI');

xlabel('DamageVal');
ylabel('Composite Damage Index (0-1)');
title('Composite Damage Index');
legend('Location','SouthEast');
grid on;
end 

%% simulateMasterSignal
function [t_out, masterSignal, c_layers_out] = simulateMasterSignal(dVal, config)

    %% Unpack & apply damage
    E_lyr     = config.E_layers_orig;
    rho_lyr   = config.rho_layers_orig;
    nu_lyr    = config.nu_layers;
    thick     = config.thick;
    numLayers = config.numLayers;
    if dVal>0
        stiffReduction = 0.5;
        E_lyr(config.damageLayer) = E_lyr(config.damageLayer) * (1 - stiffReduction * dVal);
    end

    %% Time‐axis & Hertz pulse
    N_points          = config.N_points;
    time_pad          = config.time_pad;
    waveEventsMax_val = config.waveEventsMax;

    % pull from GUI selection
    E_proj    = config.projectile.E;
    nu_proj   = config.projectile.nu;
    rho_proj  = config.projectile.rho;
    vel       = config.projectile.vel;
    radius    = config.projectile.radius;
    E_star    = 1/((1-nu_proj^2)/E_proj + (1-nu_lyr(1)^2)/E_lyr(1));
    mass_proj = 4/3 * pi * rho_proj * radius^3;
    t_c       = 2.94 * (mass_proj^2/(radius * E_star^2 * vel))^(1/5);
    F_max     = 2.41 .* radius.^2 .* (rho_proj.^3 .* E_star.^2 .* vel.^6).^(1/5);

    t_out     = linspace(0, time_pad * t_c, N_points);
    Ft_local  = F_max * sin(pi * t_out / t_c) .* (t_out <= t_c);

    %% Layer properties
    c_layers = sqrt(E_lyr ./ (rho_lyr .* (1 - nu_lyr.^2)));
    Z_layers = rho_lyr .* c_layers;
    Rcoeff   = [(Z_layers(2:end) - Z_layers(1:end-1)) ./ (Z_layers(2:end) + Z_layers(1:end-1)), -1];
    Tcoeff   = [2 * Z_layers(2:end) ./ (Z_layers(2:end) + Z_layers(1:end-1)), 0];
    tau      = thick ./ c_layers;

    %% Prepare event‐queue
    h_multi_local = zeros(size(t_out));
    dt_local      = t_out(2) - t_out(1);

    initEvent = struct('time',0,'layer',1,'dir',+1,'amplitude',1,'pathLen',0);
    lastEvent  = 1;
    if isinf(waveEventsMax_val)
        waveEvents = initEvent;
    else
        blankEvent = struct('time',0,'layer',0,'dir',0,'amplitude',0,'pathLen',0);
        waveEvents  = repmat(blankEvent, waveEventsMax_val, 1);
        waveEvents(1) = initEvent;
    end

    %% Damping and Dispersion 
    function A = applyDamping(a, L)
        if config.damageLayer == L
            alpha_eff = config.alpha_damp * (1 + config.damageVal);
        else
            alpha_eff = config.alpha_damp;
        end
        A = a * exp(-alpha_eff * thick(L));
    end
    function A = applyDispersion(a, dist)
        k_loc = config.k_factor;
        A     = a ./ sqrt(1 + (k_loc * dist).^2);
    end
    function placeSpike(tt, a)
        idx = round(tt / dt_local) + 1;
        if idx >= 1 && idx <= numel(h_multi_local)
            h_multi_local(idx) = h_multi_local(idx) + a;
        end
    end
    function enqueue(e)
        if isinf(waveEventsMax_val) || lastEvent < waveEventsMax_val
            lastEvent = lastEvent + 1;
            waveEvents(lastEvent) = e;
        end
    end

    %% Main loop (WaveEvent propagation)
    iEvent = 1;
    hWait  = waitbar(0, 'Processing wave events: 0 of ???');
    while iEvent <= lastEvent
        if mod(iEvent,100000) == 0 || iEvent == lastEvent
            waitbar(iEvent / lastEvent, hWait, ...
                   sprintf('Processing wave events: %d of %d', iEvent, lastEvent));
        end
        eData = waveEvents(iEvent); 
        iEvent = iEvent + 1;
        t0       = eData.time; 
        L        = eData.layer; 
        dirn     = eData.dir;
        amp0     = eData.amplitude; 
        pathSoFar = eData.pathLen;
        if t0 > t_out(end) || abs(amp0) < config.minAmp
            continue;
        end

        placeSpike(t0, amp0);

        if dirn == +1
            if L < numLayers
                R = Rcoeff(L);
                T = Tcoeff(L);
                dampedAmp = applyDamping(amp0, L);
                dt        = tau(L);
                nt        = t0 + dt;
                p2        = pathSoFar + thick(L);
                amp_disp  = applyDispersion(dampedAmp, p2);
                enqueue(struct('time',nt,'layer',L,  'dir',-1,'amplitude',amp_disp*R,'pathLen',p2));
                enqueue(struct('time',nt,'layer',L+1,'dir',+1,'amplitude',amp_disp*T,'pathLen',p2));
            else
                R         = Rcoeff(end);
                dampedAmp = applyDamping(amp0, numLayers);
                dt        = tau(end);
                nt        = t0 + dt;
                p2        = pathSoFar + thick(end);
                amp_disp  = applyDispersion(dampedAmp, p2);
                enqueue(struct('time',nt,'layer',numLayers,'dir',-1,'amplitude',amp_disp*R,'pathLen',p2));
            end
        else
            if L > 1
                Zup       = Z_layers(L-1); 
                Zdn       = Z_layers(L);
                R0        = (Zup - Zdn) / (Zup + Zdn);
                T0        = 2 * Zup / (Zup + Zdn);
                R         = R0;
                T         = T0;
                dampedAmp = applyDamping(amp0, L-1);
                dt        = tau(L-1);
                nt        = t0 + dt;
                p2        = pathSoFar + thick(L-1);
                amp_disp  = applyDispersion(dampedAmp, p2);
                enqueue(struct('time',nt,'layer',L,  'dir',+1,'amplitude',amp_disp*R,'pathLen',p2));
                enqueue(struct('time',nt,'layer',L-1,'dir',-1,'amplitude',amp_disp*T,'pathLen',p2));
            else
                nt       = t0;
                amp_disp = applyDispersion(amp0, pathSoFar);
                enqueue(struct('time',nt,'layer',1,'dir',+1,'amplitude',-amp_disp,'pathLen',pathSoFar));
            end
        end
    end
    delete(hWait);

    %% Final event count 
    if ~isinf(waveEventsMax_val)
        fprintf('Fixed mode: Processed %d events (cap = %d).\n', lastEvent, waveEventsMax_val);
    else
        fprintf('Adaptive mode: Processed %d events.\n', lastEvent);
    end

    %% Finalize
    c_layers_out = c_layers;
    fullConv     = conv(Ft_local, h_multi_local, 'full');
    masterSignal = fullConv(1:numel(t_out));
end

