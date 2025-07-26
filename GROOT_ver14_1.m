%General notes: This version (14.1) contains:
%1.     Auto-grouping of wells.
%2.     Gradual Grey Shade for Selected wells group.
%3.     Additional label info in graph tooltip.
%4.     Fixed blinking led
%5.     Fixed trailing zeros trimming problems when switching between horizontal and vertical time vectors
%6.     Added Mean and STD to all data groups.
%6.     More then 24 hours time vector bug is fixed
%7.     Includes custom groups creation (by user) and delta values analysis.
%8.     Including auto data gathering from multiple matrices in onw sheet
%9.     Including calculating ratios between mean and std of two chosen groups.
%10.    Including scenarios in which time vectors are mismatched - detecting,
%       notifying and skipping crash.
%11.    Including 2D Groups relation Plot (No Time Vector)
%12.    Enabling plot layering for ratio graphs. Enabling cleaning
%       graph button.
%13.    Added an option to clear the custom group tree
%14.    More features were added including minimize graph legends, editing
%option for axes and title labels, color coding by well in 2d graph, more matrces format
%supported, allowence of time discrapencies.
%15.    Added One-way ANOVA tab with the option to select specific time points for analysis.

%Legend: TBC - To bechecked |

%clear all;
clc;
clear moveToGroups;

createGUI();

function createGUI()

persistent grfrstime
if isempty(grfrstime)
    grfrstime = true;
end

persistent tdfrstime
if isempty(tdfrstime)
    tdfrstime = true;
end

persistent colorIdx
if isempty(colorIdx)
    colorIdx = 1;
end

% Create the main figure/window
Fig = uifigure('Name', 'Plate Reader Analyzer', 'Position', [100, 100, 1100, 700]);

% Create tabs at the top using uitabgroup
TabGroup = uitabgroup(Fig, 'Position', [10, 50, 1080, 640], 'TabLocation', 'top');
SetupTab = uitab(TabGroup, 'Title', 'Setup');
GraphEditingTab = uitab(TabGroup, 'Title', 'Graph Editing');

%% ─────────────────────  ADD: ANOVA TAB  ─────────────────────────────
% 1) create the tab + content panel
AnovaTab      = uitab(TabGroup,'Title','One-way ANOVA');
AnovaContent = uipanel(AnovaTab,'Position',[10 10 1060 600], ...
                       'Title','');          % <─ clear default title

% 2) build the left-side controls  (analysis pane)
uiLeft               = uipanel(AnovaContent,'Units','pixels', ...
                                'Position',[10 10 310 580], ...
                                'BorderType','none');

lblGroups            = uilabel(uiLeft,'Position',[10 540 180 20], ...
                                'Text','Select groups (multi-select):');
lstGroups            = uilistbox(uiLeft,'Position',[10 300 290 240], ...
                                'Items',{'(No groups yet)'}, ...
                                'Multiselect','on');

lblFactorName        = uilabel(uiLeft,'Position',[10 260 180 20], ...
                                'Text','Factor name (optional):');
fldFactorName        = uieditfield(uiLeft,'text', ...
                                   'Position',[10 235 290 22]);

lblPostHoc           = uilabel(uiLeft,'Position',[10 120 120 20], ...
                                'Text','Post-hoc test:');
ddPostHoc            = uidropdown(uiLeft,'Position',[130 120 170 22], ...
                                'Items',{'None','Tukey','Bonferroni'}, ...
                                'Value','None');

% --- Time-point selector ----------------------------------------------
lblTimePt = uilabel(uiLeft,'Position',[10 210 180 20], ...
                    'Text','Time point(s) for ANOVA:');

lstTimePt = uilistbox(uiLeft,'Position',[10 150 290 60], ... % ~3 rows tall
                      'Items',{'(fill after import)'}, ...
                      'Multiselect','on');
setappdata(Fig,'lstTimePt',lstTimePt);

% --- Show post-hoc plot button (disabled until an ANOVA + post-hoc exist)
btnShowPostHoc = uibutton(uiLeft,'push', ...
    'Position',[10 90 290 24], ...
    'Text','Show post-hoc plot', ...
    'Enable','off', ...
    'ButtonPushedFcn',@showPostHoc);
setappdata(Fig,'btnShowPostHoc',btnShowPostHoc);

btnRun               = uibutton(uiLeft,'push', ...
                                'Position',[10 60 140 26], ...
                                'Text','Run ANOVA', ...
                                'ButtonPushedFcn',@runAnova);

btnClearAnova        = uibutton(uiLeft,'push', ...
                                'Position',[160 60 140 26], ...
                                'Text','Clear ANOVA', ...
                                'ButtonPushedFcn',@clearAnova);

btnExport            = uibutton(uiLeft,'push', ...
                                'Position',[10 20 290 26], ...
                                'Text','Export results table', ...
                                'ButtonPushedFcn',@exportAnovaTable);

statusLabel          = uilabel(uiLeft,'Position',[10 0 290 22], ...
                                'Text','', ...
                                'FontColor',[0.8 0 0]);

% 3) build the right-side results pane
uiRight              = uipanel(AnovaContent,'Units','pixels', ...
                                'Position',[330 10 720 580], ...
                                'BorderType','none');

tblResults           = uitable(uiRight,'Position',[10 360 700 200], ...
                                'ColumnName',{'SS','df','MS','F','p'}, ...
                                'RowName',{},'Data',[]);

axAnova              = uiaxes(uiRight,'Position',[10 10 700 330]);
title(axAnova,''); xlabel(axAnova,''); ylabel(axAnova,'');

% keep handles for later use
setappdata(Fig,'lstGroups',lstGroups);
setappdata(Fig,'tblResults',tblResults);
setappdata(Fig,'axAnova',axAnova);
setappdata(Fig,'statusLabel',statusLabel);

% ─────────────────── helper: populate listbox when groups exist ─────
addlistener(Fig,'UserData','PostSet',@(s,e)populateGroupList());
populateGroupList();                 % first call at GUI init

    function populateGroupList
        groupedData = getappdata(Fig,'groupedData');
        if isempty(groupedData), return; end
        gnames = fieldnames(groupedData);
        if isempty(gnames)
            lstGroups.Items       = {'(No groups yet)'};
            lstGroups.Enable      = 'off';
            btnRun.Enable         = 'off';
        else
            lstGroups.Items       = gnames;
            lstGroups.Enable      = 'on';
            btnRun.Enable         = 'on';
        end
    end

function [commonTime,dataMatrix] = combineGroupData_A(gStruct)
    % gStruct has .plotData, .plotTime the same as the original
    nWells = numel(gStruct.plotData);
    if nWells==0, commonTime=[]; dataMatrix=[]; return; end

    commonTime = gStruct.plotTime{1}(:)';                 % assume aligned
    dataMatrix = nan(nWells, numel(commonTime));

    for k = 1:nWells
        thisT = gStruct.plotTime{k}(:)';
        if ~isequal(thisT,commonTime)                     % simple guard
            warning('Time mismatch in group; skipped well %d',k);
            continue
        end
        dataMatrix(k,:) = gStruct.plotData{k}(:)';
    end
end

% ─────────────────────────  RUN ANOVA  ───────────────────────────────
    function runAnova(~,~)
        % reset visuals
        statusLabel.Text      = '';                           %#ok
        tblResults.Data       = [];
        cla(axAnova,'reset');

        timeUnit = getappdata(Fig,'graphTimeUnit');
        switch timeUnit
            case 'Seconds', tol = 5;          % 5 s
            case 'Minutes', tol = 5/60;       % 0.0833 min
            case 'Hours',   tol = 5/3600;     % 0.00139 h
            otherwise,      tol = 5;          % fallback
        end

        % gather selection
        selGroups   = lstGroups.Value;
        if isempty(selGroups)
            statusLabel.Text = 'Select at least one group.'; return
        end

        groupedData = getappdata(Fig,'groupedData');
        if isempty(groupedData), statusLabel.Text='No grouped data.'; return; end

        % convert each group to a single value per well (mean of time-series)
        groupVectors = {}; grpLabels = {};
        for g = 1:numel(selGroups)
            gname                  = selGroups{g};
            [~,mat] = combineGroupData_A(groupedData.(gname));
            if isempty(mat), continue; end
            % ---- use only the selected time column(s) -----------------------------
            lstTimePt = getappdata(Fig,'lstTimePt');
            selTimes  = str2double(lstTimePt.Value);          % numeric vector
            if isempty(selTimes)
                statusLabel.Text = 'Select at least one time point.'; return
            end

            timeVec = groupedData.(gname).plotTime{1}(:)';   % row vector
            colIdx = zeros(size(selTimes));                  % pre-alloc

            for kSel = 1:numel(selTimes)
                % find closest time stamp in this group
                [d, idx] = min(abs(timeVec - selTimes(kSel)));
                if d <= tol
                    colIdx(kSel) = idx;                      % accept match
                end
            end

            if any(colIdx == 0)      % at least one time not found within ±5 s
                statusLabel.Text = 'One or more selected times not present in all groups.';
                return
            end

            sliceMat = mat(:, colIdx);                        % wells × chosen times
            valPerWell = mean(sliceMat, 2, 'omitnan');        % average if >1 time
            groupVectors{end+1} = valPerWell;                 % store
            % -----------------------------------------------------------------------
            grpLabels(end+1,:)     = {gname};                %#ok<AGROW>
        end
        if numel(groupVectors)<2
            statusLabel.Text = 'Need ≥2 groups for ANOVA.'; return
        end

        % run one-way anova silently
        statsFactorName = strtrim(fldFactorName.Value);
        if isempty(statsFactorName), statsFactorName = 'Factor';
        setappdata(Fig,'factorLabel',statsFactorName);end

        % flatten observations and matching labels  (handles unequal sizes)
        dataAll = [];   labCell = {};
        for g = 1:numel(groupVectors)
            v        = groupVectors{g};                    % column vector
            dataAll  = [dataAll ; v];                      %#ok<AGROW>
            labCell  = [labCell ; repmat({grpLabels{g}},numel(v),1)]; %#ok<AGROW>
        end

        % drop any NaNs (anova1 cannot accept them)
        keep          = ~isnan(dataAll);
        dataClean     = dataAll(keep);
        grpClean      = labCell(keep);

        % run ANOVA on unequal group sizes
        [p,tbl,stats] = anova1(dataClean, grpClean, 'off');

        % fill table
        tblResults.RowName = tbl(2:end-1,1);    % skip header & Total rows

        raw   = tbl(2:end-1,2:6);               % df, SS, MS, F, p
        [nR,nC] = size(raw);
        numMat  = NaN(nR,nC);                   % pre-fill with NaNs

        for r = 1:nR
            for c = 1:nC
                v = raw{r,c};
                if isnumeric(v) && isscalar(v) && ~isempty(v)
                    numMat(r,c) = v;            % copy numeric value
                end                             % leave NaN otherwise
            end
        end
        tblResults.Data = numMat;               % populate table

        % post-hoc?
        ph = ddPostHoc.Value;
        btnShowPostHoc = getappdata(Fig,'btnShowPostHoc');   % handle

        if ~strcmp(ph,'None')
            % map GUI label → multcompare ctype
            switch ph
                case 'Tukey',      ctype = 'tukey-kramer';   % or 'hsd'
                case 'Bonferroni', ctype = 'bonferroni';
                otherwise,         ctype = lower(ph);
            end
            mcmp  = multcompare(stats,'ctype',ctype,'display','off'); %#ok<NASGU>
            btnShowPostHoc.Enable = 'on';          % allow user to open the figure
            setappdata(Fig,'anovaStats',stats);    % stash for later
        else
            btnShowPostHoc.Enable = 'off';
            if isappdata(Fig,'anovaStats')
                rmappdata(Fig,'anovaStats');
            end
        end

        % --- visual: box-and-whisker inside the panel -------------------
        cla(axAnova,'reset');   hold(axAnova,'on');

        % build a single data vector and matching category vector
        dataAll = [];
        labCell = {};                               % plain cell of char-vecs
        for g = 1:numel(groupVectors)
            v          = groupVectors{g};
            dataAll    = [dataAll ; v];             %#ok<AGROW>
            labCell    = [labCell ; ...
                         repmat({grpLabels{g}}, numel(v), 1)]; %#ok<AGROW>
        end
        catAll = categorical(labCell);              % convert once, at the end

        % draw to the embedded axes, colour by group
        boxchart(axAnova, catAll, dataAll, 'GroupByColor', catAll);

        ylabel(axAnova, statsFactorName);
        lstTimePt   = getappdata(Fig,'lstTimePt');          % grab selector
        selTimesStr = strjoin(lstTimePt.Value, ', ');
        title(axAnova, sprintf('One-way ANOVA  (t = %s)  p = %.3g', selTimesStr, p));
        grid  (axAnova,'on');
        hold  (axAnova,'off');

        % store results so export can grab them
        setappdata(Fig,'anovaTbl',tbl);
        setappdata(Fig,'anovaP',p);
    end

% ─────────────────────────  CLEAR ANOVA  ─────────────────────────────
    function clearAnova(~,~)
        lstGroups.Value = {};
        fldFactorName.Value = '';
        tblResults.Data  = [];
        cla(axAnova,'reset');
        btnShowPostHoc = getappdata(Fig,'btnShowPostHoc');
        btnShowPostHoc.Enable = 'off';
        if isappdata(Fig,'anovaStats')
            rmappdata(Fig,'anovaStats');
        end
        statusLabel.Text = '';
        rmappdata(Fig,'anovaTbl'); rmappdata(Fig,'anovaP');
    end

% ─────────────────────────  EXPORT TABLE  ────────────────────────────
    function exportAnovaTable(~,~)
        if ~isappdata(Fig,'anovaTbl')
            statusLabel.Text = 'Run ANOVA first.'; return
        end
        [file, path] = uiputfile({'*.csv';'*.xlsx'},'Save ANOVA table');
        if isequal(file,0), return; end
        outPath = fullfile(path,file);
        anovaCell  = getappdata(Fig,'anovaTbl');    % full ANOVA cell array
        anovaCell  = anovaCell(2:end-1,:);          % skip header & Total
        tblOut     = cell2table(anovaCell, ...
                      'VariableNames',{'Source','df','SS','MS','F','p'});
        try
            if endsWith(outPath,'.csv'),  writetable(tblOut,outPath);
            else,                         writetable(tblOut,outPath,'FileType','spreadsheet');
            end
            statusLabel.Text = 'Results exported successfully.';
            statusLabel.FontColor = [0 0.5 0];
        catch ME
            statusLabel.Text = ['Export failed: ' ME.message];
            statusLabel.FontColor = [0.8 0 0];
        end
    end

% ───────── display MATLAB's multiple-comparison figure ─────────
    function showPostHoc(~,~)
        stats = getappdata(Fig,'anovaStats');
        if isempty(stats), return; end          % safety

        ph = ddPostHoc.Value;
        switch ph
            case 'Tukey',      ctype = 'tukey-kramer';
            case 'Bonferroni', ctype = 'bonferroni';
            otherwise,         return;          % 'None' or unknown
        end

        % multcompare → third output is the figure handle
        [~,~,hFig] = multcompare(stats,'ctype',ctype,'alpha',0.05);

        % add the user-supplied factor label to the x-axis
        if isgraphics(hFig)
            ax = findobj(hFig,'Type','axes','-not','Tag','Colorbar');
            if ~isempty(ax)
                xlabel(ax(1), getappdata(Fig,'factorLabel'));
            end
        end
    end
%% ───────────────────── END ANOVA TAB BLOCK ─────────────────────────


% Style tabs to have rounded corners
set(SetupTab, 'BackgroundColor', [0.94, 0.94, 0.94]);
set(GraphEditingTab, 'BackgroundColor', [0.94, 0.94, 0.94]);

% Create the Main Setup panel
SetupContent = uipanel(SetupTab, 'Position', [10, 10, 1060, 600], 'Title', 'Excel Data Analysis');

% Create the Graph Editing panel
GraphEditingContent = uipanel(GraphEditingTab, 'Position', [10, 10, 1060, 600], 'Title', 'Graph Editing');

% Integrate Excel Data Analysis GUI into Setup Tab
integrateExcelDataAnalysisGUI(SetupContent);

% Gui of Graph Editing Tab
graphEditingGUI(GraphEditingContent);

% Create the toolbar panel at the bottom
Panel = uipanel(Fig, 'Position', [10, 0, 1080, 40]);

% Create "log" button on the left side of the toolbar panel
LogButton = uibutton(Panel, 'Position', [10, 10, 80, 20], 'Text', 'Log', 'ButtonPushedFcn', @logCallback);

% Create a panel to display log messages to the right of the "log" button
LogPanel = uipanel(Panel, 'Position', [100, 8, 730, 25]);
LogMessage = uilabel(LogPanel, 'Position', [5, 5, 730, 15], 'Text', '');

% Create LED indicator (uilamp) on the right side
LED = uilamp(Panel, 'Position', [930, 10, 20, 20], 'Color', 'green');

% Create a button to simulate processing (blinking LED) - TBC: THE BUTTON IS A TEST CASE, but it seems the function blinkLED(LED) isn't called at the designated process.
BlinkButton = uibutton(Panel, 'Position', [840, 10, 80, 20], 'Text', 'Blink LED', 'ButtonPushedFcn', @(~, ~) blinkLED(LED));

% Create time and date display on the bottom right
TimeDate = uilabel(Panel, 'Position', [960, 10, 120, 20]);

% Create a timer to update the time and date display every second
t = timer('ExecutionMode', 'fixedRate', 'Period', 1, 'TimerFcn', @(~, ~) updateTimeDate(TimeDate));
start(t);

% Initialize the log messages
setappdata(Fig, 'txtMessages', []);
setappdata(Fig, 'allMessages', {});

% Global variables to store data from analysis process
data = {};
time = {};
normalizedWellPlate = {};
parentCounter = 1;

% Nested functions for callbacks
    function updateTimeDate(TimeDate)
        currentTime = datestr(now, 'dd/mm/yyyy HH:MM:SS');
        TimeDate.Text = currentTime;
    end

    function blinkLED(LED)
        % Toggle LED color
        if isequal(LED.Color, [0.8, 0.8, 0.8]) % If gray
            LED.Color = 'yellow'; % Turn yellow
        else
            LED.Color = [0.8, 0.8, 0.8]; % Turn gray
        end
        drawnow; % Ensure GUI updates are processed
    end

    function logCallback(~, ~)
        % Retrieve 'txtMessages' from appdata
        txtMessages = getappdata(Fig, 'txtMessages');

        % If it's empty OR the handle is invalid, make a new log window
        if isempty(txtMessages) || ~isgraphics(txtMessages)
            LogFig = uifigure('Name', 'Log Messages', 'Position', [300, 300, 600, 400]);
            txtMessages = uitextarea(LogFig, ...
                'Position', [10, 10, 580, 380], ...
                'Editable', 'off');
            setappdata(Fig, 'txtMessages', txtMessages);
        end

        % Get the current log messages
        allMessages = getappdata(Fig, 'allMessages');

        % If no messages yet, set a friendly placeholder
        if isempty(allMessages)
            allMessages = {'No messages logged yet.'};
        end

        % Finally, update the text area if it's still valid
        if isgraphics(txtMessages)
            txtMessages.Value = allMessages;
            drawnow; % Refresh GUI
        else
            warning('txtMessages is invalid or was closed. Cannot update log.');
        end
    end



    function graphEditingGUI(parentPanel)


        % Create UI components for graph selection and display in Graph Editing panel
        lblGraphType = uilabel(parentPanel, 'Position', [20, 540, 100, 20], 'Text', 'Select Graph:');
        bgGraphType = uibuttongroup(parentPanel, 'Position', [20, 430, 250, 100]);
        rbSelectedWellPlates = uiradiobutton(bgGraphType, 'Position', [10, 70, 200, 20], 'Text', 'Selected Well Plates');
        rbAvgStd = uiradiobutton(bgGraphType, 'Position', [10, 50, 200, 20], 'Text', 'Average and STD');
        rbDelta = uiradiobutton(bgGraphType, 'Position', [10, 30, 200, 20], 'Text', 'Delta');
        rbDeltaMeanSTD = uiradiobutton(bgGraphType, 'Position', [10, 10, 200, 20], 'Text', 'Average and STD Deltas');
        lblStdLegend = uilabel(parentPanel, 'Position', [350, 10, 400, 20], 'Text', 'Filled areas in the plot represent the Standard Deviation (STD)');

        % Create "Select Group" label and dropdown
        lblSelectGroup = uilabel(parentPanel, ...
            'Position', [20, 380, 100, 20], ...
            'Text', 'Select Group:');

        groupDropdown = uidropdown(parentPanel, ...
            'Position', [20, 360, 100, 22], ...
            'Items', {'(No groups yet)'}, ...
            'ValueChangedFcn', @(dd, ~) onGroupSelected(dd), 'Enable', 'off'); %groupDropdown is set to disable at start

        bgGraphType.SelectionChangedFcn = @(bg,event) onGroupSelected(groupDropdown); %Calling the group selection function when changing graph type to trigger the latest graph drawing.

        setappdata(Fig, 'groupDropdown', groupDropdown);

        % Create axes for graph display
        axGraph = uiaxes(parentPanel, 'Position', [350, 40, 700, 500]);
        axGraph.Visible = 'on';

        % Store axGraph in appdata to be accessed later
        setappdata(Fig, 'axGraph', axGraph);

        % Disable Graph Editing UI components initially
        bgGraphType.Enable = 'off';
        setappdata(Fig, 'bgGraphType', bgGraphType);
        axGraph.Visible = 'on';

        function onGroupSelected(dd)
            % Get which group was selected
            selectedGroup = dd.Value;

            % Retrieve the figure handle
            Fig = ancestor(dd, 'figure');
            groupedData = getappdata(Fig, 'groupedData');

            % Get axes and groupedData from appdata
            axGraph = getappdata(Fig, 'axGraph');

            % Get the selected radio button's text
            bgGraphType = getappdata(Fig, 'bgGraphType');
            if isvalid(bgGraphType) && ~isempty(bgGraphType.SelectedObject)
                selectedRadio = bgGraphType.SelectedObject.Text;
            else
                selectedRadio = 'Selected Well Plates'; % fallback
            end

            % Safety checks
            if isempty(selectedGroup) || ~isfield(groupedData, selectedGroup) && ~strcmp(selectedGroup, 'All Groups')
                return;
            end

            % Clear and plot only the selected group
            cla(axGraph, 'reset');          % Reset axes properties
            delete(axGraph.Children);       % Remove any lingering plots
            hold(axGraph, 'on');
            
            switch selectedRadio
                %================================================
                % CASE 1: "Selected Well Plates"
                %================================================
                case 'Selected Well Plates'
                    grfrstime = true; tdfrstime = true;
                    %-----------------------------------------------
                    % Selected Well Plates logic for plotting
                    %-----------------------------------------------
                    if strcmp(selectedGroup, 'All Groups')
                        % Plot all groups
                        prefixes = fieldnames(groupedData);
                        colorMap = lines(numel(prefixes));  % Distinct colors for each group

                        for p = 1:numel(prefixes)
                            prefix = prefixes{p};
                            color_ = colorMap(p, :);
                            % For each well in this prefix
                            for k = 1:length(groupedData.(prefix).plotData)
                                if k == 1
                                    % First plot in this group -> visible in legend
                                    plot(axGraph, ...
                                        groupedData.(prefix).plotTime{k}, ...
                                        groupedData.(prefix).plotData{k}, ...
                                        'DisplayName', prefix, ...      % Display group name only
                                        'Color', color_);
                                else
                                    % Subsequent plots in this group -> hidden from legend
                                    plot(axGraph, ...
                                        groupedData.(prefix).plotTime{k}, ...
                                        groupedData.(prefix).plotData{k}, ...
                                        'HandleVisibility','off', ...   % Hide from legend
                                        'Color', color_);
                                end
                            end
                        end

                    else
                        % Plot only the selected group

                        numPlots = length(groupedData.(selectedGroup).plotData);
                        colorMap = linspace(0, 0.75, numPlots); % Create grayscale values from black (0) to light gray (1)

                        % Plot each well in the group
                        for k = 1:length(groupedData.(selectedGroup).plotData)
                            color_ = [colorMap(k), colorMap(k), colorMap(k)]; % Convert grayscale value to RGB
                            hPlot = plot(axGraph, ...
                                groupedData.(selectedGroup).plotTime{k}, ...
                                groupedData.(selectedGroup).plotData{k}, ...
                                'DisplayName', groupedData.(selectedGroup).plotLabels{k}, ...
                                'Color', color_);

                            % Attach series name as UserData (to show later on in
                            % tooltip info)
                            hPlot.UserData = groupedData.(selectedGroup).plotLabels{k};

                            % Defining axes names for tooltip
                            hPlot.DataTipTemplate.DataTipRows(1).Label = 'X'; % Keep the x-value label
                            hPlot.DataTipTemplate.DataTipRows(2).Label = 'Y'; % Keep the y-value label

                            % Add a custom data tip row for the series name
                            hPlot.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Series', ...
                                @(~, ~) hPlot.UserData); % Use UserData via a function handle

                        end
                    end

                    %================================================
                    % CASE 2: "Average and STD"
                    %================================================
                case 'Average and STD'
                    grfrstime = true; tdfrstime = true;
                    % For "All Groups": show mean & STD
                    if strcmp(selectedGroup, 'All Groups')
                        prefixes = fieldnames(groupedData);
                        colorMap = lines(numel(prefixes));

                        for p = 1:numel(prefixes)
                            prefix = prefixes{p};
                            color_ = colorMap(p, :);

                            % Gather all data from that prefix
                            [commonTime, dataMatrix] = combineGroupData(groupedData.(prefix));

                            % Compute mean and std
                            meanData = mean(dataMatrix, 1, 'omitnan');
                            stdData  = std(dataMatrix, 0, 1, 'omitnan'); % sample std

                            % Compute upper and lower bounds
                            lowerBound = meanData - stdData;
                            upperBound = meanData + stdData;

                            % Plot the fill area for std
                            fill(axGraph, ...
                                [commonTime, fliplr(commonTime)], ...
                                [lowerBound, fliplr(upperBound)], ...
                                color_, 'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility', 'off')

                            % Plot only mean
                            plot(axGraph, commonTime, meanData, ...
                                'DisplayName', ['Mean of ' prefix], ...
                                'LineWidth', 2, ...
                                'Color', color_);
                        end

                    else
                        % If a single group: show mean + std
                        prefix = selectedGroup;

                        % Combine all data from that group
                        [commonTime, dataMatrix] = combineGroupData(groupedData.(prefix));

                        if isempty(commonTime) || isempty(dataMatrix)
                            hold(axGraph, 'off');
                            legend(axGraph, 'hide');
                            return;
                        end

                        % Compute mean and std
                        meanData = mean(dataMatrix, 1, 'omitnan');
                        stdData  = std(dataMatrix, 0, 1, 'omitnan'); % sample std

                        % Compute upper and lower bounds
                        lowerBound = meanData - stdData;
                        upperBound = meanData + stdData;

                        % Plot the fill area for std
                        fill(axGraph, ...
                            [commonTime, fliplr(commonTime)], ...
                            [lowerBound, fliplr(upperBound)], ...
                            'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility', 'off');

                        % Plot the mean line
                        plot(axGraph, commonTime, meanData, ...
                            'b', 'LineWidth', 2, ...
                            'DisplayName', [prefix ' Mean']);
                    end

                    %================================================
                    % CASE 3: "Delta"
                    %================================================
                case 'Delta'
                    grfrstime = true; tdfrstime = true;
                    %-----------------------------------------------
                    % Logic for plotting the DELTAS (differences)
                    % between consecutive data points.
                    % Example:
                    % If data = [3, 5, 8], then diff(data) = [2, 3].
                    %-----------------------------------------------
                    if strcmp(selectedGroup, 'All Groups')
                        % Plot all groups, each group gets its own color
                        prefixes = fieldnames(groupedData);
                        colorMap = lines(numel(prefixes));  % Distinct colors for each group

                        for p = 1:numel(prefixes)
                            prefix = prefixes{p};
                            color_ = colorMap(p, :);
                            % For each well in this prefix
                            for k = 1:length(groupedData.(prefix).plotData)
                                originalData = groupedData.(prefix).plotData{k};
                                if length(originalData) < 2
                                    % Not enough data points to take a difference
                                    continue;
                                end

                                % Compute the difference
                                deltaData = diff(originalData);

                                % Decide on the X values for the differences
                                % One common approach: shift time to "middle" or simply
                                % use time(2:end) to align with each consecutive difference
                                timeVals = groupedData.(prefix).plotTime{k};
                                deltaTime = timeVals(2:end);

                                % Create the plot
                                if k == 1
                                    % First plot in this group -> visible in legend
                                    plot(axGraph, ...
                                        deltaTime, ...
                                        deltaData, ...
                                        'o-', ...   % marker + line (adjust as needed)
                                        'DisplayName', prefix, ...
                                        'Color', color_);
                                else
                                    % Subsequent plots in this group -> hidden from legend
                                    plot(axGraph, ...
                                        deltaTime, ...
                                        deltaData, ...
                                        'o-', ...
                                        'HandleVisibility','off', ...
                                        'Color', color_);
                                end
                            end
                        end

                    else
                        % Plot only the selected group

                        % Extract the data from the selected group
                        numPlots = length(groupedData.(selectedGroup).plotData);
                        % Create grayscale values from black (0) to lighter gray (0.75)
                        colorMap = linspace(0, 0.75, numPlots);

                        for k = 1:numPlots
                            originalData = groupedData.(selectedGroup).plotData{k};
                            if length(originalData) < 2
                                % Not enough data to compute a difference
                                continue;
                            end

                            % Compute the difference
                            deltaData = diff(originalData);

                            % Construct the time vector for these differences
                            timeVals = groupedData.(selectedGroup).plotTime{k};
                            deltaTime = timeVals(2:end);

                            % Convert grayscale value to RGB
                            color_ = [colorMap(k), colorMap(k), colorMap(k)];

                            % Plot the data (using circles for clarity)
                            hPlot = plot(axGraph, ...
                                deltaTime, ...
                                deltaData, ...
                                'o-', ...  % Use a line plus circle marker
                                'DisplayName', groupedData.(selectedGroup).plotLabels{k}, ...
                                'Color', color_);

                            % Attach series name as UserData (for tooltip info)
                            hPlot.UserData = groupedData.(selectedGroup).plotLabels{k};

                            % Defining axes names for tooltip
                            hPlot.DataTipTemplate.DataTipRows(1).Label = 'X';
                            hPlot.DataTipTemplate.DataTipRows(2).Label = 'Δ (diff)';

                            % Add a custom data tip row for the series name
                            hPlot.DataTipTemplate.DataTipRows(end+1) = ...
                                dataTipTextRow('Series', @(~, ~) hPlot.UserData);
                        end
                    end

                    %================================================
                    % CASE 4: "TBD"
                    %================================================

                case 'Average and STD Deltas'
                    grfrstime = true; tdfrstime = true;
                    %----------------------------------------------------------
                    % Compute and plot the average + std of the delta (diff)
                    % data instead of the raw data.
                    %
                    % If a group has data = [3,5,8] for a well,
                    % its delta data = [2,3].
                    % We'll combine across wells, compute the mean & std
                    % of those deltas, and plot.
                    %----------------------------------------------------------
                    if strcmp(selectedGroup, 'All Groups')
                        prefixes = fieldnames(groupedData);
                        colorMap = lines(numel(prefixes));

                        for p = 1:numel(prefixes)
                            prefix = prefixes{p};
                            color_ = colorMap(p, :);

                            % 1) Gather all data from that prefix (interpolates onto a common time)
                            [commonTime, dataMatrix] = combineGroupData(groupedData.(prefix));

                            % Skip if there's not enough data
                            if isempty(commonTime) || size(dataMatrix, 2) < 2
                                continue;
                            end

                            % 2) Compute delta for each row (well) across time
                            %    'diff' along dimension 2 (time dimension)
                            deltaMatrix = diff(dataMatrix, 1, 2);

                            % 3) Adjust the commonTime (one fewer point after diff)
                            deltaTime = commonTime(2:end);

                            % 4) Compute mean & std of delta data (omit NaNs)
                            meanData = mean(deltaMatrix, 1, 'omitnan');
                            stdData  = std(deltaMatrix, 0, 1, 'omitnan'); % sample std by default

                            % 5) Compute upper and lower bounds
                            lowerBound = meanData - stdData;
                            upperBound = meanData + stdData;

                            % 6) Plot the fill area for std
                            fill(axGraph, ...
                                [deltaTime, fliplr(deltaTime)], ...
                                [lowerBound, fliplr(upperBound)], ...
                                color_, 'FaceAlpha', 0.1, 'EdgeColor', 'none','HandleVisibility', 'off');

                            % 7) Plot the mean
                            plot(axGraph, deltaTime, meanData, ...
                                'DisplayName', ['Mean of delta ' prefix], ...
                                'LineWidth', 2, ...
                                'Color', color_);
                        end

                    else
                        %-----------------------------------------------------
                        % If a single group is selected
                        %-----------------------------------------------------
                        prefix = selectedGroup;

                        % 1) Combine all data from that group
                        [commonTime, dataMatrix] = combineGroupData(groupedData.(prefix));

                        % Skip if there's not enough data
                        if isempty(commonTime) || size(dataMatrix, 2) < 2
                            hold(axGraph, 'off');
                            legend(axGraph, 'hide');
                            return;
                        end

                        % 2) Compute delta for each row (well) across time
                        deltaMatrix = diff(dataMatrix, 1, 2);

                        % 3) Adjust the time array
                        deltaTime = commonTime(2:end);

                        % 4) Compute mean & std
                        meanData = mean(deltaMatrix, 1, 'omitnan');
                        stdData  = std(deltaMatrix, 0, 1, 'omitnan');

                        % 5) Compute upper and lower bounds
                        lowerBound = meanData - stdData;
                        upperBound = meanData + stdData;

                        % 6) Plot the fill area for std
                        fill(axGraph, ...
                            [deltaTime, fliplr(deltaTime)], ...
                            [lowerBound, fliplr(upperBound)], ...
                            'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility', 'off');

                        % 7) Plot the mean line
                        plot(axGraph, deltaTime, meanData, ...
                            'b', 'LineWidth', 2, ...
                            'DisplayName', [prefix ' Mean of deltas']);
                    end

                case 'TBD'
                    text(axGraph, 0.5, 0.5, ...
                        'TBD logic TBC', ...
                        'Units','normalized', ...
                        'HorizontalAlignment','center');
            end

            hold(axGraph, 'off');
            lgd = legend(axGraph, 'show'); % Create the legend and store its handle
            lgd.Interpreter = 'none';     % Disable TeX interpretation for the legend labels
            xlabel(axGraph, 'Time');
            ylabel(axGraph, 'Data Value');
        end

        %------------------------------------------
        % Helper function to unify a group's data
        %------------------------------------------
        function [commonTime, dataMatrix] = combineGroupData(groupDataStruct)
            % groupDataStruct has fields:
            %     .plotData   = cell of vectors
            %     .plotTime   = cell of time vectors
            %     .plotLabels = cell of strings

            numWells = length(groupDataStruct.plotData);
            if numWells == 0
                commonTime = [];
                dataMatrix = [];
                return;
            end

            % Assume all time vectors have the same length and alignment.
            % If your data might differ in length, you'd need interpolation.
            commonTime = groupDataStruct.plotTime{1}(:)';
            dataMatrix = zeros(numWells, length(commonTime));

            for k = 1:numWells

                % Check if the current time vector matches commonTime
                currentTime = groupDataStruct.plotTime{k}(:)';
                if length(currentTime)~=length(commonTime) ...
                        || any(abs(currentTime - commonTime) > 180)
                    custom_fprintf('Warning: time mismatch >180 s in dataset %d; skipping.', k);
                    commonTime = []; dataMatrix = []; return;
                end
                % Make sure data is row-oriented for easy stacking
                rowData = groupDataStruct.plotData{k}(:)';
                dataMatrix(k, :) = rowData;
            end
            custom_fprintf('Group data combined succesfully');
        end

        %====================================================
        % 2) Create the "Groups ratio" UI below the dropdown
        %====================================================
        % Label

        lblRatioTitle = uilabel(parentPanel, ...
            'Position', [20, 310, 120, 22], ...
            'Text', 'Groups ratio:');

        % First group dropdown
        groupDropdown1 = uidropdown(parentPanel, ...
            'Position', [20, 280, 100, 22], ...
            'Items', {'(No groups yet)'});

        % Slash label
        lblSlash = uilabel(parentPanel, ...
            'Position', [130, 280, 10, 22], ...
            'Text', '/');

        % Second group dropdown
        groupDropdown2 = uidropdown(parentPanel, ...
            'Position', [140, 280, 100, 22], ...
            'Items', {'(No groups yet)'});

        % Calculate button
        btnCalculate = uibutton(parentPanel, ...
            'Position', [20, 230, 220, 22], ...  % left=20, width=220, height=22
            'Text', 'Plot groups ratio over time vector', ...
            'ButtonPushedFcn', @(~, ~) onCalculateRatio());

        % Second identical-sized button below the first button
        btnNoTimePlot = uibutton(parentPanel, ...
            'Position', [20, 200, 220, 22], ...
            'Text', '2D Relation Plot (No Time Vector)', ...
            'ButtonPushedFcn', @(~, ~) on2DRatio());

        % Error / status label
        errorLabel = uilabel(parentPanel, ...
            'Position', [20, 170, 400, 22], ...
            'Text', '', ...
            'FontColor', 'r');

        % Clear graph button
        btnClearPlot = uibutton(parentPanel, ...
            'Position', [20, 140, 220, 22], ...
            'Text', 'Clear Graph', ...
            'ButtonPushedFcn', @(~, ~) clearPlot());

        btnEditLabels = uibutton(parentPanel, ...
        'Position', [20, 110, 220, 22], ...
        'Text', 'Edit Graph Labels', ...
        'ButtonPushedFcn', @(~, ~) editLabelsCallback());

        % Store these in appdata if desired
        setappdata(Fig, 'groupDropdown1', groupDropdown1);
        setappdata(Fig, 'groupDropdown2', groupDropdown2);
        setappdata(Fig, 'ratioErrorLabel', errorLabel);

        %-------------------------------------------
        % Callback: onCalculateRatio
        %-------------------------------------------
        function onCalculateRatio()
            % Clear any old error text
            errorLabel.Text = '';
            

            nColors = 10; % total distinct colors in the gradient
            cmap = hsv(nColors) * 0.8; % generates rainbow-like colors

            tdfrstime = true;

            groupDropdown1 = getappdata(Fig, 'groupDropdown1');
            groupDropdown2 = getappdata(Fig, 'groupDropdown2');

            % Retrieve the selected group names from the two dropdowns
            gd1 = groupDropdown1.Value;
            gd2 = groupDropdown2.Value;

            % Basic validations
            if any(strcmp(gd1, '(No groups yet)')) || any(strcmp(gd2, '(No groups yet)'))
                errorLabel.Text = 'Please choose valid groups before calculating.';
                errorLabel.FontColor = 'red';  % Set text color to red
                return;
            end

            if strcmp(gd1, gd2)
                errorLabel.Text = 'Cannot compute ratio of the same group.';
                errorLabel.FontColor = 'red';  % Set text color to red
                return;
            end

            % Get the grouped data from appdata
            groupedData = getappdata(Fig, 'groupedData');
            if ~isfield(groupedData, gd1) || ~isfield(groupedData, gd2)
                errorLabel.Text = 'One or both selected groups do not exist.';
                errorLabel.FontColor = 'red';  % Set text color to red
                return;
            end

            %---------------------------
            % 1) Combine data of each group
            %---------------------------
            [commonTime1, dataMatrix1] = combineGroupData(groupedData.(gd1));
            [commonTime2, dataMatrix2] = combineGroupData(groupedData.(gd2));

            % If either group is empty, show an error
            if isempty(commonTime1) || isempty(dataMatrix1) ...
                    || isempty(commonTime2) || isempty(dataMatrix2)
                errorLabel.Text = 'One or both selected groups have no data.';
                errorLabel.FontColor = 'red';  % Set text color to red
                return;
            end

            %---------------------------
            % 2) (Optional) Check if times match
            %---------------------------
            % If your data come from the same underlying time vectors, they should match.
            % If they don't match, you need interpolation or an error message.
            if length(commonTime1)~=length(commonTime2) ...
                    || any(abs(commonTime1 - commonTime2) > 180)
                errorLabel.Text = 'Time vectors differ by more than 180 s; cannot compute ratio.';
                errorLabel.FontColor = 'red';  % Set text color to red
                return;
            end

            %---------------------------
            % 3) Compute the mean & std for each group across wells
            %    dataMatrix is [numWells x timePoints]
            %---------------------------
            mean1 = mean(dataMatrix1, 1, 'omitnan');
            std1  = std(dataMatrix1, 0, 1, 'omitnan');   % sample std

            mean2 = mean(dataMatrix2, 1, 'omitnan');
            std2  = std(dataMatrix2, 0, 1, 'omitnan');   % sample std

            %---------------------------
            % 4) Compute ratio and its approximate std via error propagation
            %---------------------------
            ratioData = mean1 ./ mean2;  % element-wise division

            % Avoid division by zero or invalid points
            badIdx = (mean2 == 0 | isnan(mean2));
            ratioData(badIdx) = 0;

            % ratioSTD(t) = ratio * sqrt( (std1/mean1)^2 + (std2/mean2)^2 )
            ratioStd = ratioData .* sqrt( (std1./mean1).^2 + (std2./mean2).^2 );
            ratioStd(badIdx) = 0;   % Set to NaN where ratio is invalid

            %---------------------------
            % 5) Plot the ratio ± std on the existing axes
            %---------------------------
            axGraph = getappdata(Fig, 'axGraph');
            if isempty(axGraph) || ~isvalid(axGraph)
                errorLabel.Text = 'Could not find valid axes to plot.';
                errorLabel.FontColor = 'red';
                return;
            end

            if grfrstime
                cla(axGraph, 'reset');
                yline(axGraph, 0, 'k--', 'LineWidth', 1, 'DisplayName', 'Reference (0)');
            end

            hold(axGraph, 'on');

            % Fill the area for ratio ± std
            upperBound = ratioData + ratioStd;
            lowerBound = ratioData - ratioStd;

            fill(axGraph, ...
                [commonTime1, fliplr(commonTime1)], ...
                [lowerBound,  fliplr(upperBound)], ...
                cmap(colorIdx, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility', 'off');

            % Plot the ratio line
            plot(axGraph, commonTime1, ratioData, ...
                'LineWidth', 2, ...
                'Color', cmap(colorIdx, :), ...
                'DisplayName', [gd1 ' / ' gd2 ' Ratio']);

            grid(axGraph, 'on');
            axGraph.XGrid = 'on';
            axGraph.YGrid = 'on';
            axGraph.GridColor = [0.7, 0.7, 0.7]; % Light gray for subtle effect
            axGraph.GridAlpha = 0.5; % Semi-transparent

            xlabel(axGraph, 'Time');
            ylabel(axGraph, ['Ratio: ' gd1 ' / ' gd2]);
            title(axGraph, 'Ratio of Group Means ± Std');
            legend(axGraph, 'show');
            hold(axGraph, 'off');
            
            disp(['Calculated ratio: ' gd1 ' / ' gd2]);
            errorLabel.Text = 'Ratio (± std) plotted successfully.';
            errorLabel.FontColor = [0, 0.5, 0];  % Dark green
            
            grfrstime = false;
            colorIdx = mod(colorIdx, nColors) + 1;

        end

        function on2DRatio()
            % Clear any old error text
            errorLabel.Text = '';
            % 2-D scatter coloured by well number (1-12)

            % ── persistent flag ─────────────────────────────────────────────
            grfrstime = true;

            % ── handles & shared data ──────────────────────────────────────
            Fig            = gcbf;                       % parent figure
            gd1            = getappdata(Fig,'groupDropdown1').Value;
            gd2            = getappdata(Fig,'groupDropdown2').Value;
            groupedData    = getappdata(Fig,'groupedData');
            axGraph        = getappdata(Fig,'axGraph');
            errH           = getappdata(Fig,'ratioErrorLabel');        % status label

            % ── validation ────────────────────────────────────────────────
            if any(strcmp({gd1 gd2},'(No groups yet)'))
                errH.Text='Please choose valid groups.'; errH.FontColor='r'; return
                errorLabel.FontColor = 'red';
            end
            if strcmp(gd1,gd2)
                errH.Text='Cannot plot the same group.'; errH.FontColor='r'; return
            end
            if ~isfield(groupedData,gd1) || ~isfield(groupedData,gd2)
                errH.Text='One or both groups do not exist.'; errH.FontColor='r'; return
            end

            % ── data & labels ─────────────────────────────────────────────
            [~,data1] = combineGroupData(groupedData.(gd1));   % rows = wells
            [~,data2] = combineGroupData(groupedData.(gd2));

            if isempty(data1) || isempty(data2)
                errH.Text='One or both groups have no data.'; errH.FontColor='r'; return
            end

            % replicate labels so they align with flattened data
            labels      = groupedData.(gd1).plotLabels(:);        % one per row
            [nRows,nCol]= size(data1);
            labMat      = repmat(labels,1,nCol);                   % rows×cols
            lab         = labMat(:);                              % column vector

            % flatten data and trim equally
            xData = data2(:);  yData = data1(:);
            nPts  = min(numel(xData),numel(yData));
            xData = xData(1:nPts);  yData = yData(1:nPts);  lab = lab(1:nPts);

            % Plot in the existing axes
            if tdfrstime
                cla(axGraph, 'reset');
            end

            % ── colour by well number ─────────────────────────────────────
            nWells  = 12;
            cmap    = hsv(nWells)*0.8;
            wellNum = ones(nPts,1);
            for k = 1:nPts
                m = regexp(lab{k},'\d+$','match','once');          % trailing digits
                if ~isempty(m)
                    w = str2double(m);
                    if w>=1 && w<=nWells, wellNum(k)=w; end
                end
            end
            ptCol = cmap(wellNum,:);

            % ── plotting ───────────────────────────────────────────────────────────

            if tdfrstime, cla(axGraph,'reset'); end
            hold(axGraph,'on');                       % <-- keep everything that follows

            % 1) main cloud
            s = scatter(axGraph, xData, yData, 30, ptCol, 'filled', 'HandleVisibility','off');

            % 1.1) drop the default "Color" row
            rows = s.DataTipTemplate.DataTipRows;
            isColor = arrayfun(@(r) strcmp(r.Label,'Color'), rows);
            s.DataTipTemplate.DataTipRows(isColor) = [];

            % 1.2) add a “Well” row with one string per point
            wellLabels = arrayfun(@(w) sprintf('%d',w), wellNum, 'UniformOutput', false);
            s.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Well', wellLabels);

            % 2) dummy points for the legend
            wellsShown = unique(wellNum);
            for k = 1:numel(wellsShown)
                w = wellsShown(k);
                scatter(axGraph,nan,nan,30,cmap(w,:),'filled', ...
                    'DisplayName',['Well ' num2str(w)]);
            end

            hold(axGraph,'off');                      % <-- finish plotting
            legend(axGraph,'show');
            xlabel(axGraph,gd2); ylabel(axGraph,gd1);
            title(axGraph,[gd1 ' vs ' gd2 ' (2-D relation)']); grid(axGraph,'on');
        end



        function clearPlot()
            cla(axGraph, 'reset');
        end

        function editLabelsCallback()
            axGraph = getappdata(Fig, 'axGraph');
            prompt = {'Enter Title:', 'Enter X Label:', 'Enter Y Label:'};
            dlgtitle = 'Edit Graph Labels';
            dims = [1 50];
            definput = {axGraph.Title.String, axGraph.XLabel.String, axGraph.YLabel.String};
            answer = inputdlg(prompt, dlgtitle, dims, definput);
        
            if ~isempty(answer)
                axGraph.Title.String = answer{1};
                axGraph.XLabel.String = answer{2};
                axGraph.YLabel.String = answer{3};
            end
        end

    end




    function integrateExcelDataAnalysisGUI(parentPanel)
        % Create a button to select Excel file
        btnSelectFile = uibutton(parentPanel, 'push', 'Position', [20, 540, 100, 20], 'Text', 'Select File', 'ButtonPushedFcn', @(btn, event) selectExcelFile());

        % Create a label to display the selected file path
        lblFilePath = uilabel(parentPanel, 'Position', [20, 510, 500, 20], 'Text', 'Selected File: None');

        % Create a label for sheet selection
        lblSheetSelect = uilabel(parentPanel, 'Position', [20, 470, 300, 20], 'Text', 'Select sheets:');

        % Create a listbox for sheet selection
        listboxSheets = uilistbox(parentPanel, 'Position', [20, 260, 200, 200], 'Items', {}, 'Multiselect', 'on');
        listboxSheets.Items = {'Select a file first'};
        listboxSheets.Multiselect = 'on';
        listboxSheets.Enable = 'off';

        % Create labels for time format dropdowns
        lblTimeFormatExcel = uilabel(parentPanel, 'Position', [20, 200, 150, 20], 'Text', 'Excel Time Format');
        lblTimeFormatGraph = uilabel(parentPanel, 'Position', [20, 135, 150, 20], 'Text', 'Graph Time Format');

        % Create dropdowns for time format selection
        ddlTimeFormatExcel = uidropdown(parentPanel, 'Position', [20, 175, 150, 20], 'Items', {'Seconds', 'Minutes', 'Hours', 'HH:MM:SS'}, 'Value', 'Seconds');
        ddlTimeFormatGraph = uidropdown(parentPanel, 'Position', [20, 110, 150, 20], 'Items', {'Seconds', 'Minutes', 'Hours'}, 'Value', 'Seconds');
        ddlTimeFormatExcel.Enable = 'off';
        ddlTimeFormatGraph.Enable = 'off';

        % Create a button to start analysis
        btnStartAnalysis = uibutton(parentPanel, 'push', 'Position', [20, 50, 100, 30], 'Text', 'Start Analysis', 'ButtonPushedFcn', @(btn, event) startAnalysis());
        btnStartAnalysis.Enable = 'off';

        % Create labels for well plates tree list
        lblWellPlatesTree = uilabel(parentPanel, 'Position', [350, 470, 100, 20], 'Text', 'Select well plates:');

        % Create listbox with preliminary tree
        tree = uitree(parentPanel, ...
            'checkbox', ...
            'Position', [350, 140, 200, 320] ...
            );
        tree.Enable = 'off';

        % Add a listbox to display Custom Groups items
        custom_tree = uitree(parentPanel, ...
            'checkbox', ...
            'Position', [700, 140, 200, 320] ...
            );

        % Add a button for moving selected checkboxes to Custom Groups
        btnMoveToGroups = uibutton(parentPanel, 'push', ...
            'Position', [600, 350, 50, 50], ... % Adjust position as per layout
            'Text', '→', ...
            'FontSize', 20, ...
            'FontWeight', 'bold', ...
            'ButtonPushedFcn', @(btn, event) moveToGroups(tree, custom_tree));
        btnMoveToGroups.Enable = 'off'; % Initially disabled

        % Add a label for the Custom Groups window
        lblCustomGroups = uilabel(parentPanel, ...
            'Position', [700, 470, 150, 20], ...
            'Text', 'Custom Groups:');

        % Moving selected nodes to Custom Groups window while creating the
        % ZXX Custom Groups
        function moveToGroups(tree, custom_tree)

            % Persistent variable to track the ZXX increment
            if isempty(parentCounter)
                parentCounter = 1; % Initialize the counter
            end

            % Ensure the counter doesn't exceed 99
            if parentCounter > 99
                error('Parent counter exceeded limit of 99.');
            end

            % Generate the parent node name
            parentName = sprintf('Z%02d', parentCounter);
            parentCounter = parentCounter + 1; % Increment the counter

            % Create a new parent node in the custom tree
            parentNode = uitreenode(custom_tree, 'Text', parentName);

            selectedNodes = tree.CheckedNodes;
            selected = getFullPath(selectedNodes);

            % If no nodes are selected, return early
            if isempty(selected)
                return;
            end

            % Add selected nodes to the Custom Groups list
            for i = 1:numel(selected)
                newNode = uitreenode(parentNode, 'Text', string(selected(i)));
            end

            tree.CheckedNodes = [];
        end

        % Create a button to auto group & start graph plotting
        btnPlotGraph = uibutton(parentPanel, 'push', 'Position', [350, 50, 200, 30], 'Text', 'Auto Group & Graph Plotting', 'ButtonPushedFcn', @(btn, event) checkSelection(tree));

        % Create a button to custom group & start graph plotting
        btnPlotGraph = uibutton(parentPanel, 'push', 'Position', [700, 50, 200, 30], 'Text', 'Custom Group & Graph Plotting', 'ButtonPushedFcn', @(btn, event) checkSelection(custom_tree));

        btnClearCustomGroups = uibutton(parentPanel, 'push', 'Position', [700, 10, 200, 30], 'Text', 'Clear Custom Groups', 'BackgroundColor', [1, 0.9, 0.9], 'ButtonPushedFcn', @(btn, event) clearCustomGroups());

        function clearCustomGroups
            custom_tree.Children.delete;
        end

        btnPlotGraph.Enable = 'off';

        % Variable to store selected file and sheets
        selectedFile = '';
        selectedSheets = {};

        function selectExcelFile() %Excel File selection

            % Prompt the user to select an Excel file
            [filename, filepath] = uigetfile({'*.xlsx;*.xls', 'Excel Files (*.xlsx, *.xls)'}, 'Select Excel File');

            % Check if the user clicked Cancel
            if isequal(filename, 0) || isequal(filepath, 0)
                custom_fprintf('User canceled the operation.');
                return;
            end

            % Read from excel file
            % Full path to the selected Excel file
            selectedFile = fullfile(filepath, filename);

            % Update the label with the selected file path
            lblFilePath.Text = ['Selected File: ', selectedFile];
            custom_fprintf('User selected: %s', selectedFile);

            % Read the sheet names from the selected file
            [~, sheetNames] = xlsfinfo(selectedFile);

            % Update the listbox with the sheet names
            listboxSheets.Items = sheetNames;
            listboxSheets.Enable = 'on';
            ddlTimeFormatExcel.Enable = 'on';
            ddlTimeFormatGraph.Enable = 'on';
            btnStartAnalysis.Enable = 'on';
        end

        function startAnalysis() %Initializing analysis
            % Read selected sheets and time formats
            selectedSheets = listboxSheets.Value;
            % Get time formats
            timeFormatExcel = ddlTimeFormatExcel.Value;
            timeFormatGraph = ddlTimeFormatGraph.Value;
            setappdata(Fig,'graphTimeUnit', timeFormatGraph);   % Seconds | Minutes | Hours

            % Check if the user made a selection
            if isempty(selectedSheets)
                uialert(Fig, 'Please select at least one sheet.', 'Error', 'Icon', 'error');
                return;
            end

            % Call the analysis function
            analyzeData(selectedFile, selectedSheets, timeFormatExcel, timeFormatGraph);
        end

        % ------Here there should be an option for custom grouping------%
        % ==============================================================%
        function selected_checkbox = getFullPath(checknodes) %Getting full paths for selected nodes by the user (from tree GUI)
            selected_checkbox = {};
            for i = 1:length(checknodes)
                node = checknodes(i);
                % Initialize the full path with the node text
                fullPath = node.Text;
                % Traverse up the parent chain to get the full path
                while ~isempty(node.Parent) && ~isa(node.Parent, 'matlab.ui.container.CheckBoxTree')
                    node = node.Parent;
                    fullText = [node.Text '_' fullPath];
                    selected_checkbox = [selected_checkbox, fullText];
                end
            end
        end

        %Here should be a CUSTOM grouping function (allowing multiple groups) that its outpout will go to
        %checkSelection(tree). In this function, for each iteration where
        %the user select nodes - a custom group will be created with a designeted ZXX prefix, where XX spans from 01 to 99%

        function checkSelection(tree) %Breaking down data from selected sheets including Wells auto group creation, and drawing the first default graph

            % Retrieve axGraph from appdata
            axGraph = getappdata(Fig, 'axGraph');

            selectedNodes = tree.CheckedNodes;
            selected = getFullPath(selectedNodes);

            % Initialize variables to hold the data to be plotted
            plotData = {};
            plotTime = {};
            plotLabels = {};

            % Loop through selected sheets to get data and time for each well plate
            %
            % In more details: Breaking down the selected sheets to sheets & well plates,
            % and extracting all the data from the sheet as distinct
            % values
            %

            for idx = 1:length(selected)
                fullPath = selected{idx};
                splitStr = strsplit(fullPath, '_');
                if ~(length(splitStr) == 2 || length(splitStr) == 3)
                    continue;
                end
                if length(splitStr) == 2
                    sheetName = splitStr{1};
                    wellPlateName = splitStr{2};
                end
                if length(splitStr) == 3
                    sheetName = splitStr{2};
                    wellPlateName = splitStr{3};
                end


                % Find the sheet index in the selected sheets
                sheetIdx = find(strcmp(selectedSheets, sheetName));
                if isempty(sheetIdx)
                    continue;
                end

                % Find well plate index within the sheet
                wellPlateNames = normalizedWellPlate{sheetIdx};
                wellPlateIdx = find(strcmp(wellPlateNames, wellPlateName));

                if isempty(wellPlateIdx) || wellPlateIdx > size(data{sheetIdx},2)
                    % Create a dummy NaN column rather than skipping
                    dataVec = nan(size(data{sheetIdx},1),1);
                else
                    dataVec = data{sheetIdx}(:, wellPlateIdx);
                end

                timeVec = time{sheetIdx};

                % Add data to the list for plotting
                plotData{end+1} = dataVec;
                plotTime{end+1} = timeVec;
                plotLabels{end+1} = fullPath;

            end

            % ===== Grouping Well Plates Process =====

            % For further exploration, examine to possibility to create a ZXX
            % artificial prefix for custom grouping and then use 'prefix =
            % wellPlateName(1:3)'

            groupedData = struct();

            for i = 1:length(plotLabels)
                % Each plotLabels{i} is something like "SheetName_WellPlateName"
                % Extract "WellPlateName" (second part)
                parts = strsplit(plotLabels{i}, '_');
                if ~(length(parts) == 2 || length(parts) == 3)
                    continue;  % Skip anything unexpected
                end
                if length(splitStr) == 2
                    customGroupName = parts{2};
                    prefix = customGroupName(1:4);
                end
                if length(splitStr) == 3
                    customGroupName = parts{1};
                    prefix = customGroupName;
                end

                % Initialize subfields if this prefix does not exist yet
                if ~isfield(groupedData, prefix)
                    groupedData.(prefix).plotData   = {};
                    groupedData.(prefix).plotTime   = {};
                    groupedData.(prefix).plotLabels = {};
                end

                % Append the data/time/labels to the correct prefix
                groupedData.(prefix).plotData{end+1}   = plotData{i};
                groupedData.(prefix).plotTime{end+1}   = plotTime{i};
                groupedData.(prefix).plotLabels{end+1} = plotLabels{i};

                % Store groupedData in appdata
                setappdata(Fig, 'groupedData', groupedData);

                % Update dropdown items with the group names (prefixes)
                prefixes = fieldnames(groupedData);
                groupDropdown = getappdata(Fig, 'groupDropdown');

                if ~isempty(groupDropdown) && isvalid(groupDropdown)
                    if isempty(prefixes)
                        groupDropdown.Items = {'All Groups', '(No groups available)'};
                        groupDropdown.Value = 'All Groups';
                        groupDropdown.Enable = 'off';  % Disable if no groups
                    else
                        % Prepend "All Groups" to the list
                        groupDropdown.Items = ['All Groups'; prefixes];
                        groupDropdown.Value = 'All Groups';  % Set default to "All Groups"
                        groupDropdown.Enable = 'on';
                    end
                end
            end

            % --- populate time-point selector ----------------------
            lstTimePt = getappdata(Fig,'lstTimePt');
            % Use the first group's time vector as reference
            anyField  = fieldnames(groupedData);
            refTime   = groupedData.(anyField{1}).plotTime{1};
            lstTimePt.Items = cellstr(num2str(refTime(:)));
            lstTimePt.Value = lstTimePt.Items(end);    % default = last time
            % ------------------------------------------------------------

            populateRatioDropdowns();

            populateGroupList();

            % Plot first default data on axGraph (All nodes for all groups)
            cla(axGraph, 'reset');          % Reset axes properties
            delete(axGraph.Children);       % Remove any lingering plots
            hold(axGraph, 'on');  % Hold on to allow multiple plots

            prefixes = fieldnames(groupedData);
            colorMap = lines(numel(prefixes));  % e.g., distinct colors for each prefix

            for p = 1:numel(prefixes)
                prefix = prefixes{p};
                color_ = colorMap(p, :);
                % For each well in this prefix
                for k = 1:length(groupedData.(prefix).plotData)
                    if k == 1
                        % First plot in this group -> visible in legend
                        plot(axGraph, ...
                            groupedData.(prefix).plotTime{k}, ...
                            groupedData.(prefix).plotData{k}, ...
                            'DisplayName', prefix, ...      % Display group name only
                            'Color', color_);
                    else
                        % Subsequent plots in this group -> hidden from legend
                        plot(axGraph, ...
                            groupedData.(prefix).plotTime{k}, ...
                            groupedData.(prefix).plotData{k}, ...
                            'HandleVisibility','off', ...   % Hide from legend
                            'Color', color_);
                    end
                end
            end

            hold(axGraph, 'off');  % Release hold
            legend(axGraph, 'show');
            xlabel(axGraph, 'Time');
            ylabel(axGraph, 'Data Value');
        end

        % Analysis and plotting
        function analyzeData(selectedFile, selectedSheets, timeFormatExcel, timeFormatGraph)

            % Read data from the selected sheets
            dataFromExcel = {};
            textFromExcel = {};

            % Create and configure the blinking timer
            LED.Color = [0.8, 0.8, 0.8]; % Initialize LED color (gray)
            t = timer('ExecutionMode', 'fixedRate', ... % Run periodically
                'Period', 0.1, ...               % Interval in seconds
                'TimerFcn', @(~,~) blinkLED(LED)); % Function to run

            start(t); % Start the blinking timer




            for i = 1:length(selectedSheets) % Extracting Numeric & Textual data
                sheetName = selectedSheets{i};
                custom_fprintf(['Reading data from sheet: ', sheetName]);
                drawnow; % Update GUI

                try
                    % -- MANDATORY LINES (unchanged) -------------------------
                    [dataMatrix, textData] = xlsread(selectedFile, sheetName);
                    dataFromExcel{i} = dataMatrix;
                    textFromExcel{i} = textData;
                    % --------------------------------------------------------

                    % Read the full raw data (numbers + text in a single cell array)
                    [~, ~, rawData] = xlsread(selectedFile, sheetName);

                    % Use rawData's dimensions for alignment
                    [numRows, numCols] = size(rawData);

                    % Initialize a NaN matrix of the same size
                    newDataMatrix = NaN(numRows, numCols);

                    % Fill newDataMatrix where rawData contains numeric values
                    for row = 1:numRows
                        for col = 1:numCols
                            if isnumeric(rawData{row, col}) && ~isempty(rawData{row, col})
                                newDataMatrix(row, col) = rawData{row, col};
                            end
                        end
                    end

                    % Overwrite the original dataMatrix with the aligned version
                    dataFromExcel{i} = newDataMatrix;

                catch ME
                    custom_fprintf('Error reading the Excel file:');
                    custom_fprintf(ME.message);
                    stop(t);
                    delete(t); % Stop and clean up the timer on error
                    LED.Color = 'red'; % Indicate error
                    drawnow;
                    return;
                end
            end



            function matrixPositions = FindMatrixPositions(textFromExcel)
                matrixPositions = {}; % Store row & column indices of matrices
                [numRows, numCols] = size(textFromExcel);

                for row = 1:numRows
                    for col = 1:numCols
                        if contains(lower(string(textFromExcel{row, col})), 'time') % Look for "time"
                            %Check within 2 cells (right & below) if "A1" exists
                            if (col+2 <= numCols && contains(lower(string(textFromExcel{row, col+2})), 'a1')) || ...
                                    (col+1 <= numCols && contains(lower(string(textFromExcel{row, col+1})), 'a1')) || ...
                                    (row+2 <= numRows && contains(lower(string(textFromExcel{row+2, col})), 'a1')) || ...
                                    (row+1 <= numRows && contains(lower(string(textFromExcel{row+1, col})), 'a1'))

                                % Store matrix start position
                                matrixPositions{end+1} = [row, col];
                            end

                        end
                    end
                end
            end

            iCounter = 1; % Counter for storing matrices sequentially

            for sheetIdx = 1:length(selectedSheets)
                custom_fprintf(['Processing data from sheet: ', selectedSheets{sheetIdx}]);
                drawnow; % Update GUI

                % Locate all matrices in the sheet
                matrixPositions = FindMatrixPositions(textFromExcel{sheetIdx});

                for posIdx = 1:length(matrixPositions)  % Process each detected matrix sequentially
                    startRow = matrixPositions{posIdx}(1);
                    startCol = matrixPositions{posIdx}(2);

                    % 0) Grab the “everything from here to end” block
                    rawTextBlock = textFromExcel{sheetIdx}(startRow:end,   startCol:end);
                    rawDataBlock = dataFromExcel{sheetIdx}(startRow:end, startCol:end);

                    % ---------------------------------------------------------------
                    % 1) Identify how far down/right the matrix really goes
                    textData_type{iCounter} = DetermineVectorType(rawTextBlock);
                    endRowLocal = size(rawDataBlock, 1);
                    endColLocal = size(rawTextBlock, 2);

                    if textData_type{iCounter} == 0
                        % Horizontal time: stop at first empty header‐cell in column 1
                        for r = 2:endRowLocal
                            if isempty(rawTextBlock{r, 1})
                            
                                endRowLocal = r - 1;
                                break;
                            end
                        end
                    else
                        % Vertical time: stop at first empty header‐cell in row 1
                        for c = 2:endColLocal
                            if isempty(rawTextBlock{1, c})
                                endColLocal = c - 1;
                                break;
                            end
                        end

                        for r = 2:endRowLocal               % scan downward
                            if isnan(rawDataBlock(r,1))     % first value in column 1 is NaN ⇒ blank row
                                endRowLocal = r-1;          % last valid data row
                                break;
                            end
                        end
                    end

                    % 2) Truncate *only* the relevant dimension
                    if textData_type{iCounter} == 0
                        % Horizontal-time: keep rows 1…endRowLocal, all columns
                        matrixText = rawTextBlock(1:endRowLocal,  :);
                        matrixData = rawDataBlock(1:endRowLocal, :);
                    else
                        % Vertical-time: keep all rows, columns 1…endColLocal
                        matrixText = rawTextBlock(1, 1:endColLocal);
                        matrixData = rawDataBlock(1:endRowLocal, 1:endColLocal);
                    end

                    % 3) Remove any fully empty rows or columns of text
                    matrixText(all(cellfun(@isempty, matrixText), 2), :) = [];
                    matrixText(:, all(cellfun(@isempty, matrixText), 1)) = [];
                    % ---------------------------------------------------------------

                    % Reduce matrixText to just the header row/column
                    if textData_type{iCounter}  % vertical time
                        matrixText = matrixText(1,:);
                    else                       % horizontal time
                        matrixText = matrixText(:,1);
                    end

                    % 4) Re-determine orientation now that headers are clean
                    textData_type{iCounter} = checkVectorType(matrixText);

                    % 5) Drop rows/columns of pure NaN from the numeric block
                    matrixData(all(isnan(matrixData), 2), :) = [];
                    matrixData(:, all(isnan(matrixData), 1)) = [];

                    matrixData(isnan(matrixData)) = 0;

                    % 6) Extract & normalize
                    timeVector{iCounter}           = TimeVecExt(matrixData,   textData_type{iCounter});
                    normalizedDataMatrix{iCounter} = NormalizedData(matrixData, textData_type{iCounter});
                    normalizedWellPlate{iCounter}  = NormalizedText(matrixText, textData_type{iCounter}, iCounter);

                    sheetIndexForMatrix(iCounter) = sheetIdx;

                    % 7) Trim trailing zeros
                    [timeVectorTrailed{iCounter}, dataMatrixTrailed{iCounter}] = ...
                        TrailZero(timeVector{iCounter}, normalizedDataMatrix{iCounter}, textData_type{iCounter});

                    % 8) Convert time format
                    time{iCounter} = TimeVectorConvert(timeFormatExcel, timeFormatGraph, timeVectorTrailed{iCounter});

                    custom_fprintf(['Finished processing matrix ', num2str(iCounter), ...
                        ' from sheet: ', selectedSheets{sheetIdx}]);
                    drawnow;

                    % Increment the matrix counter
                    iCounter = iCounter + 1;
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 1) Find all matrix indices belonging to this sheet
                idxForThisSheet = find(sheetIndexForMatrix == sheetIdx);

                % 2) Combine all well labels into one cell array
                allLabels = {};
                for m = 1:length(idxForThisSheet)
                    matrixNumber = idxForThisSheet(m);
                    % Append the well labels from that matrix
                    allLabels = [allLabels, normalizedWellPlate{matrixNumber}];
                end

                % 3) Overwrite normalizedWellPlate{sheetIdx} so it has *all* well labels from all matrices in this sheet
                normalizedWellPlate{sheetIdx} = allLabels;

                % 4) If you also need to combine numeric data so that wellPlateIdx lines up, do the same for dataMatrixTrailed:
                % Ensure every block has the same number of rows by padding to the longest
                rows = cellfun(@(x) size(x,1), dataMatrixTrailed(idxForThisSheet));
                M    = max(rows);
                allData = [];
                for m = idxForThisSheet
                    mat = dataMatrixTrailed{m};
                    if size(mat,1) < M
                        mat(end+1:M, :) = 0;    % pad shorter series with 0
                    end
                    allData = [allData, mat];
                end
                dataMatrixTrailed{sheetIdx} = allData;
                idxForThisSheet = find(sheetIndexForMatrix==sheetIdx);
                time{sheetIdx}    = time{ idxForThisSheet(1) };
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Enable Graph Editing UI components
                bgGraphType = getappdata(Fig, 'bgGraphType');
                bgGraphType.Enable = 'on';
                axGraph.Visible = 'on';
            end

            % Stop and delete the blinking timer
            stop(t);
            delete(t);
            LED.Color = 'green';
            drawnow;

            % Final output from data analysis
            data = dataMatrixTrailed;

            % Create checkbox tree for well plates selection

            % Delete existing nodes
            delete(tree.Children);

            treeSheetNodes = {};
            wellPlateNodes = {};
            for s = 1:length(selectedSheets)
                % Create a node for this sheet
                treeSheetNodes{s} = uitreenode(tree, 'Text', selectedSheets{s});

                % Use allLabels directly since it already contains all well labels from all matrices in this sheet
                wellsForSheet = normalizedWellPlate{s};

                for w = 1:length(wellsForSheet)
                    wellPlateNodes{s}{w} = uitreenode(treeSheetNodes{s}, 'Text', wellsForSheet{w});
                end
            end

            % Expand the tree
            expand(tree);

            tree.Enable = 'on';
            btnPlotGraph.Enable = 'on';
            btnMoveToGroups.Enable = 'on';

            function vec_type = checkVectorType(vec)
                % Get the size of the vector
                [rows, cols] = size(vec);

                % Check if the vector is a column vector
                if cols == 1
                    vec_type = true;
                    custom_fprintf('The text vector is a column vector, therefore time is row vector.');
                    % Check if the vector is a row vector
                elseif rows == 1
                    vec_type = false;
                    custom_fprintf('The text vector is a row vector, therefore time is column vector and each well plates samples are vertical.');
                    % If neither, then it's not a 1-dimensional vector
                else
                     % Fallback: count text vs numeric in first row/column
                     vec_type = DetermineVectorType(vec);
                    % custom_fprintf('The input is neither a column vector nor a row vector.');
                    % stop(t);
                    % delete(t); % Stop and clean up the timer on error
                    % LED.Color = 'red'; % Indicate error
                end
            end

            function gvecType = DetermineVectorType(matrixText)
                % Get number of rows and columns
                [numRows, numCols] = size(matrixText);

                % Extract the first row and first column
                firstRow = matrixText(1, :);
                firstCol = matrixText(:, 1);

                % Count how many are text (ignoring empty cells)
                textCountRow = sum(cellfun(@(x) ischar(x) && ~isempty(x), firstRow));
                textCountCol = sum(cellfun(@(x) ischar(x) && ~isempty(x), firstCol));

                % Count how many are numbers
                numCountRow = sum(cellfun(@(x) isnumeric(x) && ~isnan(x), firstRow));
                numCountCol = sum(cellfun(@(x) isnumeric(x) && ~isnan(x), firstCol));

                % Determine vector type based on majority
                if textCountRow == 1
                    gvecType = false;  % Horizontal time
                elseif textCountRow > 1
                    gvecType = true; % Vertical time
                else
                    % Default to horizontal if there's a tie (this can be adjusted)
                    gvecType = false;
                end
            end

            function TimeVec = TimeVecExt(Data, TextType)
                if TextType
                    TimeVec = Data(1, :); % Extract time as row vector
                else
                    TimeVec = Data(:, 1)'; % Extract time as column vector and transpose to row vector
                end
                custom_fprintf('Time vector extracted successfully.');
            end

            function NormVec = NormalizedData(Data, TextType)
                if TextType % Normalize to vertical data
                    NormVec = Data(2:end, :)'; % Remove the first row and normalize matrix using transpose operation
                else % keep as vertical data
                    NormVec = Data(:, 2:end); % Remove the first column
                end
                custom_fprintf('Well plates matrix extracted without time vector and normalized to vertical well plates samples data.');
            end

            function NormVec = NormalizedText(Vec, TextType, matrixIndex)
                if TextType % Normalize to Horizontal data
                    NormVec = (Vec(2:end))'; % Remove the first row and normalize using transpose operation
                else % Keep as Horizontal data
                    NormVec = Vec(2:end); % Remove the first column

                end

                % Add "MXY" prefix to each text label
                for k = 1:length(NormVec)
                    if ischar(NormVec{k}) % Only modify text labels
                        NormVec{k} = ['MX', num2str(matrixIndex), '', NormVec{k}];
                    end
                end

                custom_fprintf(['Well plates matrix extracted and prefixed with MXY: MX', num2str(matrixIndex)]);
            end



            %Trailing Zero
            function [TimeVecOut, DataVecOut] = TrailZero(TimeVecIn, DataVecIn, TextType)
                % Find the last non-zero element in the time vector
                idx = find(TimeVecIn, 1, 'last');

                if isempty(idx)
                    % No valid time points
                    TimeVecOut = [];
                    DataVecOut = [];
                    return;
                end

                % Trim the time vector up to idx
                TimeVecOut = TimeVecIn(1:idx);

                % If the time vector is vertical (TextType == false),
                % each row of DataVecIn corresponds to a time point => trim rows.
                % If the time vector is horizontal (TextType == true),
                % leave DataVecIn untouched (no row- or column-trimming).
                if ~TextType
                    % Vertical time => trim rows
                    DataVecOut = DataVecIn(1:idx, :);
                else
                    % Horizontal time => do not trim data matrix
                    DataVecOut = DataVecIn;
                end
            end


            function TimeVecOut = TimeVectorConvert(TimeInFormat, TimeOutFormat, TimeVecIn)
                % Convert the input time vector to seconds
                switch TimeInFormat
                    case 'Seconds'
                        TimeInSeconds = TimeVecIn;
                    case 'Minutes'
                        TimeInSeconds = TimeVecIn * 60;
                    case 'Hours'
                        TimeInSeconds = TimeVecIn * 3600;
                    case 'HH:MM:SS'
                        % 1.25 days = 1 day + 6 hours = 30 hours total => 108,000 seconds
                        TimeInSeconds = TimeVecIn * 86400;
                    otherwise
                        error('Invalid TimeInFormat');
                end

                % Convert from seconds to the output format
                switch TimeOutFormat
                    case 'Seconds'
                        TimeVecOut = TimeInSeconds;
                    case 'Minutes'
                        TimeVecOut = TimeInSeconds / 60;
                    case 'Hours'
                        TimeVecOut = TimeInSeconds / 3600;
                    otherwise
                        error('Invalid TimeOutFormat');
                end
                custom_fprintf(['Time Vector Converted from ', TimeInFormat, ' to ', TimeOutFormat]);
            end
        end


    end

%-------------------------------------------
% Helper: fill the ratio dropdowns
%-------------------------------------------
    function populateRatioDropdowns()
        groupedData = getappdata(Fig, 'groupedData');
        if isempty(groupedData)
            return;
        end

        groupDropdown1 = getappdata(Fig, 'groupDropdown1');
        groupDropdown2 = getappdata(Fig, 'groupDropdown2');

        groupNames = fieldnames(groupedData);
        % Exclude "All Groups" from ratio combos
        maskAll = strcmp(groupNames, 'All Groups');
        groupNames(maskAll) = [];

        if isempty(groupNames)
            groupDropdown1.Items = {'(No groups yet)'};
            groupDropdown2.Items = {'(No groups yet)'};
            groupDropdown1.Enable = 'off';
            groupDropdown2.Enable = 'off';
        else
            groupDropdown1.Items = groupNames;
            groupDropdown2.Items = groupNames;
            groupDropdown1.Value = groupNames{1};
            groupDropdown2.Value = groupNames{1};
            groupDropdown1.Enable = 'on';
            groupDropdown2.Enable = 'on';
        end
    end

% Custom fprintf function to include a timestamp and update the message box
    function custom_fprintf(varargin)
        % --- Use a global variable for txtMessages ---
        global gTxtMessages

        % Build the message with a timestamp
        timestamp = datestr(now, 'mm/dd/yy HH:MM:SS.FFF');
        message = ['[' timestamp '] ' sprintf(varargin{1}, varargin{2:end})];

        % Retrieve or initialize allMessages
        allMessages = getappdata(Fig, 'allMessages');
        if isempty(allMessages)
            allMessages = {};
        end

        % Append the new message
        allMessages{end+1} = message;
        setappdata(Fig, 'allMessages', allMessages);

        % Update LogMessage text
        LogMessage.Text = message;

        % Optionally change color if it's an error
        if contains(message, 'Error', 'IgnoreCase', true)
            % For App Designer uilabel
            LogMessage.FontColor = [1, 0, 0];
        else
            LogMessage.FontColor = [0, 0, 0];
        end

        % --- Update txtMessages using the global handle ---
        % Only if it's not empty and still a valid UI control
        if ~isempty(gTxtMessages) && isgraphics(gTxtMessages)
            gTxtMessages.Value = allMessages;
            drawnow; % Refresh GUI
        end
    end

end
