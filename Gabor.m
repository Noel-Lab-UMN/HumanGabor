%% -------------------------------------------------------------------------------------Description----------------------------------------------------------------------------------------------
% This experiment focuses on how humans update their expectations for sensory information 
% Gabor stimuli with different contrast values are used, the background also has different contrast levels
% The Gabor display contrast is integrated by the Gabor contrast and background contrast 
% The possibilities of the location of the Gabor on the left and right are 50:50, 80:20, 20:80
% function createGaborSegStreamB is used to generate three unique streams that contain different segments of trials that have different probabilities
% For the first 6 blocks, stream 1 and stream 2 will be randomly assigned to each block
% For the last 3 blocks, stream 3 will be assigned 
% The probability of the locations of first 20 trials of each stream is 50:50 
% Participants should take a break after each block, and we encourage them to blink, close and rest their eyes during the break
% Be sure to check the distance between the monitor and participants' eyes
% Change the Gabor size, intervals, N trials, N blocks, Gabor contrast, background contrast as needed
% Eye tracking and EEG is used for the experiment, check the code for the EEG and Tobii system to make sure they work successfully 

%% ------------------------------------------------------------------------------------Edit History-----------------------------------------------------------------------------------------------


%% --------------------------------------------------------------------Participant Information and Setup---------------------------------------------------------------------------------

clear;
close all;
sca;

% Get the subject file
[protfile, pname] = uigetfile('Z:\6-human ASD Gabor\3-documentation\*.xl*', 'Pick a protocol file');
[xlst, xlshts, xlsfmt] = xlsfinfo([pname, protfile]);

% Select subject to create a folder for the data
SubjMenu = figure('Name','', ...
                  'NumberTitle','off', ...
                  'MenuBar','none', 'Toolbar','none', ...
                  'Position', [20,  500,  200, 200 ], ...
                  'Resize', 'off'); 

subjtxt = uicontrol('Parent', SubjMenu, ...
                    'Style', 'text', ...
                    'String', 'Type Subject Number', 'FontSize', 10, ...
                    'BackgroundColor', [.8 .8 .8], ...
                    'Units','normalized', ...
                    'Position', [.1 .6 .8 .3]);

% Loop for subject selection
sdupflg = 2;
while sdupflg == 2
    subjflg = 1;
    while subjflg == 1
        subjstr = '';
        subjedit = uicontrol('Parent', SubjMenu, ...
                             'Style', 'edit', ...
                             'String', subjstr, 'FontSize', 14, ...
                             'HorizontalAlignment', 'left', ...
                             'Units','normalized', 'Position', [.1 .1 .8 .4], ...
                             'BackgroundColor', [1 1 1], ...
                             'Callback', 'uiresume(gcbf)');
        uicontrol(subjedit);
        uiwait(gcf);
        subjstr = get(subjedit, 'String');
        subjid = str2double(subjstr);

        % Validate subject ID
        if isempty(subjid) || isnan(subjid) || floor(subjid) ~= subjid
            subjflg = 1;
            errordlg('Please enter a valid integer for the Subject ID.', 'Invalid Input');
        else
            nshts = false;
            for i = 1:length(xlshts)
                nshts = nshts || strcmpi(xlshts{i}, ['S', num2str(subjid)]);
            end
            if subjid < 1 || ~nshts
                subjflg = 1;
                errordlg('Subject ID not found in the protocol file.', 'Invalid Subject ID');
            else
                subjflg = 0;
            end
        end
    end

    % Check for data folder for subjid and create participant folder
    participantFolder = [pname, 'S', num2str(subjid)];  % Define participant folder
    if isempty(dir(participantFolder))
        mkdir(participantFolder);
        sdupflg = 1;
    else
        sdupflg = menu(['A folder already exists for subject ID# ' num2str(subjid)], 'OK', 'Select again');
    end
end
delete(SubjMenu);

% Load subject protocol sheet
sheetname = ['S', num2str(subjid)];
[~, ~, protocolData] = readtable([pname, protfile], sheetname);

% Setup Psychtoolbox defaults
PsychDefaultSetup(2);

% Set up parallel card as trigger port
ioObj = io64;
status = io64(ioObj);
address = hex2dec('3FF8');

% Set up Tobii eye tracker
Tobii = EyeTrackingOperations();

eyetracker_address = 'tet-tcp://169.254.6.40';

eyetracker = Tobii.get_eyetracker(eyetracker_address);

if isa(eyetracker,'EyeTracker')
    disp(['Address:',eyetracker.Address]);
    disp(['Name:',eyetracker.Name]);
    disp(['Serial Number:',eyetracker.SerialNumber]);
    disp(['Model:',eyetracker.Model]);
    disp(['Firmware Version:',eyetracker.FirmwareVersion]);
    disp(['Runtime Version:',eyetracker.RuntimeVersion]);
else
    disp('Eye tracker not found!');
end

% Screen set up
screens = Screen('Screens');
screenNumber = max(screens);
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, 0.5);
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
screen_pixels = [screenXpixels screenYpixels];
[xCenter, yCenter] = RectCenter(windowRect);
Screen('TextSize', window, 36);

% Define colors as RGB triplets in the [0, 1] range
white = [1, 1, 1];
black = [0, 0, 0];
gray = [0.5, 0.5, 0.5];
red = [1, 0, 0];
green = [0, 1, 0];

% Define keys
KbName('UnifyKeyNames');
left_key = KbName('LeftArrow');
right_key = KbName('RightArrow');
escape_key = KbName('ESCAPE');

HideCursor;
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);

spaceKey = KbName('Space');
RKey = KbName('R');
% Calibration setup
dotSizePix = 30;

dotColor = [[1 0 0];[1 1 1]]; % Red and white

leftColor = [1 0 0]; % Red
rightColor = [0 0 1]; % Blue

% Calibration points
lb = 0.1;  % left bound
xc = 0.5;  % horizontal center
rb = 0.9;  % right bound
ub = 0.1;  % upper bound
yc = 0.5;  % vertical center
bb = 0.9;  % bottom bound

points_to_calibrate = [[lb,ub];[rb,ub];[xc,yc];[lb,bb];[rb,bb]];

% Create calibration object
calib = ScreenBasedCalibration(eyetracker);

calibrating = true;

while calibrating
    % Enter calibration mode
    calib.enter_calibration_mode();

    for i=1:length(points_to_calibrate)

        Screen('DrawDots', window, points_to_calibrate(i,:).*screen_pixels, dotSizePix, dotColor(1,:), [], 2);
        Screen('DrawDots', window, points_to_calibrate(i,:).*screen_pixels, dotSizePix*0.5, dotColor(2,:), [], 2);

        Screen('Flip', window);

        % Wait a moment to allow the user to focus on the point
        pause(1);

        if calib.collect_data(points_to_calibrate(i,:)) ~= CalibrationStatus.Success
            % Try again if it didn't go well the first time.
            % Not all eye tracker models will fail at this point, but instead fail on ComputeAndApply.
            calib.collect_data(points_to_calibrate(i,:));
        end

    end

    DrawFormattedText(window, 'Calculating calibration result....', 'center', 'center', white);

    Screen('Flip', window);

    % Blocking call that returns the calibration result
    calibration_result = calib.compute_and_apply();

    calib.leave_calibration_mode();

    if calibration_result.Status ~= CalibrationStatus.Success
        break
    end

    % Calibration Result

    points = calibration_result.CalibrationPoints;

    for i=1:length(points)
        Screen('DrawDots', window, points(i).PositionOnDisplayArea.*screen_pixels, dotSizePix*0.5, dotColor(2,:), [], 2);
        for j=1:length(points(i).RightEye)
            if points(i).LeftEye(j).Validity == CalibrationEyeValidity.ValidAndUsed
                Screen('DrawDots', window, points(i).LeftEye(j).PositionOnDisplayArea.*screen_pixels, dotSizePix*0.3, leftColor, [], 2);
                Screen('DrawLines', window, ([points(i).LeftEye(j).PositionOnDisplayArea; points(i).PositionOnDisplayArea].*screen_pixels)', 2, leftColor, [0 0], 2);
            end
            if points(i).RightEye(j).Validity == CalibrationEyeValidity.ValidAndUsed
                Screen('DrawDots', window, points(i).RightEye(j).PositionOnDisplayArea.*screen_pixels, dotSizePix*0.3, rightColor, [], 2);
                Screen('DrawLines', window, ([points(i).RightEye(j).PositionOnDisplayArea; points(i).PositionOnDisplayArea].*screen_pixels)', 2, rightColor, [0 0], 2);
            end
        end

    end

    DrawFormattedText(window, 'Press the ''R'' key to recalibrate or ''Space'' to continue....', 'center', screenYpixels * 0.95, white)

    Screen('Flip', window);

    while 1.
        [ keyIsDown, seconds, keyCode ] = KbCheck;
        keyCode = find(keyCode, 1);

        if keyIsDown
            if keyCode == spaceKey
                calibrating = false;
                break;
            elseif keyCode == RKey
                break;
            end
            KbReleaseWait;
        end
    end
end

% Set random seed 
rng('shuffle');

% Fixation dot
dotSize = 10;  

% Filename for testing responses in the participant's folder
filename = fullfile(participantFolder, ['S', num2str(subjid), '_experiment_responses.xlsx']);

%% ----------------------------------------------------------------------------------Experiment Parameters---------------------------------------------------------------------------------
% Number of blocks and trials
blocks_n = 9;                  
trialsPerBlock = 300;        
totalTrials = blocks_n * trialsPerBlock;  
stimOnsetTimes = nan(totalTrials, 1); 
responseTimes = nan(totalTrials, 1);

% Background contrasts
bc_options = [0, 0.25, 0.5]; 

% Gabor contrasts
gc0 = [0, 0.005, 0.01, 0.015, 0.02, 0.025];
gc = [0, 0.01, 0.02, 0.04, 0.06, 0.1];
gcN = length(gc);

% Spatial frequency in cycles per degree
spatial_freq = 0.5; 

% Screen and viewing parameters for degrees to pixels conversion
screenWidthCm = 53;
viewingDistanceCm = 51;

% Calculate pixels per cm
pixelsPerCm = screenXpixels / screenWidthCm;

% Calculate cm per degree
cmPerDegree = 2 * viewingDistanceCm * tan(deg2rad(1)); 

% Calculate pixels per degree
pixelsPerDegree = pixelsPerCm * cmPerDegree;

% Gabor positions in degrees (-6 for left and +6 for right)
gabor_positions = [-6, 6];

% Convert degrees to pixels for both positions
gp_pix = gabor_positions * pixelsPerDegree;

% Random delay between trials (700 ms to 1000 ms)
interval = 0.7 + (0.3 * rand());

% Write data in excel
allTrialData = cell(totalTrials, 15);  
headers = {'Trial', 'Block', 'StreamNumber', 'Background Contrast', 'Gabor Contrast', 'Gabor Position', 'Response', 'Correct', 'Reaction Time (S)', 'Stimulus Onset Time', 'Response Time', 'Probability Left', 'Probability Right', 'Trial Start Time', 'Trial End Time'};

try
    headerTable = cell2table(cell(0, length(headers)), 'VariableNames', headers);
    writetable(headerTable, filename, 'WriteVariableNames', true);
catch ME
    warning('Could not write headers to Excel file. Ensure the file is not open elsewhere.');
    disp(ME.message);
end

% Trial counter
trialCounter = 0;  

% Initialize cell array to store eye tracking data for all trials
allEyeTrackingData = cell(totalTrials, 1);  % Each cell will hold data for one trial

%% --------------------------------------------------------------------------Stream Generation and Assignment-----------------------------------------------------------------------
% Use the createGaborSegStreamB function
streams = createGaborSegStreamB();

% Create all unique combinations
stream12 = [1; 2];
combos12 = [];
for s = 1:length(stream12)
    for bc = 1:length(bc_options)
        combos12 = [combos12; stream12(s), bc_options(bc)];
    end
end

combos12 = combos12(randperm(size(combos12,1)), :);

combos3 = [];
for bc = 1:length(bc_options)
    combos3 = [combos3; 3, bc_options(bc)];
end
% Shuffle the last 3 combos
combos3 = combos3(randperm(size(combos3,1)), :);

% Concatenate first 6 combos and last 3 combos
blockAssignments = [combos12; combos3];

% Sanity-check that blockAssignments has size 9 x 2 (matching blocks_n=9)
if size(blockAssignments,1) ~= blocks_n
    error('Number of blocks does not match the required 9.');
end

% Now blockAssignments(b,1) is the stream, and blockAssignments(b,2) is the background contrast for block b.
% We can proceed to create the trials for each block as before.

trialBlocks = cell(blocks_n, 1);

% Define number of Gabor contrasts and repetitions
repsPerGabor = trialsPerBlock / gcN;

% Generate trials for each block
for b = 1:blocks_n
    currentStream = blockAssignments(b, 1);
    current_bc = blockAssignments(b, 2);   
    
    % Choose the appropriate Gabor contrast list based on background contrast
    if current_bc == 0
        current_gc_list = gc0;
    elseif current_bc == 0.25 || current_bc == 0.5
        current_gc_list = gc;
    else
        error('Invalid background contrast value.');
    end

    gcN = length(current_gc_list);
    repsPerGabor = trialsPerBlock / gcN;

    % Check if trialsPerBlock is divisible by the number of Gabor contrasts
    if mod(trialsPerBlock, gcN) ~= 0
        error('Number of trials per block (%d) is not divisible by the number of Gabor contrasts (%d).', trialsPerBlock, gcN);
    end

    % Create a vector of Gabor contrasts repeated appropriately
    gc_vector = repmat(current_gc_list', repsPerGabor, 1);

    % Ensure the first trial has a non-zero Gabor contrast
    nonZeroContrasts = current_gc_list(current_gc_list > 0);
    firstContrast = nonZeroContrasts(randi(length(nonZeroContrasts)));

    % Shuffle the entire vector
    shuffledContrasts = gc_vector(randperm(length(gc_vector)));

    % Assign the first trial to have a non-zero contrast
    shuffledContrasts(1) = firstContrast;

    % This ensures the first trial is not zero without affecting the overall distribution
    if shuffledContrasts(1) == 0
        nonZeroIdx = find(shuffledContrasts > 0, 1);
        if ~isempty(nonZeroIdx)
            shuffledContrasts(1) = shuffledContrasts(nonZeroIdx);
            shuffledContrasts(nonZeroIdx) = firstContrast;
        else
            error('No non-zero contrasts available to assign to the first trial.');
        end
    end

    gc_vector = shuffledContrasts;

    % Assign to blockTrials
    blockTrials = [gc_vector, repmat(current_bc, trialsPerBlock, 1), repmat(currentStream, trialsPerBlock, 1)];

    % Assign to trialBlocks
    trialBlocks{b} = blockTrials;
end

%% -------------------------------Main Experiment Loop------------------------------------
% Start gaze data stream by calling get_gaze_data() the first time
eyetracker_address = 'tet-tcp://169.254.6.40';
eyetracker = Tobii.get_eyetracker(eyetracker_address);
eyetracker.get_gaze_data();
gazeData = eyetracker.get_gaze_data();

% Wait briefly to allow gaze data to accumulate
pause(1); 

% Retrieve gaze data
gazeData = eyetracker.get_gaze_data();

% Check if gaze data is available
if isempty(gazeData)
    error('No gaze data received from the eye tracker.');
end

% Get the first gaze data sample
firstGazeData = gazeData(1);

% Assume SystemTimeStamp field exists and is in microseconds
startReference = double(firstGazeData(1).SystemTimeStamp) / 1e6;

% Eye-tracking with trial start time:
currentGazeData = eyetracker.get_gaze_data();
currentTimestamp = double(currentGazeData(1).SystemTimeStamp) / 1e6;

trialStartTime = currentTimestamp - startReference;

% Display the fields of the first gaze data sample
disp('Fields of firstGazeData:');
disp(fieldnames(firstGazeData));

eyeTrackerStartTime = double(firstGazeData(1).SystemTimeStamp) / 1e6;

matlabStartTime = GetSecs;

% Calculate the time offset
timeOffset = matlabStartTime - eyeTrackerStartTime;  % Positive if MATLAB clock is ahead

try
    % Display initial instructions
    Screen('FillRect', window, gray);
    DrawFormattedText(window, ['In each trial, you will see a circle with striped texture appear either on the left or right side of the fixation dot.\n\n' ...
        'Press the Left or Right Arrow key to indicate the position of the circle.\n\n' ...
        '← = On the left  → = On the right \n\n' ...
        'Press any key to start.'], 'center', 'center', black);
    Screen('Flip', window);
    KbStrokeWait;  
    
    % Start the block 
    for block = 1:blocks_n
        currentBlockTrials = trialBlocks{block};
        currentStream = blockAssignments(block, 1);
        current_bc = bc_options(blockAssignments(block, 2));
        
        stream_prob = streams{currentStream}.probabilities;
        stream_choices = streams{currentStream}.choices;
        
        % Display block initiation message
        Screen('FillRect', window, gray);
        blockMessage = sprintf('Block %d\n\nWait for the experimenter to start.', block);
        DrawFormattedText(window, blockMessage, 'center', 'center', black);
        Screen('Flip', window);
        GetClicks(window);  % Experimenter initiates the start of the block

        for trial = 1:trialsPerBlock
            trialCounter = trialCounter + 1;
            
            % Get trial parameters
            gaborContrast = currentBlockTrials(trial, 1);
            bgContrast = currentBlockTrials(trial, 2);
            streamNumber = currentBlockTrials(trial, 3);
            
            % Get the trial-specific probabilities and choice within the current stream
            probLeft = stream_prob(trial, 1); 
            probRight = stream_prob(trial, 2);
            trialChoice = stream_choices(trial);
            
            if trialChoice == 1
                % Left
                gaborPosX = xCenter + gp_pix(1); 
                gaborPosition = 'Left';
                correctKey = left_key;
            else
                % Right
                gaborPosX = xCenter + gp_pix(2); 
                gaborPosition = 'Right';
                correctKey = right_key;
            end
            
            % Determine the current Gabor contrast list based on background contrast
            if bgContrast == 0
                current_gc_list = gc0;
            elseif bgContrast == 0.25 || bgContrast == 0.5
                current_gc_list = gc;
            else
                error('Invalid background contrast value.');
            end
            
           % Determine triggering code for each condition
            tolerance = 1e-6;
            
            % Determine the index for the gabor position 
            if strcmp(gaborPosition, 'Left')
                pos_index = 0;  
                p_val = probLeft;  
            else
                pos_index = 1;  
                p_val = probRight;  
            end
            
            % Determine the index for the probability condition
            if abs(p_val - 0.5) < tolerance
                p_index = 1;
            elseif abs(p_val - 0.2) < tolerance
                p_index = 2;
            elseif abs(p_val - 0.8) < tolerance
                p_index = 3;
            else
                error('Invalid probability value for trigger code calculation.');
            end
            
            % Determine the index for the background contrast
            if abs(bgContrast - 0) < tolerance
                bg_index = 1;
            elseif abs(bgContrast - 0.25) < tolerance
                bg_index = 2;
            elseif abs(bgContrast - 0.5) < tolerance
                bg_index = 3;
            else
                error('Invalid background contrast for trigger code calculation.');
            end
            
            % Determine the gabor contrast index based on the current contrast list
            gaborContrastIndex = find(abs(current_gc_list - gaborContrast) < tolerance);
            if isempty(gaborContrastIndex)
                error('Gabor contrast not found in the current contrast list.');
            end
            
            % Calculate the unique stimulus trigger code
            % Multipliers: 54 for position, 18 for probability, 6 for background contrast
            base = 10;
            stim_trigger =base + pos_index * 54 + (p_index - 1) * 18 + (bg_index - 1) * 6 + (gaborContrastIndex - 1) + 1;

            % Generate background grating for the entire screen for this trial
            [xBg, yBg] = meshgrid(-screenXpixels / 2 : screenXpixels / 2 - 1, -screenYpixels / 2 : screenYpixels / 2 - 1);
            
            % Compute the background grating
            stripeFreq = spatial_freq / pixelsPerDegree;  
            sinWaveBg = cos(2 * pi * stripeFreq * xBg);
            
            % Multiply background grating by contrast
            gratingBg = bgContrast * sinWaveBg;
            
            % Normalize to [0,1] range
            gratingBgNormalized = 0.5 + 0.5 * gratingBg;
            
            % Create RGB image
            gratingBgRGB = repmat(gratingBgNormalized, [1, 1, 3]);
            
            % Create a texture for the background
            bgTexture = Screen('MakeTexture', window, gratingBgRGB);

            % Draw the background grating texture
            Screen('DrawTexture', window, bgTexture);

            % Draw Fixation dot (Red before stimulus)
            Screen('DrawDots', window, [xCenter yCenter], dotSize, red, [], 2);
            Screen('Flip', window);

            if trial == 1
                WaitSecs(0.85);  
            end

            % Record trial start time
            trialStartTime = GetSecs;

            % Clear previous gaze data
            eyetracker.get_gaze_data(); 

            % Display Gabor if gaborContrast is greater than zero
            if gaborContrast > 0
                % Define the size of the Gabor (pixels)
                gaborDimPix = 800;  
                [xGabor, yGabor] = meshgrid(-gaborDimPix / 2 : gaborDimPix / 2 - 1, -gaborDimPix / 2 : gaborDimPix / 2 - 1);
                sigma = gaborDimPix / 7;  % Standard deviation for the Gaussian
                gaussianMask = exp(-((xGabor).^2 + (yGabor).^2) / (2 * sigma^2));

                % Create the sinusoidal grating
                stripeFreq = spatial_freq / pixelsPerDegree;  
                gaborGrating = cos(2 * pi * stripeFreq * xGabor);

                % Combine the grating with the Gaussian mask
                gaborPatch = gaborGrating .* gaussianMask;
                % Adjust the contrast of the Gabor
                gaborPatch = gaborContrast * gaborPatch;
                % Create the alpha channel
                alphaChannel = gaussianMask;
                % Create the background grating portion
                bgGratingPortion = bgContrast * cos(2 * pi * stripeFreq * xGabor);
                % Combine background and Gabor
                combinedPatch = bgGratingPortion + gaborPatch;
                % Normalize to [0,1]
                combinedPatchNormalized = 0.5 + 0.5 * combinedPatch;
                % Create RGB image and include alpha channel
                gaborImage = repmat(combinedPatchNormalized, [1, 1, 3]);
                gaborImage(:, :, 4) = alphaChannel;

                % Create the texture with alpha channel
                gaborTexture = Screen('MakeTexture', window, gaborImage);
                % Destination rectangle
                dstRect = CenterRectOnPointd([0 0 gaborDimPix gaborDimPix], gaborPosX, yCenter);

                % Record stimulus onset time
                stimOnsetTime = GetSecs();
                Stimulus_Onset_Time = stimOnsetTime;

                % Draw background, Gabor texture with alpha channel, fixation dot
                Screen('DrawTexture', window, bgTexture);
                Screen('DrawTexture', window, gaborTexture, [], dstRect);
                Screen('DrawDots', window, [xCenter yCenter], dotSize, green, [], 2);
                stimOnset = Screen('Flip', window);
            
                % Send the trigger code to Biosemi
                io64(ioObj, address, stim_trigger);
            
                % Wait for stimulus duration
                WaitSecs(0.05);
            
                % Close texture
                Screen('Close', gaborTexture);
            else
                % Only show fixation dot
                Screen('DrawTexture', window, bgTexture);
                Screen('DrawDots', window, [xCenter yCenter], dotSize, green, [], 2);
                stimOnset = Screen('Flip', window);
            
                % Send the trigger code even if Gabor contrast is zero
                io64(ioObj, address, stim_trigger);
            
                WaitSecs(0.05);
            end

            % Keep fixation dot on screen
            Screen('DrawTexture', window, bgTexture);
            Screen('DrawDots', window, [xCenter yCenter], dotSize, green, [], 2);
            Screen('Flip', window);

            %% Response Collection
            responded = false;
            response = '';
            rt = NaN;
            correct = false;
            userTerminated = false;
            while ~responded
                [keyIsDown, keyTime, keyCode] = KbCheck;
                if keyIsDown
                    if keyCode(escape_key)
                        userTerminated = true;
                        break;
                    elseif keyCode(left_key) || keyCode(right_key)
                        % Record response time
                        responseTime = keyTime;
                        rt = keyTime - stimOnset;
                        % Response collection
                        if keyCode(left_key)
                            response = 'Left';
                            io64(ioObj, address, 1);
                        elseif keyCode(right_key)
                            response = 'Right';
                            io64(ioObj, address, 2);
                        end

                        responded = true;
                        correct = strcmpi(response, gaborPosition);
                    end
                end
            end

            if userTerminated
                error('Experiment terminated by user.');
            end

            % Collect gaze data for this trial
            gaze_data = eyetracker.get_gaze_data();  % Retrieve data collected during the trial

            % Store the gaze data for this trial
            allEyeTrackingData{trialCounter} = gaze_data;

            % Record trial end time
            trialEndTime = GetSecs;

            stimOnsetTimes(trialCounter) = stimOnsetTime;
            responseTimes(trialCounter) = responseTime;

            % Inter-trial interval with fixation dot
            Screen('DrawTexture', window, bgTexture);
            Screen('DrawDots', window, [xCenter yCenter], dotSize, red, [], 2);
            Screen('Flip', window);
            WaitSecs(interval);
            Screen('Close', bgTexture);

            % Store trial data in memory
            allTrialData{trialCounter, 1} = trialCounter;
            allTrialData{trialCounter, 2} = block;
            allTrialData{trialCounter, 3} = streamNumber;
            allTrialData{trialCounter, 4} = bgContrast;
            allTrialData{trialCounter, 5} = gaborContrast;
            allTrialData{trialCounter, 6} = gaborPosition;
            allTrialData{trialCounter, 7} = response;
            allTrialData{trialCounter, 8} = correct;
            allTrialData{trialCounter, 9} = rt;
            allTrialData{trialCounter, 10} = probLeft;
            allTrialData{trialCounter, 11} = probRight;
            allTrialData{trialCounter, 12} = trialStartTime;
            allTrialData{trialCounter, 13} = trialEndTime;
        end

        % After block is completed, save data
%         try
%             % Save behavioral data after each block
%             dataTable = cell2table(allTrialData(1:trialCounter, :), 'VariableNames', headers);
%             filenameBlock = fullfile(participantFolder, ['S', num2str(subjid), '_experiment_responses_block_', num2str(block), '.xlsx']);
%             writetable(dataTable, filenameBlock, 'WriteVariableNames', true);
% 
%             % Also save behavioral data as .mat and .tsv
%             behDataStruct = table2struct(dataTable);
%             behMatFilenameBlock = fullfile(participantFolder, ['S', num2str(subjid), '_experiment_responses_block_', num2str(block), '.mat']);
%             behTsvFilenameBlock = fullfile(participantFolder, ['S', num2str(subjid), '_experiment_responses_block_', num2str(block), '.tsv']);
%             save(behMatFilenameBlock, 'behDataStruct');
%             writetable(dataTable, behTsvFilenameBlock, 'FileType', 'text', 'Delimiter', '\t', 'WriteVariableNames', true);
% 
%             % Process and save eye tracking data collected up to this point
%             processedEyeData = processEyeTrackingData(allEyeTrackingData, trialCounter, timeOffset, stimOnsetTimes, responseTimes);
%             eyeDataTable = struct2table(processedEyeData);
% 
%             eyeDataFilenameBlock = fullfile(participantFolder, ['S', num2str(subjid), '_eye_tracking_data_block_', num2str(block), '.mat']);
%             eyeDataExcelFilenameBlock = fullfile(participantFolder, ['S', num2str(subjid), '_eye_tracking_data_block_', num2str(block), '.xlsx']);
%             eyeDataTSVFilenameBlock = fullfile(participantFolder, ['S', num2str(subjid), '_eye_tracking_data_block_', num2str(block), '.tsv']);
% 
%             save(eyeDataFilenameBlock, 'processedEyeData');
%             writetable(eyeDataTable, eyeDataExcelFilenameBlock, 'WriteVariableNames', true);
%             writetable(eyeDataTable, eyeDataTSVFilenameBlock, 'FileType', 'text', 'Delimiter', '\t', 'WriteVariableNames', true);
% 
%         catch dataSaveError
%             warning('Failed to save data after block %d: %s', block, dataSaveError.message);
%         end

            if block == ceil(blocks_n / 2)
            try
                % Save behavioral data at midpoint
                dataTableMid = cell2table(allTrialData(1:trialCounter, :), 'VariableNames', headers);
                filenameMid = fullfile(participantFolder, ['S', num2str(subjid), '_experiment_responses_midpoint.xlsx']);
                writetable(dataTableMid, filenameMid, 'WriteVariableNames', true);
                
                % Also save as .mat and .tsv if desired
                behDataStructMid = table2struct(dataTableMid);
                behMatFilenameMid = fullfile(participantFolder, ['S', num2str(subjid), '_experiment_responses_midpoint.mat']);
                % behTsvFilenameMid = fullfile(participantFolder, ['S', num2str(subjid), '_experiment_responses_midpoint.tsv']);
                save(behMatFilenameMid, 'behDataStructMid');
                writetable(dataTableMid, behTsvFilenameMid, 'FileType', 'text', 'Delimiter', '\t', 'WriteVariableNames', true);
                
                % Save eye tracking data at midpoint
                processedEyeDataMid = processEyeTrackingData(allEyeTrackingData, trialCounter, timeOffset, stimOnsetTimes, responseTimes);
                eyeDataTableMid = struct2table(processedEyeDataMid);
                eyeDataFilenameMid = fullfile(participantFolder, ['S', num2str(subjid), '_eye_tracking_data_midpoint.mat']);
                eyeDataExcelFilenameMid = fullfile(participantFolder, ['S', num2str(subjid), '_eye_tracking_data_midpoint.xlsx']);
                eyeDataTSVFilenameMid = fullfile(participantFolder, ['S', num2str(subjid), '_eye_tracking_data_midpoint.tsv']);
                save(eyeDataFilenameMid, 'processedEyeDataMid');
                writetable(eyeDataTableMid, eyeDataExcelFilenameMid, 'WriteVariableNames', true);
                writetable(eyeDataTableMid, eyeDataTSVFilenameMid, 'FileType', 'text', 'Delimiter', '\t', 'WriteVariableNames', true);
                
                % Provide feedback on save success
                fprintf('Midpoint data saved successfully at block %d.\n', block);
            catch midDataSaveError
                warning('Failed to save midpoint data: %s', Error.message);
            end
            end

        if block < blocks_n
            Screen('FillRect', window, gray);
            blockMessage = sprintf('Block %d completed.\n\nWait for the experimenter to proceed to the next block.', block);
            DrawFormattedText(window, blockMessage, 'center', 'center', black);  
            Screen('Flip', window);
            GetClicks(window);
        end
    end

    Screen('FillRect', window, gray);
    DrawFormattedText(window, 'End of Experiment.\n\nThank you!', 'center', 'center', white);
    Screen('Flip', window);
    KbStrokeWait;

    %% --------------------Stop Eye Tracking and Save Final Data----------------------
    eyetracker.stop_gaze_data();

    % Process and save final eye tracking data
    processedEyeData = processEyeTrackingData(allEyeTrackingData, trialCounter, timeOffset, stimOnsetTimes, responseTimes);
    eyeDataTable = struct2table(processedEyeData);

    eyeDataFilename = fullfile(participantFolder, ['S', num2str(subjid), '_eye_tracking_data.mat']);
    eyeDataExcelFilename = fullfile(participantFolder, ['S', num2str(subjid), '_eye_tracking_data.xlsx']);
    eyeDataTSVFilename = fullfile(participantFolder, ['S', num2str(subjid), '_eye_tracking_data.tsv']);

    save(eyeDataFilename, 'processedEyeData');
    writetable(eyeDataTable, eyeDataExcelFilename, 'WriteVariableNames', true);
    writetable(eyeDataTable, eyeDataTSVFilename, 'FileType', 'text', 'Delimiter', '\t', 'WriteVariableNames', true);

    % Save final behavioral data
    try
        dataTable = cell2table(allTrialData, 'VariableNames', headers);
        writetable(dataTable, filename, 'WriteVariableNames', true);

        % Also save final behavioral data as .mat and .tsv
        behDataStruct = table2struct(dataTable);
        behMatFilename = fullfile(participantFolder, ['S', num2str(subjid), '_experiment_responses.mat']);
        behTsvFilename = fullfile(participantFolder, ['S', num2str(subjid), '_experiment_responses.tsv']);
        save(behMatFilename, 'behDataStruct');
        writetable(dataTable, behTsvFilename, 'FileType', 'text', 'Delimiter', '\t', 'WriteVariableNames', true);

    catch ME
        warning('Could not write full data to Excel file or other formats.');
        disp(ME.message);
    end

    % Cleanup actions
    ShowCursor;
    sca;
    Priority(0);

catch ME
    % Display error message
    fprintf('An error occurred: %s\n', ME.message);

    % Stop gaze data stream
    eyetracker.stop_gaze_data();

    % Process and save partial eye tracking data collected up to this point
    processedEyeData = processEyeTrackingData(allEyeTrackingData, trialCounter, timeOffset, stimOnsetTimes, responseTimes);
    eyeDataTable = struct2table(processedEyeData);

    partialMatFilename = fullfile(participantFolder, ['S', num2str(subjid), '_eye_tracking_data_partial.mat']);
    partialTsvFilename = fullfile(participantFolder, ['S', num2str(subjid), '_eye_tracking_data_partial.tsv']);
    partialXlsxFilename = fullfile(participantFolder, ['S', num2str(subjid), '_eye_tracking_data_partial.xlsx']);

    save(partialMatFilename, 'processedEyeData');
    writetable(eyeDataTable, partialTsvFilename, 'FileType', 'text', 'Delimiter', '\t', 'WriteVariableNames', true);
    writetable(eyeDataTable, partialXlsxFilename, 'WriteVariableNames', true);

    % Save partial behavioral data
    if trialCounter > 0
        partialData = allTrialData(1:trialCounter, :);
        partialDataTable = cell2table(partialData, 'VariableNames', headers);

        partialBehMatFilename = fullfile(participantFolder, ['S', num2str(subjid), '_experiment_responses_partial.mat']);
        partialBehTsvFilename = fullfile(participantFolder, ['S', num2str(subjid), '_experiment_responses_partial.tsv']);
        partialBehXlsxFilename = fullfile(participantFolder, ['S', num2str(subjid), '_experiment_responses_partial.xlsx']);

        behDataStruct = table2struct(partialDataTable);
        save(partialBehMatFilename, 'behDataStruct');
        writetable(partialDataTable, partialBehTsvFilename, 'FileType', 'text', 'Delimiter', '\t', 'WriteVariableNames', true);
        writetable(partialDataTable, partialBehXlsxFilename, 'WriteVariableNames', true);
    end

    % Cleanup actions
    ShowCursor;
    sca;
    Priority(0);

    % Rethrow the error to see the stack trace
    rethrow(ME);
end

%% ---------------------------Function to create streams--------------------------
function streams = createGaborSegStreamB()
    rng(0);  % Set random seed

    % Define Streams and Trials
    streams_n = 3;
    trialsPerStream = 300;
    probOptions = [0.2 0.8; 0.8 0.2];
    totalLeft = trialsPerStream / 2;
    totalRight = trialsPerStream / 2;

    % Probabilities for the first 20 trials within each stream should be 50:50
    first_t = 20;
    firstLeft = first_t / 2;
    firstRight = first_t / 2;
    streams = cell(1, streams_n);

    for s = 1:streams_n
        probabilities = [];
        choices = [];
        cumulativeLeft = 0;
        cumulativeRight = 0;

        trialIdx = 1;
        previousProb = [0.5, 0.5];

        % First 20 trials with 50:50 probability
        firstProb = [0.5, 0.5];
        probabilities = [probabilities; repmat(firstProb, first_t, 1)];
        firstChoices = [ones(firstLeft, 1); zeros(firstRight, 1)];
        firstChoices = firstChoices(randperm(first_t));
        choices = [choices; firstChoices];
        cumulativeLeft = cumulativeLeft + firstLeft;
        cumulativeRight = cumulativeRight + firstRight;

        trialIdx = trialIdx + first_t;

        while trialIdx <= trialsPerStream
            segmentLength = randi([10, 20]);
            if trialIdx + segmentLength - 1 > trialsPerStream
                segmentLength = trialsPerStream - trialIdx + 1;
            end

            % Choose a new probability different from the previous one
            availableProbs = probOptions(~ismember(probOptions, previousProb, 'rows'), :);
            selectedProb = availableProbs(randi(size(availableProbs, 1)), :);

            expectedLeft = selectedProb(1) * segmentLength;
            expectedRight = selectedProb(2) * segmentLength;

            leftCount = round(expectedLeft);
            rightCount = segmentLength - leftCount;

            % Adjust counts if exceeding total allowed trials
            if cumulativeLeft + leftCount > totalLeft
                leftCount = totalLeft - cumulativeLeft;
                rightCount = segmentLength - leftCount;
            end
            if cumulativeRight + rightCount > totalRight
                rightCount = totalRight - cumulativeRight;
                leftCount = segmentLength - rightCount;
            end

            % Ensure counts are non-negative
            leftCount = max(leftCount, 0);
            rightCount = max(rightCount, 0);

            % Generate choices for the segment
            segmentChoices = [ones(leftCount, 1); zeros(rightCount, 1)];
            if length(segmentChoices) < segmentLength
                missingTrials = segmentLength - length(segmentChoices);
                additionalChoices = zeros(missingTrials, 1);
                segmentChoices = [segmentChoices; additionalChoices];
            end
            segmentChoices = segmentChoices(randperm(length(segmentChoices)));

            choices = [choices; segmentChoices];
            probabilities = [probabilities; repmat(selectedProb, segmentLength, 1)];

            trialIdx = trialIdx + segmentLength;

            cumulativeLeft = cumulativeLeft + sum(segmentChoices);
            cumulativeRight = cumulativeRight + (length(segmentChoices) - sum(segmentChoices));
            previousProb = selectedProb;
        end

        numLeftChoices = sum(choices);
        numRightChoices = length(choices) - numLeftChoices;
        if numLeftChoices ~= totalLeft || numRightChoices ~= totalRight
            error('Total left and right trials do not match the required counts.');
        end

        % Assign probabilities and choices to the streams structure
        streams{s}.probabilities = probabilities;
        streams{s}.choices = choices;
    end
end

%% -------------------- Function to Process Eye Tracking Data ----------------------
function processedEyeData = processEyeTrackingData(allEyeTrackingData, trialCounter, timeOffset, stimOnsetTimes, responseTimes)
   processedEyeData = struct();
    dataIdx = 0;

    for trialNum = 1:trialCounter
        gazeDataArray = allEyeTrackingData{trialNum};  % raw gaze data for this trial
        
        % Retrieve the corresponding times for this trial
        thisStimOnset = stimOnsetTimes(trialNum);
        thisResponse  = responseTimes(trialNum);

        for pointIdx = 1:length(gazeDataArray)
            dataIdx = dataIdx + 1;
            gaze_data = gazeDataArray(pointIdx);

            % Extract timestamp and adjust with timeOffset
            if isfield(gaze_data, 'SystemTimeStamp')
                timestamp = gaze_data.SystemTimeStamp / 1e6;  % microseconds to seconds
            elseif isfield(gaze_data, 'DeviceTimeStamp')
                timestamp = gaze_data.DeviceTimeStamp / 1e6;
            elseif isfield(gaze_data, 'Timestamp')
                timestamp = gaze_data.Timestamp / 1e6;
            else
                timestamp = NaN;
            end

            adjustedTimestamp = timestamp + timeOffset;

            % Store in the processed data
            processedEyeData(dataIdx).Trial = trialNum;
            processedEyeData(dataIdx).Timestamp = adjustedTimestamp;
            processedEyeData(dataIdx).Stimulus_Onset_Time = thisStimOnset;
            processedEyeData(dataIdx).Response_Time = thisResponse;
            processedEyeData(dataIdx).Reaction_Time = thisResponse - thisStimOnset;

            % Left eye data
            leftEye = gaze_data.LeftEye;

            % Right eye data
            rightEye = gaze_data.RightEye;

            % Store data in the structure array
            processedEyeData(dataIdx).Trial = trialNum;
            processedEyeData(dataIdx).Timestamp = adjustedTimestamp;
            % For each trial in 1:trialCounter
            processedEyeData(dataIdx).Stimulus_Onset_Time = stimOnsetTimes(trialNum);
            processedEyeData(dataIdx).Response_Time = responseTimes(trialNum);
            processedEyeData(dataIdx).Reaction_Time = responseTimes(trialNum) - stimOnsetTimes(trialNum);

            % Left eye data
            processedEyeData(dataIdx).Left_GazePoint_OnDisplayArea = leftEye.GazePoint.OnDisplayArea;
            processedEyeData(dataIdx).Left_GazePoint_InUserCoordinateSystem = leftEye.GazePoint.InUserCoordinateSystem;
            processedEyeData(dataIdx).Left_GazePoint_Validity = leftEye.GazePoint.Validity.value;
            processedEyeData(dataIdx).Left_Pupil_Diameter = leftEye.Pupil.Diameter;
            processedEyeData(dataIdx).Left_Pupil_Validity = leftEye.Pupil.Validity.value;
            processedEyeData(dataIdx).Left_GazeOrigin_InUserCoordinateSystem = leftEye.GazeOrigin.InUserCoordinateSystem;
            processedEyeData(dataIdx).Left_GazeOrigin_Validity = leftEye.GazeOrigin.Validity.value;

            % Right eye data
            processedEyeData(dataIdx).Right_GazePoint_OnDisplayArea = rightEye.GazePoint.OnDisplayArea;
            processedEyeData(dataIdx).Right_GazePoint_InUserCoordinateSystem = rightEye.GazePoint.InUserCoordinateSystem;
            processedEyeData(dataIdx).Right_GazePoint_Validity = rightEye.GazePoint.Validity.value;
            processedEyeData(dataIdx).Right_Pupil_Diameter = rightEye.Pupil.Diameter;
            processedEyeData(dataIdx).Right_Pupil_Validity = rightEye.Pupil.Validity.value;
            processedEyeData(dataIdx).Right_GazeOrigin_InUserCoordinateSystem = rightEye.GazeOrigin.InUserCoordinateSystem;
            processedEyeData(dataIdx).Right_GazeOrigin_Validity = rightEye.GazeOrigin.Validity.value;
        end
    end
end
