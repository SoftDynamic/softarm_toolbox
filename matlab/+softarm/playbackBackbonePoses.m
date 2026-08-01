function viewer = playbackBackbonePoses(logData,options)
%PLAYBACKBACKBONEPOSES Interactively replay logged soft-arm node poses.
arguments
    logData = []
    options.AxisLength (1,1) double {mustBePositive} = 0.08
    options.PlaybackSpeed (1,1) double {mustBePositive} = 0.5
    options.Bundle (1,1) string = ""
end

if isempty(logData)
    if evalin("base","exist('softarm_q_log','var')") == 1
        logData = evalin("base","softarm_q_log");
    else
        error("softarm:MissingQLog", ...
            "softarm_q_log is not present in the base workspace.");
    end
end
[time,qFrames] = normalizeLog(logData);
plant = loadPlaybackPlant(options.Bundle);
poseFrames = softarm.reconstructBackbonePoses(qFrames,plant);
frameCount = numel(time);
axisLength = options.AxisLength;

viewer = uifigure("Name","SoftArm Pose Playback", ...
    "Position",[100 100 940 680],"Color",[0.97 0.97 0.97]);
layout = uigridlayout(viewer,[2 4]);
layout.RowHeight = {"1x",34};
layout.ColumnWidth = {78,"1x",135,90};
layout.Padding = [8 8 8 8];

axesHandle = uiaxes(layout,"Box","on","Color",[1 1 1], ...
    "DataAspectRatio",[1 1 1],"Projection","perspective", ...
    "ZDir","reverse","NextPlot","add");
axesHandle.Layout.Row = 1;
axesHandle.Layout.Column = [1 4];
grid(axesHandle,"on");
view(axesHandle,135,22);
xlabel(axesHandle,"NED x / North (m)");
ylabel(axesHandle,"NED y / East (m)");
zlabel(axesHandle,"NED z / Down (m)");
titleHandle = title(axesHandle,"Soft-arm key-node poses", ...
    "Interpreter","none");
enableDefaultInteractivity(axesHandle);

backbone = line(axesHandle,nan,nan,nan, ...
    "Color",[0.12 0.12 0.12],"LineWidth",2.0, ...
    "Marker","o","MarkerSize",5,"MarkerFaceColor",[0.85 0.85 0.85], ...
    "Tag","softarm_backbone");
line(axesHandle,0,0,0,"LineStyle","none","Marker","+", ...
    "MarkerSize",8,"Color",[0.35 0.35 0.35],"Tag","softarm_world_origin");
xAxis = quiver3(axesHandle,nan,nan,nan,nan,nan,nan,0, ...
    "Color",[0.88 0.10 0.10],"LineWidth",1.4,"MaxHeadSize",0.45, ...
    "Tag","softarm_x_axes");
yAxis = quiver3(axesHandle,nan,nan,nan,nan,nan,nan,0, ...
    "Color",[0.10 0.62 0.18],"LineWidth",1.4,"MaxHeadSize",0.45, ...
    "Tag","softarm_y_axes");
zAxis = quiver3(axesHandle,nan,nan,nan,nan,nan,nan,0, ...
    "Color",[0.12 0.30 0.90],"LineWidth",1.4,"MaxHeadSize",0.45, ...
    "Tag","softarm_z_axes");

playButton = uibutton(layout,"push","Text","Play");
playButton.Layout.Row = 2;
playButton.Layout.Column = 1;
sliderLimits = [time(1),time(end)];
if sliderLimits(1) == sliderLimits(2)
    sliderLimits(2) = sliderLimits(1)+1;
end
timeSlider = uislider(layout,"Limits",sliderLimits, ...
    "Value",time(1));
timeSlider.Layout.Row = 2;
timeSlider.Layout.Column = 2;
timeSlider.MajorTicks = sliderTicks(time);
timeLabel = uilabel(layout,"HorizontalAlignment","center");
timeLabel.Layout.Row = 2;
timeLabel.Layout.Column = 3;
speeds = unique([0.1 0.25 0.5 1 2 options.PlaybackSpeed]);
speedControl = uidropdown(layout, ...
    "Items",compose("%gx",speeds), ...
    "ItemsData",speeds, ...
    "Value",options.PlaybackSpeed);
speedControl.Layout.Row = 2;
speedControl.Layout.Column = 4;

setFixedLimits(axesHandle,poseFrames,axisLength);
playTimer = timer("ExecutionMode","fixedSpacing","Period",0.04, ...
    "BusyMode","drop","ObjectVisibility","off","TimerFcn",@advancePlayback);
wallStart = [];
simulationStart = time(1);

timeSlider.ValueChangingFcn = @dragSlider;
timeSlider.ValueChangedFcn = @commitSlider;
playButton.ButtonPushedFcn = @togglePlayback;
viewer.CloseRequestFcn = @closePlayer;
renderFrame(1);

    function dragSlider(~,event)
        pausePlayback();
        renderAtTime(event.Value);
    end

    function commitSlider(source,~)
        renderAtTime(source.Value);
    end

    function togglePlayback(~,~)
        if strcmp(playTimer.Running,"on")
            pausePlayback();
            return
        end
        if timeSlider.Value >= time(end)
            timeSlider.Value = time(1);
        end
        simulationStart = timeSlider.Value;
        wallStart = tic;
        playButton.Text = "Pause";
        start(playTimer);
    end

    function pausePlayback()
        if isvalid(playTimer) && strcmp(playTimer.Running,"on")
            stop(playTimer);
        end
        if isvalid(playButton)
            playButton.Text = "Play";
        end
    end

    function advancePlayback(~,~)
        requestedTime = simulationStart + toc(wallStart)*speedControl.Value;
        if requestedTime >= time(end)
            requestedTime = time(end);
            pausePlayback();
        end
        timeSlider.Value = requestedTime;
        renderAtTime(requestedTime);
    end

    function renderAtTime(requestedTime)
        [~,frame] = min(abs(time-requestedTime));
        renderFrame(frame);
    end

    function renderFrame(frame)
        transformStack = poseFrames(frame,:).';
        nodeCount = numel(transformStack)/16;
        transforms = reshape(transformStack,4,4,nodeCount);
        positions = reshape(transforms(1:3,4,:),3,nodeCount).';
        xDirections = reshape(transforms(1:3,1,:),3,nodeCount).';
        yDirections = reshape(transforms(1:3,2,:),3,nodeCount).';
        zDirections = reshape(transforms(1:3,3,:),3,nodeCount).';
        set(backbone,"XData",positions(:,1),"YData",positions(:,2), ...
            "ZData",positions(:,3));
        setQuiver(xAxis,positions,axisLength*xDirections);
        setQuiver(yAxis,positions,axisLength*yDirections);
        setQuiver(zAxis,positions,axisLength*zDirections);
        timeSlider.Value = time(frame);
        timeLabel.Text = sprintf("t = %.3f / %.3f s",time(frame),time(end));
        titleHandle.String = sprintf("Soft-arm key-node poses   frame %d / %d", ...
            frame,frameCount);
        drawnow limitrate
    end

    function closePlayer(~,~)
        if isvalid(playTimer)
            stop(playTimer);
            delete(playTimer);
        end
        delete(viewer);
    end
end

function [time,poseFrames] = normalizeLog(logData)
if isa(logData,"Simulink.SimulationOutput")
    try
        logData = logData.get("softarm_q_log");
    catch
        error("softarm:MissingQLog", ...
            "SimulationOutput does not contain softarm_q_log.");
    end
end

if isa(logData,"timeseries")
    time = double(logData.Time(:));
    poseFrames = timeFirst(logData.Data,numel(time));
elseif istimetable(logData)
    rowTimes = logData.Properties.RowTimes;
    time = seconds(rowTimes-rowTimes(1));
    poseFrames = double(table2array(logData));
elseif isstruct(logData) && isfield(logData,"time") && isfield(logData,"signals")
    time = double(logData.time(:));
    poseFrames = timeFirst(logData.signals.values,numel(time));
elseif isnumeric(logData) && ismatrix(logData) && size(logData,2) > 1
    time = double(logData(:,1));
    poseFrames = double(logData(:,2:end));
else
    error("softarm:InvalidQLog","Unsupported q log format.");
end

if isempty(time) || any(~isfinite(time)) || any(diff(time) < 0)
    error("softarm:InvalidQLog","Q log time must be finite and nondecreasing.");
end
if size(poseFrames,1) ~= numel(time) || isempty(poseFrames)
    error("softarm:InvalidQLog", ...
        "Logged data dimensions do not match the time vector.");
end
[time,indices] = unique(time,"stable");
poseFrames = poseFrames(indices,:);
end

function plant = loadPlaybackPlant(bundle)
if bundle == ""
    if evalin("base","exist('softarm_bundle','var')") ~= 1
        error("softarm:MissingBundle", ...
            "A q log requires softarm_bundle in the base workspace or the Bundle option.");
    end
    bundle = string(evalin("base","softarm_bundle"));
end
plant = softarm.loadModel(bundle);
end

function frames = timeFirst(data,timeCount)
data = double(data);
if size(data,1) == timeCount
    frames = reshape(data,timeCount,[]);
elseif size(data,ndims(data)) == timeCount
    order = [ndims(data),1:ndims(data)-1];
    frames = reshape(permute(data,order),timeCount,[]);
elseif rem(numel(data),timeCount) == 0
    frames = reshape(data,[],timeCount).';
else
    error("softarm:InvalidQLog","Q data dimensions do not match its time vector.");
end
end

function ticks = sliderTicks(time)
if time(1) == time(end)
    ticks = time(1);
else
    ticks = linspace(time(1),time(end),5);
end
end

function setFixedLimits(axesHandle,poseFrames,axisLength)
frameCount = size(poseFrames,1);
nodeCount = size(poseFrames,2)/16;
positions = zeros(frameCount*nodeCount,3);
for frame = 1:frameCount
    transforms = reshape(poseFrames(frame,:),4,4,nodeCount);
    rows = (frame-1)*nodeCount+(1:nodeCount);
    positions(rows,:) = reshape(transforms(1:3,4,:),3,nodeCount).';
end
minimum = min(positions,[],1);
maximum = max(positions,[],1);
span = max(maximum-minimum);
span = max(span,4*axisLength);
margin = max(0.12*span,axisLength);
center = 0.5*(minimum+maximum);
halfSpan = 0.5*span+margin;
axesHandle.XLim = center(1)+[-halfSpan halfSpan];
axesHandle.YLim = center(2)+[-halfSpan halfSpan];
axesHandle.ZLim = center(3)+[-halfSpan halfSpan];
end

function setQuiver(quiverHandle,bases,directions)
set(quiverHandle,"XData",bases(:,1),"YData",bases(:,2),"ZData",bases(:,3), ...
    "UData",directions(:,1),"VData",directions(:,2),"WData",directions(:,3));
end
