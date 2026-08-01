function viewer = softarm_pose_playback(varargin)
%SOFTARM_POSE_PLAYBACK Open the offline soft-arm pose player.
%   SOFTARM_POSE_PLAYBACK uses softarm_q_log and softarm_bundle from the
%   base workspace, then reconstructs the arm root and node poses offline.
%   SOFTARM_POSE_PLAYBACK(LOGDATA) also accepts an existing q log.
%   Name-value options include Bundle, AxisLength, PlaybackSpeed, and
%   SamplesPerSegment (default 16, valid range 2 through 128).

viewer = softarm.playbackBackbonePoses(varargin{:});
end
