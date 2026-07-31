function viewer = softarm_pose_playback(varargin)
%SOFTARM_POSE_PLAYBACK Open the offline soft-arm pose player.
%   SOFTARM_POSE_PLAYBACK uses softarm_q_log and softarm_bundle from the
%   base workspace, then reconstructs node poses offline.
%   SOFTARM_POSE_PLAYBACK(LOGDATA) also accepts an existing pose or q log.
%   Name-value options include Bundle, AxisLength, and PlaybackSpeed.

viewer = softarm.playbackBackbonePoses(varargin{:});
end
