classdef ParforProgress < handle
    %PARFORPROGRESS Client-side progress for parfor via DataQueue.
    %   pp = ParforProgress(total, title, updateEvery, useWaitbar)
    %   - total:        total number of iterations
    %   - title:        waitbar title (string) [optional]
    %   - updateEvery:  only refresh UI every N ticks (default 1)
    %   - useWaitbar:   true -> waitbar, false -> console prints (default true)
    %
    %   Methods:
    %     pp.tick()   % advance by 1
    %     pp.close()  % close waitbar (also auto-closes on deletion)

    properties (SetAccess = private)
        Total       (1,1) double {mustBePositive}    = 1
        Count       (1,1) double {mustBeNonnegative}  = 0
        Title       (1,1) string                      = "Progress"
        UpdateEvery (1,1) double {mustBePositive}     = 1
        UseWaitbar  (1,1) logical                      = true
        StartTime   datetime                           = datetime.empty
    end

    properties (Access = private)
        H % waitbar handle (if any)
        LastPrintedPercent double = -inf
    end

    methods
        function obj = ParforProgress(total, title, updateEvery, useWaitbar)
            if nargin < 2 || isempty(title),       title = "Progress"; end
            if nargin < 3 || isempty(updateEvery), updateEvery = 1;   end
            if nargin < 4 || isempty(useWaitbar),  useWaitbar = true; end

            obj.Total       = total;
            obj.Title       = string(title);
            obj.UpdateEvery = updateEvery;
            obj.UseWaitbar  = useWaitbar;
            obj.StartTime   = datetime('now');

            if obj.UseWaitbar
                obj.H = waitbar(0, char(obj.Title));
            end
        end

        function tick(obj)
            % Advance by 1 and update UI/console occasionally.
            obj.Count = obj.Count + 1;
            if mod(obj.Count, obj.UpdateEvery) ~= 0 && obj.Count < obj.Total
                return;
            end
            p = obj.Count / obj.Total;

            if obj.UseWaitbar
                if isgraphics(obj.H)
                    waitbar(p, obj.H, sprintf('%s — %d/%d (%.1f%%)', ...
                        obj.Title, obj.Count, obj.Total, 100*p));
                end
            else
                percent = floor(100*p);
                if percent >= obj.LastPrintedPercent + 1 || obj.Count == obj.Total
                    obj.LastPrintedPercent = percent;
                    fprintf('[%s] %s: %d/%d (%.0f%%)\n', ...
                        datestr(now,'HH:MM:SS'), obj.Title, obj.Count, obj.Total, 100*p);
                end
            end
        end

        function close(obj)
            if ~isempty(obj.H) && isgraphics(obj.H)
                close(obj.H);
            end
        end

        function delete(obj)
            % Auto-cleanup on clear or scope exit
            obj.close();
        end
    end
end
