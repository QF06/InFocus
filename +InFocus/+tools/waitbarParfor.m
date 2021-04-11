function [triggerUpdateFcn] = waitbarParfor(totalLoops, varargin)
%waitbarParfor Waitbar implementation for parfor loops.
%   [triggerUpdateFcn] = waitbarParfor(totalLoops) creates a waitbar for use in parfor loops. totalLoops indicates the number of loops until
%       completion. Use the triggerUpdateFcn to increment the waitbar count. Make sure to call this function in the outer parfor loop.
%   waitbarParfor(___, message) displays the specified message on the waitbar.
%   waitbarParfor(___, Name, Value) accepts all Name-Value pairs for the waitbar function. Use the "Name" option to set the window title.
%
%     Example
%         nLoops = 100;
% 
%         updateWaitbar = waitbarParfor(nLoops, "Calculation in progress...");
%         parfor loopCnt = 1:nLoops
%             A = rand(5000);
%             updateWaitbar(); %#ok<PFBNS>
%         end
%
% Author: Girmi Schouten (girmi.schouten@uantwerpen.be), 2019. 
% Written in MATLAB 2019b, tested on Ubuntu 18.04 & Windows 10.


    %% Parse input arguments

    argParser = inputParser();
    argParser.KeepUnmatched = true;

    addRequired(argParser, "totalLoops", @isscalar);
    defaultWaitbarMessage = "Please wait ...";
    addOptional(argParser, "waitbarMessage", defaultWaitbarMessage, @(str) isstring(str) || ischar(str));

    parse(argParser, totalLoops, varargin{:});
    totalLoops = argParser.Results.totalLoops;
    waitbarMessage = argParser.Results.waitbarMessage;

    
    %% Initialize waitbar requirements
    
    parellelDataQueue = parallel.pool.DataQueue;
    afterEach(parellelDataQueue, @updateWaitbar);
    
    waitbarHandle = waitbar(0, waitbarMessage);
    % Pass unmatched parameters to the waitbar
    if ~isempty(fieldnames(argParser.Unmatched))
        waitbarOptions = namedargs2cell(argParser.Unmatched);
        set(waitbarHandle, waitbarOptions{:});
    end
    
    triggerUpdateFcn = @updateProxy;

    loopCnt = 1;
    
    
    %% Helper functions
    
    function updateWaitbar(~)
        if ~isvalid(waitbarHandle)
            warning("waitbarParfor:waitbarExpired", "The waitbar has expired, please create a new one.");
            return;
        elseif loopCnt == totalLoops
            close(waitbarHandle);
            return;
        end
        
        waitbar(loopCnt/totalLoops, waitbarHandle);
        loopCnt = loopCnt + 1;
    end

    function updateProxy()
        send(parellelDataQueue, []);
    end


end

