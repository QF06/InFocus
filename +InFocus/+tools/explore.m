% explore.m
%
% makes system call to windows explorer
%
% takes a directory,full file path,function name, or model name as input
%
% USAGE:
%   explore
%      opens explorer window to current file open in matlab editor. (if
%      no file is open, then the current explorer window at current directory
%      will be explored
%
%   explore(file)
%      opens explorer window to directory containing "file".  This input
%      can be a full-file-path, m-function name, or model name. 
%
%   explore(dir)
%     opens explorer window to specified direcotry "dir". This input can
%     be a full-path, relative-path, or just the name of any directory
%     on the matlab path.  Will open the first folder on the matlab-path
%     matching the specified name.
%
%
% NOTE: You may change the executable called by explore.m.  The default is
%       "explorer".  To change, use
%       explore SetExplorerCommand <command>
%
%       For example, to use UltraExplorer:
%       explore SetExplorerComand 'UltraExplorer /p'
%
%       Note that the path containing UltraExplorer.exe must be included
%       in your Windows path environment variable.  
%
%       This setting is stored as a persistent variable, and will not be
%       remembered between Matlab sessions, so this line should reside in
%       your startup.m file to permanently change your explorer preference.
%       
%
% ---CVS Keywords----
% $Author: jhopkin $
% $Date: 2009/11/05 17:03:18 $
% $Name: work $
% $Revision: 1.5 $

% $Log: explore.m,v $
% Revision 1.5  2009/11/05 17:03:18  jhopkin
% now can configure which explorer command to call
%
% Revision 1.4  2009/03/16 15:23:21  jhopkin
% removed dependency on isdir.mex32
%
% Revision 1.3  2009/03/16 15:21:22  jhopkin
% Now works with dirname of any directory on the matlab path.
%
% Revision 1.2  2008/11/26 16:10:58  jhopkin
% changed so that now if no arguments are given, this will open an explorer window the file currently opened in the matlab editor. If no file is open, then the Matlab current directory will be opened.
%


function explore(varargin)
	persistent explorerCommand;
	if isempty(explorerCommand)
		%set default explorer command
		explorerCommand = 'explorer';
	end
	
	bContinue = true;
	if(nargin == 0)
		string = char(com.mathworks.mlservices.MLEditorServices.builtinGetActiveDocument);
		if isempty(string)
			string = pwd;
		end
	else
		%check if user is changing the explorer command
		if isequal(lower(varargin{1}),'setexplorercommand')
			explorerCommand = varargin{2};
			bContinue = false;
		else
			%patch all the inputs together.  If there are multiple inputs, put the 
			%string back together assuming it was separated by a space.
			string = varargin{1};
			for i = 2:nargin
				string = strcat(string,' ',varargin{i});
			end
		end
	end
	if bContinue
		%see if string expands to an m-file
		if(which(string))
			string = which(string);
		end
		

		switch exist(string)
			case {2,3,4,6} %a filename
				slashes = find(string == '\' | string == '/');
				lastSlash = slashes(length(slashes));
				string = string(1:lastSlash);
				
				
			case 7 % a directory
				%see if string is a dir
				x = what(string);
				if ~isempty(x) &&  ~isempty(x(1).path)
					string = x(1).path;
				else
					disp('ERROR: not a valid filename or directory');
					bContinue = false;
				end
			otherwise
				disp('ERROR: not a valid filename or directory');
				bContinue = false;
		end
	end
	
	if bContinue
		exStr = sprintf('start %s "%s"',explorerCommand,string);
		disp(exStr);
		dos(exStr);
	end
end

