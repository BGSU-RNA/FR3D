% Find the path to a directory if it exists, starting with pwd

function pathToDirectory = zPathToDirectory(directoryName)

    pathToDirectory = '';

    % Check if the directory exists in the current working directory
    if exist([pwd filesep directoryName], 'dir') == 7
        pathToDirectory = pwd;
    elseif exist(directoryName, 'dir')
        % exists somewhere on the path, but Matlab which generally won't tell you where
        % Get the list of all directories on the MATLAB path
        pathList = strsplit(path, pathsep);

        % check for exact match to the path list
        for i = 1:length(pathList)
            % if directory name ends with directoryName, return that directory minus directoryName
            if ~isempty(strfind(pathList{i}, [filesep directoryName]))
                pathToDirectory = pathList{i}(1:end-length(directoryName)-1);
                break
            end
        end

        if length(pathToDirectory) == 0
            % Check for being a subdirectory of a directory on the path
            for i = 1:length(pathList)
                if exist([pathList{i} filesep directoryName], 'dir')
                    pathToDirectory = pathList{i};
                    break
                end
            end
        end
    end
end