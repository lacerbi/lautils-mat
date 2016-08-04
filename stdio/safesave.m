function safesave(filename)
%SAFESAVE Safeguarded save workspace variables to file.
%   SAFESAVE(FILENAME) stores all variables from the current workspace in a
%   MATLAB formatted binary file (MAT-file) called FILENAME. SAFESAVE
%   creates a temporary copy of an existing file so that if MATLAB fails 
%   while saving it avoids to ruin the old file.

if ~exist(filename,'file')
    save(filename);
else

    % File already exists, move it to temporary copy
    [pathstr,name,ext] = fileparts(filename);
    tempfile = fullfile(pathstr,[name '.old_']);
    movefile(filename,tempfile);
    for iAttempt = 1:3
        temp = [];
        try
            save(filename);
            temp = load(filename); % See if it worked
            if ~isempty(temp); break; end
        catch
            % Did not work
        end
    end
    if isempty(temp)
        movefile(tempfile,filename);
        error(['Cannot save file ''' filename '''.']);
    end

    % Remove temporary old file
    if exist(tempfile,'file'); delete(tempfile); end
end

end

