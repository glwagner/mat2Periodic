function [p, sol] = runSpectralModel(p, sol, nSteps)

    % Initialize computation time
    p.tc = tic;

    % Initialize iterate and time if not already initialized
    if ~exist('p.iit'), p.iit = 0; end
    if ~exist('p.t'), p.t = 0; end

    % Initialize saving 
    if p.dnSave ~= 0 && ~exist(p.saveCount)
        p.saveCount = 0;
        p.saveDir = fprintf('%s/data', pwd);

        if saveDir ~= 7, mkdir(p.saveDir), end

        saveModelState(p, sol)
    end

    % Loop over time-steps.
    fprintf('Time-stepping with %s...\n', p.timeStepper), tc=tic; t1=tic;
    for iStep = 1:nSteps
        
        sol = takeTimeStep(p, p.mu, sol);

        % Update iteration
        p.t = p.t + p.dt;
        p.iit = p.iit + 1;

        % Calculate diagnostics and display message
        if mod(p.iit, p.dnPrint) == 0 
            diags = getDiagnostics(p, sol);
            writeMessage(p, diags)
            % Restart display clock.
            p.tc = tic;
        end

        % Quick plot
        if mod(p.iit, p.dnPlot) == 0 
            quickPlot(p, sol);
        end

        % Save
        if mod(p.iit, p.dnSave) == 0 
            p.saveCount = p.saveCount + 1;
            saveModelState(p.saveDir, p.saveCount)
        end

    end
end

function writeMessage(p, diags)

    % Base message
    msg = sprintf('ii = %05d, t = %06.2f s, tComp = %0.2f s', ...
                    p.iit, p.t, toc(p.tc));

    % Put diagnostics into print string 
    if ~isempty(diags)
        msg = strcat(msg, ', ');
        % Get diagnostics
        diagNames = fieldnames(diags);
        for ii = 1:length(diagNames)
            if diags.(diagNames{ii}).print
                if ii < length(diagNames)
                    newMsg = sprintf(' %s = %0.3e %s,',diagNames{ii}, ...
                                        diags.(diagNames{ii}).value, ...
                                        diags.(diagNames{ii}).units );
                else
                    newMsg = sprintf(' %s = %0.3e %s',diagNames{ii}, ...
                                        diags.(diagNames{ii}).value, ...
                                        diags.(diagNames{ii}).units );
                end
            end
            % Add diagnostic to message
            msg = strcat(msg, newMsg);
        end
    end

    fprintf('%s\n', msg)
end

function saveModelState(p, sol)
    saveName = fprintf('%s/p.name_%05d.mat', p.saveDir, p.saveCount);
    save(saveName, 'sol', 'p');
end
