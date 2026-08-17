function subjectMap = build_subject_ic_map(SUBJECTS_ICS)
% BUILD_SUBJECT_IC_MAP  Invert region-organized IC table to subject-organized map.
%
% Subject IDs in SUBJECTS_ICS are offset by +4 (applied internally).

    SUBJECT_OFFSET = 4;   % subject IDs in SUBJECTS_ICS need +4

    nRegions = size(SUBJECTS_ICS, 1);

    allSubj   = [];
    allICs    = [];
    allRegion = {};

    for r = 1:nRegions
        regName = SUBJECTS_ICS{r,1};
        s       = SUBJECTS_ICS{r,2};

        if isempty(s) || ~isfield(s,'Subjects') || ~isfield(s,'ICs')
            continue;
        end

        subj = s.Subjects(:)' + SUBJECT_OFFSET;
        ics  = s.ICs(:)';

        if numel(subj) ~= numel(ics)
            error('Region "%s": Subjects and ICs length mismatch (%d vs %d).', ...
                  regName, numel(subj), numel(ics));
        end

        allSubj   = [allSubj,   subj];
        allICs    = [allICs,    ics];
        allRegion = [allRegion, repmat({regName}, 1, numel(subj))];
    end

    uniqueSubj = unique(allSubj);
    subjectMap = struct('subjectID',{},'icIdx',{},'regionNames',{});

    for k = 1:numel(uniqueSubj)
        sid  = uniqueSubj(k);
        mask = (allSubj == sid);

        subjectMap(k).subjectID   = sid;
        subjectMap(k).icIdx       = allICs(mask);
        subjectMap(k).regionNames = allRegion(mask);
    end



    % --- Detect ICs assigned to multiple regions for the same subject ---
    for k = 1:numel(subjectMap)
        sid    = subjectMap(k).subjectID;
        ics    = subjectMap(k).icIdx;
        regs   = subjectMap(k).regionNames;

        [uniqueICs, ~, grpIdx] = unique(ics, 'stable');

        if numel(uniqueICs) == numel(ics)
            continue;   % no duplicates for this subject
        end

        % There are duplicates — resolve them one at a time
        keepMask = true(1, numel(ics));

        for u = 1:numel(uniqueICs)
            positions = find(grpIdx == u);
            if numel(positions) < 2, continue; end

            ic = uniqueICs(u);
            conflictRegions = regs(positions);

            fprintf('\n*** Subject %d: IC%d is assigned to %d regions:\n', ...
                    sid, ic, numel(conflictRegions));
            for j = 1:numel(conflictRegions)
                fprintf('     [%d] %s\n', j, conflictRegions{j});
            end
            fprintf('     [0] Drop this IC entirely\n');

            choice = -1;
            while choice < 0 || choice > numel(conflictRegions)
                resp = input(sprintf('Choose region for IC%d (0-%d): ', ...
                                     ic, numel(conflictRegions)), 's');
                choice = str2double(resp);
                if isnan(choice), choice = -1; end
            end

            if choice == 0
                keepMask(positions) = false;     % drop all copies
                fprintf('  -> IC%d dropped from subject %d.\n', ic, sid);
            else
                drop = positions;
                drop(choice) = [];                % keep the chosen, drop the rest
                keepMask(drop) = false;
                fprintf('  -> IC%d kept under "%s".\n', ic, conflictRegions{choice});
            end
        end

        subjectMap(k).icIdx       = ics(keepMask);
        subjectMap(k).regionNames = regs(keepMask);
    end


    fprintf('Built IC map for %d subjects:\n', numel(subjectMap));
    for k = 1:numel(subjectMap)
        fprintf('  Subject %2d: %d ICs  [%s]\n', ...
                subjectMap(k).subjectID, ...
                numel(subjectMap(k).icIdx), ...
                strjoin(arrayfun(@num2str, subjectMap(k).icIdx, 'uni', 0), ' '));
    end
end