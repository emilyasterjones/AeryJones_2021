% newpos = POSINTERP(posdata, options)
%          posinterp returns interpolated position and direction data
%          where posdata is a position structure from readpos() and
%               options take the form 'option', value, 'option', value, ... 
%                   'diode' - 1 if only front diode is used for pos, and 0
%                               if only back is used.  'diode' can be any
%                               value in between 0 and 1, in which case the
%                               final pos will be calculated as the specified
%                               ratio between the front and back diodes
%                               (default 0.5, i.e. average of the diodes)
%                   'reversex' 1 for reversing flipping the pos image in x
%                   'reversey' 1 for reversing flipping the pos image in y

%                   'bounds' - a two by two array giving the bounds of the
%                              environment (all points outside the bounds are
%                              discarded) (default - all points included)
%                   'maxv'   - the maximum velocity allowed (default infinity)
%                   'maxdevpoints' - the maximum number of points
%                              between two high velocity pairs that are to be 
%                              zeroed out (default 10)
%                   'pinterp' 1 for position interpoloation  (default 1)
%                   'maxpinterp' is the maximum number of missing points to be
%                                interpolated over (default inf)
%                   'maxdinterp' - the maximum distance (in pixels) to be
%                                 interpolated over (default infinity)
%
%newpos is a structure with arglist, description, field, and data elements

function [newpos] = ag_posinterp(posdata, varargin)

arglist = {'diode', 'reversex', 'bounds', 'maxv', 'maxdevpoints', 'pinterp',...
    'maxpinterp', 'maxdinterp', 'reversey'}';

% assign default values to varargin options
diode = 1;            % use average between front and back diodes
reversex = 0; % don't flip pos along x
reversey = 0;
bounds = [-1e100 -1e100 ; 1e100 1e100];     % no bounds
maxv = 1e100;           % allow any velocity
maxdevpoints = 10;
pinterp = 1;            % do not interpolate
maxpinterp = 1e100;     % allow interpolation of any number of points
maxdinterp = 1e100;     % allow interpolation over any distance

if ~isempty(posdata.data)
    % process varargin and overwrite default values
    if (~isempty(varargin))
        assign(varargin{:});
    end

    % put all of the arguments in the newpos.arg field
    newpos.arg = cell(size(arglist,1), 2);
    for i = 1:size(arglist,1)
        newpos.arg{i,1} = arglist{i};
        newpos.arg{i,2} = eval(arglist{i});
    end

    % write out the command parameters to a string
    descript1 = sprintf('interpolated position data from %s structure', inputname(1));
    newpos.descript = descript1;
    newpos.fields = 'time x y';

    % get the times of all pos's
    time = posdata.data(:,1);

    % now interpolate the diode position(s).  when i = 1, the for loop
    % interpolates the front diode position.  when i = 2, the for loop
    % interpolates the back diode position.
    pos = [];
    for i = 1 %2]
        % save current diode positions into tmppos
        tmppos = posdata.data(:,(2*i):(2*i)+1);

        % find all valid (non-zero in both x and y) positions
        indvalid = find(prod(tmppos,2));

        % identify all points representing "unreasonable" velocities and set
        % them to zero
        vel = dist(tmppos(indvalid(1:end-1),:), tmppos(indvalid(2:end),:)) ./ ...
            (time(indvalid(2:end)) - time(indvalid(1:end-1)));
        indtoofast = find(vel > maxv);

        % zero out both the x and y positions from the toofast points
        indtoofastpos = indvalid(indtoofast+1);
        tmppos(indtoofastpos,:) = 0;

        % find all pairs of "toofast" points less than maxdevpoints apart
        % any two points close together and both associated with a high velocity
        % come from tracker artifacts where the wrong point is tracked for a few samples
        % first figure out how many points separate each pair of toofast points
        % then identify the ones that are closer than maxdevpoints
        indpairs = find(diff(indtoofastpos) < maxdevpoints);
        % zero out the pos values between these pairs
        for k = 1:length(indpairs)
            tmppos(indtoofastpos(indpairs(k)):indtoofastpos(indpairs(k)+1),:) = 0;
        end

        if (pinterp)
            % after the changes made above, find the new set of valid and invalid
            % pos points
            tmp = prod(tmppos,2);
            indvalid = find(tmp);
            indzeros = find(~tmp);

            % only look at those positions between first valid and last valid
            indzeros = indzeros(find((indzeros > indvalid(1)) & (indzeros < indvalid(end))));
            % go through the list and calculate the length of each interpolated segment
            % and the distances to be interpolated over
            j = 1;
            while (j < length(indzeros))
                % the first valid position is pos(zeropos(j) - 1)
                startpos = tmppos(indzeros(j)-1, :);
                startindex = j;
                % go though and find the whole list of consecutive zeros
                while ((j<length(indzeros)) & ((indzeros(j+1) - indzeros(j)) == 1))
                    j = j + 1;
                end
                % the end point is the next point after the current zero
                endpos = tmppos(indzeros(j)+1,:);
                % if either the distance to be interpolated over or the number of
                % points is too large, don't interpolated over them
                d = dist(startpos, endpos);
                if ((d > maxdinterp) | ((j-startindex+1) > maxpinterp))
                    indzeros(startindex:j) = -1;
                else
                    % put d into the list of interpdists and the difference between the
                    % indeces into interppoints
                end
                j = j + 1;
            end

            % cut out the negative ones from the list of zeropos
            indzeros = indzeros(find(indzeros ~= -1));

            % get the list of times where interpolation should be done
            interptimes = time(indzeros);

            % interpolate over the x and y coordinates
            tmppos(indzeros,1) = round(interp1(time(indvalid), ...
                tmppos(indvalid, 1), interptimes));
            tmppos(indzeros,2) = round(interp1(time(indvalid), ...
                tmppos(indvalid, 2), interptimes));
        end

        % save this interpolated pos
        pos(:,(2*i)-1:(2*i)) = tmppos;
    end

    indvalidpos = find(prod(pos(:,1:2),2));

    % get the indeces of the first and last valid pos
    pvalidfirst = indvalidpos(1);
    pvalidlast = indvalidpos(end);

    % truncate pos if there are zeros at the beginning and end
    pos = pos(pvalidfirst:pvalidlast, :);
    time = time(pvalidfirst:pvalidlast,1);

    % take all of the valid position samples and the corresponding set of direction
    % samples and save for the output
    newpos.data = [time pos];   % was [time pos dir]
else
    newpos.data = [];
end
