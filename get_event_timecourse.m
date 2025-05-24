function [counts, vals, cnts] = get_event_timecourse(metrics, tv, binw)
tmin = floor(min(tv) / binw) * binw;
tmax = ceil(max(tv) / binw) * binw;
edges = (tmin:binw:tmax)';
cnts = (edges(1:end-1) + edges(2:end))/2;
[vals, ~, counts] = groupsummary(metrics, tv, edges, ...
    'mean', 'IncludeEmptyGroups', true);
end