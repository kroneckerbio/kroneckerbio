function obs = fixObservation(con, obs)
if isnumeric(obs)
    obs = observationAll(obs);
end

if numel(obs) == 1
    obs = repmat(obs, 1,numel(con));
end