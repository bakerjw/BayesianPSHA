function [numEventsVect, totEvents] = getNumEventsVect(dataVect, totDays)
%getDatesVect returns a number of events vector containing nu,ber of events
%that have occured till every day in the time range
%specified in dataVect and starting on startDate to be used in changePoint
%calculations

totEvents = length(dataVect) + 1;   % 1 is added to account for first event on day 0

numEventsVect = ones(totDays - 1, 1);
numEvents = 1;  %Number of events till time tau. First event occured on day 0


numEventsVect(1:dataVect(1)-1) = 1; % first event is actally at day = 0, so end index has a '-1'
startDay = dataVect(1); % first day for the next set of assignments

for j = 2:length(dataVect)
    numEventsVect(startDay:startDay+dataVect(j)) = j; 
    startDay = startDay+dataVect(j);
end

numEventsVect(end) = numEventsVect(end) + 1; % one event happened on the last day

