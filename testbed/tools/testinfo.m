function [runid,testid,datadir]=testinfo
% Define some strings for plotting purposes and locate data directory

PN = fliplr(pwd);
[runid,R] = strtok(PN,'\'); runid = fliplr(runid);
[testid,R] = strtok(R,'\'); testid = fliplr(testid);
datadir=['..\..\data\',runid,'\']
