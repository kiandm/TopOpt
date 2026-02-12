function [BCs_xy,loads_xy] = loadConditions()
% Defines two nodesets from Cubit ABAQUS file, containing nodes for BCs,
% and nodes for loads.

BCs_xy = [
       2,     182,     218,     219,     220,     221,...
     222,     223,     224,     225,     226,     227,...
     228,     229,     230,     231,     232,     233,...
     234,     235,     236,     237,     238,     239,...
     240
     ];

loads_xy = [134,    147,    148];

end