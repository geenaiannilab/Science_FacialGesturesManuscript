function struct2vars(sSTMA)

names = fieldnames(sSTMA);
for i = 1:numel(names)
    assignin('caller', names{i}, sSTMA.(names{i}));
end
end