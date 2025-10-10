function error = permloss(id, label)
    idset = unique(id);
    n = length(label);
    labelset = unique(label);
    k = length(labelset);
    P = perms(labelset);
    N = size(P, 1);
    confusion = zeros(1, N);
    for i = 1 : N
        for node = 1 : n
            if (label(node) ~= P(i, id(node)))
                   confusion(i) = confusion(i) + 1;
            end
        end
    end
    error = min(confusion);





end