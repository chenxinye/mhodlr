function sortIdx = sort_by_u(u_chain)
    callCellFunc = cellfun(@(x)x.u, u_chain);
    [~, sortIdx] = sort(callCellFunc);
end