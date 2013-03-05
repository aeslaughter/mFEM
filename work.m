function work

elem_map = [10,11,11,12; 12,15,17,16]
no_id = [12,13,14,15,16,17,10,11];

v = reshape(elem_map,1,numel(elem_map))
[C,ia,ib] = intersect(no_id,v,'stable')

[Lia,Locb] = ismember(v,no_id);
elem_map = reshape(Locb,size(elem_map))

% elem_map(no_id) = 1:length(no_id)


