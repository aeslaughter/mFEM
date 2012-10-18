% Returns the row index of b that matches any of the rows in a

function idx = intersection(a,b)

[R,C] = size(a);
idx = [];
index = false(size(b));
for r = 1:R;
    for c = 1:C;
        index(:,c) = a(r,c) == b(:,c); 
    end

    idx = [idx; find(sum(index,2) == C)];
    
end









