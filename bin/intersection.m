% Returns the row index of b that matches any of the rows in a

function idx = intersection(a, b)

[R,C] = size(a);
index = zeros(length(b),R);
for r = 1:R;
        index(:,r) = a(r,1) == b(:,1) & a(r,2) == b(:,2);
end

index
idx = []










