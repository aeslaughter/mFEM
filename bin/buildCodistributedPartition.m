function out = buildCodistributedPartition(value)
    %BUILDCODISTRIBUTEDPARTITION
    %
    % See Also init getDof

       labs = 1:numlabs;
       other = labs(labs ~= labindex);
       labSend(value, other);

       out = zeros(1,numlabs);
       for i = 1:length(other)
            out(other(i)) = labReceive(other(i));
       end
       out(labindex) = value;