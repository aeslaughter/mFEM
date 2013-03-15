function delete(obj)
    delete([obj.quad]);
    delete([obj.nodes]);
    clear classes obj.quad obj.nodes;
end