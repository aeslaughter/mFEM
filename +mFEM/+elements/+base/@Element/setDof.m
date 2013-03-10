function setDof(obj)
    for i = 1:length(obj);
        obj(i).dof = [obj(i).nodes.dof];
        obj(i).n_dof = sum([obj(i).nodes.n_dof]);
    end
end