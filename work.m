function work


no(4) = mFEM.elements.base.Node();
no.init(1:4,rand(4,2));

elem = mFEM.elements.Quad4(1,no);

