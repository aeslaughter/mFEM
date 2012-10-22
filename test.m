function test

import mFEM.*

K = Matrix(6,7);

K.add_matrix([22,24; 42,44], [2,4], [2,5]);

K = K.build();

full(K)
