function [A_x,A_y, bin] = Mapping_64_QAM(m,d)

mapping_16_QAM=[3   3  0 0 0 0;
	     3   1  0 0 0 1;
         3  -3  0 0 1 0;
	     3  -1  0 0 1 1;
	     1   3  0 1 0 0;
	     1   1  0 1 0 1;
	     1  -3  0 1 1 0;
	     1  -1  0 1 1 1;
 	    -3   3  1 0 0 0;
	    -3   1  1 0 0 1;
	    -3  -3  1 0 1 0;
        -3  -1  1 0 1 1;
	    -1   3  1 1 0 0;
	    -1   1  1 1 0 1
	    -1  -3  1 1 1 0
	    -1  -1  1 1 1 1];
   
   
A_x = repmat(mapping_16_QAM(:,1)',1,ceil(m/16)).*d/2;
A_y = repmat(mapping_16_QAM(:,2)',1,ceil(m/16)).*d/2;
bin = repmat(mapping_16_QAM(:,3:6),ceil(m/16),1);
end

