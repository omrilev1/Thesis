function [ demodulated_bits_stream ] = qam_demodulate(modulated_vector, bits_per_tone)
% gray code for QPSK, QAM16 and QAM64:
% QPSK:
%    10 | 11
%    -------
%    00 | 01

% QAM16:
%    0010 | 0110 | 0111 | 0011
%    -------------------------
%    1010 | 1110 | 1111 | 1011
%    -------------------------
%    1000 | 1100 | 1101 | 1001
%    -------------------------
%    0000 | 0100 | 0101 | 0001

% QAM64
%    000010 | 010010 | 010110 | 000110 | 000111 | 010111 | 010011 | 000011
%    ---------------------------------------------------------------------
%    100010 | 110010 | 110110 | 100110 | 100111 | 110111 | 110011 | 100011
%    ---------------------------------------------------------------------
%    101010 | 111010 | 111110 | 101110 | 101111 | 111111 | 111011 | 101011
%    ---------------------------------------------------------------------
%    001010 | 011010 | 011110 | 001110 | 001111 | 011111 | 011011 | 001011
%    ---------------------------------------------------------------------
%    001000 | 011000 | 011100 | 001100 | 001101 | 011101 | 011001 | 001001
%    ---------------------------------------------------------------------
%    101000 | 111000 | 111100 | 101100 | 101101 | 111101 | 111001 | 101001
%    ---------------------------------------------------------------------
%    100000 | 110000 | 110100 | 100100 | 100101 | 110101 | 110001 | 100001
%    ---------------------------------------------------------------------
%    000000 | 010000 | 010100 | 000100 | 000101 | 010101 | 010001 | 000001

switch bits_per_tone
	case 2
		a = sqrt(1/2);
		decision_boundries = 0;
	case 4
		a = sqrt(1/10);
		decision_boundries = [-2*a, 0, 2*a];
	case 6
		a = sqrt(1/42);
		decision_boundries = [-6*a, -4*a, -2*a, 0, 2*a, 4*a, 6*a];
	otherwise
		error('bits_per_symbol out of range');
end

demodulated_bits_stream = 2*ones(length(modulated_vector)*bits_per_tone,1);

for sc=1:length(modulated_vector)
	qam_sc = modulated_vector(sc);
	real_idx = 1;
	while (real(qam_sc) > decision_boundries(real_idx))
		real_idx = real_idx+1;
		if(real_idx>length(decision_boundries));break;end
	end
	imag_idx = 1;
	while (imag(qam_sc) > decision_boundries(imag_idx));
		imag_idx = imag_idx+1;
		if(imag_idx>length(decision_boundries));break;end
	end
	
	switch(bits_per_tone)
		case 2
			% bit 0:
			if (real_idx  == 1);demodulated_bits_stream((sc-1)*bits_per_tone+1) = 0;else demodulated_bits_stream((sc-1)*bits_per_tone+1) = 1;end
			% bit 1:
			if (imag_idx  == 1);demodulated_bits_stream((sc-1)*bits_per_tone+2) = 0;else demodulated_bits_stream((sc-1)*bits_per_tone+2) = 1;end
		case 4
			% bit 0:
			if (ismember(real_idx,[1,2]));demodulated_bits_stream((sc-1)*bits_per_tone+1) = 0;else demodulated_bits_stream((sc-1)*bits_per_tone+1) = 1;end
			% bit 1:
			if (ismember(imag_idx,[1,2]));demodulated_bits_stream((sc-1)*bits_per_tone+2) = 0;else demodulated_bits_stream((sc-1)*bits_per_tone+2) = 1;end
			% bit 2:
			if (ismember(real_idx,[1,4]));demodulated_bits_stream((sc-1)*bits_per_tone+3) = 0;else demodulated_bits_stream((sc-1)*bits_per_tone+3) = 1;end
			% bit 3:
			if (ismember(imag_idx,[1,4]));demodulated_bits_stream((sc-1)*bits_per_tone+4) = 0;else demodulated_bits_stream((sc-1)*bits_per_tone+4) = 1;end
		case 6
			% bit 0:
			if (ismember(real_idx,[1,2,3,4]));demodulated_bits_stream((sc-1)*bits_per_tone+1) = 0;else demodulated_bits_stream((sc-1)*bits_per_tone+1) = 1;end
			% bit 1:
			if (ismember(imag_idx,[1,2,3,4]));demodulated_bits_stream((sc-1)*bits_per_tone+2) = 0;else demodulated_bits_stream((sc-1)*bits_per_tone+2) = 1;end
			% bit 2:
			if (ismember(real_idx,[1,2,7,8]));demodulated_bits_stream((sc-1)*bits_per_tone+3) = 0;else demodulated_bits_stream((sc-1)*bits_per_tone+3) = 1;end
			% bit 3:
			if (ismember(imag_idx,[1,2,7,8]));demodulated_bits_stream((sc-1)*bits_per_tone+4) = 0;else demodulated_bits_stream((sc-1)*bits_per_tone+4) = 1;end
			% bit 4:
			if (ismember(real_idx,[1,4,5,8]));demodulated_bits_stream((sc-1)*bits_per_tone+5) = 0;else demodulated_bits_stream((sc-1)*bits_per_tone+5) = 1;end
			% bit 5:
			if (ismember(imag_idx,[1,4,5,8]));demodulated_bits_stream((sc-1)*bits_per_tone+6) = 0;else demodulated_bits_stream((sc-1)*bits_per_tone+6) = 1;end
			
	end
	
end

end