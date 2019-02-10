
function mod = qam_modulate(data, bits_per_symbol, num_output_streams)
num_bits = size(data,1);
num_symbols = ceil(num_bits / bits_per_symbol);
num_symbols_per_stream = ceil(num_symbols / num_output_streams);
mod = zeros(num_symbols_per_stream, num_output_streams);

% extend input data to allow for ceil operations above
nbe = num_symbols_per_stream * num_output_streams * bits_per_symbol;
data = [data;rand(nbe-num_bits,1)>0.5];
%data = [data;rand_simple(nbe-num_bits,1)>0.5];

% set QAM levels for unit signal power
switch bits_per_symbol
	case 1
		a = sqrt(1/2);
	case 2
		a = sqrt(1/2);
	case 4
		a = sqrt(1/10);
	case 6
		a = sqrt(1/42);
	case 8
		a = sqrt(1/170);
	otherwise
		error('bits_per_symbol out of range');
end

% generate a gray code weighted lookup table
bt = [-1;1];
for i=2:bits_per_symbol/2
	nc = -1*ones(size(bt,1),1)*2^(i-1);
	nc = [nc; nc*-1];
	bt = [bt; (bt)*-1];
	bt = [bt, nc];
end
gray = sum(bt*a, 2);

% map input bits to symbols
b = 1;
for i=1:num_symbols_per_stream
	for j=1:num_output_streams
		vi = 0;
		vq = 0;
		if(bits_per_symbol > 1)
			for k=1:bits_per_symbol/2
				vi = 2 * vi + data(b);
				vq = 2 * vq + data(b+1);
				b = b + 2;
			end
		else		% Map BPSK equally to I and Q
			vi = data(b);
			vq = data(b);
			b = b + 1;
		end
		txi = gray(vi+1);
		txq = gray(vq+1);
		mod(i, j) = complex(txi, txq);
	end
end
end
