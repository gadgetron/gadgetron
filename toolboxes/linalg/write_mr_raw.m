function write_mr_raw(input,filename)

msize = size(input);

fid = fopen(filename,'wb');
fwrite(fid,length(msize),'int32');
fwrite(fid,msize,'int32');

out = zeros(prod(msize)*2,1);

out(1:2:end) = real(input);
out(2:2:end) = imag(input);

fwrite(fid,out,'float32');

fclose(fid);

return