% WRITEDATA writes data stream on FILENAME with N kbytes.
function writedata(filename, n)

fout = fopen(filename, 'a+');
n = round(n);

display(['Generating random stream of ' num2str(n) ' bytes.']);

% Megabytes
for i = 1:floor(n/(1024*1024))
    str = randi(256, [1, 1024*1024]) - 1; fwrite(fout, str);
end
n = n - 1024*1024*floor(n/(1024*1024));

% Kilobytes
for i = 1:floor(n/1024)
    str = randi(256, [1, 1024]) - 1; fwrite(fout, str);
end
n = n - 1024*floor(n/1024);

% Bytes
str = randi(256, [1, n]) - 1; fwrite(fout, str);

fclose(fout);

end