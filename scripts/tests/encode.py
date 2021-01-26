import matplotlib.pyplot as plt
from gcpy.decode import read_data
from gcpy.encode import encode
import struct

def expose_decompress_double_delta(f, offset):
    end = f.seek(0, 2)
    prefix = 32767  # 2^15 - 1
    multiplier = 4294967296  # 2^33
    f.seek(offset)
    read_val = diff_acc = write_val = 0

    count_a = count_b = 0
    diff_acc = []
    write_val = []
    read_val = []
    while f.tell() < end:
        read_val.append(struct.unpack('>h', f.read(2))[0])
        if read_val[-1] != prefix:  # 2^15 - 1
            diff_acc.append(diff_acc[-1] + read_val[-1])
            write_val.append(diff_acc[-1] + write_val[-1])
            count_a += 1
        else:
            wv = struct.unpack('>h', f.read(2))[0] * multiplier
            write_val.append(wv + struct.unpack('>I', f.read(4))[0])
            diff_acc.append(0)
            count_b += 1

    return count_a, count_b, read_val, diff_acc, write_val

# show typical values for start position of data, slope and intercept of linear transform
with open('../data/sample.CH', 'rb') as f:
    f.seek(264)
    offset = (struct.unpack('>i', f.read(4))[0] - 1) * 512  # i is by default 4 byte
    f.seek(4724)
    intercept = struct.unpack('>d', f.read(8))
    f.seek(4732)
    slope = struct.unpack('>d', f.read(8))
    print(offset, intercept, slope)

# show the history of the buffer in the decoder (used to make the encoder)
with open('../data/sample.CH', 'rb') as f:
    f.seek(264)  # by default the offset is from beginning of file
    offset = (struct.unpack('>i', f.read(4))[
              0] - 1) * 512  # i is by default 4 byte
    ca, cb, rv, da, wv = expose_decompress_double_delta(f, offset)
fig, axs = plt.subplots(nrows = 3, ncols = 1)
axs[0].plot(rv)
axs[1].plot(da)
axs[2].plot(wv)
plt.show()

# test the encoder
data1 = read_data('../data/sample.CH')
with open('../data/sample.CH', 'rb') as f:
    f.seek(264)
    offset = (struct.unpack('>i', f.read(4))[
              0] - 1) * 512  # i is by default 4 byte

encode('../data/psample.CH', data1['tic'], start_data=offset)

data2 = read_data('../data/psample.CH')

plt.figure()
plt.plot(data1['tic'], 'o-', label='orig')
plt.plot(data2['tic'], 'o-', label='reencoded')
plt.legend()
plt.figure()
plt.plot(100 * (data2['tic'] - data1['tic']) /
         data1['tic'], 'o-', label='residual')
# fractional residual < 1 % at all points (should be zero)
plt.ylabel('fractional residual in %')
plt.show()
