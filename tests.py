import matplotlib.pyplot as plt
from py_chemplexity import expose_decompress_double_delta, read_data
from encode import encode
import struct

# show typical values for start position of data, slope and intercept of linear transform
# with open('sample.CH', 'rb') as f:
#    f.seek(264)
#    offset = (struct.unpack('>i', f.read(4))[0] - 1) * 512  # i is by default 4 byte
#    f.seek(4724)
#    intercept = struct.unpack('>d', f.read(8))
#    f.seek(4732)
#    slope = struct.unpack('>d', f.read(8))
#    print(offset, intercept, slope)

# show the history of the buffer in the decoder (used to make the encoder)
# with open('sample.CH', 'rb') as f:
#        f.seek(264)  # by default the offset is from beginning of file
#        offset = (struct.unpack('>i', f.read(4))[
#                  0] - 1) * 512  # i is by default 4 byte
#        ca, cb, rv, da, wv = expose_decompress_double_delta(f, offset)
#fig, axs = plt.subplots(nrows = 3, ncols = 1)
# axs[0].plot(rv)
# axs[1].plot(da)
# axs[2].plot(wv)
# plt.show()

# test the encoder
data1 = read_data('data/sample.CH')
with open('data/sample.CH', 'rb') as f:
    f.seek(264)
    offset = (struct.unpack('>i', f.read(4))[
              0] - 1) * 512  # i is by default 4 byte

encode('data/psample.CH', data1['tic'], start_data=offset)

data2 = read_data('data/psample.CH')

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
