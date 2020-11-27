import struct
import numpy as np
from datetime import datetime


def simple_compress_delta(y):
    """

    Since nearby points should be closely valued, this should cause all
    data points to be represented as 0s. This may allow compression
    algorithms to get a better compression ratio. This doesn't
    actually achieve compression (same number of bytes as the original
    data). However if the original data is 4-byte precision, sometimes
    2-byte precision deltas can be used, so that there is just under a
    compression ratio of 50%.

    """
    d = y - np.roll(y, 1)
    d = d[1:]
    for i in d:
        yield struct.pack(i, '>h')


def simple_compress_double_delta(y):
    """
    Given the sequence

    y0, y1, ..., yn

    take the differences between adjacent elements (ignore edge effects)

    <>, y1 - y0, ..., yn - y(n-1)

    the cumulative sum of the these sequential differences is just
    the most recent term in the original sequence minus the base term
    y0. That is, the differences make a telescoping series.

    For a double delta, you need something for which the cumulative
    sum of the cumulative sum is the original sequence. The following
    sequence satisfies this condition:

    <>, <>, (y2-y1) - (y1-y0), (y3-y2) - (y2-y1),  ...

    This gives linear artifacts when read by the Agilent reader. In
    analogy to a second derivative, the second finite difference may
    tend to have a linear term if integrated directly (integration
    constants are non-zero).

    You might think of doing a Legendre transform compression. Define
    y to be y - x*dy/dx. But that analytical transform is effectively
    what is done here already, except subtracting dy rather than
    dy/dx*x. Also a double Legendre transform isn't very sensible
    since the independent variable becomes the derivative on the first
    transform.

    """
    d = y - np.roll(y, 1)  # backwards difference
    dd = d - np.roll(d, 1)  # second backwards difference
    dd = dd[2:]
    dd *= 32768 / abs(dd).max()  
    for i in dd:
        yield struct.pack('>h', int(i))


def compress_double_delta(f, y, limit=1000):
    """

    Compresses y to double delta representation. Provides the linear
    transform parameters. Most data points should be exact. Fractional
    error in data points should be less than 0.1%.

    """

    header = 32767  # 2^15 - 1

    intercept = y.min()
    y = y - intercept
    slope = y.max() / 2 ** 31
    y /= slope
    d = y - np.roll(y, 1)
    dd = d - np.roll(d, 1)
    dd /= slope
    dd = dd.astype(int)
    y = y.astype(int)

    f.write(struct.pack('>h', header))
    f.write(struct.pack('>h', 0))
    f.write(struct.pack('>I', y[0]))

    f.write(struct.pack('>h', header))
    f.write(struct.pack('>h', 0))
    f.write(struct.pack('>I', y[1]))

    limit_counter = 0

    for x in range(2, len(y)):
        curr = y[x]
        d2 = dd[x]

        # the difference must be less than the maximum value of a signed short
        condition = np.log2(abs(d2) + 1) < 15 and limit_counter < limit

        if condition:
            f.write(struct.pack('>h', d2))
        else:  # start a new head (reference)
            f.write(struct.pack('>h', header))
            # for whatver (possible compression) reason, the start value is represented as 2^33*h + I
            # where h is a short and I is an unsigned long
            f.write(struct.pack('>h', 0))  # assign 0 arbitrarily
            f.write(struct.pack('>I', curr))
            limit_counter = 0

        limit_counter += 1

    return slope, intercept


def encode(filename, ydata, meta=dict(), doe=dict(), xmin=0, xmax=1,
           slope=1, intercept=0, start_data=6144, options_encoding='H'):

    f = open(filename, 'wb')

    # write meta data
    # most characters are ASCII and so the H option just amounts to padding with zeroes
    # I don't know why unsigned shorts rather than "char" are used unles they
    # decided to support unicode at some point

    if options_encoding == 'H':
        m = 2
        echar = 'H'
    else:
        m = 1
        echar = 'B'

    poses = {'sample': 858,
             'description': 1369,
             'method': 2574,
             'operator': 1880,
             'date': 2391,
             'instrument': 2533,
             'inlet': 2492,
             'units': 4172}

    for key in poses:
        try:
            value = meta[key]
            position = poses[key]
        except KeyError:
            continue
        f.seek(position)
        l = len(value)
        f.write(struct.pack('B', l))
        for char in value:
            f.write(struct.pack(f'<{echar}', ord(char)))

    fmt_str = '%d-%b-%y, %H:%M:%S'
    date_str = datetime.now().strftime(fmt_str)
    f.seek(poses['date'])
    l = len(date_str)
    f.write(struct.pack('B', l))
    for char in date_str:
        f.write(struct.pack(f'<{echar}', ord(char)))

    f.seek(252)
    for key in 'sequence', 'vial', 'replicate':
        try:
            f.write(struct.pack('>h', doe[key]))
        except KeyError:
            f.write(struct.pack('>h', 0))

    f.seek(264)
    f.write(struct.pack('>i', start_data // 512 + 1))
    f.seek(282)
    f.write(struct.pack('>f', int(xmin * 60000)))
    f.write(struct.pack('>f', int(xmax * 60000)))

    f.seek(start_data)
    slope, intercept = compress_double_delta(f, ydata)

    f.seek(4724)
    f.write(struct.pack('>d', intercept))
    f.seek(4732)
    f.write(struct.pack('>d', slope))

    f.close()
