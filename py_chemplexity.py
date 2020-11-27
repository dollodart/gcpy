from datetime import datetime
import struct
import numpy as np

def read_double_array(f, offset):
    end = f.seek(0, 2)
    f.seek(offset)
    signal = []
    while f.tell() < end:
        signal.append(struct.unpack('>d', f.read(8))[0])
    return signal


def decompress_double_delta(f, offset):
    end = f.seek(0, 2)
    prefix = 32767  # 2^15 - 1
    multiplier = 4294967296  # 2^33
    f.seek(offset)
    signal = []
    read_val = diff_acc = write_val = 0

    while f.tell() < end:
        read_val = struct.unpack('>h', f.read(2))[0]
        if read_val != prefix:
            diff_acc += read_val
            write_val += diff_acc
        else:
            write_val = struct.unpack('>h', f.read(2))[0] * multiplier
            write_val += struct.unpack('>I', f.read(4))[0]
            diff_acc = 0

        signal.append(write_val)

    return signal


def expose_decompress_double_delta(f, offset):
    end = f.seek(0, 2)
    prefix = 32767  # 2^15 - 1
    multiplier = 4294967296  # 2^33
    f.seek(offset)
    signal = []
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


def decompress_delta(f, offset):
    """Untested."""
    end = f.seek(0, 2)
    start = f.seek(0, 0)

    inner_prefix = -32678
    # in twos complement for a 16 bit
    # -32678 = ~2**15 + 1
    # 1011111111111111 + 1 = 1100000000000000
    comparison = 4095
    # 4095 = 2**12 - 1
    # 0000111111111111

    signal = []
    read_outer = write_val = read_inner = back_store = 0

    while f.tell() < end:
        read_outer = struct.unpack('>h', f.read(2))[0]
        write_val = back_store
        if read_outer == comparison:  # read_outer >> 12 == 0
            break
        for _ in range(
                read_outer & comparison):  # not sure if range is correct here, no test data file
            read_inner = struct.unpack('>h', f.read(2))[0]
            if read_inner != inner_prefix:
                write_val += read_inner
            else:
                write_val = struct.unpack('>i', f.read(4))
            signal.append(write_val)
        back_store = read_inner

    return signal


def file_info(f, data, options):
    """

    Seeks to absolute file positions based on the known location of data
    fields in the file type and version, reads the data length, and then
    unpacks the data.

    """

    if options['encoding'] == 'H':
        m = 2
        char = 'H'
    else:
        m = 1
        char = 'B'

    f.seek(options['offset']['sample'])
    # one byte cannot have (byte) endianness
    val = struct.unpack('B', f.read(1))[0]
    data['sample']['name'] = struct.unpack('<' + char * val, f.read(m * val))

    f.seek(options['offset']['description'])
    val = struct.unpack('B', f.read(1))[0]
    data['sample']['description'] = struct.unpack(
        '<' + char * val, f.read(m * val))

    f.seek(options['offset']['method'])
    val = struct.unpack('B', f.read(1))[0]
    data['method']['name'] = struct.unpack('<' + char * val, f.read(m * val))

    f.seek(options['offset']['operator'])
    val = struct.unpack('B', f.read(1))[0]
    data['method']['operator'] = struct.unpack(
        '<' + char * val, f.read(m * val))

    f.seek(252)
    data['sample']['sequence'] = struct.unpack('>h', f.read(2))
    data['sample']['vial'] = struct.unpack('>h', f.read(2))
    data['sample']['replicate'] = struct.unpack('>h', f.read(2))

    f.seek(options['offset']['date'])
    val = struct.unpack('B', f.read(1))[0]
    date = struct.unpack('<' + char * val, f.read(m * val))

    data['method']['date'] = date
    data['method']['time'] = date

    f.seek(options['offset']['instrument'])
    val = struct.unpack('B', f.read(1))[0]
    data['instrument']['name'] = struct.unpack(
        '<' + char * val, f.read(m * val))

    f.seek(options['offset']['units'])
    val = struct.unpack('B', f.read(1))[0]
    data['instrument']['name'] = struct.unpack(
        '<' + char * val, f.read(m * val))

    return data


def read_data(f, version='181'):
    """

    Reads the data for an FID chromatogram in the given filename and
    returns a dictionary containing the time, signal (tics), and method
    information.

    Nominally all versions of Agilent CH are supported, only version 181
    has been tested, see the original pychemplexity MATLAB package if
    errors are encountered or for other file extensions.

    """

    # there is a more elegant way of constructing an empty dictionary structure
    # it's also not necessary since the dictionaries can be dynamically made
    f = open(f, 'rb')
    sample = ['name', 'description', 'sequence', 'vial', 'replicate']
    sample = dict(zip(sample, [None] * len(sample)))
    method = ['operator', 'name', 'date', 'time']
    method = dict(zip(method, [None] * len(method)))
    instrument = ['name', 'units']
    instrument = dict(zip(instrument, [None] * len(instrument)))
    data = {'sample': sample, 'method': method, 'instrument': instrument}
    offset = [
        'sample',
        'description',
        'method',
        'operator',
        'date',
        'instrument',
        'units',
        'inlet',
        'tic',
        'xic',
    ]
    offset = dict(zip(offset, [None] * len(offset)))
    options = {'offset': offset, 'scans': None, 'version': version}

    if options['version'] in ['8', '81']:
        options['encoding'] = 'B'
    else:
        options['encoding'] = 'H'

    if options['version'] in ['179', '181']:

        options['offset']['sample'] = 858
        options['offset']['description'] = 1369
        options['offset']['method'] = 2574
        options['offset']['operator'] = 1880
        options['offset']['date'] = 2391
        options['offset']['instrument'] = 2533
        options['offset']['inlet'] = 2492
        options['offset']['units'] = 4172

        f.seek(264)  # by default the offset is from beginning of file
        offset = (struct.unpack('>i', f.read(4))[
                  0] - 1) * 512  # i is by default 4 byte

        f.seek(282)
        xmin = struct.unpack('>f', f.read(4))[
            0] / 60000.  # f in struct is by default 4 byte, need to specify big endian
        xmax = struct.unpack('>f', f.read(4))[0] / 60000.

        f.seek(4724)
        intercept = struct.unpack('>d', f.read(8))
        f.seek(4732)
        slope = struct.unpack('>d', f.read(8))

        data = file_info(f, data, options)
        if options['version'] == '181':
            # default dtype of numpy array is float64, and the output is 64 bit
            data['tic'] = np.array(decompress_double_delta(f, offset))
        else:
            data['tic'] = np.array(read_double_array(f, offset))
        data['tic'] = data['tic'] * slope + intercept
        data['time'] = np.linspace(xmin, xmax, len(data['tic']))

    elif options['version'] == '81':

        options['offset']['sample'] = 24
        options['offset']['description'] = 86
        options['offset']['method'] = 228
        options['offset']['operator'] = 148
        options['offset']['date'] = 178
        options['offset']['instrument'] = 218
        options['offset']['inlet'] = 208
        options['offset']['units'] = 580

        f.seek(264)

        offset = (struct.unpack('>i', f.read(4))[
                  0] - 1) * 512  # i is by default 4 byte

        data = file_info(f, data, options)
        data['tic'] = decompress_double_delta(f, offset)

        f.seek(282)
        xmin = struct.unpack('>f', f.read(4))[
            0] / 60000.  # f in struct is by default 4 byte, need to specify big endian
        xmax = struct.unpack('>f', f.read(4))[0] / 60000.
        data['time'] = np.linspace(xmin, xmax, len(data['tic']))

        f.seek(636)
        intercept = struct.unpack('>d', f.read(8))
        f.seek(644)
        slope = struct.unpack('>d', f.read(8))

        data['tic'] = data['tic'] * slope + intercept

    elif options['version'] == '8':

        options['offset']['sample'] = 24
        options['offset']['description'] = 86
        options['offset']['method'] = 228
        options['offset']['operator'] = 148
        options['offset']['date'] = 178
        options['offset']['instrument'] = 218
        options['offset']['inlet'] = 208
        options['offset']['units'] = 580

        f.seek(264)

        offset = (struct.unpack('>i', f.read(4))[
                  0] - 1) * 512  # i is by default 4 byte

        data = file_info(f, data, options)
        data['tic'] = decompress_delta(f, offset)

        f.seek(282)
        xmin = struct.unpack('>i', f.read(4))[0] / 60000.
        xmax = struct.unpack('>i', f.read(4))[0] / 60000.
        data['time'] = np.linspace(xmin, xmax, len(data['tic']))

        f.seek(542)
        header = struct.unpack('>i', f.read(4))[0]

        if header in [1, 2, 3]:
            data['tic'] = np.array(data['tic']) * 1.3321110047553
        else:
            f.seek(636)
            intercept = struct.unpack('>d', f.read(8))[0]
            f.seek(644)
            slope = struct.unpack('>d', f.read(8))[0]
            data['tic'] = np.array(data['tic']) * slope + intercept

    data['method']['name'] = ''.join(
        [chr(x) for x in data['method']['name']])
    data['sample']['name'] = ''.join(
        [chr(x) for x in data['sample']['name']])
    data['method']['operator'] = ''.join(
        [chr(x) for x in data['method']['operator']])
    data['method']['date'] = datetime.strptime(
        ''.join([chr(x) for x in data['method']['date']]), '%d-%b-%y, %H:%M:%S').date()
    data['method']['time'] = datetime.strptime(
        ''.join([chr(x) for x in data['method']['time']]), '%d-%b-%y, %H:%M:%S').time()

    f.close()
    return data
