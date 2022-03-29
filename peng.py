#!/usr/bin/env python3
# Copyright (C) 2016  Niklas Rosenstein
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from __future__ import division, print_function

__author__ = 'Niklas Rosenstein <rosensteinniklas@gmail.com>'
__version__ = '1.0.0-dev'

import argparse
import gzip, zlib, bz2
import math
import png
import random
import sys
import textwrap

try: import brotli
except ImportError: brotli = None

HEADER_SIZE = 9
MAX_ENCODE_SIZE = 256 * 256 * 3 - HEADER_SIZE


def decompose(num, bytes=3):
  """
  Decompose the value *num* into a list of *bytes* 8-bit values.
  """

  return [(num >> (i * 8)) & 0xFF for i in range(bytes)]


def compose(*values):
  """
  Compose the specified 8-bit *values* into a single integer.
  """

  result = 0
  for i, value in enumerate(values):
    result |= value << (i * 8)
  return result


def compress(data, level=9, method='gz'):
  """
  Compress *data* with the specified *method*. It can be ``'zip'``,
  ``'gz'`` or ``'bz'``. If the :mod:`brotli` module is installed,
  ``'brt'`` is also an accepted method.

  :raise ValueError: If an invalid *method* was supplied.
  """

  if method == 'gz':
    return gzip.compress(data, level)
  elif method == 'zip':
    return zlib.compress(data, level)
  elif method == 'bz':
    return bz2.compress(data, level)
  elif method == 'brt':
    if not brotli:
      raise ImportError('brotli')
    return brotli.compress(data)
  else:
    raise ValueError('invalid method: {0!r}'.format(method))


def decompress(data, method='gz'):
  """
  Decompress *data* previously compressed with the specified
  compression *method*.
  """

  if method == 'gz':
    return gzip.decompress(data)
  elif method == 'zip':
    return zlib.decompress(data)
  elif method == 'bz':
    return zlib.decompress(data)
  elif method == 'brt':
    if not brotli:
      raise ImportError('brotli')
    return brotli.decompress(bytes(data))
  else:
    raise ValueError('invalid method: {0!r}'.format(method))


def encode(data, width, height, scale=1, private=0, fill='repeat', seed=None):
  """
  Encode the binary *data* into an RGB image of the specified *width*
  and *height*. An additional 24-bit unsigned integer can be saved with
  the *private* parameter. The *fill* parameter describes the method that
  is used to generate the padding bytes. Available modes are

  ``'repeat'``
    Repeat the *data* until the padding is complete.
  ``'zero'``
    Fill the padding with zero value components.
  ``'random'``
    Generate random color values. If this mode is used, the *seed*
    parameter will be taken into account.

  :raise ValueError:
    * If *data* can not be fit into an image of the specified width and
      height. Note that the image requires three pixels for header data.
    * If the specified *width*, *height* or *scale* are zero or smaller
      or are larger than 255 (maximum image size).
    * If *private* exceeds the range of a 24-bit unsigned integer.
  :raise TypeError:
    * If *data* is not a :class:`bytes` object.
    * If any of *width*, *height*, *scale* or *private* is not an integer.
  :return:
    A flat-row flat-pixel representation of the RGB image as a
    :class:`bytearray` object.
  """

  if not isinstance(data, bytes):
    raise TypeError('data must be bytes')
  if not isinstance(width, int):
    raise TypeError('width must be int')
  if not isinstance(height, int):
    raise TypeError('height must be int')
  if not isinstance(scale, int):
    raise TypeError('scale must be int')
  if not isinstance(private, int):
    raise TypeError('private must be int')

  if not (0 < width < 256):
    raise ValueError('invalid width: {0}'.format(width))
  if not (0 < height < 256):
    raise ValueError('invalid height: {0}'.format(height))
  if not (0 < scale < 256):
    raise ValueError('invalid scale: {0}'.format(scale))
  if not (0 <= private < 0xFFFFFF):
    raise ValueError('invalid private: {0}'.format(private))
  if fill not in ('repeat', 'zero', 'random'):
    raise ValueError('invalid fill: {0!r}'.format(fill))

  result_size = width * height * 3
  max_data_size = result_size - HEADER_SIZE
  if len(data) > max_data_size:
    raise ValueError('data too large: {0} > {1} bytes'.format(len(data), max_data_size))

  result = bytearray()
  result.extend([width, height, scale])
  result.extend(decompose(private))
  result.extend(decompose(len(data)))
  result.extend(data)

  if fill == 'repeat':
    while len(result) < result_size:
      fill_size = result_size - len(result)
      result.extend(data[:fill_size])
  elif fill == 'zero':
    result.extend(0 for __ in range(max_data_size - len(data)))
  elif fill == 'random':
    if seed is None:
      randint = random.randint
    else:
      randint = random.Random(seed).randint
    result.extend(randint(0, 255) for __ in range(max_data_size - len(data)))
  else:
    assert False

  # Scale the image.
  if scale > 1:
    temp = bytearray()
    for i in range(height * scale):
      i //= scale
      for j in range(width * scale):
        j //= scale
        offset = i * width * 3 + j * 3
        temp.append(result[offset + 0])
        temp.append(result[offset + 1])
        temp.append(result[offset + 2])
    result = temp

  return result


def decode(data):
  """
  Decode the flat-row flat-pixel *data* that was previously encoded
  with :meth:`encoded` back to its original data.

  :raise ValueError:
    * If the encoded data size is larger than the actual available data.
  :return:
    A tuple of ``(data, private)`` where *data* is the decoded
    :class:`bytes` object and *private* is a 32-bit integer that
    was passed to the respective argument of :class:`encode`.
  """

  width, height, scale = data[:3]

  # De-scale the image by using the pixel of the top-left corner
  # of each enlarged pixel.
  result = bytearray()
  for i in range(height):
    i *= scale
    for j in range(width):
      j *= scale
      offset = i * 3 * width * scale + j * 3
      result.extend(data[offset:offset+3])

  private = compose(*result[3:6])
  data_size = compose(*result[6:9])
  if len(result) - HEADER_SIZE < data_size:
    raise ValueError('too few data for encoded data size')

  # Strip the header.
  result = result[HEADER_SIZE:]
  if len(result) < data_size:
    raise ValueError('invalid data size')

  return (result[:data_size], private)


def perfect_fit(data_size, constraint=None):
  """
  Returns a number of pixels for which data with the specified *data_size*
  can be encoded into an image using :meth:`encode`. If *constraint* is
  specified, the size of the encodeded image is constraint in one direction.
  Otherwise, an equirectangular image is assumed. This calculation includes
  the header size of the encoding.
  """

  data_size = math.ceil((data_size + HEADER_SIZE) / 3.0)
  if constraint is None:
    return int(math.ceil(data_size ** 0.5))
  else:
    return int(math.ceil(data_size / float(constraint)))


def main():
  parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent("""
      With `peng` you can encode arbitrary data to a PNG image file. The
      data can optionally be compressed with zip, gzip or bzip2 to reduce
      data size and generate a more colorful result (especially for plain
      text data). The resulting PNG image file has only a slight overhead
      over the original data, that being the PNG header plus 9 bytes of
      the `peng` header information.

      The amount of data that can be encoded with `peng` is limited to
      196,599 bytes.
      """))
  parser.add_argument('input', help="input file or '-' for stdin")
  parser.add_argument('output', help="output file or '-' for stdout")
  parser.add_argument(
    '-d', '--decode',
    action='store_true',
    help="decode input to output")
  parser.add_argument(
    '--compress',
    default=None,
    choices=('gz', 'zip', 'bz', 'brt'),
    help="compression method, omit for no compression")
  parser.add_argument(
    '--scale', type=int, default=1,
    help="post-encode scaling, defaults to 1")
  parser.add_argument('--width', type=int, help="fixed width (< 256)")
  parser.add_argument('--height', type=int, help="fixed height (< 256)")
  parser.add_argument(
    '--fill',
    choices=('repeat', 'zero', 'random'),
    default='repeat',
    help="padding fill method, defaults to repeat")
  parser.add_argument('--seed', type=int, help="seed for the --fill=random mode")
  args = parser.parse_args()

  if args.input == '-':
    infile = sys.__stdin__.buffer
  else:
    infile = open(args.input, 'rb')
  if args.output == '-':
    outfile = sys.__stdout__.buffer
  else:
    outfile = open(args.output, 'wb')

  def check_method(method):
    if method == 'brt' and not brotli:
      parser.error('brotli compression module not installed')

  if args.decode:
    w, h, encoded_data, info = png.Reader(file=infile).read_flat()
    if info['alpha'] or info['bitdepth'] != 8:
      parser.error('input must be 8-bit RGB PNG file')
    data, private = decode(encoded_data)
    if private != 0:
      method = ''.join(chr(c) for c in decompose(private)).strip()
      check_method(method)
      data = decompress(data, method)
    outfile.write(data)
  else:
    check_method(args.compress)
    data = infile.read()
    if len(data) > MAX_ENCODE_SIZE:
      parser.error('maximum encodable data size exceeded: {0} > {1} bytes'.format(len(data), MAX_ENCODE_SIZE))
    if args.compress is None:
      private = 0
    else:
      private = compose(*map(ord, args.compress.ljust(3)))
      data = compress(data, method=args.compress)

    if args.width:
      w = args.width
      if args.height:
        h = args.height
      else:
        h = perfect_fit(len(data), constraint=w)
    elif args.height:
      h = args.height
      w = perfect_fit(len(data), constraint=h)
    else:
      w = h = perfect_fit(len(data))

    try:
      data = encode(data, w, h, args.scale, private, fill=args.fill, seed=args.seed)
    except ValueError as exc:
      parser.error(str(exc))
    writer = png.Writer(w * args.scale, h * args.scale)
    writer.write_array(outfile, data)

  infile.close()
  outfile.close()


if __name__ == '__main__':
  main()