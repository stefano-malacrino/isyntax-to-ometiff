# !/usr/bin/python

import struct
import sys
from ctypes import *
from ctypes.util import find_library
from enum import IntEnum

libtiff_name = find_library('tiff')
if libtiff_name is not None:
    libtiff = cdll.LoadLibrary(libtiff_name)
else:
    try:
        if 'win' in sys.platform:
            libtiff = cdll.LoadLibrary('libtiff.dll')
        elif sys.platform == 'darwin':
            libtiff = cdll.LoadLibrary('libtiff.dylib')
        else:
            libtiff = cdll.LoadLibrary('libtiff.so')
    except OSError:
        raise ImportError('Failed to load libtiff')

class TIFFTAG(IntEnum):
    SUBFILETYPE = 254
    IMAGEWIDTH = 256
    IMAGELENGTH = 257
    BITSPERSAMPLE = 258
    COMPRESSION = 259
    PHOTOMETRIC = 262
    THRESHHOLDING = 263
    FILLORDER = 266
    IMAGEDESCRIPTION = 270
    MAKE = 271
    MODEL = 272
    STRIPOFFSETS = 273
    ORIENTATION = 274
    SAMPLESPERPIXEL = 277
    ROWSPERSTRIP = 278
    STRIPBYTECOUNTS = 279
    MINSAMPLEVALUE = 280
    MAXSAMPLEVALUE = 281
    XRESOLUTION = 282
    YRESOLUTION = 283
    PLANARCONFIG = 284
    RESOLUTIONUNIT = 296
    SOFTWARE = 305
    DATETIME = 306
    ARTIST = 315
    HOSTCOMPUTER = 316
    COLORMAP = 320
    TILEWIDTH = 322
    TILELENGTH = 323
    TILEOFFSETS = 324
    TILEBYTECOUNTS = 325
    SUBIFD = 330
    EXTRASAMPLES = 338
    SAMPLEFORMAT = 339
    SMINSAMPLEVALUE = 340
    SMAXSAMPLEVALUE = 341
    JPEGTABLES = 347
    YCBCRCOEFFICIENTS = 529
    YCBCRSUBSAMPLING = 530
    YCBCRPOSITIONING = 531
    COPYRIGHT = 33432

class TIFFTAG_PSEUDO(IntEnum):
    JPEGQUALITY = 65537
    JPEGCOLORMODE = 65538
    
class TIFFDATATYPE(IntEnum):
    NOTYPE = 0      # placeholder
    BYTE = 1        # 8-bit unsigned integer
    ASCII = 2       # 8-bit bytes w/ last byte null
    SHORT = 3       # 16-bit unsigned integer
    LONG = 4        # 32-bit unsigned integer
    RATIONAL = 5    # 64-bit unsigned fraction
    SBYTE = 6       # !8-bit signed integer
    UNDEFINED = 7   # !8-bit untyped data
    SSHORT = 8      # !16-bit signed integer
    SLONG = 9       # !32-bit signed integer
    SRATIONAL = 10  # !64-bit signed fraction
    FLOAT = 11      # !32-bit IEEE floating point
    DOUBLE = 12     # !64-bit IEEE floating point
    IFD = 13        # %32-bit unsigned integer (offset)
    LONG8 = 16      # BigTIFF 64-bit unsigned integer
    SLONG8 = 17     # BigTIFF 64-bit signed integer
    IFD8 = 18       # BigTIFF 64-bit unsigned integer (offset)

class COMPRESSION(IntEnum):
    NONE = 1
    CCITTRLE = 2
    CCITT_T4 = 3
    CCITT_T6 = 4
    LZW = 5
    OJPEG = 6
    JPEG = 7
    NEXT = 32766
    CCITTRLEW = 32771
    PACKBITS = 32773
    THUNDERSCAN = 32809
    IT8CTPAD = 32895
    IT8LW = 32896
    IT8MP = 32897
    IT8BL = 32898
    PIXARFILM = 32908
    PIXARLOG = 32909
    DEFLATE = 32946
    ADOBE_DEFLATE = 8
    DCS = 32947
    JBIG = 34661
    SGILOG = 34676
    SGILOG24 = 34677
    JP2000 = 34712
    
class PHOTOMETRIC(IntEnum):
    MINISWHITE = 0
    MINISBLACK = 1
    RGB = 2
    PALETTE = 3
    MASK = 4
    SEPARATED = 5
    YCBCR = 6
    CIELAB = 8
    ICCLAB = 9
    ITULAB = 10
    LOGL = 32844
    LOGLUV = 32845

class EXTRASAMPLE(IntEnum):
    UNSPECIFIED = 0
    ASSOCALPHA = 1
    UNASSALPHA = 2
    
class FILETYPE(IntEnum):
    REDUCEDIMAGE = 1
    PAGE = 2
    MASK = 4
    
class PLANARCONFIG(IntEnum):
    CONTIG = 1
    SEPARATE = 2

class ORIENTATION(IntEnum):
    TOPLEFT = 1
    TOPRIGHT = 2
    BOTRIGHT = 3
    BOTLEFT = 4
    LEFTTOP = 5
    RIGHTTOP = 6
    RIGHTBOT = 7
    LEFTBOT = 8
    
class RESUNIT(IntEnum):
    NONE = 1
    INCH = 2
    CENTIMETER = 3

class JPEGCOLORMODE(IntEnum):
    RAW = 0
    RGB = 1
    
get_ptr = lambda d, t: cast((t*len(d))(*d), POINTER(t))
get_count = lambda d, t: t(len(d))

TIFFTAGS = {
    TIFFTAG.SUBFILETYPE: (lambda d: (c_uint32(d),), lambda d: d.value, lambda: (c_uint32(),)),
    TIFFTAG.IMAGEWIDTH: (lambda d: (c_uint32(d),), lambda d: d.value, lambda: (c_uint32(),)),
    TIFFTAG.IMAGELENGTH: (lambda d: (c_uint32(d),), lambda d: d.value, lambda: (c_uint32(),)),
    TIFFTAG.BITSPERSAMPLE: (lambda d: (c_uint16(d),), lambda d: d.value, lambda: (c_uint16(),)),
    TIFFTAG.COMPRESSION: (lambda d: (c_uint16(d),), lambda d: d.value, lambda: (c_uint16(),)),
    TIFFTAG.PHOTOMETRIC: (lambda d: (c_uint16(d),), lambda d: d.value, lambda: (c_uint16(),)),
    TIFFTAG.THRESHHOLDING: (lambda d: (c_uint16(d),), lambda d: d.value, lambda: (c_uint16(),)),
    TIFFTAG.FILLORDER: (lambda d: (c_uint16(d),), lambda d: d.value, lambda: (c_uint16(),)),
    TIFFTAG.IMAGEDESCRIPTION: (lambda d: (c_char_p(d),), lambda d: d.value, lambda: (c_char_p(),)),
    TIFFTAG.MAKE: (lambda d: (c_char_p(d),), lambda d: d.value, lambda: (c_char_p(),)),
    TIFFTAG.MODEL: (lambda d: (c_char_p(d),), lambda d: d.value, lambda: (c_char_p(),)),
    TIFFTAG.STRIPOFFSETS: (None, lambda d, c: d[:c], lambda: (POINTER(c_ulong)(),)),
    TIFFTAG.ORIENTATION: (lambda d: (c_uint16(d),), lambda d: d.value, lambda: (c_uint16(),)),
    TIFFTAG.SAMPLESPERPIXEL: (lambda d: (c_uint16(d),), lambda d: d.value, lambda: (c_uint16(),)),
    TIFFTAG.ROWSPERSTRIP: (lambda d: (c_uint32(d),), lambda d: d.value, lambda: (c_uint32(),)),
    TIFFTAG.STRIPBYTECOUNTS: (None, lambda d, c: d[:c], lambda: (POINTER(c_ulong)(),)),
    TIFFTAG.MINSAMPLEVALUE: (lambda d: (c_uint16(d),), lambda d: d.value, lambda: (c_uint16(),)),
    TIFFTAG.MAXSAMPLEVALUE: (lambda d: (c_uint16(d),), lambda d: d.value, lambda: (c_uint16(),)),
    TIFFTAG.XRESOLUTION: (lambda d: (c_float(d),), lambda d: d.value, lambda: (c_float(),)),
    TIFFTAG.YRESOLUTION: (lambda d: (c_float(d),), lambda d: d.value, lambda: (c_float(),)),
    TIFFTAG.PLANARCONFIG: (lambda d: (c_uint16(d),), lambda d: d.value, lambda: (c_uint16(),)),
    TIFFTAG.RESOLUTIONUNIT: (lambda d: (c_uint16(d),), lambda d: d.value, lambda: (c_uint16(),)),
    TIFFTAG.SOFTWARE: (lambda d: (c_char_p(d),), lambda d: d.value, lambda: (c_char_p(),)),
    TIFFTAG.DATETIME: (lambda d: (c_char_p(d),), lambda d: d.value, lambda: (c_char_p(),)),
    TIFFTAG.ARTIST: (lambda d: (c_char_p(d),), lambda d: d.value, lambda: (c_char_p(),)),
    TIFFTAG.HOSTCOMPUTER: (lambda d: (c_char_p(d),), lambda d: d.value, lambda: (c_char_p(),)),
    TIFFTAG.COLORMAP: (lambda d: (get_ptr(d, c_uint16),), lambda d: d[:3], lambda: (POINTER(c_uint16)(),)),
    TIFFTAG.EXTRASAMPLES: (lambda d: (get_count(d, c_uint16), get_ptr(d, c_uint16)), lambda c, d: d[:c.value], lambda: (c_uint16(), POINTER(c_uint16)())),
    TIFFTAG.TILELENGTH: (lambda d: (c_uint32(d),), lambda d: d.value, lambda: (c_uint32(),)),
    TIFFTAG.TILEWIDTH: (lambda d: (c_uint32(d),), lambda d: d.value, lambda: (c_uint32(),)),
    TIFFTAG.TILEOFFSETS: (None, lambda d, c: d[:c], lambda: (POINTER(c_ulong)(),)),
    TIFFTAG.TILEBYTECOUNTS: (None, lambda d, c: d[:c], lambda: (POINTER(c_ulong)(),)),
    TIFFTAG.SUBIFD: (lambda d: (get_count(d, c_uint16), get_ptr(d, c_ulong)), lambda c, d: d[:c.value], lambda: (c_uint16(), POINTER(c_ulong)())),
    TIFFTAG.SAMPLEFORMAT: (lambda d: (c_uint16(d),), lambda d: d.value, lambda: (c_uint16(),)),
    TIFFTAG.SMINSAMPLEVALUE: (lambda d: (c_double(d),), lambda d: d.value, lambda: (c_double(),)),
    TIFFTAG.SMAXSAMPLEVALUE: (lambda d: (c_double(d),), lambda d: d.value, lambda: (c_double(),)),
    TIFFTAG.JPEGTABLES: (lambda d: (get_count(d, c_ushort), get_ptr(d, c_void_p)), lambda c, d: d[:c.value], lambda: (c_ushort(), c_void_p())),   
    TIFFTAG.YCBCRCOEFFICIENTS: (lambda d: (get_ptr(d, c_float),), lambda d: d[:3], lambda: (POINTER(c_float)(),)),
    TIFFTAG.YCBCRSUBSAMPLING: (lambda d: (c_uint16(d[0]), c_uint16(d[1])), lambda d: d[:2], lambda: (c_uint16(), c_uint16())),
    TIFFTAG.YCBCRPOSITIONING: (lambda d: (c_uint16(d),), lambda d: d.value, lambda: (c_uint16(),)),
    TIFFTAG.COPYRIGHT: (lambda d: (c_char_p(d),), lambda d: d.value, lambda: (c_char_p(),)),
    TIFFTAG_PSEUDO.JPEGQUALITY: (lambda d: (c_int(d),), None, None),
    TIFFTAG_PSEUDO.JPEGCOLORMODE: (lambda d: (c_int(d),), None, None),
}

class Tiff(c_void_p):     
    
    @classmethod
    def open(cls, filename, mode='r', bigtiff=False):
        mode = mode + '8' if mode == 'w' and bigtiff else mode
        tiff = libtiff.TIFFOpen(str(filename).encode('ascii'), mode.encode('ascii'))
        if tiff.value is None:
            raise TypeError(f'Error opening {filename}')
        return tiff
        
    @property
    def image_width(self):
        return self.get_field(TIFFTAG.IMAGEWIDTH)
        
    @property
    def image_length(self):
        return self.get_field(TIFFTAG.IMAGELENGTH)
        
    @property
    def tile_width(self):
        return self.get_field(TIFFTAG.TILEWIDTH)

    @property
    def tile_length(self):
        return self.get_field(TIFFTAG.TILELENGTH)
    
    @property
    def tiles_across(self):
        return (self.image_width + self.tile_width - 1) // self.tile_width
    
    @property
    def tiles_down(self):
        return (self.image_length + self.tile_length - 1) // self.tile_length
        
    @property
    def tiles_per_image(self):
        return self.tiles_across * self.tiles_down
        
    def close(self):
        libtiff.TIFFClose(self)
        
    def compute_tile(self, x, y, z, sample):
        return libtiff.TIFFComputeTile(self, x, y, z, sample)
        
    def current_directory(self):
        return libtiff.TIFFCurrentDirectory(self)
        
    def fileno(self):
        return libtiff.TIFFFileno(self)
        
    def get_field(self, tag, count=None):
        try:
            _, convert_fn_get, buffers_factory = TIFFTAGS[tag]
        except KeyError:
            raise ValueError(f'Tag {tag} not defined')
        
        buffers = buffers_factory()
        res = libtiff.TIFFGetField(self, ttag_t(tag), *[byref(buff) for buff in buffers])
        
        args = buffers + ((count,) if count is not None else ())

        if res == 1:
            return convert_fn_get(*args)
        else:
            return None
        
    def get_mode(self):
        return libtiff.TIFFGetMode(self)
        
    def is_bigtiff(self):
        with open(self.fileno(), 'rb', closefd=False) as f:
            ptr = f.tell()
            f.seek(0)
            header = f.read(4)
            f.seek(ptr)
            byteorder = {b'II': '<', b'MM': '>', b'EP': '<'}[header[:2]]
            version = struct.unpack(byteorder + 'H', header[2:4])[0]
        return version == 43
        
    def is_tiled(self):
        return bool(libtiff.TIFFIsTiled(self))
    
    def read_directory(self):
        return libtiff.TIFFReadDirectory(self)
        
    def read_raw_tile(self, tile, buff, size):
        res = libtiff.TIFFReadRawTile(self, tile, buff, size)
        if res < 0:
            raise RuntimeError('Error reading raw tile {tile}')
        return res
        
    def set_directory(self, dirnum):
        res = libtiff.TIFFSetDirectory(self, dirnum)
        if res < 0:
            raise RuntimeError('Error setting directory {dirnum}')
        return res

    def set_field(self, tag, value):
        try:
            convert_fn_set, _, _ = TIFFTAGS[tag]
        except KeyError:
            raise ValueError(f'Tag {tag} not defined')

        res = libtiff.TIFFSetField(self, ttag_t(tag), *convert_fn_set(value))
        
        if res != 1:
            raise RuntimeError(f'Could not set tag {tag}')
            
    def tile_size(self):
        return libtiff.TIFFTileSize(self)

    def write_directory(self):
        res = libtiff.TIFFWriteDirectory(self)
        if res != 1:
            raise RuntimeError('Error writing directory')
            
    def write_encoded_tile(self, tile, data, size):
        res = libtiff.TIFFWriteEncodedTile(self, tile, data, size)
        if res < 0:
            raise RuntimeError('Error writing encoded tile {tile}')
        return res
        
    def write_raw_strip(self, strip, data, size):
        res = libtiff.TIFFWriteRawStrip(self, strip, data, size)
        if res < 0:
            raise RuntimeError('Error writing raw strip {strip}')
        return res
            
    def write_raw_tile(self, tile, data, size):
        res = libtiff.TIFFWriteRawTile(self, tile, data, size)
        if res < 0:
            raise RuntimeError('Error writing raw tile {tile}')
        return res

tdir_t = c_uint16
tsample_t = c_uint16
tsize_t = c_int32
ttag_t = c_uint32
ttile_t = c_uint32

libtiff.TIFFComputeTile.restype = ttile_t
libtiff.TIFFComputeTile.argtypes = [Tiff, c_uint32, c_uint32, c_uint32, tsample_t]
libtiff.TIFFClose.restype = None
libtiff.TIFFClose.argtypes = [Tiff]
libtiff.TIFFCurrentDirectory.restype = tdir_t
libtiff.TIFFCurrentDirectory.argtypes = [Tiff]
libtiff.TIFFFieldDataType.restype = c_uint
libtiff.TIFFFieldDataType.argtypes = [c_void_p]
libtiff.TIFFFieldReadCount.restype = c_int
libtiff.TIFFFieldReadCount.argtypes = [c_void_p]
libtiff.TIFFFieldWithTag.restype = c_void_p
libtiff.TIFFFieldWithTag.argtypes = [Tiff, ttag_t]
libtiff.TIFFFileno.restype = c_int
libtiff.TIFFFileno.argtypes = [Tiff]
libtiff.TIFFGetField.restype = c_int
libtiff.TIFFGetMode.restype = c_int
libtiff.TIFFGetMode.argtypes = [Tiff]
libtiff.TIFFIsTiled.restype = c_int
libtiff.TIFFIsTiled.argtypes = [Tiff]
libtiff.TIFFOpen.restype = Tiff
libtiff.TIFFOpen.argtypes = [c_char_p, c_char_p]
libtiff.TIFFReadDirectory.restype = c_int
libtiff.TIFFReadDirectory.argtypes = [Tiff]
libtiff.TIFFReadRawTile.restype = tsize_t
libtiff.TIFFReadRawTile.argtypes = [Tiff, c_uint32, c_void_p, c_int32]
libtiff.TIFFSetDirectory.restype = c_int
libtiff.TIFFSetDirectory.argtypes = [Tiff, tdir_t]
libtiff.TIFFSetField.restype = c_int
libtiff.TIFFTileSize.restype = tsize_t
libtiff.TIFFTileSize.argtypes = [Tiff]
libtiff.TIFFWriteDirectory.restype = c_int
libtiff.TIFFWriteDirectory.argtypes = [Tiff]
libtiff.TIFFWriteRawStrip.restype = tsize_t
libtiff.TIFFWriteRawStrip.argtypes = [Tiff, c_uint32, c_void_p, c_int32]
libtiff.TIFFWriteRawTile.restype = tsize_t
libtiff.TIFFWriteRawTile.argtypes = [Tiff, c_uint32, c_void_p, c_int32]
libtiff.TIFFWriteEncodedTile.restype = tsize_t
libtiff.TIFFWriteEncodedTile.argtypes = [Tiff, c_uint32, c_void_p, c_int32]
