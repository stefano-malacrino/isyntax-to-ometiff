import argparse
import glymur
import math
import numpy as np
import sys
import tempfile
import traceback
import uuid

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from queue import Queue
from threading import Thread, Lock
from xml.etree.ElementTree import ElementTree, Element, SubElement, Comment, tostring

from libtiff import *
from pixelengine import PixelEngine
from softwarerendercontext import SoftwareRenderContext
from softwarerenderbackend import SoftwareRenderBackend


class Jpeg2000TileWriter:

    def __init__(self, lock, tiff_handle, tile_dims, psnr):
        self.lock = lock
        self.tiff_handle = tiff_handle
        self.tile_width, self.tile_height = tile_dims
        self.psnr = psnr
        self.tmp = tempfile.NamedTemporaryFile()
        self.buff = glymur.Jp2k(self.tmp.name, shape=(self.tile_height, self.tile_width, 3))
        self.buff_reader = open(self.tmp.name, 'rb')
        
    def __call__(self, tile_index, data):
        data.shape = (self.tile_width, self.tile_height, 3)
        self.buff._write(data, numres=6, psnr=[self.psnr])
        self.buff_reader.seek(0)
        data = self.buff_reader.read()
        data_size = len(data)
        with self.lock:
            res = self.tiff_handle.write_raw_tile(tile_index, data, data_size)
        return res   

class JpegTileWriter:

    def __init__(self, lock, tiff_handle):
        self.lock = lock
        self.tiff_handle = tiff_handle
        
    def __call__(self, tile_index, data):
        data_size = len(data)
        with self.lock:
            res = self.tiff_handle.write_encoded_tile(tile_index, data.ctypes.data, data_size)
        return res
        
def worker_write_tiles(queue_in, lock, tiff_handle, tile_dims, compression, psnr):
    if compression == 'jpeg2000':
        tile_writer = Jpeg2000TileWriter(lock, tiff_handle, tile_dims, psnr)
    else:
        tile_writer = JpegTileWriter(lock, tiff_handle)
            
    exc_list = []

    while True:
        data_in = queue_in.get()
        if data_in is None:
            break
        try:
            x, y, tile = data_in
            tile_index = tiff_handle.compute_tile(x, y, 0, 0)
            tile_writer(tile_index, tile)
        except Exception as exc:
            traceback.print_exc()
            exc_list.append(exc)

    if len(exc_list) > 0:
        raise Exception(exc_list)
        
class TiffWriter:
    
    def __init__(self, conf, pixel_engine):
        self.conf = conf
        self.pixel_engine = pixel_engine

    def __call__(self):
        output_file = self.conf.output_dir / (self.conf.input_file.stem + '.tiff')
        if output_file.exists():
            raise FileExistsError(f'{output_file} already exists')
            
        self.tiff_handle = Tiff.open(output_file, mode='w', bigtiff=self.conf.bigtiff)
        self.pixel_engine['in'].open(str(self.conf.input_file))
        n_levels = self.pixel_engine['in']['WSI'].source_view.num_derived_levels + 1
        if self.conf.start_level >= n_levels:
            self.pixel_engine['in'].close()
            raise ValueError('Invalid start_level Input')
            
        view = self.pixel_engine['in']['WSI'].source_view
        
        pixel_res = view.scale[:2]

        for level in range(self.conf.start_level, n_levels):
            print(f'Level {level}')
            
            tiff_dims = self.get_tiff_dimensions(view, level)
            
            self.set_tiff_tags(tiff_dims, pixel_res, level, n_levels)
            patches = self.get_patches(view, level, tiff_dims)
            self.write_tiles(level, patches)
            self.tiff_handle.write_directory()

        images = [self.pixel_engine['in'][i].image_type for i in range(self.pixel_engine['in'].num_images)]
        
        if self.conf.macro and 'MACROIMAGE' in images:
            self.write_additional_image('MACROIMAGE')
        
        if self.conf.label and 'LABELIMAGE' in images:
            self.write_additional_image('LABELIMAGE')
        
        self.tiff_handle.close()
        self.pixel_engine['in'].close()
        
    def round_to_mul(self, x, mul):
        return mul * round(x / mul)
        
    def get_tiff_dimensions(self, view, level):
        x_step, x_end = view.dimension_ranges(level)[0][1:]
        y_step, y_end = view.dimension_ranges(level)[1][1:]
        tiff_width = self.round_to_mul(x_end + 1, x_step) // x_step
        tiff_height = self.round_to_mul(y_end + 1, y_step) // y_step
        return tiff_width, tiff_height
     
    def patch_within_data_envelopes(self, patch, data_envelopes):
        
        x_min, x_max, y_min, y_max = patch[:4]
        
        for bounding_box in data_envelopes:
            bound_x_min, bound_x_max, bound_y_min, bound_y_max = bounding_box
            outside = (x_max < bound_x_min) or (bound_x_max < x_min) \
                or (y_max < bound_y_min) or (bound_y_max < y_min)
            if not outside:
                return True
        return False

    def write_tiles(self, level, patches):
    
        n_workers = min(self.conf.n_workers, len(patches))
        queue_in = Queue(self.conf.queue_size)
        lock = Lock()

        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            futures = [
                executor.submit(
                    worker_write_tiles,
                    queue_in, lock,
                    self.tiff_handle,
                    (self.conf.tiff_tile_width, self.conf.tiff_tile_height),
                    self.conf.compression, self.conf.psnr) for _ in range(n_workers)]
        
            view = self.pixel_engine['in']['WSI'].source_view
            trunc_level = {0: [0, 0, 0]}
            view.truncation(False, False, trunc_level)
            
            x_step = view.dimension_ranges(level)[0][1]
            y_step = view.dimension_ranges(level)[1][1]
            
            data_envelopes = view.data_envelopes(level)
            regions = view.request_regions(patches, data_envelopes, True, [255,255,255],
                        self.pixel_engine.BufferType(0))
            buff_size = self.conf.tiff_tile_width * self.conf.tiff_tile_height * 3
            try:
                while regions:
                    for region in self.pixel_engine.wait_any():
                        regions.remove(region)
                        patch = np.empty(buff_size, dtype=np.uint8)
                        region.get(patch)
                        x, y = self.round_to_mul(region.range[0], x_step) // x_step, \
                            self.round_to_mul(region.range[2], y_step) // y_step
                        queue_in.put((x, y, patch))
            except:
                raise
            finally:
                for _ in range(n_workers):
                    queue_in.put(None)
                [future.result() for future in futures]

    def get_additional_image_metadata(self, image_type, root, ifd):
        image_name = {'MACROIMAGE': 'MACRO', 'LABELIMAGE': 'LABEL'}[image_type]
        view = self.pixel_engine['in'][image_type].source_view
        x_end = view.dimension_ranges(0)[0][2]
        y_end = view.dimension_ranges(0)[1][2]
        width = int(x_end + 1)
        height = int(y_end + 1)

        image = SubElement(root, 'Image', {'ID':f'Image:{ifd}', 'Name':image_name})
        pixels = SubElement(image, 'Pixels', {
            'BigEndian':'false',
            'DimensionOrder':'XYCZT',
            'ID':'Pixels:0',
            'Interleaved':'true',
            'SignificantBits':'8',
            'SizeC':'3',
            'SizeZ':'1',
            'SizeT':'1',
            'SizeX':str(width),
            'SizeY':str(height),
            'Type':'uint8'
            })
        channel = SubElement(pixels, 'Channel', {
            'ID':'Channel:0:0',
            'SamplesPerPixel':'3'
            })
        tiff_data = SubElement(pixels, 'TiffData', {'IFD':str(ifd), 'PlaneCount':'1'})

    def get_ome_metadata(self, tiff_dims, pixel_res):
        root = Element('OME', {
            'xmlns':'http://www.openmicroscopy.org/Schemas/OME/2016-06',
            'xmlns:xsi':'http://www.w3.org/2001/XMLSchema-instance',
            'Creator':'isyntax2tiff',
            'UUID':'urn:uuid:{}'.format(uuid.uuid1()),
            'xsi:schemaLocation':'http://www.openmicroscopy.org/Schemas/OME/2016-06 http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd'
            })
        image = SubElement(root, 'Image', {'ID':'Image:0', 'Name':'WSI'})
        pixels = SubElement(image, 'Pixels', {
            'BigEndian':'false',
            'DimensionOrder':'XYCZT',
            'ID':'Pixels:0',
            'PhysicalSizeX':str(pixel_res[0]),
            'PhysicalSizeY':str(pixel_res[1]),
            'PhysicalSizeXUnit':'\u00B5m',
            'PhysicalSizeYUnit':'\u00B5m',    
            'Interleaved':'true',
            'SignificantBits':'8',
            'SizeC':'3',
            'SizeZ':'1',
            'SizeT':'1',
            'SizeX':str(tiff_dims[0]),
            'SizeY':str(tiff_dims[1]),
            'Type':'uint8'
            })
        channel = SubElement(pixels, 'Channel', {
            'ID':'Channel:0:0',
            'SamplesPerPixel':'3'
            })
        tiff_data = SubElement(pixels, 'TiffData', {'IFD':'0', 'PlaneCount':'1'})
        ifd = 1
        
        if self.conf.macro:
            self.get_additional_image_metadata('MACROIMAGE', root, ifd)
            ifd += 1
            
        if self.conf.label:
            self.get_additional_image_metadata('LABELIMAGE', root, ifd)
    
        xml_declaration = b'<?xml version="1.0" encoding="UTF-8"?>\n'
        xml_root = tostring(root, encoding='UTF-8', method='xml')
        metadata = xml_declaration + xml_root

        return metadata

    def set_tiff_tags(self, tiff_dims, pixel_res, level, n_levels):
        if level > self.conf.start_level:
            self.tiff_handle.set_field(TIFFTAG.SUBFILETYPE, FILETYPE.REDUCEDIMAGE)
        else:
            self.tiff_handle.set_field(TIFFTAG.SOFTWARE, b'isyntax2tiff')
            self.tiff_handle.set_field(TIFFTAG.XRESOLUTION, 1e4 / pixel_res[0])
            self.tiff_handle.set_field(TIFFTAG.YRESOLUTION, 1e4 / pixel_res[1])
            if not self.conf.no_ome:
                # Use subIFDs to store sub-resolutions following OME-TIFF specification
                offsets = [0 for _ in range(level+1, n_levels)]
                self.tiff_handle.set_field(TIFFTAG.SUBIFD, offsets)
                metadata = self.get_ome_metadata(tiff_dims, pixel_res)
                self.tiff_handle.set_field(TIFFTAG.IMAGEDESCRIPTION, metadata)
            

        self.tiff_handle.set_field(TIFFTAG.IMAGEWIDTH, tiff_dims[0])
        self.tiff_handle.set_field(TIFFTAG.IMAGELENGTH, tiff_dims[1])
        self.tiff_handle.set_field(TIFFTAG.RESOLUTIONUNIT, RESUNIT.CENTIMETER)
        self.tiff_handle.set_field(TIFFTAG.TILEWIDTH, self.conf.tiff_tile_width)
        self.tiff_handle.set_field(TIFFTAG.TILELENGTH, self.conf.tiff_tile_height)
        self.tiff_handle.set_field(TIFFTAG.BITSPERSAMPLE, 8)
        self.tiff_handle.set_field(TIFFTAG.SAMPLESPERPIXEL, 3)
        self.tiff_handle.set_field(TIFFTAG.PLANARCONFIG, PLANARCONFIG.CONTIG)
        self.tiff_handle.set_field(TIFFTAG.ORIENTATION, ORIENTATION.TOPLEFT)
        if self.conf.compression == 'jpeg2000':
            self.tiff_handle.set_field(TIFFTAG.COMPRESSION, COMPRESSION.JP2000)
            self.tiff_handle.set_field(TIFFTAG.PHOTOMETRIC, PHOTOMETRIC.RGB)
        elif self.conf.compression == 'jpeg':
            self.tiff_handle.set_field(TIFFTAG.COMPRESSION, COMPRESSION.JPEG)
            self.tiff_handle.set_field(TIFFTAG_PSEUDO.JPEGQUALITY, self.conf.jpeg_quality)
            if self.conf.crsub is None:
                self.tiff_handle.set_field(TIFFTAG.PHOTOMETRIC, PHOTOMETRIC.RGB)
            else:
                self.tiff_handle.set_field(TIFFTAG.PHOTOMETRIC, PHOTOMETRIC.YCBCR)
                self.tiff_handle.set_field(TIFFTAG.YCBCRSUBSAMPLING, tuple(self.conf.crsub))
                self.tiff_handle.set_field(TIFFTAG_PSEUDO.JPEGCOLORMODE, JPEGCOLORMODE.RGB)

    def get_patches(self, view, level, tiff_dims):
        x_step, x_end = view.dimension_ranges(level)[0][1:]
        y_step, y_end = view.dimension_ranges(level)[1][1:]
        
        patch_width = self.conf.tiff_tile_width * x_step
        patch_height = self.conf.tiff_tile_height * y_step
        
        if self.conf.sparse:
            data_envelopes = view.data_envelopes(level).as_rectangles()
            
        patches = []

        for y_patch_start in range(0, y_end, patch_height):
            y_patch_end = y_patch_start + patch_height
            for x_patch_start in range(0, x_end, patch_width):
                x_patch_end = x_patch_start + patch_width
                patch = [x_patch_start, x_patch_end - x_step,
                         y_patch_start, y_patch_end - y_step, level]

                if (not self.conf.sparse) or (self.patch_within_data_envelopes(patch, data_envelopes)):
                    patches.append(patch)

        return patches
        
    def write_additional_image(self, image_type):
        view = self.pixel_engine['in'][image_type].source_view
        width, height = self.get_tiff_dimensions(view, 0)

        self.tiff_handle.set_field(TIFFTAG.IMAGEWIDTH, width)
        self.tiff_handle.set_field(TIFFTAG.IMAGELENGTH, height)
        self.tiff_handle.set_field(TIFFTAG.BITSPERSAMPLE, 8)
        self.tiff_handle.set_field(TIFFTAG.SAMPLESPERPIXEL, 3)
        self.tiff_handle.set_field(TIFFTAG.PLANARCONFIG, PLANARCONFIG.CONTIG)
        self.tiff_handle.set_field(TIFFTAG.ORIENTATION, ORIENTATION.TOPLEFT)
        self.tiff_handle.set_field(TIFFTAG.COMPRESSION, COMPRESSION.JPEG)
        self.tiff_handle.set_field(TIFFTAG.PHOTOMETRIC, PHOTOMETRIC.RGB)
        
        data = np.asarray(self.pixel_engine['in'][image_type].image_data)
        data_size = len(data)
        
        self.tiff_handle.write_raw_strip(0, data.ctypes.data, data_size)
        self.tiff_handle.write_directory()
        

def main(args):
    render_context = SoftwareRenderContext()
    render_backend = SoftwareRenderBackend()
    pixel_engine = PixelEngine(render_backend, render_context)    
    tiff_writer = TiffWriter(args, pixel_engine)()

def non_neg_int(value):
    try:
        ivalue = int(value)
        if ivalue < 0:
            raise argparse.ArgumentError()
    except:
        raise argparse.ArgumentError('Value must be non negative integer')
    return ivalue
    
def mul_16(value):
    try:
        value = non_neg_int(value)
        if value == 0 or value % 16 != 0:
            raise argparse.ArgumentError()
    except:
        raise argparse.ArgumentError('Value must be non negative multiple of 16')
    return value
    
def jpeg_q(value):
    try:
        value = non_neg_int(value)
        if value < 1 or value > 100:
            raise argparse.ArgumentError()
    except:
        raise argparse.ArgumentError('Value must be integer between 1-100')
    return value 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_file', type=Path, help='Input image file')
    parser.add_argument('output_dir', type=Path, help='Output directory')
    parser.add_argument('--tile-w', dest='tiff_tile_width', type=mul_16, default=512, help='Tile width')
    parser.add_argument('--tile-h', dest='tiff_tile_height', type=mul_16, default=512, help='Tile height')
    parser.add_argument('--bigtiff', action='store_true', help='Use BigTIFF format')
    parser.add_argument('--sparse', action='store_true', help='Use sparse format')
    parser.add_argument('--start-level', dest='start_level', type=non_neg_int, default=0, help='Starting resolution level')
    parser.add_argument('--compression', choices=['jpeg', 'jpeg2000'], default='jpeg', help='Compression')
    parser.add_argument('--Q', dest='jpeg_quality', type=jpeg_q, default=80, help='Compression quality (JPEG only)')
    parser.add_argument('--crsub', type=int, nargs=2, choices=[1,2,4], help='Chroma subsampling factors (JPEG only)')
    parser.add_argument('--psnr', type=non_neg_int, default=40, help='PSNR (JPEG2000 only)')
    parser.add_argument('--macro', action='store_true', help='Include macro image')
    parser.add_argument('--label', action='store_true', help='Include label image')
    parser.add_argument('--n-workers', dest='n_workers', type=non_neg_int, default=8, help='Number of workers')
    parser.add_argument('--queue-size', dest='queue_size', type=non_neg_int, default=512, help='Workers queue size')
    parser.add_argument('--no-ome', dest='no_ome', action='store_true', help='Convert as regular tiff')
    args = parser.parse_args()
    main(args)
