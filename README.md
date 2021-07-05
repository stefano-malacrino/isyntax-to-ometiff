# isyntax-to-ometiff
Python script to convert Philips iSyntax files to [OME-TIFF](https://docs.openmicroscopy.org/ome-model/6.2.0/ome-tiff/specification.html).

## Requirements
* Python 3.6+
* numpy
* glymur
* libtiff 4.0.10+
* OpenJPEG 2.3.1+
* [Philips iSyntax SDK 2.0](https://www.openpathology.philips.com)

## Usage
Basic usage is
```
python isyntax_to_tiff.py /path/to/input/file.isyntax /path/to/output/dir
```
See `python isyntax_to_tiff.py --help` for the list of optional parameters and their default values.
