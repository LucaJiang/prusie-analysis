"""Render the published architecture SVG using librsvg and Cairo."""
import ctypes
from pathlib import Path
import xml.etree.ElementTree as ET

from ctypes.util import find_library


root_dir = Path(__file__).resolve().parents[1]
svg = (root_dir / "03_results/architecture.svg").read_bytes()
root = ET.fromstring(svg)
width = float(root.attrib["width"]) * 0.75
height = float(root.attrib["height"]) * 0.75
output = root_dir / "04_figures_tables/figures/architecture.pdf"
output.parent.mkdir(parents=True, exist_ok=True)
def library(name):
    location = find_library(name)
    if not location:
        raise RuntimeError(f"Missing system library {name}; see environment/README.md")
    return ctypes.CDLL(location)


cairo = library("cairo")
rsvg = library("rsvg-2")
gobject = library("gobject-2.0")


class Rectangle(ctypes.Structure):
    _fields_ = [(name, ctypes.c_double) for name in ("x", "y", "width", "height")]


def signature(library, name, arguments, result):
    function = getattr(library, name)
    function.argtypes = arguments
    function.restype = result
    return function


pointer = ctypes.c_void_p
create_handle = signature(rsvg, "rsvg_handle_new_from_data", [ctypes.c_char_p, ctypes.c_size_t, pointer], pointer)
render = signature(rsvg, "rsvg_handle_render_document", [pointer, pointer, ctypes.POINTER(Rectangle), pointer], ctypes.c_int)
create_surface = signature(cairo, "cairo_pdf_surface_create", [ctypes.c_char_p, ctypes.c_double, ctypes.c_double], pointer)
create_context = signature(cairo, "cairo_create", [pointer], pointer)
show_page = signature(cairo, "cairo_show_page", [pointer], None)
destroy_context = signature(cairo, "cairo_destroy", [pointer], None)
finish_surface = signature(cairo, "cairo_surface_finish", [pointer], None)
surface_status = signature(cairo, "cairo_surface_status", [pointer], ctypes.c_int)
destroy_surface = signature(cairo, "cairo_surface_destroy", [pointer], None)
release_handle = signature(gobject, "g_object_unref", [pointer], None)
handle = create_handle(svg, len(svg), None)
assert handle
surface = create_surface(str(output).encode(), width, height)
context = create_context(surface)
viewport = Rectangle(0, 0, width, height)
assert render(handle, context, ctypes.byref(viewport), None)
show_page(context)
destroy_context(context)
finish_surface(surface)
assert surface_status(surface) == 0
destroy_surface(surface)
release_handle(handle)

print("Rendered architecture.pdf")
