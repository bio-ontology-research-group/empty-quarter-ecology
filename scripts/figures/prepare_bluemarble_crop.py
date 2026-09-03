#!/usr/bin/env python3
"""Crop the NASA Blue Marble Next Generation tile to the Figure 1 map extent.

Source
------
NASA Earth Observatory, Blue Marble Next Generation with topography and
bathymetry, July 2004, tile C1 (0-90 degrees E, 0-90 degrees N, 21,600 x
21,600 pixels, 240 pixels per degree, equirectangular):

    https://eoimages.gsfc.nasa.gov/images/imagerecords/73000/73751/
        world.topo.bathy.200407.3x21600x21600.C1.jpg

The tile is a NASA product in the public domain (Reto Stoeckli, NASA Earth
Observatory; MODIS data).  It is not redistributed in this repository; this
script records its SHA-256 and reproduces the committed crop
``metadata/geodata/bluemarble_arabia_200407_120ppd.png`` from it.

Procedure
---------
1. Verify the tile checksum.
2. Crop longitude 44-57 E and latitude 16-25 N (3,120 x 2,160 pixels).
3. Reduce by an exact factor of two with a Lanczos filter to 120 pixels per
   degree (1,560 x 1,080 pixels), so that the embedded raster stays below
   2 MB while exceeding 400 dots per inch at the printed panel width.
4. Write an optimised PNG and a JSON sidecar with the geographic extent.

Run with the pinned environment (Pillow is a Matplotlib dependency).
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from PIL import Image

SOURCE_URL = (
    "https://eoimages.gsfc.nasa.gov/images/imagerecords/73000/73751/"
    "world.topo.bathy.200407.3x21600x21600.C1.jpg"
)
SOURCE_SHA256 = "ee8490ab1eb35d620d8d1ad8e69b3234c0b050e4eddb80e7232a2d165e475aa0"
TILE_LON0, TILE_LAT1 = 0.0, 90.0  # upper-left corner of tile C1
PIXELS_PER_DEGREE = 240
EXTENT = {"lon_min": 44.0, "lon_max": 57.0, "lat_min": 16.0, "lat_max": 25.0}
REDUCTION = 2


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tile", type=Path, help="downloaded Blue Marble tile C1 (JPEG)")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[2]
        / "metadata/geodata/bluemarble_arabia_200407_120ppd.png",
    )
    args = parser.parse_args()

    observed = sha256(args.tile)
    if observed != SOURCE_SHA256:
        raise SystemExit(f"Tile checksum {observed} does not match {SOURCE_SHA256}")

    Image.MAX_IMAGE_PIXELS = None
    tile = Image.open(args.tile)
    if tile.size != (21600, 21600):
        raise SystemExit(f"Unexpected tile size {tile.size}")
    left = int((EXTENT["lon_min"] - TILE_LON0) * PIXELS_PER_DEGREE)
    right = int((EXTENT["lon_max"] - TILE_LON0) * PIXELS_PER_DEGREE)
    top = int((TILE_LAT1 - EXTENT["lat_max"]) * PIXELS_PER_DEGREE)
    bottom = int((TILE_LAT1 - EXTENT["lat_min"]) * PIXELS_PER_DEGREE)
    crop = tile.crop((left, top, right, bottom))
    reduced = crop.resize(
        (crop.size[0] // REDUCTION, crop.size[1] // REDUCTION), Image.LANCZOS
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    reduced.save(args.output, format="PNG", optimize=True)
    sidecar = {
        "source_url": SOURCE_URL,
        "source_sha256": SOURCE_SHA256,
        "source_product": (
            "NASA Blue Marble Next Generation with topography and bathymetry, "
            "July 2004, tile C1; public domain (NASA Earth Observatory)"
        ),
        "crop_pixels_in_tile": [left, top, right, bottom],
        "extent_degrees": EXTENT,
        "pixels_per_degree": PIXELS_PER_DEGREE // REDUCTION,
        "size_pixels": list(reduced.size),
        "resample": "Lanczos, exact factor 2",
        "output_sha256": sha256(args.output),
    }
    args.output.with_suffix(".json").write_text(
        json.dumps(sidecar, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(sidecar, indent=2))


if __name__ == "__main__":
    main()
