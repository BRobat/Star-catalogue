# Gaia Star Atlas

A browser-based 3D map of 118,658 stars around the Sun, positioned with Gaia DR3 distances. It runs in a single HTML page using three.js; there is no build step.

## Features

- **3D star field**: galactic coordinates with the Sun at the origin. Each star's size and brightness are its apparent magnitude as seen from the camera's current position.
- **Stellar drift**: play or scrub ±5 Myr along each star's measured space velocity, with motion trails. Choose a Sun-at-rest or galactic (LSR) frame. Drift can be turned off.
- **Gliese 710 flyby**: a preset that jumps to its pass about 0.05 pc from the Sun, roughly 1.29 Myr from now.
- **Spotlight**: shows the neighbours around the selected star within an adjustable radius (10 ly to about 3,200 ly), labelled with their distance from that star.
- **Star lanes**: a nearest-neighbour graph (k links per star, capped by a maximum jump length) inside a region around the Sun, with Dijkstra route-finding between stars. **Copy graph as JSON** exports the current graph.
- **HR diagram**: drag a box on the diagram to isolate those stars in 3D.
- **Your own Gaia data**: drop in a CSV exported from the [Gaia Archive](https://gea.esac.esa.int/archive/). The ADQL query to use is shown in the page.

## Running locally

The page loads its data with `fetch`, so it needs a local web server rather than opening the file directly:

```sh
python3 -m http.server 8000
# open http://localhost:8000
```

## Files

| Path | What it is |
|---|---|
| `index.html` | The whole viewer: markup, styles and script. three.js r128 is loaded from cdnjs/jsDelivr. |
| `stars.bin` | Star table: 13 float32 columns, column-major (`x y z mag absmag ci hip gaia spect vx vy vz hasrv`). |
| `stars.b64.txt` | The same bytes as base64. It's a fallback for hosts that can't serve binary files, such as the claude.ai artifact. |
| `meta.json` | Row count, column names, spectral-type lookup table and star names. |
| `tools/prep.py` | Builds `stars.bin`, `stars.b64.txt` and `meta.json` from the source catalogue. |
| `tools/lanes.py` | Exports a lane graph: `python3 tools/lanes.py [radius_pc] [max_jump_pc] [links_per_star]`. |
| `exports/star_lanes_25pc.json` | Example lane graph: 25 pc radius, 3.5 pc maximum jump, 3 links per star. |

### Coordinates and units

- Positions are in parsecs, heliocentric galactic: +x toward the galactic centre, +y toward galactic rotation (l = 90°), +z toward the north galactic pole.
- Velocities are stored as pc/Myr in `stars.bin` and as km/s in the JSON exports.
- `ci` is the B−V colour index. `gaia` is 1 when the distance comes from a Gaia parallax. `hasrv` is 1 when a radial velocity was measured; stars without one only move across the sky.

### Lane graph format

```json
{
  "meta": { "units": { "position": "parsec", "velocity": "km/s" }, "regionRadiusPc": 25, "maxJumpPc": 3.5, "linksPerStar": 3 },
  "stars": [{ "id": 0, "name": "Sol", "x": 0, "y": 0, "z": 0, "vx": 0, "vy": 0, "vz": 0, "absmag": 4.85, "tempK": 5760, "cluster": 0, "spect": "G2 V" }],
  "lanes": [[0, 1347, 1.302]]
}
```

`lanes` entries are `[star id, star id, length in pc]`. `cluster` is the connected-component ID.

## Rebuilding the data

1. Download `hyg/athyg_v3/hyglike_from_athyg_v32.csv.gz` from [astronexus/HYG-Database](https://github.com/astronexus/HYG-Database).
2. From the repo root, run:

```sh
pip install pandas numpy scipy
python3 tools/prep.py path/to/hyglike_from_athyg_v32.csv.gz
python3 tools/lanes.py 25 3.5 3
```

`prep.py` renames Rigil Kentaurus and Toliman to Alpha Centauri A and B, keeping the old names as aliases. It also fills in Gliese designations for stars that have no other name.

## Data source, limits and licence

The star data comes from the HYGLike subset of **AT-HYG v3.2** by David Nash (astronexus), which merges Tycho-2, Hipparcos and Gaia DR3. 98% of these stars take their distance from Gaia DR3; the rest, mostly very bright stars that Gaia saturates on, keep Hipparcos distances.

The catalogue is limited to stars Tycho-2 and Hipparcos could see, so most faint red dwarfs beyond about 20 pc are missing. Lanes far from the Sun are therefore longer than a complete map would give.

AT-HYG is licensed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/). The files derived from it (`stars.bin`, `stars.b64.txt`, `meta.json` and `exports/`) are shared under the same licence, with attribution to the HYG/AT-HYG project.
