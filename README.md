# detector_geometry
#### A tool to project resolution cones at different detector geometries (tilt, rotation, offset) and X-ray energies.
 - Main application is to visualize the maximum achievable resolution at a given geometry.
 - The math used is not meant to bring people to the moon but to provide a quick and simple preview.
 - The module building code is designed for [Dectris](https://www.dectris.com) [Pilatus3](https://www.dectris.com/detectors/x-ray-detectors/pilatus3/) and [Eiger2](https://www.dectris.com/detectors/x-ray-detectors/eiger2/) Detectors but one-module systems like the [Bruker](https://www.bruker.com/en.html) [Photon II](https://www.bruker.com/en/products-and-solutions/diffractometers-and-scattering-systems/single-crystal-x-ray-diffractometers/sc-xrd-components/detectors.html) are possible as well.
 - It uses [python3](https://www.python.org), [numpy](https://numpy.org) and [matplotlib 3.5.1](https://matplotlib.org).

## Short how-to:
 - Edit the geometry section in the *.py* file
 - run it
#### Or
 - run it (plo.interactive = True)
 - Use the (highly unsatisfying) text input to change geometry:
   - Detector [cm]
   - Tilt [˚]
   - Offset [cm]
   - Rotation [˚]
 - And the (slightly better) radio buttons (top left) to change the unit:
   - d: d-space
   - t: 2-theta
   - q: q-space
   - s: sin(theta)/lambda

## Here's an example showing a Pilatus3 2M offset vertically:

 |   Geometry   |        Value      | Hint |
 |--------------|-------------------|------|
 | geo.det_type | 'Pilatus3 X CdTe' | [str]  Pilatus3 / Eiger2
 | geo.det_size | '2M'              | [str]  300K 1M 2M 6M / 1M 4M 9M 16M
 | geo.dist     | 17.5              | [cm]   Detector distance
 | geo.tilt     | 0.0               | [deg]  Detector tilt
 | geo.rota     | 0.0               | [deg]  Detector rotation
 | geo.yoff     | 13.0              | [cm]   Detector offset (vertical)
 | geo.energy   | 35.0              | [keV]  Beam energy
 | geo.unit     | 'q'               | [tdqs] Contour legend (t: 2-Theta, d: d-spacing, q: q-space, s: sin(theta)/lambda)
 | geo.origin   | True              | [bool] Plot contour lines for original geometry?

![Image](../main/Pilatus3_X_CdTe_2M.png)

## And a rotated Eiger2 4M:

 |   Geometry   |        Value      | Hint |
 |--------------|-------------------|------|
 | geo.det_type | 'Eiger2 CdTe'     | [str]  Pilatus3 / Eiger2
 | geo.det_size | '4M'              | [str]  300K 1M 2M 6M / 1M 4M 9M 16M
 | geo.dist     | 7.5               | [cm]   Detector distance
 | geo.tilt     | 0.0               | [deg]  Detector tilt
 | geo.rota     | 20.0              | [deg]  Detector rotation
 | geo.yoff     | 0.0               | [cm]   Detector offset (vertical)
 | geo.energy   | 22.0              | [keV]  Beam energy
 | geo.unit     | 'd'               | [tdqs] Contour legend (t: 2-Theta, d: d-spacing, q: q-space, s: sin(theta)/lambda)
 | geo.origin   | True              | [bool] Plot contour lines for original geometry?
 
![Image](../main/Eiger2_CdTe_4M.png)

## To add a detector:
 - Choose 'Name' and 'Version' (will be used in the figure title).
 - 'Name' must not start with 'Pilatus' or 'Eiger', those are pre-set.
 - Adjust the "ADD CUSTOM DETECTOR SPECIFICATIONS HERE" section.

 |   Geometry   |   Value   | Hint |
 |--------------|-----------|------|
 | geo.det_type | 'Name'    | [str]  Pilatus3 / Eiger2
 | geo.det_size | 'Version' | [str]  300K 1M 2M 6M / 1M 4M 9M 16M
 | geo.dist     | 6.0       | [cm]   Detector distance
 | geo.tilt     | 0.0       | [deg]  Detector tilt
 | geo.rota     | 0.0       | [deg]  Detector rotation
 | geo.yoff     | 0.0       | [cm]   Detector offset (vertical)
 | geo.energy   | 17.0      | [keV]  Beam energy
 | geo.unit     | 'd'       | [tdqs] Contour legend (t: 2-Theta, d: d-spacing, q: q-space, s: sin(theta)/lambda)
 | geo.origin   | True      | [bool] Plot contour lines for original geometry?

 | Detector |       Value       | Hint |
 |----------|-------------------|------|
 | det.name | 'Name Version'    | [auto] Generated
 | det.hms  | 10.0              | [cm]   Module size (horizontal)
 | det.vms  | 14.0              | [cm]   Module size (vertical)
 | det.pxs  | 50e-4             | [cm]   Pixel size
 | det.hgp  | 0                 | [pix]  Gap between modules (horizontal)
 | det.vgp  | 0                 | [pix]  Gap between modules (vertical)
 | det.hmn  | 1                 | [int]  Number of modules (horizontal)
 | det.vmn  | 1                 | [int]  Number of modules (vertical)
 
