# detector_geometry
#### A tool to project resolution cones at different detector geometries (tilt, rotation, offset) at given X-ray energies.
 - Main application is to visualize the maximum achievable resolution at a given geometry.
 - The math used is not meant to bring people to the moon but to provide a quick and simple preview.
 - It uses python3, numpy and matplotlib.

## Short how-to:
### Edit the following sections in the *.py* file:

 |     Geometry      |   |
 |-------------------|---|
 | geo_dist = 17.5   | [cm]   Detector distance
 | geo_tilt = 0.0    | [deg]  Detector tilt
 | geo_rota = 0.0    | [deg]  detector rotation
 | geo_yoff = 13.0   | [cm]   Detector offset (vertical)
 | geo_energy = 35.0 | [keV]  Beam energy
 | plt_unit = 'q'    | [tdqs] Contour legend (t: 2-Theta, d: d-spacing, q: q-space, s: sin(theta)/lambda)
 | plt_origin = True | [bool] plot contour lines for original geometry?

## To add a detector:
 - add/change the following lines.
 - comment/uncomment to change between

 |  Specifications   |   |
 |-------------------|---|
 | det_name = 'Pilatus3 X CdTe 2M' |
 | det_hms = 8.38     | [cm]  module size (horizontal)
 | det_vms = 3.35     | [cm]  module size (vertical)
 | det_pxs = 172e-4   | [cm]  pixel size
 | det_hgp = 7        | [pix] gap between modules (horizontal)
 | det_vgp = 17       | [pix] gap between modules (vertical)
 | det_hmn = 3        | [int] number of modules (horizontal)
 | det_vmn = 8        | [int] number of modules (vertical)

## Here's an example showing an Pilatus3 2M offset vertically:
![Image](../main/Pilatus3_X_CdTe_2M.png)

## And a rotated Eiger2 4M:
![Image](../main/Eiger2_CdTe_4M.png)
