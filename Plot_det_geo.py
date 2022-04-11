import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm, colors, patches

# Geometry
geo_det_family = 'Pilatus3 X CdTe' # [--] Pilatus3 / Eiger2
geo_det_size = '2M'         # [--] 300K 1M 2M 6M / 1M 4M 9M 16M
geo_dist = 17.6             # [cm]   Detector distance
geo_tilt = 0.0              # [deg]  Detector tilt
geo_rota = 0.0              # [deg]  detector rotation
geo_yoff = 13.0             # [cm]   Detector offset (vertical)
geo_energy = 35.0           # [keV]  Beam energy
plt_unit = 'q'              # [tdqs] Contour legend (t: 2-Theta, d: d-spacing, q: q-space, s: sin(theta)/lambda)
plt_origin = True           # [bool] plot contour lines for original geometry?

# Detector Specifications
if geo_det_family.startswith('Pilatus'):
    det_name = f'{geo_det_family} {geo_det_size}'
    det_hms = 8.38    # [cm]  module size (horizontal)
    det_vms = 3.35    # [cm]  module size (vertical)
    det_pxs = 172e-4  # [cm]  pixel size
    det_hgp = 7       # [pix] gap between modules (horizontal)
    det_vgp = 17      # [pix] gap between modules (vertical)
    det_sizes = {'300K':(1,3),'1M':(2,5),'2M':(3,8),'6M':(5,12)}
    det_hmn, det_vmn = det_sizes[geo_det_size]
elif geo_det_family.startswith('Eiger'):
    det_name = f'{geo_det_family} {geo_det_size}'
    det_hms = 7.71    # [cm]  module size (horizontal)
    det_vms = 3.84    # [cm]  module size (vertical)
    det_pxs = 75e-4   # [cm]  pixel size
    det_hgp = 38      # [pix] gap between modules (horizontal)
    det_vgp = 12      # [pix] gap between modules (vertical)
    det_sizes = {'1M':(1,2),'4M':(2,4),'9M':(3,6),'16M':(4,8)}
    det_hmn, det_vmn = det_sizes[geo_det_size]
else:
    # Add custom detector specs here
    det_name = f'{geo_det_family} {geo_det_size}'
    det_hms = 10.0    # [cm]  module size (horizontal)
    det_vms = 14.0    # [cm]  module size (vertical)
    det_pxs = 10e-4   # [cm]  pixel size
    det_hgp = 0       # [pix] gap between modules (horizontal)
    det_vgp = 0       # [pix] gap between modules (vertical)
    det_hmn = 1       # [int] number of modules (horizontal)
    det_vmn = 1       # [int] number of modules (vertical)

# Plot Details
plt_lines = np.logspace(-2, 1, num=25, base=10)/2
plt_cmap = cm.get_cmap('viridis_r')
plt_contour_fsize = 8
plt_module_alpha = 0.25
plt_cont_geom_alpha = 1.00
plt_cont_orig_alpha = 0.10
plt_top_margin = 0.93
plt_plot_size = 7

# internals
fix_conv = 12.398 # Conversion keV to Angstrom
fix_reso = 500    # contour/cone steps
fix_xmax = 50     # max x/y for drawing contours/cones

# debug
plot_3d = False   # [bool] DEBUG plot 3D cones?

def build_detector():
    # build detector modules
    for i in range(-det_hmn//2+det_hmn%2, det_hmn-det_hmn//2):
        for j in range(-det_vmn//2+det_vmn%2, det_vmn-det_vmn//2):
            origin_x = i*(det_hms+det_hgp*det_pxs) - ((det_hms+det_hgp*det_pxs)/2)*(det_hmn%2) + (det_hgp*det_pxs)/2
            origin_y = j*(det_vms+det_vgp*det_pxs) - ((det_vms+det_vgp*det_pxs)/2)*(det_vmn%2) + (det_vgp*det_pxs)/2
            ax.add_patch(patches.Rectangle((origin_x, origin_y),  det_hms, det_vms, color='black', alpha=plt_module_alpha))
    # limit axes
    xdim = (det_hms * det_hmn + det_pxs * det_hgp * det_hmn) / 2
    ydim = (det_vms * det_vmn + det_pxs * det_vgp * det_vmn) / 2
    ax.set_xlim(-xdim, xdim)
    ax.set_ylim(-ydim, ydim)
    # draw contour lines
    for n,i in enumerate(plt_lines):
        # calculate resolution rings
        # 2-Theta: np.arctan(dist/(dist/i))
        thr = np.arctan(1/i)/2
        # sin(t)/l: np.sin(Theta) / lambda -> (fix_conv/geo_energy)
        stl = np.sin(thr)/(fix_conv/geo_energy)
        # d-spacing: l = 2 d sin(t) -> 1/2(sin(t)/l)
        dsp = 1/(2*stl)
        # figure out the labels
        unit_dict = {'n':None, 't':np.rad2deg(2*thr),
                     'd':dsp, 'q':stl*4*np.pi, 's':stl}
        # draw contours for the original geometry
        if plt_origin:
            X,Y,Z = create_cone(i, 0, 0, 0)
            # don't draw contour lines that are out of bounds
            # make sure Z is large enough to draw the contour
            if np.max(Z) >= geo_dist:
                c0 = ax.contour(X, Y, Z, [geo_dist], colors='black', alpha=plt_cont_orig_alpha)
                # label original geometry contours
                fmt = {c0.levels[0]:f'{np.round(unit_dict[plt_unit],2):.2f}'}
                ax.clabel(c0, c0.levels, inline=True, fontsize=plt_contour_fsize, fmt=fmt, manual=[(xdim,ydim)])
        # draw contours for the tilted/rotated/moved geometry
        X,Y,Z = create_cone(i, geo_rota, geo_tilt, geo_yoff)
        # make sure Z is large enough to draw the contour
        if np.max(Z) >= geo_dist:
            c1 = ax.contour(X, Y, Z, [geo_dist], colors=colors.to_hex(plt_cmap((n+1)/len(plt_lines))), alpha=plt_cont_geom_alpha)
            # label moved geometry contours
            fmt = {c1.levels[0]:f'{np.round(unit_dict[plt_unit],2):.2f}'}
            ax.clabel(c1, c1.levels, inline=True, fontsize=plt_contour_fsize, fmt=fmt, manual=[(0,ydim)])

def create_cone(dim, rota, tilt, yoff):
    # creating grid
    x = np.linspace(-fix_xmax,fix_xmax,fix_reso)
    y = np.linspace(-fix_xmax,fix_xmax,fix_reso)
    X,Y = np.meshgrid(x,y)
    # set z values
    Z = np.sqrt(X**2+Y**2)*dim
    # rotate the sample around y
    a = np.deg2rad(tilt) + np.deg2rad(rota)
    t = np.transpose(np.array([X,Y,Z]), (1,2,0))
    m = [[np.cos(a), 0, np.sin(a)],[0,1,0],[-np.sin(a), 0, np.cos(a)]]
    X,Y,Z = np.transpose(np.dot(t, m), (2,0,1))
    # compensate for tilt
    comp = np.deg2rad(tilt) * geo_dist
    return Y,X+comp-yoff,Z

# translate unit for plot title
unit_names = {'t':'2-Theta',
              'd':'d-spacing',
              'q':'q-space',
              's':'sin(Theta)/lambda'}
if plt_unit not in unit_names.keys():
    raise SystemExit
# figure out proper plot dimensions
xdim = (det_hms * det_hmn + det_pxs * det_hgp * det_hmn)
ydim = (det_vms * det_vmn + det_pxs * det_vgp * det_vmn)
plt_fig_ratio = xdim / ydim
# init the plot
fig = plt.figure(figsize=(plt_plot_size*plt_top_margin*plt_fig_ratio,plt_plot_size))
ax = fig.add_subplot(111)
# setup detector and geometry
build_detector()
# add title / information
plt.suptitle(f'{det_name} | Energy: {geo_energy} keV | Distance: {geo_dist} cm\nRotation: {geo_rota}° | Tilt: {geo_tilt}° | Offset: {geo_yoff} cm | Units: {unit_names[plt_unit]}', size=10)
# adjust the margins
plt.subplots_adjust(top=plt_top_margin, bottom=0, right=1, left=0, hspace=0, wspace=0)
ax.set_aspect('equal')
plt.axis('off')
plt.show()

#####################################################
# - debug - debug - debug - debug - debug - debug - #
# - to check geometry, offset, tilt and rotation  - #
#####################################################
if plot_3d:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #plt_lines = [0.01, 0.10, 0.20, 0.30, 0.39, 0.42, 0.45, 0.49, 0.53, 0.58, 0.65, 0.75, 0.87, 1.00, 1.16, 1.50, 2.00, 3.00, 5.00, 10.00]
    for n,i in enumerate(plt_lines[-3:]):
        #ax.plot_surface(*create_cone(i, geo_tilt, geo_yoff), alpha=0.25)
        ax.plot_wireframe(*create_cone(i, 0, 0, 0), alpha=0.1, color='gray')
        ax.contour(*create_cone(i, 0, 0, 0), [geo_dist], alpha=0.1, color='gray')

        ax.plot_wireframe(*create_cone(i, geo_rota, geo_tilt, geo_yoff), alpha=0.25, color='red')
        ax.contour(*create_cone(i, geo_rota, geo_tilt, geo_yoff), [geo_dist], color='red')

    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    ax.set_zlim(-15, 20)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()