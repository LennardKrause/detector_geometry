import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm, colors, patches

# Geometry
geo_det_type = 'Pilatus3 X CdTe' # [str] Pilatus3 / Eiger2
geo_det_size = '1M' # [str]  300K 1M 2M 6M / 1M 4M 9M 16M
geo_dist = 13.0     # [cm]   Detector distance
geo_tilt = 0.0      # [deg]  Detector tilt
geo_rota = 20.0     # [deg]  detector rotation
geo_yoff = 0.0      # [cm]   Detector offset (vertical)
geo_energy = 50.0   # [keV]  Beam energy
geo_unit = 'd'      # [tdqs] Contour legend (t: 2-Theta, d: d-spacing, q: q-space, s: sin(theta)/lambda)
geo_origin = True   # [bool] plot contour lines for original geometry?

# Detector Specifications
if geo_det_type.startswith('Pilatus'):
    det_name = f'{geo_det_type} {geo_det_size}'
    det_hms = 8.38    # [cm]  module size (horizontal)
    det_vms = 3.35    # [cm]  module size (vertical)
    det_pxs = 172e-4  # [cm]  pixel size
    det_hgp = 7       # [pix] gap between modules (horizontal)
    det_vgp = 17      # [pix] gap between modules (vertical)
    det_sizes = {'300K':(1,3),'1M':(2,5),'2M':(3,8),'6M':(5,12)}
    if geo_det_size not in det_sizes.keys():
        print('Unknown detector type/size combination!')
        raise SystemExit
    det_hmn, det_vmn = det_sizes[geo_det_size]
elif geo_det_type.startswith('Eiger'):
    det_name = f'{geo_det_type} {geo_det_size}'
    det_hms = 7.71    # [cm]  module size (horizontal)
    det_vms = 3.84    # [cm]  module size (vertical)
    det_pxs = 75e-4   # [cm]  pixel size
    det_hgp = 38      # [pix] gap between modules (horizontal)
    det_vgp = 12      # [pix] gap between modules (vertical)
    det_sizes = {'1M':(1,2),'4M':(2,4),'9M':(3,6),'16M':(4,8)}
    if geo_det_size not in det_sizes.keys():
        print('Unknown detector type/size combination!')
        raise SystemExit
    det_hmn, det_vmn = det_sizes[geo_det_size]
else:
    # Add custom detector specs here
    det_name = f'{geo_det_type} {geo_det_size}'
    det_hms = 10.0    # [cm]  module size (horizontal)
    det_vms = 14.0    # [cm]  module size (vertical)
    det_pxs = 10e-4   # [cm]  pixel size
    det_hgp = 0       # [pix] gap between modules (horizontal)
    det_vgp = 0       # [pix] gap between modules (vertical)
    det_hmn = 1       # [int] number of modules (horizontal)
    det_vmn = 1       # [int] number of modules (vertical)

# Plot Details
plt_cont_levels = np.logspace(-1,1,num=25)/2  # contour levels
plt_cont_fsize = 8                            # contour label size
plt_cont_geom_alpha = 1.00                    # contour alpha (geometry)
plt_cont_geom_cmap = cm.get_cmap('viridis_r') # contour colormap (geometry)
plt_cont_orig_alpha = 0.10                    # contour alpha (original)
plt_cont_orig_color = 'black'                 # contour color (original)
plt_cont_reso = 500                           # contour steps (accuracy)
plt_cont_xmax = 50                            # max x/y for drawing contours
plt_module_alpha = 0.20                       # detector module alpha
plt_module_color = 'black'                    # detector module color
plt_margin_top = 0.93                         # plot margin for title
plt_plot_size = 7                             # plot size

# Debug
debug_plt_3d = False # [bool] DEBUG plot 3D cones?

def build_detector():
    # build detector modules
    for i in range(-det_hmn//2+det_hmn%2, det_hmn-det_hmn//2):
        for j in range(-det_vmn//2+det_vmn%2, det_vmn-det_vmn//2):
            origin_x = i*(det_hms+det_hgp*det_pxs) - ((det_hms+det_hgp*det_pxs)/2)*(det_hmn%2) + (det_hgp*det_pxs)/2
            origin_y = j*(det_vms+det_vgp*det_pxs) - ((det_vms+det_vgp*det_pxs)/2)*(det_vmn%2) + (det_vgp*det_pxs)/2
            ax.add_patch(patches.Rectangle((origin_x, origin_y),  det_hms, det_vms, color=plt_module_color, alpha=plt_module_alpha))
    # draw contour lines
    for n,i in enumerate(plt_cont_levels):
        # calculate resolution rings
        # 2-Theta: np.arctan(dist/(dist/i))
        thr = np.arctan(1/i)/2
        # Conversion factor keV to Angstrom: 12.398
        # sin(t)/l: np.sin(Theta) / lambda -> (12.398/geo_energy)
        stl = np.sin(thr)/(12.398/geo_energy)
        # d-spacing: l = 2 d sin(t) -> 1/2(sin(t)/l)
        dsp = 1/(2*stl)
        # figure out the labels
        unit_dict = {'n':None, 't':np.rad2deg(2*thr),
                     'd':dsp, 'q':stl*4*np.pi, 's':stl}
        # draw contours for the original geometry
        if geo_origin:
            X,Y,Z = create_cone(i, 0, 0, 0, geo_dist, plt_cont_xmax, plt_cont_reso)
            # don't draw contour lines that are out of bounds
            # make sure Z is large enough to draw the contour
            if np.max(Z) >= geo_dist:
                c0 = ax.contour(X, Y, Z, [geo_dist], colors=plt_cont_orig_color, alpha=plt_cont_orig_alpha)
                # label original geometry contours
                fmt = {c0.levels[0]:f'{np.round(unit_dict[geo_unit],2):.2f}'}
                ax.clabel(c0, c0.levels, inline=True, fontsize=plt_cont_fsize, fmt=fmt, manual=[(xdim,ydim)])
        # draw contours for the tilted/rotated/moved geometry
        X,Y,Z = create_cone(i, geo_rota, geo_tilt, geo_yoff, geo_dist, plt_cont_xmax, plt_cont_reso)
        # make sure Z is large enough to draw the contour
        if np.max(Z) >= geo_dist:
            c1 = ax.contour(X, Y, Z, [geo_dist], colors=colors.to_hex(plt_cont_geom_cmap((n+1)/len(plt_cont_levels))), alpha=plt_cont_geom_alpha)
            # label moved geometry contours
            fmt = {c1.levels[0]:f'{np.round(unit_dict[geo_unit],2):.2f}'}
            ax.clabel(c1, c1.levels, inline=True, fontsize=plt_cont_fsize, fmt=fmt, manual=[(0,ydim)])

def create_cone(dim, rota, tilt, yoff, geo_dist, plt_cont_xmax, plt_cont_reso):
    # creating grid
    x = np.linspace(-plt_cont_xmax,plt_cont_xmax,plt_cont_reso)
    X,Y = np.meshgrid(x,x)
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
if geo_unit not in unit_names.keys():
    raise SystemExit

# figure out proper plot dimensions
xdim = (det_hms * det_hmn + det_pxs * det_hgp * det_hmn)
ydim = (det_vms * det_vmn + det_pxs * det_vgp * det_vmn)
plt_fig_ratio = xdim / ydim
# init the plot
fig = plt.figure(figsize=(plt_plot_size*plt_margin_top*plt_fig_ratio,plt_plot_size))
ax = fig.add_subplot(111)
# limit the axes
ax.set_xlim(-xdim/2, xdim/2)
ax.set_ylim(-ydim/2, ydim/2)
# setup detector and geometry
build_detector()
# add title / information
plt.suptitle(f'{det_name} | Energy: {geo_energy} keV | Distance: {geo_dist} cm\nRotation: {geo_rota}° | Tilt: {geo_tilt}° | Offset: {geo_yoff} cm | Units: {unit_names[geo_unit]}', size=10)
# adjust the margins
plt.subplots_adjust(top=plt_margin_top, bottom=0, right=1, left=0, hspace=0, wspace=0)
ax.set_aspect('equal')
plt.axis('off')
plt.show()

#####################################################
# - debug - debug - debug - debug - debug - debug - #
# - to check geometry, offset, tilt and rotation  - #
#####################################################
def plot_3d():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #plt_cont_levels = [0.01, 0.10, 0.20, 0.30, 0.39, 0.42, 0.45, 0.49, 0.53, 0.58, 0.65, 0.75, 0.87, 1.00, 1.16, 1.50, 2.00, 3.00, 5.00, 10.00]
    for n,i in enumerate(plt_cont_levels[-3:]):
        #ax.plot_surface(*create_cone(i, geo_tilt, geo_yoff), alpha=0.25)
        ax.plot_wireframe(*create_cone(i, 0, 0, 0, geo_dist, plt_cont_xmax, plt_cont_reso), alpha=0.1, colors='gray')
        ax.contour(*create_cone(i, 0, 0, 0, geo_dist, plt_cont_xmax, plt_cont_reso), [geo_dist], alpha=0.1, colors='gray')

        ax.plot_wireframe(*create_cone(i, geo_rota, geo_tilt, geo_yoff, geo_dist, plt_cont_xmax, plt_cont_reso), alpha=0.25, colors='red')
        ax.contour(*create_cone(i, geo_rota, geo_tilt, geo_yoff, geo_dist, plt_cont_xmax, plt_cont_reso), [geo_dist], colors='red')

    ax.set_xlim(-plt_cont_xmax/2, plt_cont_xmax/2)
    ax.set_ylim(-plt_cont_xmax/2, plt_cont_xmax/2)
    ax.set_zlim(0, geo_dist)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

if debug_plt_3d:
    plot_3d()
