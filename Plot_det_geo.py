import numpy as np
from matplotlib.widgets import TextBox, RadioButtons
from matplotlib import pyplot as plt
plt.rcParams['savefig.dpi'] = 300
from matplotlib import cm, colors, patches

def get_specs():
    ######################
    # Setup the geometry #
    ######################
    geo = container()
    geo.det_type = 'Eiger2 CdTe' # [str]  Pilatus3 / Eiger2
    geo.det_size = '4M'          # [str]  300K 1M 2M 6M / 1M 4M 9M 16M
    geo.dist = 8.0               # [cm]   Detector distance
    geo.tilt = 0.0               # [deg]  Detector tilt
    geo.rota = 0.0               # [deg]  Detector rotation
    geo.yoff = 0.0               # [cm]   Detector offset (vertical)
    geo.energy = 21.0            # [keV]  Beam energy
    geo.unit = 'd'               # [tdqs] Contour legend (t: 2-Theta, d: d-spacing, q: q-space, s: sin(theta)/lambda)
    geo.origin = True            # [bool] Plot contour lines for original geometry?

    ###########################
    # Detector Specifications #
    ###########################
    det = container()
    if geo.det_type.startswith('Pilatus'):
        ###############################
        # Specifications for Pilatus3 #
        ###############################
        det.hms = 8.38    # [cm]  Module size (horizontal)
        det.vms = 3.35    # [cm]  Module size (vertical)
        det.pxs = 172e-4  # [cm]  Pixel size
        det.hgp = 7       # [pix] Gap between modules (horizontal)
        det.vgp = 17      # [pix] Gap between modules (vertical)
        det.name = f'{geo.det_type} {geo.det_size}'
        det.sizes = {'300K':(1,3),'1M':(2,5),'2M':(3,8),'6M':(5,12)}
        if geo.det_size not in det.sizes.keys():
            print('Unknown detector type/size combination!')
            raise SystemExit
        det.hmn, det.vmn = det.sizes[geo.det_size]
    elif geo.det_type.startswith('Eiger'):
        #############################
        # Specifications for Eiger2 #
        #############################
        det.hms = 7.71    # [cm]  Module size (horizontal)
        det.vms = 3.84    # [cm]  Module size (vertical)
        det.pxs = 75e-4   # [cm]  Pixel size
        det.hgp = 38      # [pix] Gap between modules (horizontal)
        det.vgp = 12      # [pix] Gap between modules (vertical)
        det.name = f'{geo.det_type} {geo.det_size}'
        det.sizes = {'1M':(1,2),'4M':(2,4),'9M':(3,6),'16M':(4,8)}
        if geo.det_size not in det.sizes.keys():
            print('Unknown detector type/size combination!')
            raise SystemExit
        det.hmn, det.vmn = det.sizes[geo.det_size]
    else:
        ###########################################
        # ADD CUSTOM DETECTOR SPECIFICATIONS HERE #
        ###########################################
        det.hms = 10.0    # [cm]  Module size (horizontal)
        det.vms = 14.0    # [cm]  Module size (vertical)
        det.pxs = 10e-4   # [cm]  Pixel size
        det.hgp = 0       # [pix] Gap between modules (horizontal)
        det.vgp = 0       # [pix] Gap between modules (vertical)
        det.hmn = 1       # [int] Number of modules (horizontal)
        det.vmn = 1       # [int] Number of modules (vertical)
        det.name = f'{geo.det_type} {geo.det_size}'

    ################
    # Plot Details #
    ################
    plo = container()
    plo.cont_levels = np.logspace(-1,1,num=15)/2  # [list]  Contour levels
    plo.cont_fsize = 8                            # [int]   Contour label size
    plo.cont_geom_alpha = 1.00                    # [float] Contour alpha (geometry)
    plo.cont_geom_cmap = cm.get_cmap('viridis_r') # [cmap]  Contour colormap (geometry)
    plo.cont_orig_alpha = 0.25                    # [float] Contour alpha (original)
    plo.cont_orig_color = 'gray'                  # [color] Contour color (original)
    plo.cont_reso = 500                           # [int]   Contour steps (accuracy)
    plo.cont_xmax = 50                            # [int]   Max x/y for drawing contours
    plo.module_alpha = 0.20                       # [float] Detector module alpha
    plo.module_color = 'gray'                     # [color] Detector module color
    plo.margin_top = 0.93                         # [float] Plot margin for title
    plo.plot_size = 7                             # [int]   Plot size
    plo.interactive = True                        # [bool]  Make the plot interactive
    plo.debug_3d = False                          # [bool]  DEBUG plot 3D cones?

    ###################################
    # !!! Don't change below here !!! #
    ###################################
    return geo, det, plo

def main():
    # fetch the geometry, detector and plot specifications
    geo, det, plo = get_specs()
    # translate unit for plot title
    geo.unit_names = {'t':r'2$\theta$',
                      'd':r'$d_{space}$',
                      'q':r'$q_{space}$',
                      's':r'$sin(\theta)/\lambda$'}
    if geo.unit not in geo.unit_names.keys():
        print('Unknown contour label unit!')
        raise SystemExit
    # figure out proper plot dimensions
    plo.xdim = (det.hms * det.hmn + det.pxs * det.hgp * det.hmn)
    plo.ydim = (det.vms * det.vmn + det.pxs * det.vgp * det.vmn)
    plo.fig_ratio = plo.xdim / plo.ydim
    # init the plot
    fig = plt.figure(figsize=(plo.plot_size * plo.margin_top * plo.fig_ratio, plo.plot_size))
    # needed to avoid the following warning:
    # MatplotlibDeprecationWarning: Toggling axes navigation from the keyboard is deprecated
    # since 3.3 and will be removed two minor releases later.
    #fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)
    # add axis for the detector modules
    bg = fig.add_subplot(111, aspect='equal')
    # add axis for the contours
    ax = fig.add_subplot(111, aspect='equal')
    # remove the axis
    bg.set_axis_off()
    ax.set_axis_off()
    # limit the axes
    bg.set_xlim(-plo.xdim/2, plo.xdim/2)
    ax.set_xlim(-plo.xdim/2, plo.xdim/2)
    bg.set_ylim(-plo.ydim/2, plo.ydim/2)
    ax.set_ylim(-plo.ydim/2, plo.ydim/2)
    # setup detector and geometry
    build_detector(bg, det, plo)
    # create cones and draw contour lines
    draw_contours(ax, geo, plo)
    # add title / information
    plt.suptitle(f'{det.name} | Energy: {geo.energy} keV | Distance: {geo.dist} cm\nRotation: {geo.rota}° | Tilt: {geo.tilt}° | Offset: {geo.yoff} cm | Units: {geo.unit_names[geo.unit]}', size=10)
    # adjust the margins
    plt.subplots_adjust(top=plo.margin_top, bottom=0, right=1, left=0, hspace=0, wspace=0)
    # generate some sense of interactivity
    if plo.interactive:
        # cut the title to make room for the text boxes
        plt.suptitle(f'{det.name}', size=10)
        # add an axis for the radio buttons
        axs_unit = fig.add_axes([0.0, 0.92, 0.09, 0.09], frameon=False)
        # put the buttons in the axis
        box_unit = RadioButtons(axs_unit, ('t','d','q','s'), active=1, activecolor=u'#1f77b4')
        # update the plot on click
        box_unit.on_clicked(lambda val: update_plot('unit', val, fig, geo, plo, det, ax))
        # make a new axis for the TextBox
        axs_ener = fig.add_axes([0.20, 0.93, 0.06, 0.03], frameon=False)
        # add the box with label and initial value
        box_ener = TextBox(axs_ener, 'Energy:', initial=geo.energy)
        # update the plot on submit (press enter or leaving the box)
        box_ener.on_submit(lambda val: update_plot('energy', val, fig, geo, plo, det, ax))
        # make the TextBox respond to changing text
        box_ener.on_text_change(lambda val: update_box(box_tilt, val))
        # and distance
        axs_dist = fig.add_axes([0.40, 0.93, 0.06, 0.03], frameon=False)
        box_dist = TextBox(axs_dist, 'Distance:', initial=geo.dist)
        box_dist.on_submit(lambda val: update_plot('dist', val, fig, geo, plo, det, ax))
        box_dist.on_text_change(lambda val: update_box(box_dist, val))
        # proceed with rotation
        axs_rota = fig.add_axes([0.60, 0.93, 0.06, 0.03], frameon=False)
        box_rota = TextBox(axs_rota, 'Rotation:', initial=geo.rota)
        box_rota.on_submit(lambda val: update_plot('rota', val, fig, geo, plo, det, ax))
        box_rota.on_text_change(lambda val: update_box(box_rota, val))
        # and offset
        axs_yoff = fig.add_axes([0.75, 0.93, 0.06, 0.03], frameon=False)
        box_yoff = TextBox(axs_yoff, 'Offset:', initial=geo.yoff)
        box_yoff.on_submit(lambda val: update_plot('yoff', val, fig, geo, plo, det, ax))
        box_yoff.on_text_change(lambda val: update_box(box_yoff, val))
        # and tilt
        axs_tilt = fig.add_axes([0.90, 0.93, 0.06, 0.03], frameon=False)
        box_tilt = TextBox(axs_tilt, 'Tilt:', initial=geo.tilt)
        box_tilt.on_submit(lambda val: update_plot('tilt', val, fig, geo, plo, det, ax))
        box_tilt.on_text_change(lambda val: update_box(box_tilt, val))

    # show the plot
    plt.show()
    # plot the 3d cones?
    if plo.debug_3d:
        plot_3d(geo, plo)

def build_detector(bg, det, plo):
    # build detector modules
    # beam position is between the modules (even) or at the center module (odd)
    # determined by the "+det.hmn%2" part
    for i in range(-det.hmn//2+det.hmn%2, det.hmn-det.hmn//2):
        for j in range(-det.vmn//2+det.vmn%2, det.vmn-det.vmn//2):
            # place modules along x (i) and y (j) keeping the gaps in mind ( + (det.hgp*det.pxs)/2)
            # the " - ((det.hms+det.hgp*det.pxs)/2)" positions the origin (the beam) at the center of a module
            # and "det.hmn%2" makes sure this is only active for detectors with an odd number of modules
            origin_x = i*(det.hms+det.hgp*det.pxs) - ((det.hms+det.hgp*det.pxs)/2)*(det.hmn%2) + (det.hgp*det.pxs)/2
            origin_y = j*(det.vms+det.vgp*det.pxs) - ((det.vms+det.vgp*det.pxs)/2)*(det.vmn%2) + (det.vgp*det.pxs)/2
            bg.add_patch(patches.Rectangle((origin_x, origin_y),  det.hms, det.vms, color=plo.module_color, alpha=plo.module_alpha))

def draw_contours(ax, geo, plo):
    # draw contour lines
    for n,i in enumerate(plo.cont_levels):
        # calculate resolution rings
        # 2-Theta: np.arctan(dist/(dist/i))
        thr = np.arctan(1/i)/2
        # Conversion factor keV to Angstrom: 12.398
        # sin(t)/l: np.sin(Theta) / lambda -> (12.398/geo_energy)
        stl = np.sin(thr)/(12.398/geo.energy)
        # d-spacing: l = 2 d sin(t) -> 1/2(sin(t)/l)
        dsp = 1/(2*stl)
        # figure out the labels
        unit_dict = {'n':None, 't':np.rad2deg(2*thr),
                     'd':dsp, 'q':stl*4*np.pi, 's':stl}
        # draw contours for the original geometry
        if geo.origin:
            X,Y,Z = create_cone(i, 0, 0, 0, geo.dist, plo.cont_xmax, plo.cont_reso)
            # don't draw contour lines that are out of bounds
            # make sure Z is large enough to draw the contour
            if np.max(Z) >= geo.dist:
                c0 = ax.contour(X, Y, Z, [geo.dist], colors=plo.cont_orig_color, alpha=plo.cont_orig_alpha)
                # label original geometry contours
                fmt = {c0.levels[0]:f'{np.round(unit_dict[geo.unit],2):.2f}'}
                ax.clabel(c0, c0.levels, inline=True, fontsize=plo.cont_fsize, fmt=fmt, manual=[(plo.xdim,plo.ydim)])
        # draw contours for the tilted/rotated/moved geometry
        X,Y,Z = create_cone(i, geo.rota, geo.tilt, geo.yoff, geo.dist, plo.cont_xmax, plo.cont_reso)
        # make sure Z is large enough to draw the contour
        if np.max(Z) >= geo.dist:
            c1 = ax.contour(X, Y, Z, [geo.dist], colors=colors.to_hex(plo.cont_geom_cmap((n+1)/len(plo.cont_levels))), alpha=plo.cont_geom_alpha)
            # label moved geometry contours
            fmt = {c1.levels[0]:f'{np.round(unit_dict[geo.unit],2):.2f}'}
            ax.clabel(c1, c1.levels, inline=True, fontsize=plo.cont_fsize, fmt=fmt, manual=[(0,plo.ydim)])

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

def update_plot(nam, val, fig, geo, plo, det, ax):
    ##################################################
    # This is a sloppy and hacky way to achieve some #
    #   interactivity without building a proper GUI  #
    ##################################################
    if nam == 'dist':
        geo.dist = float(val)
    elif nam == 'rota':
        geo.rota = float(val)
    elif nam == 'tilt':
        geo.tilt = float(val)
    elif nam == 'yoff':
        geo.yoff = float(val)
    elif nam == 'unit':
        geo.unit = str(val)
    elif nam == 'energy':
        geo.energy = float(val)
    # remove contour lines
    for _ in ax.collections:
        ax.collections.pop()
    # remove contour labels
    for _ in ax.texts:
        ax.texts.pop()
    # re-calculate cones and re-draw contours
    draw_contours(ax, geo, plo)
    plt.suptitle(f'{det.name}', size=10)
    fig.canvas.blit(ax)
    fig.canvas.flush_events()

def update_box(box, val):
    box.canvas.blit()

def plot_3d(geo, plo):
    #####################################################
    # - debug - debug - debug - debug - debug - debug - #
    # - to check geometry, offset, tilt and rotation  - #
    #####################################################
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for n,i in enumerate(plo.cont_levels[-3:]):
        ax.plot_wireframe(*create_cone(i, 0, 0, 0, geo.dist, plo.cont_xmax, plo.cont_reso), alpha=0.1, colors='gray')
        ax.contour(*create_cone(i, 0, 0, 0, geo.dist, plo.cont_xmax, plo.cont_reso), [geo.dist], alpha=0.1, colors='gray')
        ax.plot_wireframe(*create_cone(i, geo.rota, geo.tilt, geo.yoff, geo.dist, plo.cont_xmax, plo.cont_reso), alpha=0.25, colors='red')
        ax.contour(*create_cone(i, geo.rota, geo.tilt, geo.yoff, geo.dist, plo.cont_xmax, plo.cont_reso), [geo.dist], colors='red')
    ax.set_xlim(-plo.cont_xmax/2, plo.cont_xmax/2)
    ax.set_ylim(-plo.cont_xmax/2, plo.cont_xmax/2)
    ax.set_zlim(0, geo.dist)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

class container(object):
    pass

if __name__ == '__main__':
    main()
