import numpy as np
from matplotlib.widgets import TextBox, RadioButtons, Slider
from matplotlib import pyplot as plt
plt.rcParams['savefig.dpi'] = 300
from matplotlib import cm, colors, patches

def get_specs():
    ######################
    # Setup the geometry #
    ######################
    geo = container()
    geo.det_type = 'Eiger2 CdTe'   # [str]  Pilatus3 / Eiger2
    geo.det_size = '4M'            # [str]  300K 1M 2M 6M / 1M 4M 9M 16M
    geo.ener = 21.0                # [keV]  Beam energy
    geo.dist = 75.0                # [mm]   Detector distance
    geo.yoff = 0.0                 # [mm]   Detector offset (vertical)
    geo.rota = 32.0                # [deg]  Detector rotation
    geo.tilt = 0.0                 # [deg]  Detector tilt
    geo.unit = 'd'                 # [tdqs] Contour legend (t: 2-Theta, d: d-spacing, q: q-space, s: sin(theta)/lambda)

    ###########################
    # Detector Specifications #
    ###########################
    det = container()
    if geo.det_type.startswith('Pilatus'):
        ###############################
        # Specifications for Pilatus3 #
        ###############################
        det.hms = 83.8    # [mm]  Module size (horizontal)
        det.vms = 33.5    # [mm]  Module size (vertical)
        det.pxs = 172e-3  # [mm]  Pixel size
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
        det.hms = 77.1    # [mm]  Module size (horizontal)
        det.vms = 38.4    # [mm]  Module size (vertical)
        det.pxs = 75e-3   # [mm]  Pixel size
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
        det.hms = 100.0   # [mm]  Module size (horizontal)
        det.vms = 140.0   # [mm]  Module size (vertical)
        det.pxs = 10e-3   # [mm]  Pixel size
        det.hgp = 0       # [pix] Gap between modules (horizontal)
        det.vgp = 0       # [pix] Gap between modules (vertical)
        det.hmn = 1       # [int] Number of modules (horizontal)
        det.vmn = 1       # [int] Number of modules (vertical)
        det.name = f'{geo.det_type} {geo.det_size}'

    ################
    # Plot Details #
    ################
    plo = container()
    plo.cont_levels = np.logspace(1,-1,num=15)/2  # [list]   Contour levels (high -> low)
    plo.cont_fsize = 9                            # [int]    Contour label size
    plo.cont_geom_cmark = 'o'                     # [marker] Beam center marker (geometry)
    plo.cont_geom_csize = 4                       # [int]    Beam center size (geometry)
    plo.cont_geom_alpha = 1.00                    # [float]  Contour alpha (geometry)
    plo.cont_geom_cmap = cm.get_cmap('viridis')   # [cmap]   Contour colormap (geometry)
    plo.origin = True                             # [bool]   Plot contour lines for original geometry?
    plo.cont_orig_cmark = 'o'                     # [marker] Beam center marker (original)
    plo.cont_orig_csize = 4                       # [int]    Beam center size (original)
    plo.cont_orig_alpha = 0.25                    # [float]  Contour alpha (original)
    plo.cont_orig_color = 'gray'                  # [color]  Contour color (original)
    plo.cont_reso = 100                           # [int]    Minimum contour steps
    plo.module_alpha = 0.20                       # [float]  Detector module alpha
    plo.module_color = 'gray'                     # [color]  Detector module color
    plo.margin_top = 0.93                         # [float]  Plot margin for title
    plo.plot_size = 8                             # [int]    Plot size
    plo.interactive = True                        # [bool]   Make the plot interactive

    ##########
    # Limits #
    ##########
    lmt = container()
    lmt.ener_min = 1.0   # [float] Energy minimum [keV]
    lmt.ener_max = 100.0 # [float] Energy maximum [keV]
    lmt.ener_stp = 0.1   # [float] Energy step size [keV]
    lmt.dist_min = 40.0  # [float] Distance minimum [mm]
    lmt.dist_max = 450.0 # [float] Distance maximum [mm]
    lmt.dist_stp = 1.0   # [float] Distance step size [mm]
    lmt.yoff_min = -0.0  # [float] Offset minimum [mm]
    lmt.yoff_max = 200.0 # [float] Offset maximum [mm]
    lmt.yoff_stp = 1.0   # [float] Offset step size [mm]
    lmt.rota_min = 0.0   # [float] Rotation minimum [deg]
    lmt.rota_max = 45.0  # [float] Rotation maximum [deg]
    lmt.rota_stp = 1.0   # [float] Rotation step size [deg]
    lmt.tilt_min = 0.0   # [float] Tilt minimum [deg]
    lmt.tilt_max = 45.0  # [float] Tilt maximum [deg]
    lmt.tilt_stp = 1.0   # [float] Tilt step size [deg]

    ###################################
    # !!! Don't change below here !!! #
    ###################################
    return geo, det, plo, lmt

def main():
    # fetch the geometry, detector and plot specifications
    geo, det, plo, lmt = get_specs()
    # translate unit for plot title
    geo.unit_names = {'t':r'2$\theta$',
                      'd':r'$d_{space}$',
                      'q':r'$q_{space}$',
                      's':r'$sin(\theta)/\lambda$'}
    if geo.unit not in geo.unit_names.keys():
        print('Unknown contour label unit!')
        raise SystemExit
    # the reversed (v:k) dict is needed for the radio buttons
    geo.unit_names_r = dict(map(reversed, geo.unit_names.items()))
    # figure out proper plot dimensions
    plo.xdim = (det.hms * det.hmn + det.pxs * det.hgp * det.hmn)/2
    plo.ydim = (det.vms * det.vmn + det.pxs * det.vgp * det.vmn)/2
    plo.fig_ratio = plo.xdim / plo.ydim
    # scale contour grid to detector size
    plo.cont_grid_max = int(np.ceil(max(plo.xdim, plo.ydim)))
    # make room for the interactive sliders
    if plo.interactive:
        plo.margin_top -= 0.08
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
    bg.set_xlim(-plo.xdim, plo.xdim)
    ax.set_xlim(-plo.xdim, plo.xdim)
    bg.set_ylim(-plo.ydim, plo.ydim)
    ax.set_ylim(-plo.ydim, plo.ydim)
    # setup detector and geometry
    build_detector(bg, det, plo)
    # create cones and draw contour lines
    draw_contours(ax, geo, plo)
    # add title / information
    plt.suptitle(f'{det.name} | Energy: {geo.ener} keV | Distance: {geo.dist} cm\nRotation: {geo.rota}° | Tilt: {geo.tilt}° | Offset: {geo.yoff} cm | Units: {geo.unit_names[geo.unit]}', size=10)
    # adjust the margins
    plt.subplots_adjust(top=plo.margin_top, bottom=0, right=1, left=0, hspace=0, wspace=0)
    # generate some sense of interactivity
    if plo.interactive:
        # cut the title to make room for the text boxes
        plt.suptitle(f'{det.name}', size=10, fontweight='bold')
        # add radio buttons and an axis for the buttons
        axs_unit = fig.add_axes([0.0, 0.85, 0.15, 0.15], frameon=False, aspect='equal')
        box_unit = RadioButtons(axs_unit, (geo.unit_names.values()), active=1, activecolor=u'#1f77b4')
        box_unit.on_clicked(lambda val: update_plot('unit', geo.unit_names_r[val], fig, geo, plo, det, ax))
        # change label size
        for l in box_unit.labels:
            l.set_size(plo.cont_fsize)
        # Define sliders
        #add_slider(label, name, left, bottom, width, height, val, vmin, vmax, step, fig, ax, geo, det, plo)
        sli_tilt = add_slider('Tilt [˚] '     , 'tilt', 0.3, 0.85, 0.6, 0.025, geo.tilt, lmt.tilt_min, lmt.tilt_max, lmt.tilt_stp, fig, ax, geo, det, plo)
        sli_rota = add_slider('Rotation [˚] ' , 'rota', 0.3, 0.87, 0.6, 0.025, geo.rota, lmt.rota_min, lmt.rota_max, lmt.rota_stp, fig, ax, geo, det, plo)
        sli_yoff = add_slider('Offset [mm] '  , 'yoff', 0.3, 0.89, 0.6, 0.025, geo.yoff, lmt.yoff_min, lmt.yoff_max, lmt.yoff_stp, fig, ax, geo, det, plo)
        sli_dist = add_slider('Distance [mm] ', 'dist', 0.3, 0.91, 0.6, 0.025, geo.dist, lmt.dist_min, lmt.dist_max, lmt.dist_stp, fig, ax, geo, det, plo)
        sli_ener = add_slider('Energy [keV] ' , 'ener', 0.3, 0.93, 0.6, 0.025, geo.ener, lmt.ener_min, lmt.ener_max, lmt.ener_stp, fig, ax, geo, det, plo)
    # show the plot
    plt.show()

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
    # calculate the y offset resulting from the tilt and rotation
    _geo_offset = -(geo.yoff + np.deg2rad(geo.rota)*geo.dist)
    # draw the cones slightly larger than the visible detector
    _geo_margin = 1.25
    # draw beam center
    ax.plot(0, 0, color=plo.cont_orig_color, marker=plo.cont_orig_cmark, ms=plo.cont_orig_csize, alpha=plo.cont_orig_alpha)
    ax.plot(0, _geo_offset, color=colors.to_hex(plo.cont_geom_cmap(1)), marker=plo.cont_geom_cmark, ms=plo.cont_geom_csize, alpha=plo.cont_geom_alpha)
    # draw contour lines
    for n,_scale in enumerate(plo.cont_levels):
        # prepare the grid for the cones/contours
        # adjust the resolution using i (-> plo.cont_levels),
        # as smaller cones/contours (large i) need higher sampling
        # but make sure the sampling rate doesn't fall below the
        # user set plo.cont_reso value
        _x0 = np.linspace(-_geo_margin*plo.cont_grid_max, _geo_margin*plo.cont_grid_max, max(int(plo.cont_reso*_scale), plo.cont_reso))
        # the grid position needs to adjusted upon change of geometry
        # the center needs to be shifted by _geo_offset to make sure 
        # sll contour lines are drawn
        _x1 = np.linspace(-_geo_margin*plo.cont_grid_max-_geo_offset, _geo_margin*plo.cont_grid_max-_geo_offset, max(int(plo.cont_reso*_scale), plo.cont_reso))
        # calculate resolution rings
        # 2-Theta: np.arctan(dist/(dist/i))
        _thr = np.arctan(1/_scale)/2
        # Conversion factor keV to Angstrom: 12.398
        # sin(t)/l: np.sin(Theta) / lambda -> (12.398/geo_energy)
        _stl = np.sin(_thr)/(12.398/geo.ener)
        # d-spacing: l = 2 d sin(t) -> 1/2(sin(t)/l)
        _dsp = 1/(2*_stl)
        # figure out the labels
        _units = {'n':None, 't':np.rad2deg(2*_thr),
                  'd':_dsp, 'q':_stl*4*np.pi, 's':_stl}
        # draw contours for the original geometry
        if plo.origin:
            X0, Y0 = np.meshgrid(_x0,_x0)
            Z0 = np.sqrt(X0**2+Y0**2)*_scale
            X,Y,Z = geo_cone(X0, Y0, Z0, 0, 0, 0, geo.dist)
            # don't draw contour lines that are out of bounds
            # make sure Z is large enough to draw the contour
            if np.max(Z) >= geo.dist:
                c0 = ax.contour(X, Y, Z, [geo.dist], colors=plo.cont_orig_color, alpha=plo.cont_orig_alpha)
                # label original geometry contours
                fmt = {c0.levels[0]:f'{np.round(_units[geo.unit],2):.2f}'}
                ax.clabel(c0, c0.levels, inline=True, fontsize=plo.cont_fsize, fmt=fmt, manual=[(plo.xdim,plo.ydim)])
        # draw contours for the tilted/rotated/moved geometry
        # use the offset adjusted value x1 to prepare the grid
        X0, Y0 = np.meshgrid(_x1,_x0)
        Z0 = np.sqrt(X0**2+Y0**2)*_scale
        X,Y,Z = geo_cone(X0, Y0, Z0, geo.rota, geo.tilt, geo.yoff, geo.dist)
        # make sure Z is large enough to draw the contour
        if np.max(Z) > geo.dist:
            # make sure the rotated geometry contour can
            # be drawn, is within the adjusted grid
            if plo.cont_grid_max -_geo_offset - np.min(Z) > 0:
                c1 = ax.contour(X, Y, Z, [geo.dist], colors=colors.to_hex(plo.cont_geom_cmap((n+1)/len(plo.cont_levels))), alpha=plo.cont_geom_alpha)
                # label moved geometry contours
                fmt = {c1.levels[0]:f'{np.round(_units[geo.unit],2):.2f}'}
                ax.clabel(c1, c1.levels, inline=True, fontsize=plo.cont_fsize, fmt=fmt, manual=[(0,plo.ydim)])
            else:
                print(plo.cont_grid_max, np.min(Z))
        else:
            # if the Z*i is too small, break as the following cycles
            # will make Z only smaller -> leaving no contours to draw
            # - only True if the contours are iterated high to low!
            break

def geo_cone(X, Y, Z, rota, tilt, yoff, dist):
    # rotate the sample around y
    a = np.deg2rad(tilt) + np.deg2rad(rota)
    t = np.transpose(np.array([X,Y,Z]), (1,2,0))
    m = [[np.cos(a), 0, np.sin(a)],[0,1,0],[-np.sin(a), 0, np.cos(a)]]
    X,Y,Z = np.transpose(np.dot(t, m), (2,0,1))
    # compensate for tilt
    comp = np.deg2rad(tilt) * dist
    return Y,X+comp-yoff,Z

def add_slider(label, name, left, bottom, width, height, val, vmin, vmax, step, fig, ax, geo, det, plo):
    axs = fig.add_axes([left, bottom, width, height])
    sli = Slider(axs, label, valmin=vmin, valmax=vmax, valinit=val, handle_style={'size':plo.cont_fsize}, valstep=step)
    sli.vline.set_alpha(0) # Remove the mark on the slider
    sli.on_changed(lambda val: update_plot(name, val, fig, geo, plo, det, ax))
    sli.label.set_size(plo.cont_fsize)
    sli.valtext.set_size(plo.cont_fsize)
    return sli

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
    elif nam == 'ener':
        geo.ener = float(val)
    # we need to clear the axis
    # deleting the contours and labels
    # individually didn't work!
    ax.clear()
    # re-adjust the axis
    ax.set_aspect('equal')
    ax.set_axis_off()
    ax.set_xlim(-plo.xdim, plo.xdim)
    ax.set_ylim(-plo.ydim, plo.ydim)
    # re-calculate cones and re-draw contours
    draw_contours(ax, geo, plo)
    fig.canvas.blit(ax)

class container(object):
    pass

if __name__ == '__main__':
    main()
