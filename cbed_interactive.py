#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

class TwoCircles:
    def __init__(self, ax, c1, c2, r):

        self.c1 = plt.Circle(c1, r, ec="k", fill=False)
        self.c2 = plt.Circle(c2, r, ec="k", fill=False, ls="--")

        self.ax = ax
        self.ax.add_artist(self.c1)
        self.ax.add_artist(self.c2)

        self.drag = None

        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.ax.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        if event.inaxes != self.ax:
            return

        dx = (event.xdata - self.c1.center[0])
        dy = (event.ydata - self.c1.center[1])
        r1 = np.sqrt(dx*dx + dy*dy)/self.c1.radius

        dx = (event.xdata - self.c2.center[0])
        dy = (event.ydata - self.c2.center[1])
        r2 = np.sqrt(dx*dx + dy*dy)/self.c2.radius
        
        f = 0.25
        if r1 < f:
            self.drag = "center", self.c1
        elif r2 < f:
            self.drag = "center", self.c2
        elif (1-f) < r1 < (1+f):
            self.drag = "radius", self.c1
        elif (1-f) < r2 < (1+f):
            self.drag = "radius", self.c2

    def on_motion(self, event):
        if event.inaxes != self.ax or self.drag is None:
            return
        
        d, w = self.drag
        if d == "center":
            w.set_center( (event.xdata, event.ydata) )
            self.ax.figure.canvas.draw_idle()

        elif d == "radius":
            dx = event.xdata - w.center[0]
            dy = event.ydata - w.center[1]
            r = np.sqrt(dx*dx + dy*dy)

            self.c1.set_radius( r ) 
            self.c2.set_radius( r )
            self.ax.figure.canvas.draw_idle()

    def on_release(self, event):
        self.drag = None


if __name__ == "__main__":
    import sys, json
    from PIL import Image

    img_fname = sys.argv[1]
    info_fname = img_fname+".json"

    #
    # Calibrate CBED size, spacing and rotation
    # 

    img = Image.open(img_fname)
    data = np.asarray(img)

    plt.figure("Image calibration")
    plt.imshow(data)
    
    ax = plt.gca()
    
    try:
        with open(info_fname, "r") as fh:
            info = json.load(fh)
    except:
        a, b = ax.get_xlim()
        c, d = ax.get_ylim()
        c1 = (2*a+b)/3, (c+d)/2
        c2 = (a+2*b)/3, (c+d)/2
        r = (b-a)/6
        info = dict(c1=c1, c2=c2, r=r)

    c = TwoCircles( ax, info["c1"], info["c2"], info["r"])

    plt.show()
    
    c1, c2, r = c.c1.center, c.c2.center, c.c2.radius
    info.update(c1=c1, c2=c2, r=r)

    #
    # Extract data after rotation and scaling
    #

    angle = np.arctan2( c2[1] - c1[1], (c2[0] - c1[0]) )

    from scipy import ndimage
    data = ndimage.rotate(data, np.rad2deg(angle), reshape=False)
    c = data.shape[1]/2, data.shape[0]/2

    c1 = ( np.cos(angle)*(c1[0]-c[0]) + np.sin(angle)*(c1[1]-c[1]) + c[0],
           np.cos(angle)*(c1[1]-c[1]) - np.sin(angle)*(c1[0]-c[0]) + c[1] )

    c2 = ( np.cos(angle)*(c2[0]-c[0]) + np.sin(angle)*(c2[1]-c[1]) + c[0],
           np.cos(angle)*(c2[1]-c[1]) - np.sin(angle)*(c2[0]-c[0]) + c[1] )

    y, x = np.mgrid[0:data.shape[0],0:data.shape[1]]
    m1 = (x-c1[0])*(x-c1[0]) + (y-c1[1])*(y-c1[1]) < r*r
    m2 = (x-c2[0])*(x-c2[0]) + (y-c2[1])*(y-c2[1]) < r*r

    if 0:
        plt.figure()
        plt.imshow(data)
        #plt.imshow(m1 | m2)
        ax = plt.gca()
        TwoCircles(ax, c1, c2, r)

    a, b = int(c1[1] - r), int(c1[1]+r)
    line1 = (data*m1).sum(0)/m1.sum(0)
    line2 = (data*m2).sum(0)/m2.sum(0)
    m2l = m2.any(0)
    x = (np.arange(data.shape[1]) - c1[0])/(c2[0]-c1[0])

    plt.figure("CBED Fit")
    ax_img = plt.axes([0.1, 0.65, 0.8, 0.30])
    ax_img.imshow(data[a:b]*(m1|m2)[a:b])

    ax_line = plt.axes([0.1, 0.30, 0.8, 0.30])
    ax_line.plot(x, line1)
    ax_line.plot(x, line2)
    line_calc, = ax_line.plot([])

    #
    # Fit calculated profile
    #

    ax_thick = plt.axes([0.1, 0.15, 0.8, 0.03])
    ax_pos = plt.axes([0.1, 0.09, 0.8, 0.03])
    ax_xi_g = plt.axes([0.1, 0.03, 0.8, 0.03])

    slider_thick = plt.Slider(ax_thick, 'thick [nm]', 0., 500., valinit=info.get("thick", 100))
    slider_pos = plt.Slider(ax_pos, 'pos', -.5, .5, valinit=info.get("pos", 0))
    slider_xi_g = plt.Slider(ax_xi_g, 'xi_g [nm]', 0., 200., valinit=info.get("xi_g", 100))

    g_len = info.setdefault("g_len", 7.07) # Al (220) 
    k_len = info.setdefault("k_len", 398.)

    def update(*args):
        z = slider_thick.val
        pos = slider_pos.val
        xi_g = slider_xi_g.val

        #s_g = np.sqrt(k_len*k_len + 2*(x-1+pos)*g_len*g_len ) - k_len
        s_g = (x-1+pos)*g_len*g_len/k_len

        s_ef = np.sqrt(s_g*s_g + 1/xi_g/xi_g)
        I = np.sin(np.pi*z*s_ef)**2/(xi_g*s_ef)**2
        
        #
        # least square fit of scale and offset
        #
        scale, offset, slope  = np.linalg.lstsq(np.column_stack([I, np.ones_like(x), x])[m2l,:], line2[m2l], rcond=None)[0]
        line_calc.set_data(x[m2l], (I*scale + offset + slope*x)[m2l])


        ax_line.figure.canvas.draw_idle()

    slider_thick.on_changed(update)
    slider_pos.on_changed(update)
    slider_xi_g.on_changed(update)
    update()

    plt.show()

    info.update(thick=slider_thick.val, pos=slider_pos.val, xi_g=slider_xi_g.val)
    with open(info_fname, "w") as fh:
        json.dump(info, fh, indent=True)

