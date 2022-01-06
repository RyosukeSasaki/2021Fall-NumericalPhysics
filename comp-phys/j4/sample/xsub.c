#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <stdio.h>
#include <unistd.h>

Display	*d;
Window	w;
Pixmap	p;
GC	gc;

#define SXLEN 640
#define SYLEN 480
#define SXORG 0
#define SYORG 0
#define SBWIDTH 2

unsigned long MyColor();
unsigned long white, black;
unsigned long fg, bg, bd;
unsigned long colors[8];
XSetWindowAttributes	a;

ginit(wname)
char *wname;
{
	int	i, j, n;
	char *c;

	d = XOpenDisplay (NULL);

	black = XBlackPixel(d, 0);
	white = XWhitePixel(d, 0);
	colors[0] = MyColor(d, "black");
	colors[1] = MyColor(d, "blue");
	colors[2] = MyColor(d, "red");
	colors[3] = MyColor(d, "purple");
	colors[4] = MyColor(d, "green");
	colors[5] = MyColor(d, "skyblue");
	colors[6] = MyColor(d, "yellow");
	colors[7] = MyColor(d, "white");
	fg = white;
	bg = black;
	bd = white;

	w = XCreateSimpleWindow(d, RootWindow (d, 0),
		SXORG, SYORG, SXLEN, SYLEN, SBWIDTH, bd, bg);

     XStoreName(d, w, wname);
     XSetIconName(d, w, wname);
/************************************************************
	a.override_redirect = 1;
	XChangeWindowAttributes (d, w , CWOverrideRedirect, &a);
************************************************************/
	XMapWindow(d, w );
	XFlush (d);

	gc = XCreateGC(d, w , 0, 0);
	XSetForeground(d, gc , fg);

	p = XCreatePixmap(d, w, SXLEN, SYLEN, DefaultDepth(d, 0));

	return;
}

gdisp()
{
	XCopyArea(d, p, w, gc, 0, 0, SXLEN, SYLEN, 0, 0);
	XFlush (d);
	return;
}

gpoint(ix1, iy1, icol)
long int *ix1, *iy1, *icol;
{
	XSetForeground (d, gc, colors[*icol]);
	XDrawPoint(d, p, gc, *ix1, *iy1);
	return;
}

gline(ix1, iy1, ix2, iy2, icol)
long int *ix1, *iy1, *ix2, *iy2, *icol;
{
	XSetForeground (d, gc, colors[*icol]);
	XDrawLine(d, p, gc, *ix1, *iy1, *ix2, *iy2);
	return;
}

gbox(ix, iy, ixw, iyw, icol)
long int *ix, *iy, *ixw, *iyw, *icol;
{
/*******
	int x, y;
	unsigned int width, height;

	x = *ix1;
	width=*ix2-*ix1;
	if(width < 0)
	{
		x = *ix2;
		width = -width;
	}

	y=*iy1;
	height=*iy2-*iy1;
	if(height < 0)
	{
		y = *iy2;
		height = -height;
	}


	XSetForeground (d, gc, colors[*icol]);
	XFillRectangle(d, p, gc, x, y, width, height);
*******/
	XSetForeground (d, gc, colors[*icol]);
	XFillRectangle(d, p, gc,
		*ix, *iy, (unsigned int)*ixw, (unsigned int)*iyw);
	return;
}

gcls()
{
	XSetForeground (d, gc, bg);
	XFillRectangle(d, p, gc, 0, 0, SXLEN, SYLEN);
	return;
}

unsigned long MyColor(display, color)
Display *display;
char    *color;
{
     Colormap cmap;
     XColor c0, c1;

     cmap = DefaultColormap(display, 0);
     XAllocNamedColor(display, cmap, color, &c1, &c0);
     return (c1.pixel);
}
