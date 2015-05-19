#include <X11/Intrinsic.h>

/* STATUS RETURNS */
#	define SG_VISUAL_SUCCESS	1
#	define SG_VISUAL_DEFAULT	(SG_VISUAL_SUCCESS + 1)

#	define SG_NO_VISUAL		(-1)
#	define SG_BAD_DISPLAY		(SG_NO_VISUAL-1)
#	define SG_NO_TYPE_AND_CLASS	(SG_BAD_DISPLAY-1)
#	define SG_NO_SUCH_VISUAL	(SG_NO_TYPE_AND_CLASS-1)


/* Visual types */
#	define SG_DEFAULT_PLANES	0	/* matches resource default */
#	define SG_UNDERLAY_PLANES	1
#	define SG_NORMAL_PLANES		2
#	define SG_OVERLAY_PLANES	3
#	define SG_POPUP_PLANES		4
#	define SG_MAX_TYPES		(SG_POPUP_PLANES + 1)


/* External declarations */

#ifdef _SG_IS_STATIC
#define _SG_STORAGE static
#else
#define _SG_STORAGE extern
#endif /* _SG_IS_STATIC */

#ifdef _NO_PROTO
_SG_STORAGE int		_SG_defaultDepthAndTypeResources ();
_SG_STORAGE Colormap	_SG_getDefaultColorMap();
_SG_STORAGE int		_SG_getDefaultDepth ();
_SG_STORAGE XVisualInfo *_SG_getMatchingVisual ();
_SG_STORAGE int		_SG_getMaxDepth ();
_SG_STORAGE int		_SG_getNormalArgs();
_SG_STORAGE int		_SG_getOverlayArgs();
_SG_STORAGE int		 iv_SG_getPopupArgs();
_SG_STORAGE int		_SG_getUnderlayArgs();
_SG_STORAGE int		_SG_getVisualArgs();
#else /* _NO_PROTO */
_SG_STORAGE int	    _SG_defaultDepthAndTypeResources
	( Display *display, int screen, int *requestedClass,
	  char *requestedType, int *requestedTypeV, int *requestedDepth,
	  Visual **requestedVisual, Colormap *requestedColormap,
	  Drawable *requestedDrawable);
_SG_STORAGE Colormap _SG_getDefaultColormap
	(Display *dpy, int scr, Visual *vsl);
_SG_STORAGE int      _SG_getDefaultDepth
	(Display *dpy, int scr, int *class, int type);
_SG_STORAGE XVisualInfo *_SG_getMatchingVisual
	(Display *dpy, int scr, VisualID vsl, int *class, int type, int depth);
_SG_STORAGE int      _SG_getMaxDepth
	(Display *dpy, int scr, int *class, int type);
_SG_STORAGE int      _SG_getNormalArgs
	(Display *dpy, int scr, ArgList args, int *n);
_SG_STORAGE int     _SG_getOverlayArgs
	(Display *dpy, int scr, ArgList args, int *n);
_SG_STORAGE int     _SG_getOverlay2Args
	(Display *dpy, int scr, ArgList args, int *n);
_SG_STORAGE int     _SG_getOverlay4Args
	(Display *dpy, int scr, ArgList args, int *n);
_SG_STORAGE int       iv_SG_getPopupArgs
	(Display *dpy, int scr, ArgList args, int *n);
_SG_STORAGE int      _SG_getUnderlayArgs
	(Display *dpy, int scr, ArgList args, int *n);
_SG_STORAGE int      _SG_getVisualArgs
	(Display *dpy, int scr, int dpth, int *class, int type, ArgList args, int *n);
#endif /* _NO_PROTO */
