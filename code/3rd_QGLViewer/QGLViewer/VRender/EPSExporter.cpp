#include <stdio.h>
#include "Primitive.h"
#include "Exporter.h"
#include "math.h"

using namespace vrender ;
using namespace std ;

const double EPSExporter::EPS_GOURAUD_THRESHOLD = 0.05 ;
const char *EPSExporter::CREATOR = "VRender library - (c) Cyril Soler 2005" ;

float EPSExporter::last_r = -1.0 ;
float EPSExporter::last_g = -1.0 ;
float EPSExporter::last_b = -1.0 ;

EPSExporter::EPSExporter()
{
	last_r = -1 ;
	last_g = -1 ;
	last_b = -1 ;
}

void EPSExporter::writeHeader(QTextStream& out) const
{
	/* Emit EPS header. */

	out << "%!PS-Adobe-2.0 EPSF-2.0\n";

	out << "%%%%HiResBoundingBox: " << _xmin << " " << _ymin << " " << _xmax << " " << _ymax << "\n";

	out << "%%%%Creator: " << CREATOR << " (using OpenGL feedback)\n";
	out << "%%EndComments\n\ngsave\n\n";

	out << "%\n";
	out << "% Contributors:\n";
	out << "%\n";
	out << "%   Frederic Delhoume (delhoume@ilog.fr):\n";
	out << "%        Gouraud triangle PostScript fragment\n";
	out << "%\n";
	out << "%   Cyril Soler       (csoler@imag.fr):\n";
	out << "%        BSP Sort,\n";
	out << "%        Topological and advanced topological Sort,\n";
	out << "%        Hidden surface removal,\n";
	out << "%        Xfig3.2 (and EPS) format\n";
	out << "%\n\n";

	out << "/threshold " << EPS_GOURAUD_THRESHOLD << " def\n";

	for(int i = 0; GOURAUD_TRIANGLE_EPS[i] != nullptr; i++)
		out << GOURAUD_TRIANGLE_EPS[i] << "\n";
#ifdef A_VOIR
	out <<  "\n" <<  << " setlinewidth\n\n", _lineWidth;
#endif
	/* Clear the background like OpenGL had it. */

	if(_clearBG)
	{
		out << _clearR << " " << _clearG << " " << _clearB << " setrgbcolor\n";
		out << _xmin << " " << _ymin << " " << _xmax << " " << _ymax << " rectfill\n\n";
	}
}

void EPSExporter::writeFooter(QTextStream& out) const
{
	out << "grestore\n\n";

	out << "% uncomment next line to be able to print to a printer.\n";
	out << "% showpage\n";
}

void PSExporter::writeFooter(QTextStream& out) const
{
	out << "showpage\n";
}

const char *EPSExporter::GOURAUD_TRIANGLE_EPS[] =
{
	"/bd{bind def}bind def /triangle { aload pop   ",
	"setrgbcolor aload pop 5 3 roll 4 2 roll 3 2 roll exch moveto lineto ",
	"lineto closepath fill } bd /computediff1 { 2 copy sub abs threshold ",
	"ge {pop pop pop true} { exch 2 index sub abs threshold ge { pop pop ",
	"true} { sub abs threshold ge } ifelse } ifelse } bd /computediff3 { 3 ",
	"copy 0 get 3 1 roll 0 get 3 1 roll 0 get computediff1 {true} { 3 copy ",
	"1 get 3 1 roll 1 get 3 1 roll 1 get computediff1 {true} { 3 copy 2 ",
	"get 3 1 roll  2 get 3 1 roll 2 get computediff1 } ifelse } ifelse } ",
	"bd /middlecolor { aload pop 4 -1 roll aload pop 4 -1 roll add 2 div 5 ",
	"1 roll 3 -1 roll add 2 div 3 1 roll add 2 div 3 1 roll exch 3 array ",
	"astore } bd /gdt { computediff3 { 4 -1 roll aload 7 1 roll ",
	"6 -1 roll pop 3 -1 roll pop add 2 div 3 1 roll add 2 div exch 3 -1 roll ",
	"aload 7 1 roll exch pop 4 -1 roll pop add 2 div 3 1 roll add 2 div ",
	"exch 3 -1 roll aload 7 1 roll pop 3 -1 roll pop add 2 div 3 1 roll add ",
	"2 div exch 7 3 roll 10 -3 roll dup 3 index middlecolor 4 1 roll 2 copy ",
	"middlecolor 4 1 roll 3 copy pop middlecolor 4 1 roll 13 -1 roll aload ",
	"pop 17 index 6 index 15 index 19 index 6 index 17 index 6 array astore ",
	"10 index 10 index 14 index gdt 17 index 5 index 17 index ",
	"19 index 5 index 19 index 6 array astore 10 index 9 index 13 index ",
	"gdt 13 index 16 index 5 index 15 index 18 index 5 index 6 ",
	"array astore 12 index 12 index 9 index gdt 17 index 16 ",
	"index 15 index 19 index 18 index 17 index 6 array astore 10 index 12 ",
	"index 14 index gdt 18 {pop} repeat } { aload pop 5 3 roll ",
	"aload pop 7 3 roll aload pop 9 3 roll 8 index 6 index 4 index add add 3 ",
	"div 10 1 roll 7 index 5 index 3 index add add 3 div 10 1 roll 6 index 4 ",
	"index 2 index add add 3 div 10 1 roll 9 {pop} repeat 3 array astore ",
	"triangle } ifelse } bd",
	nullptr
};

void EPSExporter::spewPolygone(const Polygone *P, QTextStream& out)
{
	int nvertices;
	GLfloat red, green, blue;
	bool smooth;

	nvertices = P->nbVertices() ;

    const Feedback3DColor& vertex = Feedback3DColor(P->sommet3DColor(0)) ;

	if(nvertices > 0)
	{
		red   = vertex.red();
		green = vertex.green();
		blue  = vertex.blue();

		smooth = false;

		for(int i=1;i < nvertices && !smooth; i++)
			if(fabs(red - P->sommet3DColor(i).red()) > 0.01 || fabs(green - P->sommet3DColor(i).green()) > 0.01 || fabs(blue - P->sommet3DColor(i).blue()) > 0.01)
				smooth = true;

		if(smooth && !_blackAndWhite)
		{
			/* Smooth shaded polygon; varying colors at vertices. */
			/* Break polygon into "nvertices-2" triangle fans. */

			for (int j = 0; j < nvertices - 2; j++)
			{
				out <<  "[" << P->sommet3DColor(0).x() << " " << P->sommet3DColor(j + 1).x() << " " << P->sommet3DColor(j + 2).x() 
					<< " "	<< P->sommet3DColor(0).y() << " " << P->sommet3DColor(j + 1).y() << " " << P->sommet3DColor(j + 2).y() << "]";

				out <<  " [" << P->sommet3DColor(0    ).red() << " " << P->sommet3DColor(0    ).green() << " " << P->sommet3DColor(0    ).blue() 
					<< "] [" << P->sommet3DColor(j + 1).red() << " " << P->sommet3DColor(j + 1).green() << " " << P->sommet3DColor(j + 1).blue() 
					<< "] [" << P->sommet3DColor(j + 2).red() << " " << P->sommet3DColor(j + 2).green() << " " << P->sommet3DColor(j + 2).blue() << "] gdt\n";

				last_r = last_g = last_b = -1.0 ;
			}
		}
		else
		{
			/* Flat shaded polygon and white polygons; all vertex colors the same. */

			out <<  "newpath\n";

			if(_blackAndWhite)
				setColor(out,1.0,1.0,1.0) ;
			else
				setColor(out,red,green,blue) ;

			/* Draw a filled triangle. */

			out << P->sommet3DColor(0).x() << " " << P->sommet3DColor(0).y() << " moveto\n";

			for (int i = 1; i < nvertices; i++)
				out << P->sommet3DColor(i).x() << " " << P->sommet3DColor(i).y() << " lineto\n";

			out << "closepath fill\n\n";
		}
	}
}

void EPSExporter::spewSegment(const Segment *S, QTextStream& out)
{
  GLdouble dx, dy;
  GLfloat dr, dg, db, absR, absG, absB, colormax;
  int steps;
  GLdouble xstep=0.0, ystep=0.0;
  GLfloat rstep=0.0, gstep=0.0, bstep=0.0;
  GLdouble xnext=0.0, ynext=0.0, distance=0.0;
  GLfloat rnext=0.0, gnext=0.0, bnext=0.0;

  const Feedback3DColor& P1 = Feedback3DColor(S->sommet3DColor(0)) ;
  const Feedback3DColor& P2 = Feedback3DColor(S->sommet3DColor(1)) ;

  dr = P2.red()   - P1.red();
  dg = P2.green() - P1.green();
  db = P2.blue()  - P1.blue();

  if((!_blackAndWhite)&&(dr != 0 || dg != 0 || db != 0))
  {
	  /* Smooth shaded line. */

	  dx = P2.x() - P1.x();
	  dy = P2.y() - P1.y();

	  distance = sqrt(dx*dx + dy*dy);

	  absR = fabs(dr);
	  absG = fabs(dg);
	  absB = fabs(db);

	  colormax = max(absR, max(absG, absB));
	  steps = int(0.5f + max(1.0, colormax * distance * EPS_SMOOTH_LINE_FACTOR));

	  xstep = dx / steps;
	  ystep = dy / steps;

	  rstep = dr / steps;
	  gstep = dg / steps;
	  bstep = db / steps;

	  xnext = P1.x();
	  ynext = P1.y();
	  rnext = P1.red();
	  gnext = P1.green();
	  bnext = P1.blue();

	  /* Back up half a step; we want the end points to be
		  exactly the their endpoint colors. */

	  xnext -= xstep / 2.0;
	  ynext -= ystep / 2.0;
	  rnext -= rstep / 2.0f;
	  gnext -= gstep / 2.0f;
	  bnext -= bstep / 2.0f;
  }
  else
	  steps = 0; /* Single color line. */

  if(_blackAndWhite)
	  setColor(out,0.0,0.0,0.0) ;
  else
	  setColor(out,P1.red(),P1.green(),P1.blue()) ;

  out << P1.x() << " " << P1.y() << " moveto\n";

  for(int i = 0;i < steps;i++)
  {
	  xnext += xstep;
	  ynext += ystep;
	  rnext += rstep;
	  gnext += gstep;
	  bnext += bstep;

	  out << xnext << " " << ynext << " lineto stroke\n";
	  out << rnext << " " << gnext << " " << bnext << " setrgbcolor\n";
	  out << xnext << " " << ynext << " moveto\n";

	  last_r = last_g = last_b = -1.0 ;
  }
  out << P2.x() << " " << P2.y() << " lineto stroke\n";
}

void EPSExporter::spewPoint(const Point *P, QTextStream& out)
{
	const Feedback3DColor& p = Feedback3DColor(P->sommet3DColor(0)) ;

	if(_blackAndWhite)
		setColor(out,0.0,0.0,0.0) ;
	else
		setColor(out,p.red(),p.green(),p.blue()) ;

	out << p.x() << " " << p.y() << " " << (_pointSize / 2.0) << " 0 360 arc fill\n\n";
}

void EPSExporter::setColor(QTextStream& out, float red, float green, float blue)
{
	if(last_r != red || last_g != green || last_b != blue)
		out << red << " " << green << " " << blue << " setrgbcolor\n";

	last_r = red ;
	last_g = green ;
	last_b = blue ;
}

