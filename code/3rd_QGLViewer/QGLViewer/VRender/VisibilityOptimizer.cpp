#include <vector>
#include "VRender.h"
#include "Optimizer.h"
#include "Primitive.h"
#include "gpc.h"
#include "math.h"

using namespace vrender ;
using namespace std ;

#ifdef A_FAIRE
void VisibilityOptimizer::optimize(vector<PtrPrimitive>& primitives,float& percentage_finished,string& message)
#else
void VisibilityOptimizer::optimize(vector<PtrPrimitive>& primitives,VRenderParams& vparams)
#endif
{
#ifdef DEBUG_VO
        cout << "Optimizing visibility." << endl ;
#endif
        unsigned long N = primitives.size()/200 + 1 ;

#ifdef DEBUG_EPSRENDER__SHOW1
//	cout << "Showing viewer." << endl ;
//	myViewer viewer ;
//	viewer.show();
        double minx =  FLT_MAX ;
        double miny =  FLT_MAX ;
        double maxx = -FLT_MAX ;
        double maxy = -FLT_MAX ;
        for(unsigned int i=0;i<primitives.size();++i)
                for(int j=0;j<primitives[i]->nbVertices();++j)
                {
                        if(maxx < primitives[i]->vertex(j).x()) maxx = primitives[i]->vertex(j).x() ;
                        if(maxy < primitives[i]->vertex(j).y()) maxy = primitives[i]->vertex(j).y() ;
                        if(minx > primitives[i]->vertex(j).x()) minx = primitives[i]->vertex(j).x() ;
                        if(miny > primitives[i]->vertex(j).y()) miny = primitives[i]->vertex(j).y() ;
                }

        glMatrixMode(GL_PROJECTION) ;
        glLoadIdentity() ;
        glOrtho(minx,maxx,miny,maxy,-1,1) ;
        glMatrixMode(GL_MODELVIEW) ;
        glLoadIdentity() ;

        cout << "Window set to " << minx << " " << maxx << " " << miny << " " << maxy << endl ;
        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT) ;
        glLineWidth(3.0) ;
#endif

        int nb_culled = 0 ;

        // Ca serait pas mal mieux avec une interface c++...

        gpc_polygon cumulated_union ;
        cumulated_union.num_contours = 0 ;
        cumulated_union.hole = nullptr ;
        cumulated_union.contour = nullptr ;
        size_t nboptimised = 0 ;

        for(size_t pindex = primitives.size() - 1; long(pindex) >= 0;--pindex,++nboptimised)
                if(primitives[pindex] != nullptr)
                {
#ifdef A_FAIRE
                        percentage_finished = pindex / (float)primitives.size() ;
#endif

                        if(primitives[pindex]->nbVertices() > 1)
                        {
#ifdef DEBUG_VO
                                if(pindex%50==0)
                                {
                                        char buff[500] ;
                                        sprintf(buff,"Left: % 6ld - Culled: % 6ld", pindex,(long)nb_culled) ;
                                        fprintf(stdout,buff);

                                        for(unsigned int j=0;j<strlen(buff);++j)
                                                fprintf(stdout,"\b") ;

                                        fflush(stdout) ;
                                }
#endif

                                try
                                {
                                        PtrPrimitive p(primitives[pindex]) ;
                                        gpc_polygon difference ;
                                        gpc_polygon new_poly ;
                                        gpc_polygon new_poly_reduced ;
                                        new_poly.num_contours = 0 ;
                                        new_poly.hole = nullptr ;
                                        new_poly.contour = nullptr ;
                                        new_poly_reduced.num_contours = 0 ;
                                        new_poly_reduced.hole = nullptr ;
                                        new_poly_reduced.contour = nullptr ;

                                        // 1 - creates a gpc_polygon corresponding to the current primitive

                                        gpc_vertex_list *new_poly_verts = new gpc_vertex_list ;
                                        gpc_vertex_list *new_poly_reduced_verts = new gpc_vertex_list ;

                                        double mx = 0.0 ;
                                        double my = 0.0 ;

                                        if(p->nbVertices() == 2)
                                        {
                                                new_poly_verts->num_vertices = 4 ;
                                                new_poly_verts->vertex = new gpc_vertex[4] ;
                                                new_poly_reduced_verts->num_vertices = 4 ;
                                                new_poly_reduced_verts->vertex = new gpc_vertex[4] ;

                                                double deps = 0.001 ;
                                                double du = p->vertex(1).y()-p->vertex(0).y() ;
                                                double dv = p->vertex(1).x()-p->vertex(0).x() ;
                                                double n = sqrt(du*du+dv*dv) ;
                                                du *= deps/n ;
                                                dv *= deps/n ;
                                                new_poly_verts->vertex[0].x = p->vertex(0).x() + du ;
                                                new_poly_verts->vertex[0].y = p->vertex(0).y() + dv ;
                                                new_poly_verts->vertex[1].x = p->vertex(1).x() + du ;
                                                new_poly_verts->vertex[1].y = p->vertex(1).y() + dv ;
                                                new_poly_verts->vertex[2].x = p->vertex(1).x() - du ;
                                                new_poly_verts->vertex[2].y = p->vertex(1).y() - dv ;
                                                new_poly_verts->vertex[3].x = p->vertex(0).x() - du ;
                                                new_poly_verts->vertex[3].y = p->vertex(0).y() - dv ;

                                                new_poly_reduced_verts->vertex[0].x = p->vertex(0).x() + du ;
                                                new_poly_reduced_verts->vertex[0].y = p->vertex(0).y() + dv ;
                                                new_poly_reduced_verts->vertex[1].x = p->vertex(1).x() + du ;
                                                new_poly_reduced_verts->vertex[1].y = p->vertex(1).y() + dv ;
                                                new_poly_reduced_verts->vertex[2].x = p->vertex(1).x() - du ;
                                                new_poly_reduced_verts->vertex[2].y = p->vertex(1).y() - dv ;
                                                new_poly_reduced_verts->vertex[3].x = p->vertex(0).x() - du ;
                                                new_poly_reduced_verts->vertex[3].y = p->vertex(0).y() - dv ;
                                        }
                                        else
                                        {
                                                new_poly_verts->num_vertices = p->nbVertices() ;
                                                new_poly_verts->vertex = new gpc_vertex[p->nbVertices()] ;

                                                for(size_t i=0;i<p->nbVertices();++i)
                                                {
                                                        new_poly_verts->vertex[i].x = p->vertex(i).x() ;
                                                        new_poly_verts->vertex[i].y = p->vertex(i).y() ;
                                                        mx += p->vertex(i).x() ;
                                                        my += p->vertex(i).y() ;
                                                }
                                                mx /= p->nbVertices() ;
                                                my /= p->nbVertices() ;

                                                new_poly_reduced_verts->num_vertices = p->nbVertices() ;
                                                new_poly_reduced_verts->vertex = new gpc_vertex[p->nbVertices()] ;

                                                for(size_t j=0;j<p->nbVertices();++j)
                                                {
                                                        new_poly_reduced_verts->vertex[j].x = mx + (p->vertex(j).x() - mx)*0.999 ;
                                                        new_poly_reduced_verts->vertex[j].y = my + (p->vertex(j).y() - my)*0.999 ;
                                                }
                                        }
                                        gpc_add_contour(&new_poly,new_poly_verts,false) ;
                                        gpc_add_contour(&new_poly_reduced,new_poly_reduced_verts,false) ;

                                        // 2 - computes the difference between this polygon, and the union of the
                                        // 	preceeding ones.

                                        gpc_polygon_clip(GPC_DIFF,&new_poly_reduced,&cumulated_union,&difference) ;

                                        // 3 - checks the difference. If void, the primitive is not visible: skip it
                                        // 	and go to next primitive.

                                        if(difference.num_contours == 0)
                                        {
                                                ++nb_culled ;
                                                delete p ;
                                                primitives[pindex] = nullptr ;
                                                continue ;
                                        }

                                        // 4 - The primitive is visible. Let's add it to the cumulated union of
                                        // 	primitives.

                                        if(p->nbVertices() > 2)
                                        {
                                                gpc_polygon cumulated_union_tmp ;
                                                cumulated_union_tmp.num_contours = 0 ;
                                                cumulated_union_tmp.hole = nullptr ;
                                                cumulated_union_tmp.contour = nullptr ;

                                                gpc_polygon_clip(GPC_UNION,&new_poly,&cumulated_union,&cumulated_union_tmp) ;

                                                gpc_free_polygon(&cumulated_union) ;
                                                cumulated_union = cumulated_union_tmp ;
                                        }

                                        gpc_free_polygon(&new_poly) ;
                                        gpc_free_polygon(&new_poly_reduced) ;
                                        gpc_free_polygon(&difference) ;
#ifdef DEBUG_EPSRENDER__SHOW1
                                        glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT) ;

                                        glColor3f(1.0,0.0,0.0) ;

                                        for(unsigned long i=0;i<cumulated_union.num_contours;++i)
                                        {
                                                glBegin(GL_LINE_LOOP) ;
                                                for(unsigned long j=0;j<cumulated_union.contour[i].num_vertices;++j)
                                                        glVertex2f(cumulated_union.contour[i].vertex[j].x,cumulated_union.contour[i].vertex[j].y) ;
                                                glEnd() ;
                                        }

                                        glFlush() ;
                                        glXSwapBuffers(glXGetCurrentDisplay(),glXGetCurrentDrawable()) ;
#endif
                                }
                                catch(exception& )
                                {
                                        ; // std::cout << "Could not treat primitive " << pindex << ": internal gpc error." << endl ;
                                }
                        }

                        if(nboptimised%N==0)
                                vparams.progress(nboptimised/(float)primitives.size(), QGLViewer::tr("Visibility optimization")) ;
                }

#ifdef DEBUG_VO
        cout << nb_culled << " primitives culled over " << primitives.size() << "." << endl ;
#endif

        gpc_free_polygon(&cumulated_union) ;
}


