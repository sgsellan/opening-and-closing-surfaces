 #include <iostream>
#include <fstream>
#include <iomanip>
#include <queue>
#include <list>
#include <cmath>
#include <limits>
#include <set>
// #include <matlab>
// Lib IGL includes
#include <igl/adjacency_list.h>
#include <igl/adjacency_matrix.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/avg_edge_length.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/principal_curvature.h>
#include <igl/C_STR.h>
#include <igl/matlab/mexErrMsgTxt.h>
#include <igl/matlab/MexStream.h>
#include <igl/matlab/parse_rhs.h>
#include <igl/matlab/prepare_lhs.h>
#include <igl/matlab/validate_arg.h>
#include "../../src/per_vertex_signed_prin_curvature.h"


// Los que tienen el tipo de objeto definido son input, el resto output??

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    using namespace igl;
    using namespace igl::matlab;
   // using namespace igl::copyleft::cgal;
    using namespace Eigen;
    MatrixXd V,N,PD1,PD2;
    VectorXd PV1,PV2;
    MatrixXi F;
    int radius = 2;
bool useKring = true;
    
    igl::matlab::MexStream mout;
    std::streambuf *outbuf = std::cout.rdbuf(&mout);
    
    mexErrMsgTxt(nrhs>=2,"nrhs should be >= 2");
    parse_rhs_double(prhs,V); // Aquí se pasa del prhs a la matriz V de Eigen
    parse_rhs_index(prhs+1,F); // Aquí se pasa del prhs a la matriz F de Eigen
    mexErrMsgTxt(V.cols()==3,"V must be #V by 3");
    mexErrMsgTxt(F.cols()==3,"F must be #F by 3");
    
    per_vertex_signed_prin_curvature(V,F,PD1,PD2,PV1,PV2);
    

    
    switch(nlhs)
    {
        case 4:
            prepare_lhs_double(PV2,plhs+3);
        case 3:
            prepare_lhs_double(PV1,plhs+2);
        case 2:
            prepare_lhs_double(PD2,plhs+1); // Sustituir G con las direcciones?
        case 1:
            prepare_lhs_double(PD1,plhs+0); // Sustituir W con las curvaturas?
        default:break;
    }
    
    // Restore the std stream buffer Important!
    std::cout.rdbuf(outbuf);
    
}
