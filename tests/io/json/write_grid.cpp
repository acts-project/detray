/** Detray library, part of the ACTS project (R&D line)
 * 
 * (c) 2020 CERN for the benefit of the ACTS project
 * 
 * Mozilla Public License Version 2.0
 */

#include "plugins/array_defs.hpp"
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/populator.hpp"
#include "grids/serializer2.hpp"
#include "json/json_grids.hpp"

#include <fstream>

int main(int argc, char **argv)
{

    using namespace detray;
    std::ofstream output_file;

    replace_populator<guaranteed_index, std::numeric_limits<guaranteed_index>::max()> replacer;
    serializer2 serializer;

    // A rectangular grid 
    axis::closed<> xaxis{10, -5., 5.};
    axis::closed<> yaxis{10, -5., 5.};
    using grid2r = grid2<decltype(replacer), decltype(xaxis), decltype(yaxis), decltype(serializer)>;

    grid2r gr(std::move(xaxis), std::move(yaxis));

    output_file.open("grid_rect.json");
    output_file << json::write_grid(gr).dump(2);
    output_file.close();

    // A polar sectoral grid
    axis::closed<> raxis{3, 1., 5.};
    axis::closed<> phiaxis{10, -0.35, 0.35};

    using grid2ps = grid2<decltype(replacer), decltype(raxis), decltype(phiaxis), decltype(serializer)>;
    grid2ps gps(std::move(raxis), std::move(phiaxis));

    output_file.open("grid_pol_sect.json");
    output_file << json::write_grid(gps).dump(2);
    output_file.close();
    

}