
#include <assert.h>

#include "openmc/geom_util.h"

#include "openmc/constants.h"

namespace openmc {
 const bool EXIT_EARLY = false;


    /* Function to return the vertex with the lowest coordinates. To force the same
       ray-edge computation, the Plücker test needs to use consistent edge
       representation. This would be more simple with MOAB handles instead of
       coordinates... */
    inline bool first( const Position& a, const Position& b )
    {
        if( a[0] < b[0] ) { return true; }
        else if( a[0] == b[0] )
        {
            if( a[1] < b[1] ) { return true; }
            else if( a[1] == b[1] )
            {
                if( a[2] < b[2] ) { return true; }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }

    double plucker_edge_test( const Position& vertexa, const Position& vertexb, const Position& ray,
                              const Position& ray_normal )
    {

        double pip;
        const double near_zero = 10 * std::numeric_limits< double >::epsilon();

        if( first( vertexa, vertexb ) )
        {
            const Position edge        = vertexb - vertexa;
            const Position edge_normal = edge.cross(vertexa);
            pip                        = ray.dot(edge_normal) + ray_normal.dot(edge);
        }
        else
        {
            const Position edge        = vertexa - vertexb;
            const Position edge_normal = edge.cross(vertexb);
            pip                        = ray.dot(edge_normal) + ray_normal.dot(edge);
            pip                        = -pip;
        }

        if( near_zero > fabs( pip ) ) pip = 0.0;

        return pip;
    }

    bool plucker_ray_tri_intersect( const std::array<Position, 3> vertices, const Position& origin, const Position& direction,
                                    double& dist_out, const double* neg_ray_len,
                                    const int* orientation)
    {
        const double nonneg_ray_len = INFTY;
        dist_out = INFTY;

        const Position raya = direction;
        const Position rayb = direction.cross(origin);

        // Determine the value of the first Plucker coordinate from edge 0
        double plucker_coord0 = plucker_edge_test(vertices[0], vertices[1], raya, rayb);

        // If orientation is set, confirm that sign of plucker_coordinate indicate
        // correct orientation of intersection
        if( orientation && ( *orientation ) * plucker_coord0 > 0 ) { return EXIT_EARLY; }

        // Determine the value of the second Plucker coordinate from edge 1
        double plucker_coord1 = plucker_edge_test( vertices[1], vertices[2], raya, rayb );

        // If orientation is set, confirm that sign of plucker_coordinate indicate
        // correct orientation of intersection
        if( orientation )
        {
            if( ( *orientation ) * plucker_coord1 > 0 ) { return EXIT_EARLY; }
            // If the orientation is not specified, all plucker_coords must be the same sign or
            // zero.
        }
        else if( ( 0.0 < plucker_coord0 && 0.0 > plucker_coord1 ) || ( 0.0 > plucker_coord0 && 0.0 < plucker_coord1 ) )
        {
            return EXIT_EARLY;
        }

        // Determine the value of the second Plucker coordinate from edge 2
        double plucker_coord2 = plucker_edge_test( vertices[2], vertices[0], raya, rayb );

        // If orientation is set, confirm that sign of plucker_coordinate indicate
        // correct orientation of intersection
        if( orientation )
        {
            if( ( *orientation ) * plucker_coord2 > 0 ) { return EXIT_EARLY; }
            // If the orientation is not specified, all plucker_coords must be the same sign or
            // zero.
        }
        else if( ( 0.0 < plucker_coord1 && 0.0 > plucker_coord2 ) || ( 0.0 > plucker_coord1 && 0.0 < plucker_coord2 ) ||
                 ( 0.0 < plucker_coord0 && 0.0 > plucker_coord2 ) || ( 0.0 > plucker_coord0 && 0.0 < plucker_coord2 ) )
        {
            return EXIT_EARLY;
        }

        // check for coplanar case to avoid dividing by zero
        if( 0.0 == plucker_coord0 && 0.0 == plucker_coord1 && 0.0 == plucker_coord2 ) { return EXIT_EARLY; }

        // get the distance to intersection
        const double inverse_sum = 1.0 / ( plucker_coord0 + plucker_coord1 + plucker_coord2 );
        assert( 0.0 != inverse_sum );
        const Position intersection( plucker_coord0 * inverse_sum * vertices[2] +
                                     plucker_coord1 * inverse_sum * vertices[0] +
                                     plucker_coord2 * inverse_sum * vertices[1] );

        // To minimize numerical error, get index of largest magnitude direction.
        int idx            = 0;
        double max_abs_dir = 0;
        for( unsigned int i = 0; i < 3; ++i )
        {
            if( fabs( direction[i] ) > max_abs_dir )
            {
                idx         = i;
                max_abs_dir = fabs( direction[i] );
            }
        }
        const double dist = ( intersection[idx] - origin[idx] ) / direction[idx];

        // is the intersection within distance limits?
        // if( ( nonneg_ray_len && nonneg_ray_len < dist ) ||  // intersection is beyond positive limit
        //     ( neg_ray_len && *neg_ray_len >= dist ) ||       // intersection is behind negative limit
        //     ( !neg_ray_len && 0 > dist ) )
        // {  // Unless a neg_ray_len is used, don't return negative distances
        //     return EXIT_EARLY;
        // }

        dist_out = dist;

        // if( type )
        //     *type = type_list[( ( 0.0 == plucker_coord2 ) << 2 ) + ( ( 0.0 == plucker_coord1 ) << 1 ) +
        //                       ( 0.0 == plucker_coord0 )];

        return true;
    }

// bool first(const Position& v0, const Position& v1) {
//   if( v0[0] < v1[0] ) { return true; }
//     else if( v0[0] == v1[0] ) {
//       if( v0[1] < v1[1] ) { return true; }
//       else if( v0[1] == v1[1] ) {
//         if( v0[2] < v1[2] ) { return true; }
//         else { return false; }
//       } else {
//         return false;
//       }
//     } else {
//       return false;
//   }
// }

// double plucker_edge_test(const Position& v0, const Position& v1, const Position& u, const Position& ray_cross) {
//   double pip;
//   if (first(v0, v1)) {
//     const Position edge = v1 - v0;
//     const Position edge_normal = edge.cross(v0);
//     pip = u.dot(edge_normal) + ray_cross.dot(edge);
//   } else {
//     const Position edge = v0 - v1;
//     const Position edge_normal = edge.cross(v1);
//     pip = -(u.dot(edge_normal) + ray_cross.dot(edge));
//   }

//   if (NEAR_ZERO > std::abs(pip)) pip = 0.0;

//   return pip;
// }

// IntersectionType plucker_intersect(const std::array<Position, 3> coords, const Position& r, const Direction& u,
// double& dist_out) {
//   dist_out = INFTY;

//   const Position ray_cross = u.cross(r);

//   double plucker0 = plucker_edge_test(coords[0], coords[1], u, ray_cross);
//   double plucker1 = plucker_edge_test(coords[1], coords[2], u, ray_cross);
//   double plucker2 = plucker_edge_test(coords[2], coords[0], u, ray_cross);

//   // check on the correct side of all edges
//   if (plucker0 > 0.0 || plucker1 > 0.0 || plucker2 > 0.0) return IntersectionType::MISS;

//   if ( (0.0 < plucker1 && 0.0 > plucker2) || (0.0 > plucker1 && 0.0 < plucker2) ||
//        (0.0 < plucker0 && 0.0 > plucker2) || (0.0 > plucker0 && 0.0 < plucker2) ) return IntersectionType::MISS;

//   // check for a copanar intersection
//   if (plucker0 == 0.0 && plucker1 == 0.0 && plucker2 == 0.0) return IntersectionType::MISS;

//   // compute the distance to intersection
//   const double inv_sum = 1.0 / (plucker0 + plucker1 + plucker2);
//   const Position intersection {plucker0 * inv_sum * coords[2] +
//                                plucker1 * inv_sum * coords[0] +
//                                plucker2 * inv_sum * coords[1]};

//   // get index of largest magnitude in direction
//   int max_idx = 0;
//   double max_dir = std::max({std::abs(u.x), std::abs(u.y), std::abs(u.z)});
//   for (int i = 0; i < 3; i++) {
//     if (max_dir == u[i]) max_idx = i;
//   }

//   dist_out = (intersection[max_idx] - r[max_idx]) / u[max_idx];

//   if (plucker0 == 0.0) return IntersectionType::EDGE0;
//   if (plucker1 == 0.0) return IntersectionType::EDGE1;
//   if (plucker2 == 0.0) return IntersectionType::EDGE2;
//   if (plucker0 == 0.0 && plucker1 == 0.0) return IntersectionType::NODE0;
//   if (plucker0 == 0.0 && plucker1 == 0.0) return IntersectionType::NODE1;
//   if (plucker1 == 0.0 && plucker2 == 0.0) return IntersectionType::NODE2;

//   return IntersectionType::INTERIOR;
// }

} // namespace openmc