#ifndef OPENMC_POSITION_H
#define OPENMC_POSITION_H

namespace openmc {

//==============================================================================
//! Type representing a position in Cartesian coordinates
//==============================================================================

struct Position {
  // Constructors
  Position() = default;
  Position(double x_, double y_, double z_) : x{x_}, y{y_}, z{z_} { };
  Position(const double xyz[]) : x{xyz[0]}, y{xyz[1]}, z{xyz[2]} { };

  // Unary operators
  Position& operator+=(Position);
  Position& operator+=(double);
  Position& operator-=(Position);
  Position& operator-=(double);
  Position& operator*=(Position);
  Position& operator*=(double);
  const double& operator[](int i) const { return xyz[i]; }
  
        double& operator[](int i)       { return xyz[i]; }

  // Other member functions

  //! Dot product of two vectors
  //! \param[in] other Vector to take dot product with
  //! \result Resulting dot product
  inline double dot(Position other) {
    return x*other.x + y*other.y + z*other.z;
  }

  // Data members
  union{ double xyz[3]; struct{ double x,y,z; }; };
};

// Binary operators
inline Position operator+(Position a, Position b) { return a += b; }
inline Position operator+(Position a, double b)   { return a += b; }
inline Position operator+(double a, Position b)   { return b += a; }

inline Position operator-(Position a, Position b) { return a -= b; }
inline Position operator-(Position a, double b)   { return a -= b; }
inline Position operator-(double a, Position b)   { return b -= a; }

inline Position operator*(Position a, Position b) { return a *= b; }
inline Position operator*(Position a, double b)   { return a *= b; }
inline Position operator*(double a, Position b)   { return b *= a; }

//==============================================================================
//! Type representing a vector direction in Cartesian coordinates
//==============================================================================

using Direction = Position;

} // namespace openmc

#endif // OPENMC_POSITION_H
