#ifndef VECTOR3D_HPP_
#define VECTOR3D_HPP_
#include <cmath>
#include <string>
struct Vector3d {
    Vector3d(double const x_ = 0., double const y_ = 0., double const z_ = 0.)
        : x(x_), y(y_), z(z_){};

    double dot(Vector3d const &other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    double length2() const {
        return dot(*this);
    }
    double length() const {
        return std::sqrt(length2());
    }

    Vector3d operator-(Vector3d const &other) const {
        return Vector3d(x - other.x, y - other.y, z - other.z);
    }
    Vector3d operator+(Vector3d const &other) const {
        return Vector3d(x + other.x, y + other.y, z + other.z);
    }
    Vector3d operator*(Vector3d const &other) const {
        return Vector3d(x * other.x, y * other.y, z * other.z);
    }
    Vector3d operator/(Vector3d const &other) const {
        return Vector3d(x / other.x, y / other.y, z / other.z);
    }

    Vector3d operator*(double const a) const {
        return Vector3d(a * x, a * y, a * z);
    }
    Vector3d operator/(double const a) const {
        return Vector3d(x / a, y / a, z / a);
    }
    void operator+=(Vector3d const &other) {
        x += other.x;
        y += other.y;
        z += other.z;
    }
    void operator+=(double const a) {
        x += a;
        y += a;
        z += a;
    }
    void operator-=(Vector3d const &other) {
        x -= other.x;
        y -= other.y;
        z -= other.z;
    }
    void operator-=(double const a) {
        x -= a;
        y -= a;
        z -= a;
    }
    void operator*=(double const a) {
        x *= a;
        y *= a;
        z *= a;
    }
    void operator/=(double const a) {
        x /= a;
        y /= a;
        z /= a;
    }

    double &operator[](int const index) {
        if (index == 0)
            return x;
        else if (index == 1)
            return y;
        else if (index == 2)
            return z;
        else {
            printf("FATAL BAD Vector3d INDEX %d\n", index);
            exit(1);
        }
    }

    double operator[](int const index) const {
        if (index == 0)
            return x;
        else if (index == 1)
            return y;
        else if (index == 2)
            return z;
        else {
            printf("FATAL BAD Vector3d INDEX %d\n", index);
            exit(1);
        }
    }

    std::string const str() const {
        return std::string("(" + std::to_string(x) + ", " + std::to_string(y) + ", " +
                           std::to_string(z) + ")");
    };

    double x;
    double y;
    double z;
};

inline Vector3d operator*(double const a, Vector3d const &v) {
    return v * a;
}

#endif /* VECTOR3D_HPP_ */
