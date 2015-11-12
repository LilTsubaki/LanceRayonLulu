#ifndef PLANE_H
#define PLANE_H
#include <glm/vec3.hpp>
#include <glm/gtc/matrix_transform.hpp>

class Plane
{
private:
    glm::vec3 normal;
    glm::vec3 point;
public:
    Plane(const glm::vec3& n, const glm::vec3& p);
    Plane(const glm::vec3 &a, const glm::vec3& b, const glm::vec3& c);
    bool intersects(const glm::vec3& o, const glm::vec3& d, glm::vec3& out);

};

#endif // PLANE_H
