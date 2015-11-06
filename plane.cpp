#include "plane.h"

Plane::Plane(const glm::vec3 &n, const glm::vec3 &p) : normal(n), point(p)
{

}

Plane::Plane(const glm::vec3 &a, const glm::vec3 &b, const glm::vec3 &c) : point(a)
{
    // Vecteurs du plan
    glm::vec3 u(a - c);
    glm::vec3 v(b - c);
    normal = glm::cross(u,v);
}

bool Plane::intersects(const glm::vec3 &o, const glm::vec3 &d, glm::vec3 &out)
{
    // Produit scalaire, si 0 pas de solution ou une infinité
    double scalaire = glm::dot(d,normal);//Vector3D::dotProduct(direction, w); //direction.x()*w.x()+direction.y()*w.y()+direction.z()*w.z();
    if(abs(scalaire) < 0.0001)
        return false;
    double delta = glm::dot(point,normal);//Vector3D::dotProduct(a, w);

    double t = (delta - glm::dot(o,normal))/scalaire; //Vector3D::dotProduct(origine, w)


    out = o+d*(float)t;
    return true;
}

