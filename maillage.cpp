#include "maillage.h"
#include <iostream>

const QVector<QVector3D>& Maillage::getGeom() const
{
    return geom;
}

const QVector<int>& Maillage::getTopo() const
{
    return topo;
}

const QVector<QVector3D>& Maillage::getNormales() const
{
    return normales;
}

const QVector<int>& Maillage::getNormalIds() const
{
    return normalIds;
}
Maillage::Maillage()
{
}

Maillage::Maillage(const Maillage &m):geom(m.geom),topo(m.topo),normales(m.normales),normalIds(m.normalIds)
{
    std::cout<<"copie "<<geom.size()<<std::endl;
}

void Maillage::Merge(Maillage figure2)
{
    QVector<int> changement;
    int temp;
    //on parcours la liste des points de la seconde figure
    //on vérifie si il y a des points communs entre les deux figures
    for(int p1 = 0; p1 < figure2.geom.size(); p1++)
    {
//        temp = this->ComparePoint(figure2.geom.at(p1));
//        si le point n'existe pas deja alors on l'ajoute à la geométrie de la premiere figure
//        on stock également le nouvel indice de ce point
//        if(temp == -1)
//        {
            this->geom.append(figure2.geom.at(p1));
            temp = this->geom.size()-1;
        //}
        //std::cout << "l'indice " << p1 << "est maintenant l'indice " << temp << std::endl;
        changement.append(temp);
    }

    //on ajoute a la premiere topo la seconde en mettant a jour les indices
    for(int i = 0; i < figure2.topo.size(); i++)
    {
        this->topo.append(changement.at(figure2.topo.at(i)));
    }
}



void Maillage::Ecriture(std::string nom)
{
    std::ofstream fichier(nom, std::ios::out | std::ios::trunc);
    if(fichier)
    {
        fichier << "o " << nom << std::endl;
        fichier << std::endl;
        QVector3D temp;
        for (int i = 0; i < geom.size(); ++i)
        {
             temp = geom.at(i);
             fichier << "v " << temp.x() << " " << temp.y() << " " << temp.z() << std::endl;
        }
        fichier << std::endl;
        for (int i = 0; i < normales.size(); ++i)
        {
             temp = normales.at(i);
             fichier << "vn " << temp.x() << " " << temp.y() << " " << temp.z() << std::endl;
        }
        fichier << std::endl;
        for (int i = 0; i < topo.size(); i+=3)
        {
             fichier << "f " << topo.at(i)+1 << " " << topo.at(i+1)+1 << " " << topo.at(i+2)+1 << std::endl;
        }
        fichier.close();  // on referme le fichier
    }
    else
    {
        std::cerr << "Erreur à l'ouverture !" << std::endl;
    }
}

void Maillage::translate(const QVector3D &t, glm::vec3 &min, glm::vec3 &max)
{
    for(int i;i<geom.size();++i){
        geom[i]+=t;
    }
    min = glm::vec3(min.x+t.x(), min.y+t.y(),min.z+t.z());
    max = glm::vec3(max.x+t.x(), max.y+t.y(),max.z+t.z());
}

Maillage Maillage::Rotation(const double matrice[3][3])
{
    QVector<QVector3D> geom2;
    QVector3D temp;
    for (int i = 0; i < geom.size(); ++i)
    {
         temp = geom.at(i);
         temp.setX(temp.x() * matrice[0][0] + temp.y() * matrice[0][1] + temp.z() * matrice[0][2]);
         temp.setY(temp.x() * matrice[2][0] + temp.y() * matrice[2][1] + temp.z() * matrice[2][2]);
         temp.setZ(temp.x() * matrice[1][0] + temp.y() * matrice[1][1] + temp.z() * matrice[1][2]);
         geom2.append(temp);
    }

    Maillage * sphere = new Maillage(geom2, topo, normales);
    return *sphere;
}




Maillage::~Maillage()
{

}

void Maillage::geometry(const QVector3D &center, const char* obj, glm::vec3 &min, glm::vec3 &max) {
    QVector3D minVal(1E100, 1E100, 1E100), maxVal(-1E100, -1E100, -1E100);
    FILE* f = fopen(obj, "r");
    while (!feof(f)) {
        char line[255];
        fgets(line, 255, f);
        if (line[0]=='v' && line[1]==' ') {
            QVector3D vec;
            sscanf(line, "v %f %f %f\n", &vec[0], &vec[2], &vec[1]);
            vec[2] = -vec[2];
            QVector3D p = vec*50. + center;
            geom.push_back(p);
            maxVal[0] = std::max(maxVal[0], p[0]);
            maxVal[1] = std::max(maxVal[1], p[1]);
            maxVal[2] = std::max(maxVal[2], p[2]);
            minVal[0] = std::min(minVal[0], p[0]);
            minVal[1] = std::min(minVal[1], p[1]);
            minVal[2] = std::min(minVal[2], p[2]);
            min = glm::vec3(minVal[0], minVal[1], minVal[2]);
            max = glm::vec3(maxVal[0], maxVal[1], maxVal[2]);
        }
        if (line[0]=='v' && line[1]=='n') {
            QVector3D vec;
            sscanf(line, "vn %f %f %f\n", &vec[0], &vec[2], &vec[1]);
            vec[2] = -vec[2];
            normales.push_back(vec);
        }
        if (line[0]=='f') {
            int i0, i1, i2;
            int j0,j1,j2;
            int k0,k1,k2;
            sscanf(line, "f %u/%u/%u %u/%u/%u %u/%u/%u\n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2 );
            topo.push_back(i0-1);
            topo.push_back(i1-1);
            topo.push_back(i2-1);
            normalIds.push_back(k0-1);
            normalIds.push_back(k1-1);
            normalIds.push_back(k2-1);
        }

    }

    /*boundingSphere.C = 0.5*(minVal+maxVal);
    boundingSphere.R = sqrt((maxVal-minVal).sqrNorm())*0.5;*/

    fclose(f);
}
