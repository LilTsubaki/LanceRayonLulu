#ifndef MAILLAGE_H
#define MAILLAGE_H
#include <QVector>
#include <QVector3D>
#include <iostream>
#include <fstream>
#include <time.h>


class Maillage
{
    QVector<QVector3D> geom; //contient les points du maillage
    QVector<int> topo; //contient des triplets d'indice pour chaque triangle du maillage
    QVector<QVector3D> normales; //contient les normales au point
    QVector<int> normalIds;

    public:
        Maillage();
        Maillage(const Maillage& m);
        Maillage(QVector<QVector3D> geo, QVector<int> top, QVector<QVector3D> norm) : geom(geo), topo(top), normales(norm) {}
        Maillage(QVector<QVector3D> geo, QVector<int> top) : geom(geo), topo(top) {}



        // Utilitaires

        Maillage Rotation(const double matrice[3][3]);
        void Merge(Maillage figure2);
        void Ecriture(std::string nom);
        void translate(const QVector3D& t);
        void geometry(const QVector3D &center, const char* obj);

        ~Maillage();
        const QVector<QVector3D>& getGeom() const;
        void setGeom(const QVector<QVector3D> &value);
        const QVector<int>& getTopo() const;
        void setTopo(const QVector<int> &value);
        const QVector<QVector3D>& getNormales() const;
        void setNormales(const QVector<QVector3D> &value);
        const QVector<int>& getNormalIds() const;
        void setNormalIds(const QVector<int> &value);
};

static const int taille = 500;
static const float coeffDetail = 0.03f;


#endif // MAILLAGE_H
