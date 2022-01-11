# ecl-mos22-informatique-graphique
Travail réalisé pour le MOS 2.2 - Informatique Graphique à l'École Centrale de Lyon. L'objectif est d'implémenter un PathTracer en C++

# BE1

## Objectif

- Image d'une sphère
- Source de lumière
- Boîte encapsulant la scène

## Vector

La classe `Vector` permet d'effectuer des opérations classiques avec des vecteurs en dimensions 3.

```c++
class Vector
{
private:
    double vec[3];
public:
    Vector();
    Vector(double x, double y, double z);
    ~Vector();

    ...
};
```

### Opérations simples

Les opérations suivantes permettent de récupérer les valeurs du vecteurs directement à partir de l'object.

```c++
double operator[](int i) const { return vec[i]; };
double& operator[](int i) {return vec[i];};
```

Les opérations suivantes permettent les opérations `+`, `-`, `*`, et `/`.
```c++
Vector operator+(const Vector&) const;
Vector operator-(const Vector&) const;
Vector operator*(const double) const;
Vector operator/(const double) const;
```

## Résultat

![outputs](./outputs/be1-1.png)

## Améliorations

- [ ] Calcul parallèle avec `<thread>`
    - Mesure de le différence de temps de calcul
- [ ] Écrire le code dans plusieurs fichiers au lieu d'avoir un grand fichier


