# Informatique Graphique

> AdrKacz (@AdrKacz), ECL, MOS 2.2, Informatique Graphique, Nicolas Bonneel (@nbonneel)

Réalisation d'un *RayTracer* en **C++** avec les options suivantes:
- Sphère
- Source de lumière ponctuelle
    - Éclairage direct
- *Multi-threading* (réduction du temps de calcul)
- Matériaux opaque avec une couleur
- Matériaux réfléchissants
- Matériaux transparents
- Source de lumière sphérique
    - Éclairage indirect
- *Anti-aliasing*
- Profondeur de champ
- *Mesh*
    - Boîte englobante (réduction du temps de calcul)
    - ~~*Bounding Volume Hierarchy*(réduction du temps de calcul)~~
- Déplacement dans l'espace

### Comment exécuter le code ?

> Le compilateur doit supporter l'option *-fopenmp* pour pouvoir utiliser *PRAGMA*

```sh
make && ./InformatiqueGraphique && make clean
```

# Création d'une sphère

# Création d'une source de lumière ponctuelle

# Réduction du temps de calcul avec le *multi-threading*

# Ajouts des matériaux

## Matériaux opaques

## Matériaux réfléchissants

## Matériaux transparents

# Création d'une source de lumière sphérique

# Réduction du crénelage par *anti-aliasing*

# Paramétrage de la profondeur de champ de la caméra

# Création de *mesh*

## Réduction du temps de calcul avec une boîte englobante

## ~~Réduction du temps de calcul avec un *Bounding Volume Hierarchy*~~

# Déplacement de la caméra dynamique dans l'espace

![Mouvement dans l'espace](./outputs/be2-extra-movement.gif)

## Comment se déplacer

> Séquence utiliser pour la présentation : `zzzzzzddddddzzzzzzzzzqqqqqzzzzzzzzqqqqqsssss`

1. Lancer le programme comme précedement
2. Utiliser `zqsd` pour se déplacer
 - `z` : en avant
 - `s` : en arrière
 - `q` : tourner à gauche
 - `d` : tourner à droite
3. `ENTER` pour valider la suite de mouvements (tu peux entrer plus d'un déplacement à la fois pour une animation plus longue)

Par exemple, si je souhaite faire deux pas en avant puis regarder à gaude de deux crans, je peux entrer : 

```
zzqq
```

Puis `ENTER` (le programme ne s'actualise pas tant que tu ne cliques pas sur `ENTER`, c'est une limitation que je vais essayer de corriger prochainement)


---

> Les informations ci-dessous ne sont plus à jour.

# Notes
- [x] Calcul parallèle avec `<thread>`
    - Mesure de le différence de temps de calcul
- [X] Écrire le code dans plusieurs fichiers au lieu d'avoir un grand fichier

# BE1

## Objectif

- Image d'une sphère
- Source de lumière
- Boîte encapsulant la scène

## [Vector](src/Vector)

## [Ray](src/Ray)

## [Sphere](src/Sphere)

## [Environment](src/Environment)

## Main

```c++
int main(int argc, char* argv[]) {
	int W = 512;
	int H = 512;

    Environment E = Environment();
    E.add_sphere(Sphere(Vector(0, 0, 0), 10));

    E.add_sphere(Sphere(Vector(-10050, 0, 0), 10000));
    E.add_sphere(Sphere(Vector(+10050, 0, 0), 10000));
    E.add_sphere(Sphere(Vector(0, -10050, 0), 10000));
    E.add_sphere(Sphere(Vector(0, +10050, 0), 10000));
    E.add_sphere(Sphere(Vector(0, 0, -10050), 10000));
    E.add_sphere(Sphere(Vector(0, 0, +10050), 10000));

    E.add_light(Vector(30, 30, 25));
    // E.add_light(Vector(-20, -20, -15));
    

    Vector C = Vector(0, 0, 40);
    double alpha = 90 * M_PI / 180;
    double z = - W / (2 * tan(alpha / 2));
	
	std::vector<unsigned char> image(W*H * 3, 0);

    double rho = M_PI;
    double I = 5e6;

    auto start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {    
            Vector u = Vector(j - W / 2 + 0.5, H - i - H / 2 + 0.5, z);
            u.normalize();
            Ray r = Ray(C, u);

            Vector P, N;
            if (E.intersect(r, P, N))
            {
                double intensity = E.get_intensity(N, P, I, rho);
                intensity = clamp(intensity, 0., 255.);
                
                image[(i*W + j) * 3 + 0] = intensity;
			    image[(i*W + j) * 3 + 1] = intensity;
			    image[(i*W + j) * 3 + 2] = intensity;
            }
		}
	}
    auto end = std::chrono::high_resolution_clock::now();
    auto diff_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    std::cout << "Temps pour la création de l'image: " << diff_sec.count() / 1000000.0 << "ms\n";

    std::string filename = "image";
    if (argc > 1) {
        filename = argv[1];
    }
    filename = "outputs/" + filename + ".png";
    stbi_write_png(filename.c_str(), W, H, 3, &image[0], 0);
	

	return 0;
}
```

## Résultat

![BE1-1](./outputs/be1-1.png)

> Temps pour la création de l'image: 233.621ms

## BE 1 - Extra

### Parallélisation

```c++
#include <thread>
#define NB_THREAD_GRID 4 // grid n * n > NB_THREAD = NB_THREAD_GRID * NB_THREAD_GRID
```

```c++
int main(int argc, char* argv[]) {
	// ...

    // Thread
    const int step = std::ceil(512 / NB_THREAD_GRID);
    std::thread threads[NB_THREAD_GRID * NB_THREAD_GRID];

    // ...
    for (int x = 0; x < NB_THREAD_GRID; x++)    
    {
        for (int y = 0; y < NB_THREAD_GRID; y++) {
            threads[NB_THREAD_GRID * x + y] = std::thread([x, y, &step, &z, &E, &W, &H, &C, &rho, &I, &image]()
            {
                // ...
            });
        }
    }
    
    for (int x = 0; x < NB_THREAD_GRID; x++)    
    {
        for (int y = 0; y < NB_THREAD_GRID; y++) {
            threads[x * NB_THREAD_GRID + y].join();
        }
    }
    // ...
}
```

![BE1-Parallel](./outputs/be1-parallel.png)

> Temps pour la création de l'image: 54.0334ms

**Nous avons diviser le temps de calcul par environ `4.5`.**

Pour voir l'exécution des `thread`, nous annulons l'éxécution sur la diagonale.

```c++
if (x == y)
{
    continue;
}
```

![BE1-Parallel-Diag](./outputs/be1-parallel-diag.png)

> Temps pour la création de l'image: 41.6291ms

# BE 2

## Objectif

- Construire des ombres projetées
- Construire des matériaux réfléchissants
- Construire des matériaux transparents

## BE 2 - Extra

### Déplacement dans l'espace

```c++
int main(int argc, char* argv[]) {
    // ...
    std::string mvt_sequence = "";
    bool is_alive = true;
    char c = ' ';
    while (is_alive) {
            refresh(filename, threads, step, z, E, W, H, C, rho, I, I_pow_factor, image);
            std::cin >> c;
            switch (c)
            {
            case 'z':
                C.move_forward(1.);
                break;
            case 's':
                C.move_forward(-1.);
                break;
            case 'a':
                C.move_right(1.);
                break;
            case 'e':
                C.move_right(-1.);
                break;
            case 'q':
                C.rotate(- M_PI / 12);
                break;
            case 'd':
                C.rotate(+ M_PI / 12);
                break;
            default:
                is_alive = false;
                break;
            }
            mvt_sequence += c;
    }
    std::cout << "Movement Sequence:\n" << mvt_sequence;    
	
	return 0;
}
```

### Réflection totale dans les matériaux transparent (TODO)

# BE 3

> BRDF Databse: https://www.merl.com/brdf/
> Global Illumination Compendium: https://people.cs.kuleuven.be/~philip.dutre/GI/


128 rayons par pixels.

![be3-indirect-lighting-128-rays](./outputs/indirect-lighting-128-rays.png)

> Temps pour la création de l'image: 26219.4ms

Anti-aliasing avec Box-Muller, `sigma = 0.5`, `32 rayons`.

![be3-indirect-lighting-128-rays](./outputs/be3-box-muller-sigma05.png)

> Temps pour la création de l'image: 6749.88.4ms

# BE 4 (absent)

![be4-light-r5.png](./outputs/be4-light-r5.png)

> Temps pour la création de l'image: 3922.12ms

![be4-light-r20.png](./outputs/be4-light-r20.png)

> Temps pour la création de l'image: 4106.09ms

![be4-light-smart-r5.png](./outputs/be4-light-smart-r5.png)

> Temps pour la création de l'image: 8059.9ms

![be4-light-smart-r20.png](./outputs/be4-light-smart-r20.png)

> Temps pour la création de l'image: 20706.6ms

![be4-2lights-r5-1e9.png](./outputs/be4-2lights-r5-1e9.png)

> Temps pour la création de l'image: 8310.03ms

![be4-2lights-r5-1e8.png](./outputs/be4-2lights-r5-1e8.png)

> Temps pour la création de l'image: 8540.31ms

*__QUESTION__ : Avec intensity de la lumière a 1e8, les sphéres emissives sont __GRISES__... est-ce dû à la reflection partielle qui ne devrait pas être prise en compte sur une lampe? (voir exemple au dessus)*

*__QUESTION__ : Quel paramètre prendre pour une diffraction correcte ? (essaie fait sans l'indirect lighting : explosion de point)*


# BE 5 - Maillage

![be5-one-triangle.png](./outputs/be5-one-triangle.png)

> Temps pour la création de l'image: 159.27ms