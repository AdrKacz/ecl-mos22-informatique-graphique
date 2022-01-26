![Mouvement dans l'espace](./outputs/be2-extra-movement.gif)

# ecl-mos22-informatique-graphique
Travail réalisé pour le MOS 2.2 - Informatique Graphique à l'École Centrale de Lyon. L'objectif est d'implémenter un PathTracer en C++

# Comment exécuter le code ?

1. Installer [`cmake`](https://cmake.org/install/) *(cela devrait déjà être le cas)*

2. Cloner le répertoire

```
git clone https://github.com/AdrKacz/ecl-mos22-informatique-graphique.git
cd ecl-mos22-informatique-graphique
```

3. Créer un dossier `./build` et lancer le *build*

```
mkdir build
cd build
cmake ..
make
```

4. Exécuter le code

```
./InformatiqueGraphique
./InformatiqueGraphique mon-image.png
```

> Pour compiler et exécuter le code : `make && ./InformatiqueGraphique`

# Comment se déplacer

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
