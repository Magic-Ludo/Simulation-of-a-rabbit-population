/******************************************************************************
 *           ██╗   ██╗██████╗        ██╗       ██╗      ██████╗               *
 *           ██║   ██║██╔══██╗       ██║       ██║     ██╔════╝               *
 *           ██║   ██║██████╔╝    ████████╗    ██║     ██║                    *
 *           ╚██╗ ██╔╝██╔══██╗    ██╔═██╔═╝    ██║     ██║                    *
 *            ╚████╔╝ ██████╔╝    ██████║      ███████╗╚██████╗               *
 *             ╚═══╝  ╚═════╝     ╚═════╝      ╚══════╝ ╚═════╝               *
 *                                                                            *
 *                                                                            *
 *      ██████╗ ██████╗  ██████╗  ██████╗ ██████╗  █████╗ ███╗   ███╗         *
 *      ██╔══██╗██╔══██╗██╔═══██╗██╔════╝ ██╔══██╗██╔══██╗████╗ ████║         *
 *      ██████╔╝██████╔╝██║   ██║██║  ███╗██████╔╝███████║██╔████╔██║         *
 *      ██╔═══╝ ██╔══██╗██║   ██║██║   ██║██╔══██╗██╔══██║██║╚██╔╝██║         *
 *      ██║     ██║  ██║╚██████╔╝╚██████╔╝██║  ██║██║  ██║██║ ╚═╝ ██║         *
 *      ╚═╝     ╚═╝  ╚═╝ ╚═════╝  ╚═════╝ ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝         *
 *                                                                            *   
 *                                                                            *         
 *      Auteur : Boursat Vincent                                              *
 *               Corcos  Ludovic                                              *                
 *                                                                            *
 *      Université Clermont Auvergne | L2 Informatique                        *
 *                                                                            *
 *      Date : 17/03/2020                                                     *
 *                                                                            *
 *      Programme : simu_lapin.c                                              *
 *                                                                            *
 *      Description :                                                         *
 *      Ce programme permet de simuler la croissance d'une population         *
 *      de lapins en fonction d'une certaine probabilité au niveau des        *
 *      naissances, du sexe, de l'âge, de la maturité sexuelle et de          *
 *      la mortalité.                                                         *
 *      Il se compile comme suit :                                            *
 *      gcc -Wall -fopenmp simu_lapin.c -o simu_lapin                         *
 *      Puis :                                                                *
 *      ./simu_lapin                                                          *
 *                                                                            *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <omp.h>

#include "mt19937ar.c"

/* -------------------------------------------------------------------------- */
/*                          Prototypes des fonctions                          */
/* -------------------------------------------------------------------------- */

void AfficheTableau(unsigned long long ***tableau, int nb_annee_simu);

double Uniform(double borne_inf, double borne_sup);

int nbLapinPortee();

int nbPortee();

int SexeLapin();

int MortPetit();

int MortAdulte(double decroissance);

unsigned long long *NaissanceSexuee(unsigned long long ***tableau, int annee);

unsigned long long **Mortalite(unsigned long long ***tableau, unsigned long long *tab_naissances, int annee);

unsigned long long ***Evolution(unsigned long long ***tableau, int nb_annee);

unsigned long long ***AllocationTab3D(int nb_annee_simu, int ligne, int age);

/* -------------------------------------------------------------------------- */
/*                         Fonction 'main' principale                         */
/* -------------------------------------------------------------------------- */

/******************************************************************************
 *                                                                            *
 * Pour expliquer la manière dont est construit l'algorithme, on représente   *
 * un tableau sous cette forme, ce tableau représente une année. On se sert   *
 * d'un tableau en trois dimenssions pour représenter chaque années de la     *
 * même manière.                                                              *
 *                                                                            *
 * ┌───────────────────────┬─────┬─────┬─────┬─────┬─────┬─────┬─────┬─────┐  *
 * │          Age          │  0  │  1  │  2  │  3  │  4  │  5  │ ... │ 15  │  *
 * ├───────────────────────┼─────┼─────┼─────┼─────┼─────┼─────┼─────┼─────┤  *
 * │ Nb de femelles        │   0 │   0 │  25 │   3 │   2 │ ... │   0 │ ... │  *
 * │ Nb de femelles mortes │ ... │ ... │ ... │ ... │ ... │ ... │ ... │ ... │  *
 * │ Nb de mâles           │ ... │ ... │ ... │ ... │ ... │ ... │ ... │ ... │  *
 * │ Nb de mâles morts     │ ... │ ... │ ... │ ... │ ... │ ... │ ... │ ... │  *
 * └───────────────────────┴─────┴─────┴─────┴─────┴─────┴─────┴─────┴─────┘  *
 *                                                                            *
 ******************************************************************************/

int main(int argc, char *argv[])
{

    int i;
    int nombre_annee_simu = 28;
    unsigned long long ***matrix_result, ***result_fin;

    printf("Nombre d’arguments passes au programme : %d\n", argc);
    for (i = 0; i < argc; i++)
    {
        printf(" argv[%d] : '%s'\n", i, argv[i]);
    }

    //  On alloue de l'espace en mémoire pour les tableaux matrix_result et result_fin
    //  Ces deux tableaux ont la même représentation que indiqué dans les commentaires
    //  ci-dessus.

    matrix_result = AllocationTab3D(nombre_annee_simu, 4, 16);
    result_fin = AllocationTab3D(nombre_annee_simu, 4, 16);

    //  On initialise ici la tableau matrix_result avec des valeurs pour les premiers
    //  lapins. Ici en l'occurrence, on initialise avec 10 lapins mâles et femelles
    //  qui ont respectivement 4 ans.
    matrix_result[0][0][10] = 10;
    matrix_result[0][2][10] = 10;

    //  On met dans ce tableau les résultats de la simulation calculé sur le nombre
    //  d'année pris en deuxième paramètre de la fonction Evolution.
    result_fin = Evolution(matrix_result, 27);

    //  On affiche maintenant le tableau pour visualiser les résultats.
    AfficheTableau(result_fin, nombre_annee_simu);

    //desallocation(matrix_result, 20, 4);
    //desallocation(result_fin, 20, 4);

    return EXIT_SUCCESS;
}

/* -------------------------------------------------------------------------- */
/*                       Fonctions servant au programme                       */
/* -------------------------------------------------------------------------- */

/******************************************************************************
 *                                                                            *
 * Fonction : int ***Evolution (int ***tableau, int nb_annee)                 *
 *                                                                            *
 * Permet de calculer le nombre de simulation correspondant à l'année entrée  *
 * sur une population de lapins initialisé avant son appel.                   *
 *                                                                            *
 * En entrée : Un tableau en 3 dimensions avec juste les pemières valeurs     *
 *             initialisées.                                                  *
 *             Le nombre d'années sur lequel l'algorithme doit simuler la     *
 *             population de lapins.                                          *
 *                                                                            *
 * En sortie : Un tableau complet en 3 dimensions contenant les résultats     *
 *             générés.                                                       *
 *                                                                            *
 * Le tableau naissance est représenté de la manière suivante :               *
 * ┌────────────────┬─────────────┐                                           * 
 * │ Nb bb femelles │ Nb bb mâles │                                           * 
 * └────────────────┴─────────────┘                                           *           
 *                                                                            *
 * Le tableau mort est représenté de la manière suivante :                    *
 * ┌────────────────────────────┬─────┬───┬───┬───┬─────┬────┐                *
 * │                            │ Age │   │   │   │     │    │                *
 * ├────────────────────────────┼─────┼───┼───┼───┼─────┼────┤                *
 * │ Nombre de femelles mortes  │  0  │ 1 │ 2 │ 3 │ ... │ 16 │                *
 * │ Nombre de mâles mort       │  0  │ 1 │ 2 │ 3 │ ... │ 16 │                *
 * └────────────────────────────┴─────┴───┴───┴───┴─────┴────┘                *
 *                                                                            *
 ******************************************************************************/

unsigned long long ***Evolution(unsigned long long ***tableau, int nb_annee)
{

    int i, annee;
    unsigned long long *naissance;
    unsigned long long **mort;

    for (annee = 1; annee < nb_annee; annee++)
    {

        //  On rempli ici le tableau des naissances avec le nombre de bébé
        //  lapins mâles et femelles obtenue durant l'année précédente.
        naissance = NaissanceSexuee(tableau, annee - 1);

        //  On rempli ici le tableau des morts avec, le nombre de lapins mort en
        //  fonction de leur âge que l'on obtient à la fin de l'année précédente
        //  en prenant en considération le nombre de naissances.
        mort = Mortalite(tableau, naissance, annee - 1);

        printf("\rAnnées simulées : %d sur %d", annee, nb_annee);
        fflush(stdout);

        tableau[annee - 1][0][0] = naissance[0];
        tableau[annee - 1][2][0] = naissance[1];

#pragma omp parallel for
        for (i = 0; i < 16; i++)
        {
            tableau[annee - 1][1][i] = mort[0][i];
            tableau[annee - 1][3][i] = mort[1][i];
        }

        //  On calcul le nombre de lapins de l'année n - 1 à l'année n :
        //  On remplie le tableau de l'année en cours avec le nombre de lapins qui on
        //  survécue à l'année précedente en les viellisant d'un an.
        //#pragma omp parallel for
        for (i = 1; i < 16; i++)
        {
            tableau[annee][0][i] = tableau[annee - 1][0][i - 1] - tableau[annee - 1][1][i - 1];
            tableau[annee][2][i] = tableau[annee - 1][2][i - 1] - tableau[annee - 1][3][i - 1];
        }
    }

    return tableau;
}

/******************************************************************************
 *                                                                            *
 * Fonction : int **Mortalite (int ***tableau, int *tab_naissances, int annee)*
 *                                                                            *
 * Permet de calculer la mortalité des lapins en fonctions de leur âge et de  *
 * leur sexe.                                                                 *
 *                                                                            *
 * En entrée : Un tableau en 3 dimensions initialisé au fur et à mesure       *
 *             grâce à la fonction Evolution.                                 *
 *             Un tableau correspondant au nombre de naissances (mâles et     *
 *             femelles).                                                     *
 *             L'année sur laquelle ont veut calculer la mortalité.           *
 *                                                                            *
 * En sortie : Un tableau de mortalité contenant les résultats générés.       *
 *                                                                            *
 * Le tableau naissance est représenté de la manière suivante :               *
 * ┌────────────────┬─────────────┐                                           * 
 * │ Nb bb femelles │ Nb bb mâles │                                           * 
 * └────────────────┴─────────────┘                                           *           
 *                                                                            *
 * Le tableau mort est représenté de la manière suivante :                    *
 * ┌────────────────────────────┬─────┬───┬───┬───┬─────┬────┐                *
 * │                            │ Age │   │   │   │     │    │                *
 * ├────────────────────────────┼─────┼───┼───┼───┼─────┼────┤                *
 * │ Nombre de femelles mortes  │  0  │ 1 │ 2 │ 3 │ ... │ 16 │                *
 * │ Nombre de mâles mort       │  0  │ 1 │ 2 │ 3 │ ... │ 16 │                *
 * └────────────────────────────┴─────┴───┴───┴───┴─────┴────┘                *
 *                                                                            *
 ******************************************************************************/

unsigned long long **Mortalite(unsigned long long ***tableau, unsigned long long *tab_naissances, int annee)
{

    int i, j, k;
    double decroissance = 0;
    unsigned long long **tab_mort = (unsigned long long **)calloc(2, sizeof(unsigned long long *));

    //  Allocation dynamique d'un tableau en 2 dimensions.
    for (i = 0; i < 2; i++)
    {
        tab_mort[i] = (unsigned long long *)calloc(16, sizeof(unsigned long long));
    }

    //  Remplissage du tableau mort avec le nombre de bébé lapins morts
    //  mâles et femelles générés.
    for (i = 0; i < 2; i++)
    {

        for (j = 0; j < tab_naissances[i]; j++)
        {
            tab_mort[i][0] += MortPetit();
        }
    }

    //  Remplissage du tableau mort avec le nombre de lapins adultes morts
    //  mâles et femelles générés, sauf qu'à partir de 10 ans, leurs chances
    //  de survie diminue de 10 % chaque année.

    for (i = 0; i < 2; i++)
    {
        decroissance = 0;
        for (j = 1; j < 16; j++)
        {

            if (j >= 10)
            {
                decroissance += 0.1;
            }

            for (k = 0; k < tableau[annee][2 * i][j]; k++)
            {
                tab_mort[i][j] += MortAdulte(decroissance);
            }
        }
    }

    return tab_mort;
}

/******************************************************************************
 *                                                                            *
 * Fonction : int MortPetit()                                                 *
 *                                                                            *
 * Sert à calculer la mortalité des bébés. Ils ont 12 % de chance de survie.  *
 *                                                                            *
 * En entrée : Rien.                                                          *
 *                                                                            *
 * En sortie : 1 si le bébé lapin est mort                                    *
 *             0 sinon.                                                       *
 *                                                                            *
 ******************************************************************************/

int MortPetit()
{

    double val_aleatoire = genrand_real1();
    int val_retour = 0;
    if (val_aleatoire >= 0.12)
    {
        val_retour++;
    }

    return val_retour;
}

/******************************************************************************
 *                                                                            *
 * Fonction : int MortAdulte (double decroissance)                            *
 *                                                                            *
 * Sert à calculer la mortalité des lapins selon leurs âge.                   *
 * Ils ont 60 % de chance de survie de 1 à 10 ans, mais à partir de 10 ans,   *
 * leurs chances de survie diminue de 10 % tous les ans.                      *
 *                                                                            *
 * En entrée : La décroissance, elle est modifié dans la fonction Evolution   *
 *             lorsque que le lapin est âgé de plus de 10 ans.                *
 *                                                                            *
 * En sortie : 1 si le  lapin est mort                                        *
 *             0 sinon.                                                       *
 *                                                                            *
 ******************************************************************************/

int MortAdulte(double decroissance)
{

    double val_aleatoire = genrand_real1();
    int val_retour = 0;

    if (val_aleatoire >= (0.60 - decroissance))
    {
        val_retour++;
    }

    return val_retour;
}

/******************************************************************************
 *                                                                            *
 * Fonction : int *NaissanceSexuee (int ***tableau, int annee)                *
 *                                                                            *
 * Permet de calculer le nombre de bébés lapins mâles et femelles en fonction *
 * du nombre de portées et du nombre de lapins par portées.                   *
 *                                                                            *
 * En entrée : Un tableau en 3 dimensions initialisé au fur et à mesure       *
 *             grâce à la fonction Evolution.                                 *   
 *             L'année sur laquelle ont veut calculer le nombre de naissances *
 *             ainsi que le sexe des nouveaux lapins.                         *
 *                                                                            *
 * En sortie : Un tableau de naissance contenant les résultats générés.       *
 *                                                                            *
 * Le tableau naissance est représenté de la manière suivante :               *
 * ┌────────────────┬─────────────┐                                           * 
 * │ Nb bb femelles │ Nb bb mâles │                                           * 
 * └────────────────┴─────────────┘                                           *           
 *                                                                            *
 * Il y a 50 % de chances d'obtenir un mâle ou une femelle.                   *
 * Il y a 4 à 8 portées chaques année par lapines adulte, mais il est plus    *
 * probable d'en avoir 5 à 7.                                                 *
 * Il y a une équiprobabilité d'obtenir entre 3 et 6 lapins par portée.       *
 *                                                                            *
 ******************************************************************************/

unsigned long long *NaissanceSexuee(unsigned long long ***tableau, int annee)
{

    int i, j, k;
    unsigned long long nb_femelles_mature = 0,
                       nb_bb_males = 0,
                       nb_bb_femelles = 0,
                       nb_bb_portee = 0,
                       nb_bb_tot = 0;

    //  tab_result est le tableau où seront stocké les informations des
    //  naissances. C'est pour celà que l'on lui alloue de la mémoire ici.
    unsigned long long *tab_result = (unsigned long long *)malloc(2 * sizeof(unsigned long long));

    for (k = 1; k <= 15; k++)
    {
        nb_femelles_mature += tableau[annee][0][k];
    }

    unsigned long long *tab_portee = (unsigned long long *)malloc(nb_femelles_mature * sizeof(unsigned long long));
    unsigned long long *tab_naissance = (unsigned long long *)malloc(nb_femelles_mature * sizeof(unsigned long long));

    //  On défini ici le nombres de mâles et de femelles créé pour chaque
    //  femelles mature (âge supérieur à 1 an) et pour chaque portées qu'elles
    //  donneront.

    for (i = 0; i < nb_femelles_mature; i++)
    {

        tab_portee[i] = nbPortee();

        for (j = 0; j < tab_portee[i]; j++)
        {

            nb_bb_portee = nbLapinPortee();

            for (k = 0; k < nb_bb_portee; k++)
            {

                tab_naissance[k] = SexeLapin();

                if (tab_naissance[k] == 1)
                {
                    nb_bb_males++;
                }
                else
                {
                    nb_bb_femelles++;
                }
            }
            nb_bb_tot += nb_bb_portee;
        }
    }

    //  On rempli le tableau
    tab_result[0] = nb_bb_femelles;
    tab_result[1] = nb_bb_males;

    return tab_result;
}

/******************************************************************************
 *                                                                            *
 * Fonction : int nbPortee()                                                  *
 *                                                                            *
 * Sert à calculer le nombre de portées total par lapine sur une année, il y  *
 * a environ 4 à 8 portées par an, mais il y a plus de chance d'en obtenir    *
 * entre 5 et 7.                                                              *
 *                                                                            *
 * En entrée : Rien.                                                          *
 *                                                                            *
 * En sortie : Le nombre de portée.                                           *
 *                                                                            *
 * La répartition du nombre de portées se fait comme suit :                   *
 *                                                                            *
 *                          ██████  4 - 10 %                                  *
 *                          ██████████  5 - 20 %                              *
 *                          ██████████████  6 - 40 %                          *
 *                          ██████████  7 - 20 %                              *
 *                          ██████  8 - 10 %                                  *
 *                                                                            *
 ******************************************************************************/

int nbPortee()
{

    int i;
    double valGene = genrand_real1();
    double pourcentage[5] = {0.1, 0.3, 0.7, 0.9, 1.0};

    for (i = 0; i <= 4; i++)
    {

        if (valGene <= pourcentage[i])
        {

            return (4 + i);
        }
    }

    return EXIT_SUCCESS;
}

/******************************************************************************
 *                                                                            *
 * Fonction : void AfficheTableau (int ***tableau, int nb_annee_simu)         *
 *                                                                            *
 * Permet simplement d'afficher un tableau en 3 dimenssions.                  *
 *                                                                            *
 * En entrée : Un tableau à 3 dimensions                                      *
 *             Le nombre d'années sur lesquelles ont doit afficher le tableau *
 *                                                                            *
 * En sortie : Rien, cette fonction ne fait que de l'affichage.               *
 *                                                                            *
 ******************************************************************************/

void AfficheTableau(unsigned long long ***tableau, int nb_annee_simu)
{

    int i, j, k;

    for (i = 0; i < nb_annee_simu; i++)
    {

        printf("Année %d\n", i);

        for (j = 0; j <= 3; j++)
        {

            for (k = 0; k <= 15; k++)
            {

                printf("%11lld\t", tableau[i][j][k]);
            }
            printf("\n");
        }

        printf("\n\n");
    }
}

/******************************************************************************
 *                                                                            *
 * Fonction : int SexeLapin()                                                 *
 *                                                                            *
 * Permet de déterminer si un lapin est un mâle ou une femelle, il y a 50 %   *
 * de chance que se soit l'un ou l'autre.                                     *
 *                                                                            *
 * En entrée : Rien.                                                          *
 *                                                                            *
 * En sortie : 1 si le bébé lapin est mâle                                    *
 *             0 si c'est une femelle                                         *
 *                                                                            *
 ******************************************************************************/

int SexeLapin()
{

    double val = genrand_real1();
    if (val <= 0.5)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

/******************************************************************************
 *                                                                            *
 * Fonction : int nbLapinPortee()                                             *
 *                                                                            *
 * Permet de calculer le nombre de lapin par portées. Il y a environ 4 à 8    *
 * portées par an, mais il y a plus de chance d'en obtenir entre 5 et 7.      *
 *                                                                            *
 * En entrée : Rien.                                                          *
 *                                                                            *
 * En sortie : Le nombre de lapins par portées                                *
 *                                                                            *
 * La répartition du nombre de lapins se fait comme suit :                    *
 *                                                                            *
 *                              ██████  3 - 25 %                              *
 *                              ██████  4 - 25 %                              *
 *                              ██████  5 - 25 %                              *
 *                              ██████  6 - 25 %                              *
 *                                                                            *
 ******************************************************************************/

int nbLapinPortee()
{

    int x = (int)(Uniform(2.0, 6.0) + 1);

    return x;
}

/******************************************************************************
 *                                                                            *
 * Fonction : double Uniform (double borne_inf, double borne_sup)             *
 *                                                                            *
 * Permet de générer aléatoirement un nombre de type double compris entre     *
 * borne_inf et borne_sup.                                                    *
 *                                                                            *
 * En entrée : Une bonre inférieur : borne_inf                                *
 *             Une borne supérieur : borne_sup                                *
 *                                                                            *
 * En sortie : Le nombre compris entre ces bornes.                            *
 *                                                                            *
 ******************************************************************************/

double Uniform(double borne_inf, double borne_sup)
{

    return (borne_inf + (borne_sup - borne_inf) * genrand_real1());
}

/******************************************************************************
 *                                                                            *
 * Fonction : int ***AllocationTab3D(int nb_annee_simu, int ligne, int age)   *
 *                                                                            *
 * Permet d'allouer de manière dynamique la mémoire pour un tableau de 3      *
 * dimensions.                                                                *
 *                                                                            *
 * En entrée : Le nombre d'années à simuler                                   *
 *             Le nombre de ligne du tableau                                  *
 *             L'âge maximum des lapins + 1                                   *
 *                                                                            *
 * En sortie : Un tableau alloué prêt à être utilisé                          *
 *                                                                            *
 * Le tableau est de ce type : [Année][Ligne][Âge]                            *
 *                                                                            *
 ******************************************************************************/

unsigned long long ***AllocationTab3D(int nb_annee_simu, int ligne, int age)
{
    unsigned long long ***matrix_result;
    int i, j;

    matrix_result = malloc(nb_annee_simu * sizeof(unsigned long long **));
    matrix_result[0] = malloc(nb_annee_simu * ligne * sizeof(unsigned long long *));
    matrix_result[0][0] = malloc(nb_annee_simu * ligne * age * sizeof(unsigned long long));

    for (i = 0; i < nb_annee_simu; i++)
    {
        matrix_result[i] = matrix_result[0] + i * ligne;
        for (j = 0; j < ligne; j++)
        {
            matrix_result[i][j] = matrix_result[0][0] + i * ligne * age + j * age;
        }
    }
    return matrix_result;
}
