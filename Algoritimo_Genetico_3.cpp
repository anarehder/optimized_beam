#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>

#include <array>
#include <deque>
#include <fstream>
#include <iostream>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <list>
#include <cmath>
#include <unistd.h>
#include <random>
#include <set>
#include <vector>

#include <lsAdvect.hpp>
#include <lsBooleanOperation.hpp>
#include <lsDomain.hpp>
#include <lsExpand.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsGeometries.hpp>
#include <lsMakeGeometry.hpp>
#include <lsPrune.hpp>
#include <lsSmartPointer.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsWriteVisualizationMesh.hpp>
#include <lsWriter.hpp>

#include "fourRateInterpolation.hpp"

using namespace std;

const int ARRAY_SIZE = 16;
const int POPULATION_SIZE = 12;
const int MAX_GENERATIONS = 1000;
const double stopFitness = 0.1;
const int SELECTED = 4;
double bestAllFitness;
bool parada = false;

constexpr int D = 3;
using NumericType = double;
constexpr NumericType gridDelta = 0.5; // 0.2; //0.51; //0.51;//0.27;//0.92;
unsigned outputNum = 0;
using LevelSetType = lsSmartPointer<lsDomain<NumericType, D>>;

constexpr double PI = 3.1415296;

// This function creates a box starting in minCorner spanning to maxCorner
using T = NumericType;
lsSmartPointer<lsMesh<NumericType>> makeBox(hrleCoordType *minCorner, hrleCoordType *maxCorner, double angle)
{
    // draw all triangles for the surface and then import from the mesh
    std::vector<std::array<T, 3>> corners;
    corners.resize(std::pow(2, D), {0, 0, 0});

    // first corner is the minCorner
    for (unsigned i = 0; i < D; ++i)
        corners[0][i] = minCorner[i];

    // last corner is maxCorner
    for (unsigned i = 0; i < D; ++i)
        corners.back()[i] = maxCorner[i];

    // calculate all missing corners
    corners[1] = corners[0];
    corners[1][0] = corners.back()[0];

    corners[2] = corners[0];
    corners[2][1] = corners.back()[1];

    if (D == 3)
    {
        corners[3] = corners.back();
        corners[3][2] = corners[0][2];

        corners[4] = corners[0];
        corners[4][2] = corners.back()[2];

        corners[5] = corners.back();
        corners[5][1] = corners[0][1];

        corners[6] = corners.back();
        corners[6][0] = corners[0][0];
    }

    // now add all corners to mesh
    auto mesh = lsSmartPointer<lsMesh<T>>::New();
    for (unsigned i = 0; i < corners.size(); ++i)
    {
        mesh->insertNextNode(corners[i]);
    }

    if (D == 2)
    {
        std::array<unsigned, 2> lines[4] = {{0, 2}, {2, 3}, {3, 1}, {1, 0}};
        for (unsigned i = 0; i < 4; ++i)
            mesh->insertNextLine(lines[i]);
    }
    else
    {
        std::array<unsigned, 3> triangles[12] = {
            {0, 3, 1}, {0, 2, 3}, {0, 1, 5}, {0, 5, 4}, {0, 4, 2}, {4, 6, 2}, {7, 6, 4}, {7, 4, 5}, {7, 2, 6}, {7, 3, 2}, {1, 3, 5}, {3, 7, 5}};
        for (unsigned i = 0; i < 12; ++i)
            mesh->insertNextTriangle(triangles[i]);
    }

    std::array<double, D> axis;
    axis[D - 1] = 1.;
    lsTransformMesh(mesh, lsTransformEnum::ROTATION, axis, angle).apply();

    return mesh;
}

void writeSurface(LevelSetType LS, std::string fileName)
{
    auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();
    lsToSurfaceMesh<NumericType, D>(LS, mesh).apply();
    lsVTKWriter<NumericType>(mesh, fileName).apply();
}

void writeLS(LevelSetType LS, std::string fileName)
{
    auto mesh = lsSmartPointer<lsMesh<NumericType>>::New();
    lsToMesh<NumericType, D>(LS, mesh).apply();
    lsVTKWriter<NumericType>(mesh, fileName).apply();
}

class WetEtch : public lsVelocityField<NumericType>
{
    static constexpr std::array<NumericType, D> direction100{0.707106781187, 0.707106781187, 0};
    // static constexpr std::array<NumericType, D> direction100{1,0.,0};
    static constexpr std::array<NumericType, D> direction010{-0.707106781187, 0.707106781187, 0};

    NumericType r100 = 0.013283; // 0.0166666666667;
    NumericType r110 = 0.024167; // 0.0309166666667;
    NumericType r111 = 0.000083; // 0.000121666666667;
    NumericType r311 = 0.023933; // 0.0300166666667;

    std::vector<double> velocities;

public:
    WetEtch(std::vector<double> vel) : velocities(vel) {}

    double getScalarVelocity(const std::array<NumericType, 3> & /*coordinate*/,
                             int material,
                             const std::array<NumericType, 3> &normal,
                             unsigned long /* pointID */) final
    {
        if (std::abs(velocities[material]) < 1e-3)
        {
            return 0;
        }

        return velocities[material] * fourRateInterpolation<double, 3>(normal, direction100, direction010, r100, r110, r111, r311);
    }
};

// Função para gerar um array aleatório de tamanho ARRAY_SIZE
std::vector<int> generateRandomArray()
{
    std::vector<int> array;
    for (int i = 0; i < ARRAY_SIZE; i++)
    {
        int number = rand() % 2;
        array.push_back(number);
        std::cout << array[i] << ' ';
    }
    std::cout << ' ' << std::endl;
    return array;
}

// Função para calcular a aptidão de um indivíduo (quanto menor, melhor)
double calculateFitness(std::vector<int> positions)
{
    std::cout << "Fitness: " << std::endl;
    omp_set_num_threads(24);

    LevelSetType substrate;
    hrleCoordType extent = 15;
    int dimensao = 4;

    {
        hrleCoordType bounds[2 * D] = {-extent, extent, -extent, extent};
        lsDomain<NumericType, D>::BoundaryType boundaryCons[D];
        boundaryCons[0] = lsDomain<NumericType, D>::BoundaryType::PERIODIC_BOUNDARY;
        boundaryCons[1] = lsDomain<NumericType, D>::BoundaryType::PERIODIC_BOUNDARY;
        boundaryCons[2] = lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;
        substrate = LevelSetType::New(bounds, boundaryCons, gridDelta);
        NumericType origin[D] = {};
        origin[D - 1] = 0.01;
        NumericType planeNormal[D] = {};
        planeNormal[D - 1] = 1.;
        lsMakeGeometry(substrate, lsSmartPointer<lsPlane<NumericType, D>>::New(origin, planeNormal)).apply();
    }

    // make mask
    auto mask = LevelSetType::New(substrate->getGrid());
    double maskAngle = 0.;
    {
        NumericType origin[D] = {};
        origin[D - 1] = 2.;
        NumericType planeNormal[D] = {};
        planeNormal[D - 1] = 1.;
        auto plane = LevelSetType::New(mask->getGrid());
        lsMakeGeometry(plane, lsSmartPointer<lsPlane<NumericType, D>>::New(origin, planeNormal)).apply();

        origin[D - 1] = 0.;
        planeNormal[D - 1] = -1.;
        lsMakeGeometry(mask, lsSmartPointer<lsPlane<NumericType, D>>::New(origin, planeNormal)).apply();

        lsBooleanOperation(mask, plane, lsBooleanOperationEnum::INTERSECT).apply();

        // FAZENDO UMA MÁSCARA COM O TAMANHO DE S (S.SIZE)
        for (int i = 0; i < positions.size(); ++i)
        {
            auto maskHole = LevelSetType::New(mask->getGrid());
            if (positions[i] == 1)
            {
                int block = i + 1;
                hrleCoordType block_double = static_cast<hrleCoordType>(block);
                hrleCoordType dimensao_double = static_cast<hrleCoordType>(dimensao);
                hrleCoordType calculoX = static_cast<hrleCoordType>(std::ceil(block_double / dimensao_double));
                hrleCoordType calculoY = static_cast<hrleCoordType>(block % dimensao);
                hrleCoordType coordX = -extent + (calculoX * 5);
                hrleCoordType coordY = -extent + (calculoY * 5);
                if (calculoY == 0)
                {
                    coordY = (extent - 10);
                }

                hrleCoordType minCorner[D] = {coordX, coordY, -1.};
                hrleCoordType maxCorner[D] = {coordX + 5, coordY + 5, 3.};

                auto mesh = makeBox(minCorner, maxCorner, maskAngle);
                lsFromSurfaceMesh(maskHole, mesh).apply();
                lsBooleanOperation(mask, maskHole, lsBooleanOperationEnum::RELATIVE_COMPLEMENT).apply();
                writeSurface(mask, "MASK.vtp");
            }
        }
    }

    lsBooleanOperation(substrate, mask, lsBooleanOperationEnum::UNION).apply();

    std::vector<double> epiVels = {0, 16};
    auto velocity = lsSmartPointer<WetEtch>::New(epiVels);

    std::vector<LevelSetType> LSs;
    LSs.push_back(mask);
    LSs.push_back(substrate);

    lsAdvect<NumericType, D> advection(LSs, velocity);
    advection.setSingleStep(true);
    advection.setSaveAdvectionVelocities(true);
    advection.setIntegrationScheme(lsIntegrationSchemeEnum::STENCIL_LOCAL_LAX_FRIEDRICHS_1ST_ORDER);
    double time = 15; // 15; //51.; //41;
    double nextOutput = 40.5;
    unsigned counter = 1;
    while (time > 0.01)
    {
        advection.setAdvectionTime(time);
        time -= advection.getAdvectedTime();
        if (time < nextOutput)
        {
            nextOutput -= 10.;
            ++counter;
        }
        advection.apply();
    }

    NumericType origin[D] = {};
    origin[D - 1] = -1.;
    NumericType planeNormal[D] = {};
    planeNormal[D - 1] = 1.;
    auto plane1 = LevelSetType::New(mask->getGrid());
    lsMakeGeometry(plane1, lsSmartPointer<lsPlane<NumericType, D>>::New(origin, planeNormal)).apply();

    lsBooleanOperation(substrate, plane1, lsBooleanOperationEnum::INTERSECT).apply();
    writeSurface(substrate, "ARQ_FINAL.vtp");

    lsWriteVisualizationMesh<NumericType, D> visMesh;
    visMesh.insertNextLevelSet(mask);
    visMesh.insertNextLevelSet(substrate);
    if (maskAngle < 0.01)
    {
        visMesh.setFileName("aligned");
    }
    else
    {
        visMesh.setFileName("misaligned");
    }
    visMesh.apply();

    int res0;
    // ALTERANDO O NOME DO ARQUIVO GERADO
    res0 = system("mv aligned_volume.vtu DISP_TESTE.vtu");

    // ALTERANDO PARA STL
    int res1;
    res1 = system("jupyter nbconvert --to notebook --execute makeSTL.ipynb");

    // FAZENDO A COMPARAÇÃO
    int res2;
    res2 = system("openscad -o DISP_SUBTRAIDO.stl DIFFERENCE.scad");

    // CALCULANDO O VOLUME
    int res3;
    res3 = system("jupyter nbconvert --to notebook --execute volume.ipynb");

    // OBTENDO O VOLUME CALCULADO
    // FAZER O CALCULO IGUAL O LER ARQ
    char texto[100];
    char *result;
    int ip;

    // para declarar um arquivo
    FILE *fp;

    // MODELO - fp = fopen(const char filename,const char mode);
    fp = fopen("VOLUME.txt", "r");

    if (fp == NULL) // Se houve erro na abertura
    {
        printf("Problemas na abertura do arquivo\n");
    }

    while (!feof(fp)) /* Enquanto não se chegar no final do arquivo */
    {
        result = fgets(texto, 100, fp); // o 'fgets' lê até 99 caracteres ou até o '\n'
    }

    double VOLUME = atof(texto);
    printf("Volume é %f\n", VOLUME);
    fclose(fp);

    int res4 = system("mv DISP_TESTE.vtu teste.vtu");
    int res5 = system("mv DISP_TESTE.stl teste.stl");
    int res6 = system("mv DISP_SUBTRAIDO.stl subtraido.stl");
    int res7 = system("mv VOLUME.txt volume_anterior.txt");
    return VOLUME;
}

// Função para selecionar um indivíduo da população com base em suas aptidões
std::vector<int> selection(const std::vector<double> &allFitness)
{
    // std::vector<int> bestFitness = allFitness;
    //  recebo a array allFitness e escolho quais serao adicionados ao array bestFitness (os 4 menores valores)
    // std::partial_sort(bestFitness.begin(), bestFitness.begin() + 4, bestFitness.end());
    std::cout << "Seleção: " << std::endl;
    std::vector<int> indices(allFitness.size());
    for (size_t i = 0; i < indices.size(); ++i)
    {
        indices[i] = i;
    }

    std::partial_sort(indices.begin(), indices.begin() + SELECTED, indices.end(),
                      [&allFitness](int a, int b)
                      { return allFitness[a] < allFitness[b]; });

    return std::vector<int>(indices.begin(), indices.begin() + SELECTED);
}

// Função para realizar o crossover entre dois indivíduos
std::vector<int> crossover1(const std::vector<int> &parent1, const std::vector<int> &parent2)
{
    std::cout << "Cruzamento: " << std::endl;
    std::vector<int> child;
    // posso optar por quebrar sempre ao meio ou deixar aleatório
    int crossoverPoint = rand() % ARRAY_SIZE;
    // olho ponto a ponto, se estiver do lado menor que i pego de um pai1 se estiver no maios pego do pai2
    for (int i = 0; i < ARRAY_SIZE; i++)
    {
        if (i < crossoverPoint)
        {
            child.push_back(parent1[i]);
        }
        else
        {
            child.push_back(parent2[i]);
        }
    }
    return child;
}

std::vector<int> crossover2(const std::vector<int> &parent1, const std::vector<int> &parent2)
{
    std::cout << "Cruzamento: " << std::endl;
    std::vector<int> child;
    // posso optar por quebrar sempre ao meio ou deixar aleatório
    int crossoverPoint = rand() % ARRAY_SIZE;
    // olho ponto a ponto, quero pegar ametade final do pai 2 e a metade inicial do pai 1 em ordem trocada
    for (int i = 0; i < ARRAY_SIZE; i++)
    {   
        int ponto = ARRAY_SIZE - 1 - i;
        if (i < crossoverPoint)
        {
            child.push_back(parent2[ponto]);
        }
        else
        {
            child.push_back(parent1[ponto]);
        }
    }
    return child;
}

// Função para realizar uma mutação em um indivíduo - altera de 0 para 1
void mutate(std::vector<int> &individual)
{
    std::cout << "Mutação: " << std::endl;
    int max = ARRAY_SIZE - 1;
    int mutationPoint = rand() % ARRAY_SIZE;
    if (individual[mutationPoint] == 1)
    {
        individual[mutationPoint] = 0;
    }
    if (individual[mutationPoint] == 0)
    {
        individual[mutationPoint] = 1;
    }
}

int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    // Inicialização do gerador de números aleatórios
    srand(time(NULL));

    // Criação da população inicial Array com 8 vetores de 16 pontos
    std::vector<std::vector<int>> population;
    for (int i = 0; i < POPULATION_SIZE; i++)
    {
        population.push_back(generateRandomArray());
    }

    // Execução das gerações
    for (int generation = 0; generation < MAX_GENERATIONS; generation++)
    {
        // Avaliação da população atual
        std::vector<double> allFitness;

        for (int i = 0; i < POPULATION_SIZE; i++)
        {
            double fitness = calculateFitness(population[i]);
            allFitness.push_back(fitness); // o array all fitness tem os mesmos indices do array population
            // ou seja, 1o individuo é population[0] e allFitness[0]
            if (fitness < stopFitness)
            {
                bestAllFitness = fitness;
                std::cout << "Generation: " << generation << std::endl;
                cout << "Estrutura final encontrada, com diferença de volume = " << bestAllFitness << endl;
                parada = true;
                /* allFitness = population[i]; */
                break;
            }
        }
        if (parada == true)
        {
            break;
        }
        std::vector<int> bestFitness = selection(allFitness); // seleciona os 4 melhores fitness e guarda seus indices aqui.

        for (int i = 0; i < 4; i++)
        {
            std::cout << bestFitness[i] << ' ';
        }
        std::cout << ' ' << std::endl;
        // Exibir informações sobre a geração atual
        std::cout << "Generation: " << generation << /* "\tBest Fitness: " << bestFitness << */ std::endl;

        // Seleção, crossover e mutação para criar a nova população
        std::vector<std::vector<int>> newPopulation;
        for (int i = 0; i < SELECTED; i++)
        {
            int indice1 = bestFitness[i];
            int indice2 = bestFitness[i + 1];
            if (indice2 == SELECTED)
            {
                indice2 = 0;
            }
            std::vector<int> parent1 = population[indice1];
            std::vector<int> parent2 = population[indice2];
            std::vector<int> child = crossover1(parent1, parent2);
            mutate(child);
            newPopulation.push_back(child);
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                std::cout << child[i] << ' ';
            }
            std::cout << ' ' << std::endl;
        }
        for (int i = 0; i < SELECTED; i++)
        {
            int indice1 = bestFitness[i];
            int indice2 = bestFitness[i + 1];
            if (indice2 == SELECTED)
            {
                indice2 = 0;
            }
            std::vector<int> parent1 = population[indice1];
            std::vector<int> parent2 = population[indice2];
            std::vector<int> child = crossover2(parent1, parent2);
            mutate(child);
            newPopulation.push_back(child);
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                std::cout << child[i] << ' ';
            }
            std::cout << ' ' << std::endl;
        }
        for (int i = 0; i < SELECTED; i++)
        {
            int indice = bestFitness[i];
            std::vector<int> actualPopulation = population[indice];
            mutate(actualPopulation);
            mutate(actualPopulation);
            newPopulation.push_back(actualPopulation);
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                std::cout << actualPopulation[i] << ' ';
            }
            std::cout << ' ' << std::endl;
        }
        // Atualização da população
        population = newPopulation;
    }
    // Obter o tempo final e calcular a diferença
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    auto segundos = duration / 1000000;
    auto minutos = segundos / 60;
    // Imprimir o tempo de execução
    std::cout << "Tempo de execucao: " << minutos.count() << " minutos" << std::endl;
    return 0;
}