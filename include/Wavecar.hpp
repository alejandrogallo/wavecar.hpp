#ifndef _WAVECAR_HPP_DEFINED
#define _WAVECAR_HPP_DEFINED

#include <iostream>
#include <fstream>
#include <map>
#include <array>
#include <assert.h>
#include <string>
#include <numeric>
#include <algorithm>
#include <vector>
#include <math.h>
#include <complex>

using CellVector = std::array<double, 3>;

struct Cell {
  std::array<CellVector, 3> basis;
  double volume;
};

struct WavecarHeader {
  size_t recordLength;
  size_t nSpin;
  size_t version;
  size_t nKpoints;
  size_t nBands;
  size_t encut;
  double eFermi;
  size_t nPlaneWaves;
  Cell real, reciprocal;
};

using KPoint = CellVector;
using OrbitalCoefficients = std::vector<std::complex<double>>;
struct HorizontalBand {
  std::vector<double> realEnergies;
  std::vector<double> imagEnergies;
  std::vector<double> occupancies;
  std::vector<OrbitalCoefficients> coefficients;
};
using KBand = std::pair<KPoint, HorizontalBand>;

struct WaveDescriptor {
  std::vector<KBand> bands;
  inline size_t nKpoints() { return bands.size(); }
};

struct Wavecar {
  WavecarHeader header;
  std::vector<WaveDescriptor> descriptors;
};
using GGrid = std::pair<KPoint, std::vector<CellVector>>;
using HorizontalGrid = std::vector<GGrid>;
constexpr double hbarConst = 0.26246582250210965422;

CellVector crossProduct(const CellVector &a, const CellVector &b) {
  return
    { a[1] * b[2] - a[2] * b[1]
    , a[2] * b[0] - a[0] * b[2]
    , a[0] * b[1] - a[1] * b[0]
    };
}

double dotProduct(const CellVector &a, const CellVector &b) {
  return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}
double cellVolume(const Cell &c) {
  return
    dotProduct( c.basis[0]
              , crossProduct( c.basis[1]
                            , c.basis[2]
                            )
              );
}
Cell reciprocalCell(const Cell &c) {
  Cell r;
  for (size_t i=0; i<3; i++) {
    auto itB(r.basis[i].begin());
    r.basis[i] = crossProduct( c.basis[(i+1) % 3]
                             , c.basis[(i+2) % 3]);
    std::transform( itB, itB + 3, itB
                  , [&c](double i){ return i*2*M_PI/c.volume; });
  }
  r.volume = cellVolume(r);
  return r;
}

KBand
readWaveWaveDescriptor( const std::string &fileName
                      , const WavecarHeader &header
                      , const size_t &spinIndex
                      ) {

  std::vector<double> realEnergies(header.nBands)
    , imagEnergies(header.nBands)
    , occupancies(header.nBands)
    ;

  double buffer;
  std::fstream file(fileName, std::ios::binary | std::ios::in);
  size_t numberPlaneWaves;
  CellVector kpoint;
  std::vector<OrbitalCoefficients> coefficients;

  file.seekg((spinIndex + 2) * header.recordLength);

  // read numberPlaneWaves
  file.read((char*)&buffer, sizeof(double));
  numberPlaneWaves = size_t(buffer);

  //C.resize(header.nBands * numberPlaneWaves);

  file.read((char*)&kpoint, 3*sizeof(double));

  for (size_t n=0; n < header.nBands; n++) {
    file.read((char*)(realEnergies.data() + n), sizeof(double));
    file.read((char*)(imagEnergies.data() + n), sizeof(double));
    file.read((char*)(occupancies.data() + n), sizeof(double));
  }

  for (size_t n=0; n < header.nBands; n++) {
    OrbitalCoefficients C;
    C.resize(numberPlaneWaves);
    coefficients.push_back(C);
    file.read((char*)C.data(), numberPlaneWaves * 2 * sizeof(double));
  }

  return { kpoint
         , {realEnergies, imagEnergies, occupancies, coefficients}
         };

}
WavecarHeader readWavecarHeader(const std::string &fileName) {
  WavecarHeader header;
  std::fstream file(fileName, std::ios::binary | std::ios::in);
  double buffer;
  std::vector<double> vvbuffer;

  assert(sizeof(double) == 8);
  assert(sizeof(header.real.basis) == 72);
  assert(sizeof(CellVector) == 3 * sizeof(double));

  file.read((char*)&buffer, sizeof(double));
  header.recordLength = size_t(buffer);
  file.read((char*)&buffer, sizeof(double));
  header.nSpin = size_t(buffer);
  file.read((char*)&buffer, sizeof(double));
  header.version = size_t(buffer);

  if (header.version != 53300)
    throw "This program only supports VASP5 format (RTAG: 53300)";

  file.seekg(header.recordLength);

  file.read((char*)&buffer, sizeof(double));
  header.nKpoints = size_t(buffer);
  file.read((char*)&buffer, sizeof(double));
  header.nBands = size_t(buffer);
  file.read((char*)&buffer, sizeof(double));
  header.encut = size_t(buffer);

  // Setup real cell
  file.read((char*)&header.real.basis, sizeof(header.real.basis));
  header.real.volume = cellVolume(header.real);

  file.read((char*)&buffer, sizeof(double));
  header.eFermi = buffer;

  // Setup Reciprocal cell
  header.reciprocal = reciprocalCell(header.real);



  return header;
}
Wavecar readWavecar(const std::string &fileName) {
  auto header(readWavecarHeader(fileName));
  std::vector<WaveDescriptor> descriptors;

  for (size_t i=0; i < header.nSpin; i++) {
    WaveDescriptor descriptor;
    for (size_t k=0; k < header.nKpoints; k++) {
      auto kBand(readWaveWaveDescriptor(fileName, header, i));
      descriptor.bands.push_back(kBand);
    }
    descriptors.push_back(descriptor);
  }
  return {header, descriptors};

}
void writeToWavecar(std::ofstream &f, const CellVector &v) {
  f.write((char*)v.data(), sizeof(CellVector));
}
void writeToWavecar(std::ofstream &f, const Cell &c) {
  for (const auto& v: c.basis) writeToWavecar(f, v);
}
void writeToWavecar(std::ofstream &f, const WavecarHeader &h) {
  const auto writeInt
    = [&f](const size_t &i) {
        const double j(i);
        f.write((char*)&j, sizeof(double));
    };
  writeInt(h.recordLength);
  writeInt(h.nSpin);
  writeInt(h.version);

  f.seekp(h.recordLength);

  writeInt(h.nKpoints);
  writeInt(h.nBands);
  writeInt(h.encut);

  writeToWavecar(f, h.real);
  f.write((char*)&h.eFermi, sizeof(double));
}
void writeToWavecar(std::ofstream &f, const WaveDescriptor &d) {
  for (auto const& kband: d.bands)  { // spin loop
    const auto& kVector(kband.first);
    const auto& vband(kband.second);
    const double numberPlaneWaves(vband.coefficients[0].size());

    f.write((char*)&numberPlaneWaves, sizeof(double));
    writeToWavecar(f, kVector);

    for (size_t n(0); n < vband.realEnergies.size(); n++) {
      f.write((char*)&vband.realEnergies[n], sizeof(double));
      f.write((char*)&vband.imagEnergies[n], sizeof(double));
      f.write((char*)&vband.occupancies[n], sizeof(double));
    }

    for (size_t n(0); n < vband.realEnergies.size(); n++) {
      f.write((char*)vband.coefficients[n].data(),
              numberPlaneWaves * 2 * sizeof(double));
    }

  }
}
void writeToWavecar(std::ofstream &f, const Wavecar &w) {
  writeToWavecar(f, w.header);

  for (size_t ispin(0); ispin < w.descriptors.size(); ispin++) {
    f.seekp((ispin + 2) * w.header.recordLength);
    writeToWavecar(f, w.descriptors[ispin]);
  }

}
std::vector<HorizontalGrid>
mkGrid(const Wavecar &wavecar) {
  std::vector<HorizontalGrid> result;

  for (const auto& descriptor: wavecar.descriptors) {

    HorizontalGrid hGrid;

    for (const auto& kBand: descriptor.bands) {


    size_t count(0);
    const double actualNumberOfPlaneWaves = kBand.second.coefficients[0].size();

    const auto& K(kBand.first);
    std::vector<CellVector> gs(actualNumberOfPlaneWaves);

    const int numberOfPlaneWaves
      = 3 * std::ceil(std::pow(actualNumberOfPlaneWaves, 1.0/3));

    for (int z(0); z < 2 * numberOfPlaneWaves; z++) {
    for (int y(0); y < 2 * numberOfPlaneWaves; y++) {
    for (int x(0); x < 2 * numberOfPlaneWaves; x++) {

      const auto& cell(wavecar.header.reciprocal);
      const int G_i[]
        = { x > numberOfPlaneWaves ? x - 2 * numberOfPlaneWaves : x
          , y > numberOfPlaneWaves ? y - 2 * numberOfPlaneWaves : y
          , z > numberOfPlaneWaves ? z - 2 * numberOfPlaneWaves : z
          };

      double energy(0);
      CellVector g;
      for (size_t i(0); i < 3; i++) {
        double component
          = cell.basis[0][i] * ((double)G_i[0])
          + cell.basis[1][i] * ((double)G_i[1])
          + cell.basis[2][i] * ((double)G_i[2])
          ;
        component += K[i];

        g[i] = component;
        energy += component * component / hbarConst;
      }

      if (energy < wavecar.header.encut) {
        count++;
        gs[count] = g;
      }

    } // z
    } // y
    } // x

    if (count != actualNumberOfPlaneWaves)
      throw "Count and actualNumberOfPlaneWaves are different";

    hGrid.push_back({K, gs});

  } // kBand

    result.push_back(hGrid);

  } // descriptor

  return result;

}

#endif
