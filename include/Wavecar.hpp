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
struct VerticalBand {
  std::vector<double> realEnergies;
  std::vector<double> imagEnergies;
  std::vector<double> occupancies;
  std::vector<OrbitalCoefficients> coefficients;
};
using KBand = std::pair<KPoint, VerticalBand>;

struct WaveDescriptor {
  std::vector<KBand> bands;
  inline size_t nKpoints() { return bands.size(); }
};

struct Wavecar {
  WavecarHeader header;
  std::vector<WaveDescriptor> descriptors;
};

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

double norm(const CellVector &a) {
 return std::sqrt(dotProduct(a, a));
}

double angleBetween(const CellVector &a, const CellVector &b) {
  return std::acos(dotProduct(a, b) / norm(a) / norm(b));
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
  // const double hbarConst = 0.26246582250210965422; // 1/eV Ang^2

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

  for (uint8_t i=0; i < header.nSpin; i++) {
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
