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

using Band = std::vector<std::complex<double>>;
using KBand = std::map<CellVector, Band>;

struct WaveDescriptor {
  std::vector<KBand> bands;
};

struct Wavecar {
  WavecarHeader header;
  WaveDescriptor descriptor;
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

std::pair<CellVector, Band>
readWaveWaveDescriptor( const std::string &fileName
                      , const WavecarHeader &header
                      , const size_t &spinIndex
                      , const size_t &kIndex
                      ) {

  std::vector<double> realEnergies(header.nBands)
    , imagEnergies(header.nBands)
    , occupancies(header.nBands)
    ;

  double buffer, vbuffer[3];
  std::fstream file(fileName, std::ios::binary | std::ios::in);
  size_t numberPlaneWaves;
  CellVector kpoint;
  std::vector<std::complex<double>> C;

  file.seekg((spinIndex + 2) * header.recordLength);

  // read numberPlaneWaves
  file.read((char*)&buffer, sizeof(double));
  numberPlaneWaves = size_t(buffer);

  C.resize(header.nBands * numberPlaneWaves);

  file.read((char*)&kpoint, 3*sizeof(double));

  for (size_t n=0; n < header.nBands; n++) {
    file.read((char*)(realEnergies.data() + n), sizeof(double));
    file.read((char*)(imagEnergies.data() + n), sizeof(double));
    file.read((char*)(occupancies.data() + n), sizeof(double));
  }

  for (size_t n=0; n < header.nBands * numberPlaneWaves; n++) {
    file.read((char*)&C[n], 2*sizeof(double));
  }

  return {kpoint, C };

}

WavecarHeader readWavecarHeader(const std::string &fileName) {
  WavecarHeader header;
  std::fstream file(fileName, std::ios::binary | std::ios::in);
  double buffer, vbuffer[3];
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
  WaveDescriptor descriptor;

  for (uint8_t i=0; i < header.nSpin; i++) {
    KBand band;
    for (size_t k=0; k < header.nKpoints; k++) {
      auto pair(readWaveWaveDescriptor(fileName, header, i, k));
      band.emplace(pair);
    }
    descriptor.bands.push_back(band);
  }
  return {header, descriptor};

}

int main (int argc, char **argv) {
  std::cout << "WAVIS" << std::endl;
  std::cout << "=====" << std::endl;
  auto wavecar(readWavecar("WAVECAR"));
  auto& header(wavecar.header);

  std::cout << "recordLength: " << header.recordLength << "\n"
            << "nSpin: " << header.nSpin << "\n"
            << "version: " << header.version << "\n"
            << "nKpoints: " << header.nKpoints << "\n"
            << "nBands: " << header.nBands << "\n"
            << "encut: " << header.encut << "\n"
            << "eFermi: " << header.eFermi << "\n"
            << "volume: " << header.real.volume << "\n"
            << "\n";

  std::cout << "Lattice vectors: \n";
  for (const auto &b: header.real.basis)
    printf("- %f %f %f\n", b[0], b[1], b[2]);

  std::cout << "Reciprocal vectors: \n";
  for (const auto &b: header.reciprocal.basis)
    printf("- %f %f %f\n", b[0], b[1], b[2]);

}
