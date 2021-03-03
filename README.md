
# Table of Contents

-   [Wavecar.hpp](#org7bd91aa)
    -   [Introduction](#org58ffd86)
    -   [Header](#org51de487)
    -   [Main types](#org654a80e)
    -   [Mathematical functions](#org1347338)
        -   [Products](#org80ff42b)
        -   [Cell volume](#org64d43a3)
        -   [Reciprocal cell](#orgc31e112)
    -   [`WAVECAR` parsing](#org8c31a69)
        -   [Understanding the code](#orgaea8686)
        -   [Wave descriptor](#org534d65c)
        -   [Header](#org7dc2a8d)
        -   [The whole `WAVECAR`](#orgf30cbcd)
    -   [`WAVECAR` writing](#orgbac033b)
        -   [`CellVector`](#orge548fb5)
        -   [`Cell`](#orgde3601d)
        -   [`WavecarHeader`](#org59bb23f)
        -   [`WaveDescriptor`](#org2a3a593)
        -   [`Wavecar`](#orgaf414dd)
-   [Tests](#org0b53118)
    -   [Main function](#org9730958)



<a id="org7bd91aa"></a>

# Wavecar.hpp


<a id="org58ffd86"></a>

## Introduction

This package implements a C++11 header-only library
for reading `WAVECAR` files in the `VASP5` format.

It is written as a literate program <sup id="39f041f6b1d2d698620dbd1d6c83c888"><a href="#LiteratePrograKnuth1984" title="Knuth, Literate Programming, {The Computer Journal}, v(), 97--111 (1984).">LiteratePrograKnuth1984</a></sup>,
so we will try to keep the documentation and the code in the same place.

As a convenience you can find the header file in the `include` directory.


<a id="org51de487"></a>

## Header

We will only need libraries available in the standard library:

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


<a id="org654a80e"></a>

## Main types

We define the main structure for the WAVECAR header,

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


<a id="org1347338"></a>

## Mathematical functions


<a id="org80ff42b"></a>

### Products

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


<a id="org64d43a3"></a>

### Cell volume

We can use the products to calculate the cell volume given by a basis

    double cellVolume(const Cell &c) {
      return
        dotProduct( c.basis[0]
                  , crossProduct( c.basis[1]
                                , c.basis[2]
                                )
                  );
    }


<a id="orgc31e112"></a>

### Reciprocal cell

Here we calculate the reciprocal cell of a given cell

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


<a id="org8c31a69"></a>

## `WAVECAR` parsing


<a id="orgaea8686"></a>

### Understanding the code

It can be quite challenging sometimes to understand binary formats
coming from Fortran, or binary formats in general. It is quite instructive
to peek into the binary structure of the files and try to
backwards engineer the main structure of the file.

Here is an excerpt of a hexdump of a typical `VASP5` format `WAVECAR` file:

<div class="figure" id="orgaba9eba">
<pre class="example" id="orgb35ffea">
                  BYTES      
 ADDRESS  1 2  3 4  4 6  7 8  | Comments
==============================|=========
00000000: 0000 0000 0018 d940 | Fortran record length
00000008: 0000 0000 0000 f03f | number of spin channels
00000010: 0000 0000 8006 ea40 | format version (RTAG)
00000018: 0000 0000 0000 0000 |
*                               (zero padding)
00006460: 0000 0000 0000 f03f | Number of k-points
00006468: 0000 0000 0020 b940 | Number of bands
00006470: 0000 0000 00e0 8540 | ENCUT
00006478: ba49 0c02 2b07 1140 | a[1,1] (lattice vectors)
00006480: ba49 0c02 2b07 1140 | a[1,2]
00006488: 0000 0000 0000 0000 | a[1,3]
00006490: 0000 0000 0000 0000 | a[2,1]
00006498: ba49 0c02 2b07 1140 | a[2,2]
000064a0: ba49 0c02 2b07 1140 | a[2,3]
000064a8: ba49 0c02 2b07 1140 | a[3,1]
000064b0: 0000 0000 0000 0000 | a[3,2]
000064b8: ba49 0c02 2b07 1140 | a[3,3]
000064c0: 1a6f 1d53 35e9 0540 | fermi energy
000064c8: 0000 0000 0000 0000 |
*                               (zero padding)
0000c8c0: 0000 0000 0018 a940 | number of plane waves
0000c8c8: 0000 0000 0000 0000 | kpoint[0]  \
0000c8d0: 0000 0000 0000 0000 | kpoint[1]   &gt; Gamma point
0000c8d8: 0000 0000 0000 0000 | kpoint[2]  /
0000c8e0: f1e1 932b 1cee 56c0 | energy-real     x
0000c8e8: 0000 0000 0000 0000 | energy-complex  0
0000c8f0: 0000 0000 0000 f03f | occupation      1
</pre>

</div>

Here there are a couple of things we should remark,

-   `VASP5` format writes out everything using floating point numbers,
    even quantities that ought to be integers, this is done so that the
    binary format is compatible throughout machines (and most programming
    languages) since they follow IEEE standards.
-   Fortran can use an offset to read different parts of a file
    separated in records. The first quantity we get in the `WAVECAR`
    is this record length. In this particular case it is equal
    to `25696` or `0x6460` in hexadecimal notation.
    Notice that this is equal the address where the `WAVECAR`
    header begins by providing the number of \( k \)-points.
-   The second time this record length arises, is when
    the first chunk of data for the wavefunction arises, i.e.
    `2*25696 = 0xc8c0`, where we obtain for the first spin
    channel and first \( k \)-point
    -   the number of plane-waves
    -   the real part of the eigenenergies
    -   the complex part of the eigenenergies
    -   the occupation numbers
    -   the plane-wave coefficients.


<a id="org534d65c"></a>

### Wave descriptor

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


<a id="org7dc2a8d"></a>

### Header

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


<a id="orgf30cbcd"></a>

### The whole `WAVECAR`

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


<a id="orgbac033b"></a>

## `WAVECAR` writing

Our writer writes `WAVECAR` files in the `VASP5` version.


<a id="orge548fb5"></a>

### `CellVector`

    void writeToWavecar(std::ofstream &f, const CellVector &v) {
      f.write((char*)v.data(), sizeof(CellVector));
    }


<a id="orgde3601d"></a>

### `Cell`

In the case of a cell we only write the basis elements in order,

    void writeToWavecar(std::ofstream &f, const Cell &c) {
      for (const auto& v: c.basis) writeToWavecar(f, v);
    }


<a id="org59bb23f"></a>

### `WavecarHeader`

In the case of the `WavecarHeader` we have to make sure the order is the
correct one that `VASP` is expecting.

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


<a id="org2a3a593"></a>

### `WaveDescriptor`

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


<a id="orgaf414dd"></a>

### `Wavecar`

Writing a `WAVECAR` consists in writing first the header
and then the wave descriptor.

    void writeToWavecar(std::ofstream &f, const Wavecar &w) {
      writeToWavecar(f, w.header);
    
      for (size_t ispin(0); ispin < w.descriptors.size(); ispin++) {
        f.seekp((ispin + 2) * w.header.recordLength);
        writeToWavecar(f, w.descriptors[ispin]);
      }
    
    }


<a id="org0b53118"></a>

# Tests


<a id="org9730958"></a>

## Main function

    #include <Wavecar.hpp>
    
    int main (int argc, char **argv) {
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
    
      auto wavecar2(std::ofstream("WAVECAR-2", std::ios::binary));
      std::cout << "Writing WAVECAR-2" << std::endl;
      writeToWavecar(wavecar2, wavecar);
    
    }


# Bibliography
<a id="LiteratePrograKnuth1984"></a>[LiteratePrograKnuth1984] Knuth, Literate Programming, <i>The Computer Journal</i>, <b>27</b>, 97-111 (1984). <a href="http://dx.doi.org/10.1093/comjnl/27.2.97">link</a>. <a href="http://dx.doi.org/10.1093/comjnl/27.2.97">doi</a>. [↩](#39f041f6b1d2d698620dbd1d6c83c888)

