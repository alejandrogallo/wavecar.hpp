#include <Wavecar.hpp>

int main (int argc, char **argv) {

  std::cout << ">>> Parsing WAVECAR" << std::endl;
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

  std::cout << ">>> Creating grids" << std::endl;
  auto grids(mkGrid(wavecar));
  std::cout << grids.size() << " grids created" << std::endl;

  std::cout << "Lattice vectors: \n";
  for (const auto &b: header.real.basis)
    printf("  %f %f %f\n", b[0], b[1], b[2]);

  std::cout << "Reciprocal vectors: \n";
  for (const auto &b: header.reciprocal.basis)
    printf("  %f %f %f\n", b[0], b[1], b[2]);

  auto wavecar2(std::ofstream("WAVECAR-2", std::ios::binary));
  std::cout << "Writing WAVECAR-2" << std::endl;
  writeToWavecar(wavecar2, wavecar);

}
