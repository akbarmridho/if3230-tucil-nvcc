# Tugas Kecil - Paralel Inverse Matrix

Tugas Kecil IF3230 Sistem Paralel dan Terdistribusi

- Rachel Gabriela Chen (13521044)
- Jeffrey Chow (13521046)
- Akbar Maulana Ridho (13521093)

## How to Run

Jalankan perintah `make` pada root folder repository. 

### Serial

Manual

`cat ./test_cases/{test_case}.txt | ./bin/serial`

Script

`sh ./test_serial.sh`

### OpenMPI

Manual

`cat ./test_cases/{test_case}.txt | mpirun -n {n_core} ./bin/mpi`

Script

`sh ./test_mpi.sh`